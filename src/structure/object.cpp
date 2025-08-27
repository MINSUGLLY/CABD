#include "object.hpp"

namespace labd{

Object::Object(const std::string& configPath){
    std::ifstream file(configPath);
    nlohmann::json j;
    file >> j;

    // parameter
    rho = j["parameter"]["rho"];
    nv = j["parameter"]["nv"];
    kappa = j["parameter"]["kappa"];
    soft = j["parameter"]["soft"];

    // mesh
    utils::loadNpy(meshTV, j["mesh"]["pathV"]);
    utils::loadNpy(meshTT, j["mesh"]["pathT"]);
    utils::loadNpy(meshTF, j["mesh"]["pathF"]);

    // origin
    vec3 xyz = utils::jsonV3(j["origin"]["xyz"]);
    vec3 rpy = utils::jsonV3(j["origin"]["rpy"]);
    scalar scale = j["origin"]["scale"];

    meshTV = (utils::rpy2R(rpy) * scale * meshTV.transpose()).transpose().rowwise() + xyz.transpose();

    // linearMass
    generateLinearHessianInv();

    // force
    vec3 forceE = vec3::Zero();
    vec3 torqueE = vec3::Zero();

    // FEM
    if(soft){
        utils::loadNpy(fixVs, j["fixP"]);
        FEMprecompute();
    }

    // collision
    localAABBtree.init(meshTV, meshTF);
}

void Object::generateLinearHessianInv(){
    mat6 Cm = mat6::Zero();
    Cm.block<3, 3>(0, 0) = (mat3() << 1 - nv, nv, nv, nv, 1 - nv, nv, nv, nv, 1 - nv).finished();
    Cm.block<3, 3>(3, 3) = 0.5 * (1 - 2 * nv) * mat3::Identity();
    Cm *= kappa / ((1 + nv) * (1 - 2 * nv));

    Eigen::Matrix<scalar, 3, 4> nablaNlocal = (Eigen::Matrix<scalar, 3, 4>() << -1.0, 1.0, 0.0, 0.0,
                                                                                -1.0, 0.0, 1.0, 0.0,
                                                                                -1.0, 0.0, 0.0, 1.0).finished();

    Eigen::Matrix<scalar, 3, 4> nablaN;
    Eigen::Matrix<scalar, 6, 12> B;

    Mtet.setZero();

    vec3 vo = vec3(0.0, 0.0, 0.0);
    vec3 vx = vec3(1.0, 0.0, 0.0);
    vec3 vy = vec3(0.0, 1.0, 0.0);
    vec3 vz = vec3(0.0, 0.0, 1.0);
    scalar volume = std::abs((vx - vo).dot((vy - vo).cross(vz - vo))) / 6.0;

    mat3 J;
    mat3 Jinv;
    J.col(0) = vx - vo;
    J.col(1) = vy - vo;
    J.col(2) = vz - vo;
    Jinv = J.inverse();
    nablaN = Jinv.transpose() * nablaNlocal;

    B.block<3, 3>(0, 0) = nablaN.col(0).asDiagonal();
    B.block<3, 3>(0, 3) = nablaN.col(1).asDiagonal();
    B.block<3, 3>(0, 6) = nablaN.col(2).asDiagonal();
    B.block<3, 3>(0, 9) = nablaN.col(3).asDiagonal();

    B.block<3, 3>(3, 0) = BMatrixLower(nablaN.col(0));
    B.block<3, 3>(3, 3) = BMatrixLower(nablaN.col(1));
    B.block<3, 3>(3, 6) = BMatrixLower(nablaN.col(2));
    B.block<3, 3>(3, 9) = BMatrixLower(nablaN.col(3));

    Ktet = volume * B.transpose() * Cm * B;

    mat12 A;
    A.setZero();
    A.block<3, 3>(0, 0) = mat3::Identity();
    A.block<3, 3>(3, 0) = mat3::Identity();
    A.block<3, 3>(6, 0) = mat3::Identity();
    A.block<3, 3>(9, 0) = mat3::Identity();
    A.block<3, 3>(3, 3) = mat3::Identity();
    A.block<3, 3>(6, 6) = mat3::Identity();
    A.block<3, 3>(9, 9) = mat3::Identity();

    Ktet = A.transpose() * Ktet * A;

    vec3 v0temp;
    vec3 v1temp;
    vec3 v2temp;
    vec3 v3temp;
    scalar tetV;
    for(int tetIdx = 0; tetIdx < meshTT.rows(); tetIdx++){
        v0temp = meshTV.row(meshTT(tetIdx, 0));
        v1temp = meshTV.row(meshTT(tetIdx, 1));
        v2temp = meshTV.row(meshTT(tetIdx, 2));
        v3temp = meshTV.row(meshTT(tetIdx, 3));
        tetV = std::abs((v1temp - v0temp).dot((v2temp - v0temp).cross(v3temp - v0temp))) / 6.0;

        massTotal += tetV * rho;
        massCenter += v0temp * tetV / 4 * rho;
        massCenter += v1temp * tetV / 4 * rho;
        massCenter += v2temp * tetV / 4 * rho;
        massCenter += v3temp * tetV / 4 * rho;
    }
    massCenter /= massTotal;

    for(int verIdx = 0; verIdx < meshTV.rows(); verIdx++){
        meshTV.row(verIdx) -= massCenter;
    }

    for(int tetIdx = 0; tetIdx < meshTT.rows(); tetIdx++){
        v0temp = meshTV.row(meshTT(tetIdx, 0));
        v1temp = meshTV.row(meshTT(tetIdx, 1));
        v2temp = meshTV.row(meshTT(tetIdx, 2));
        v3temp = meshTV.row(meshTT(tetIdx, 3));
        tetV = std::abs((v1temp - v0temp).dot((v2temp - v0temp).cross(v3temp - v0temp))) / 6.0;

        matJ J;
        utils::getJ(v0temp, J);
        Mtet += rho * tetV / 4.0 * J.transpose() * J;
        utils::getJ(v1temp, J);
        Mtet += rho * tetV / 4.0 * J.transpose() * J;
        utils::getJ(v2temp, J);
        Mtet += rho * tetV / 4.0 * J.transpose() * J;
        utils::getJ(v3temp, J);
        Mtet += rho * tetV / 4.0 * J.transpose() * J;
    }

    // cube 
    // Mtet.setZero();
    // Mtet.block<3, 3>(0, 0).setIdentity();
    // Mtet.block<9, 9>(3, 3) = 1.0 / 12.0 * mat9::Identity();
    // spdlog::warn("use cube");

    // ball
    // Mtet.setZero();
    // Mtet.block<3, 3>(0, 0) = pi / 6.0 * mat3::Identity();;
    // Mtet.block<9, 9>(3, 3) = pi / 120.0 * mat9::Identity();
    // spdlog::warn("use ball");

    // rod
    // Mtet.setZero();
    // Mtet.block<3, 3>(0, 0) = 10.0 * mat3::Identity();
    // Mtet.block<9, 9>(3, 3) = 5.0 / 6.0 * mat9::Identity();
    // Mtet.block<3, 3>(6, 6) = 250.0 / 3.0 * mat3::Identity();
    // spdlog::warn("use rod");

    // if(soft){
    //     Mtet *= 1.0e-1;
    // }
    linearHessian = Mtet + Ktet;
    linearHessianInv = linearHessian.inverse();
}

void Object::FEMprecompute(){
    // surface points number
    std::set<int> Pset;
    for(int i = 0; i < meshTF.rows(); i++){
        Pset.insert(meshTF(i, 0));
        Pset.insert(meshTF(i, 1));
        Pset.insert(meshTF(i, 2));
    }
    nPsurface = Pset.size();

    // FEM
    meshTVtemp = meshTV;

    meshTViter.resize(meshTV.rows(), 3);
    meshTViter.setZero();
    meshTVtilde.resize(meshTV.rows(), 3);
    meshTVtilde.setZero();
    meshTVforce.resize(meshTV.rows(), 3);
    meshTVforce.setZero();
    meshTVdot.resize(meshTV.rows(), 3);
    meshTVdot.setZero();
    meshTVacc.resize(meshTV.rows(), 3);
    meshTVacc.setZero();

    FEMM.resize(meshTV.rows());
    FEMM.setZero();
    FEMK.resize(3 * meshTV.rows(), 3 * meshTV.rows());
    FEMK.setZero();

    vec3 v0temp;
    vec3 v1temp;
    vec3 v2temp;
    vec3 v3temp;
    mat3 J;
    mat3 Jinv;
    scalar volume = 0.0;

    mat6 Cm = mat6::Zero();
    Cm.block<3, 3>(0, 0) = (mat3() << 1 - nv, nv, nv, nv, 1 - nv, nv, nv, nv, 1 - nv).finished();
    Cm.block<3, 3>(3, 3) = 0.5 * (1 - 2 * nv) * mat3::Identity();
    Cm *= kappa / ((1 + nv) * (1 - 2 * nv));

    Eigen::Matrix<scalar, 3, 4> nablaNlocal = (Eigen::Matrix<scalar, 3, 4>() << -1.0, 1.0, 0.0, 0.0,
                                                                                -1.0, 0.0, 1.0, 0.0,
                                                                                -1.0, 0.0, 0.0, 1.0).finished();
    Eigen::Matrix<scalar, 3, 4> nablaN;
    Eigen::Matrix<scalar, 6, 12> B;
    mat12 Ktet;

    for(int tetIdx = 0; tetIdx < meshTT.rows(); tetIdx++){
        v0temp = meshTV.row(meshTT(tetIdx, 0));
        v1temp = meshTV.row(meshTT(tetIdx, 1));
        v2temp = meshTV.row(meshTT(tetIdx, 2));
        v3temp = meshTV.row(meshTT(tetIdx, 3));
        volume = std::abs((v1temp - v0temp).dot((v2temp - v0temp).cross(v3temp - v0temp))) / 6.0;

        FEMM[meshTT(tetIdx, 0)] += rho * volume / 4.0;
        FEMM[meshTT(tetIdx, 1)] += rho * volume / 4.0;
        FEMM[meshTT(tetIdx, 2)] += rho * volume / 4.0;
        FEMM[meshTT(tetIdx, 3)] += rho * volume / 4.0;

        J.col(0) = v1temp - v0temp;
        J.col(1) = v2temp - v0temp;
        J.col(2) = v3temp - v0temp;
        Jinv = J.inverse();
        nablaN = Jinv.transpose() * nablaNlocal;

        B.block<3, 3>(0, 0) = nablaN.col(0).asDiagonal();
        B.block<3, 3>(0, 3) = nablaN.col(1).asDiagonal();
        B.block<3, 3>(0, 6) = nablaN.col(2).asDiagonal();
        B.block<3, 3>(0, 9) = nablaN.col(3).asDiagonal();

        B.block<3, 3>(3, 0) = BMatrixLower(nablaN.col(0));
        B.block<3, 3>(3, 3) = BMatrixLower(nablaN.col(1));
        B.block<3, 3>(3, 6) = BMatrixLower(nablaN.col(2));
        B.block<3, 3>(3, 9) = BMatrixLower(nablaN.col(3));

        Ktet = volume * B.transpose() * Cm * B;
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                FEMK.block<3, 3>(3 * meshTT(tetIdx, i), 3 * meshTT(tetIdx, j)) += Ktet.block<3, 3>(3 * i, 3 * j);
            }
        } 
    }

    matX denseFEMM(3 * meshTV.rows(), 3 * meshTV.rows());
    for(int vIdx = 0; vIdx < meshTV.rows(); vIdx++){
        denseFEMM.block<3, 3>(3 * vIdx, 3 * vIdx) = FEMM[vIdx] * Eigen::Matrix3d::Identity();
    }

    FEMhessian = FEMK + denseFEMM;

    for(int fixIdx = 0; fixIdx < fixVs.rows(); fixIdx++){
        FEMhessian.block(0, 3 * fixVs(fixIdx, 0), 3 * meshTV.rows(), 3).setZero();
        FEMhessian.block(3 * fixVs(fixIdx, 0), 0, 3, 3 * meshTV.rows()).setZero();
        FEMhessian.block<3, 3>(3 * fixVs(fixIdx, 0), 3 * fixVs(fixIdx, 0)).setIdentity();
    }

    // condensation
    matX Kmm = FEMhessian.block(0, 0, 3 * nPsurface, 3 * nPsurface);
    matX Kss = FEMhessian.block(3 * nPsurface, 3 * nPsurface, 3 * meshTV.rows() - 3 * nPsurface, 3 * meshTV.rows() - 3 * nPsurface);
    matX Kms = FEMhessian.block(0, 3 * nPsurface, 3 * nPsurface, 3 * meshTV.rows() - 3 * nPsurface);
    matX Ksm = FEMhessian.block(3 * nPsurface, 0, 3 * meshTV.rows() - 3 * nPsurface, 3 * nPsurface);

    Kssinv = Kss.inverse();
    KssinvKsm = Kssinv * Ksm;
    KmsKssinv = Kms * Kssinv;
    FEMhessianCondInv = (Kmm - Kms * Kssinv * Ksm).inverse();

    FEMgradCond.resize(3 * nPsurface);
    FEMgradCond.setZero();

    spdlog::info("fix points: {} / {}   condensation points: {} / {}    tet: {}", fixVs.rows(), meshTV.rows(), nPsurface, meshTV.rows(), meshTT.rows());
}

mat3 Object::BMatrixLower(const vec3& src){
    mat3 result = (mat3() << src[1], src[0], 0.0,
                             0.0, src[2], src[1],
                             src[2], 0.0, src[0]).finished();
    return result;
}

}