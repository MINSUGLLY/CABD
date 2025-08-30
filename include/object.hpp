#pragma once

#include "utils.hpp"

namespace CABD{

class Object{
public:
    // initialize
    Object(const std::string& configPath);
    Object() = default;
    Object(std::shared_ptr<Object> src);

    // mesh
    matX meshTV;
    matXi meshTT;
    matXi meshTF;

    // idx
    int idx;
    void setIdx(int idxInput){idx = idxInput;}

    int meshIdx;
    void setMeshIdx(int meshIdxInput){meshIdx = meshIdxInput;};

    // parameter
    scalar rho;
    scalar nv;
    scalar kappa;

    // status
    vec12 q;
    vec12 qdot = vec12::Zero();
    vec12 qacc = vec12::Zero();
    vec12 qtilde = vec12::Zero();
    vec12 qtemp = vec12::Zero();
    vec12 dq = vec12::Zero();

    // linearABD
    mat12 Mtet;
    mat12 Ktet;
    mat12 linearHessian;
    mat12 linearHessianInv;
    scalar massTotal = 0.0;
    vec3 massCenter = vec3::Zero();
    void generateLinearHessianInv();

    // force
    vec3 forceE = vec3::Zero();
    vec3 torqueE = vec3::Zero();

    bool fixed = false;
    std::vector<int> begList;
    vec3 pt;
    mat3 Rt;
    vec3 pd;
    mat3 Rd;
    vec3 dp;
    mat3 dR;

    // FEM
    void FEMprecompute();

    bool soft = false;
    int nPsurface = 0;
    matX meshTVtemp;
    matX meshTViter;
    matX meshTVtilde;
    matX meshTVforce;
    matX meshTVdot;
    matX meshTVacc;

    vecX FEMM;
    matX FEMK;
    matX FEMhessian;

    matX Kssinv;
    matX KssinvKsm;
    matX KmsKssinv;
    vecX KssinvFs;

    matX FEMhessianCondInv;
    vecX FEMgradCond;

    matXi fixVs;

    // collision
    igl::AABB<Eigen::MatrixXd, 3> localAABBtree;
    
private:
    mat3 BMatrixLower(const vec3& src);
};

}