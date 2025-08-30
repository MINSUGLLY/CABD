#include "utils.hpp"

namespace CABD{
namespace utils{

void loadNpy(matX& Mat, const std::string& path){
    cnpy::NpyArray npy = cnpy::npy_load(path);
    size_t& row = npy.shape.front();
    size_t& col = npy.shape.back();
    scalar* dat = const_cast<scalar*>(npy.data<scalar>());
    Mat = Eigen::Map<matX>(dat, row, col);
}

void loadNpy(matXi& Mat, const std::string& path){
    cnpy::NpyArray npy = cnpy::npy_load(path);
    size_t& row = npy.shape.front();
    size_t& col = npy.shape.back();
    int* dat = const_cast<int*>(npy.data<int>());
    Mat = Eigen::Map<matXi>(dat, row, col);
}

void saveNpy(const matX& Mat, const std::string& path){
    std::vector<size_t> shape = {static_cast<size_t>(Mat.rows()), static_cast<size_t>(Mat.cols())};
    Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp = Mat;
    cnpy::npy_save(path, tmp.data(), shape, "w");
}

vec3 jsonV3(const nlohmann::json& src){
    return vec3(src[0], src[1], src[2]);
}

mat3 jsonM3(const nlohmann::json& src){
    return (mat3() << src[0][0], src[0][1], src[0][2],
                      src[1][0], src[1][1], src[1][2],
                      src[2][0], src[2][1], src[2][2]).finished();
}

void q2Rp(const Eigen::Ref<const vec12> q_src, mat3& R_dst, vec3& p_dst){
    p_dst = q_src.block<3, 1>(0, 0);
    R_dst.col(0) = q_src.block<3, 1>(3, 0);
    R_dst.col(1) = q_src.block<3, 1>(6, 0);
    R_dst.col(2) = q_src.block<3, 1>(9, 0);
}

void Rp2q(Eigen::Ref<vec12> q_dst, const mat3& R_src, const vec3& p_src){
    q_dst.block<3, 1>(0, 0) = p_src;
    q_dst.block<3, 1>(3, 0) = R_src.col(0);
    q_dst.block<3, 1>(6, 0) = R_src.col(1);
    q_dst.block<3, 1>(9, 0) = R_src.col(2);
}

mat3 rpy2R(const vec3& src){
    scalar roll = src[0];
    scalar pitch = src[1];
    scalar yaw = src[2];

    mat3 Rx;
    Rx << 1, 0, 0,
          0, cos(roll), -sin(roll),
          0, sin(roll), cos(roll);

    mat3 Ry;
    Ry << cos(pitch), 0, sin(pitch),
          0, 1, 0,
          -sin(pitch), 0, cos(pitch);

    mat3 Rz;
    Rz << cos(yaw), -sin(yaw), 0,
          sin(yaw), cos(yaw), 0,
          0, 0, 1;

    mat3 R = Rz * Ry * Rx;

    return R;
}

vec3 str2vec3(const std::string& src){
    std::istringstream iss(src);
    scalar x, y, z;

    if (!(iss >> x >> y >> z)) {
        throw std::runtime_error("Failed to parse string");
    }  

    return vec3(x, y, z);
}

std::string padWithZeros(const std::string& numStr, int totalLength){
    if (numStr.length() >= totalLength) {
        return numStr;
    }
    int zerosToAdd = totalLength - numStr.length();
    std::string padding(zerosToAdd, '0');
    return padding + numStr;
}

std::string path2filename(const std::string& path){
    std::filesystem::path filePath(path);
    return filePath.filename().string();
}

void q2glm(const Eigen::Ref<vec12>& q_src, glm::mat4& Tdst){
    Tdst[0][0] = q_src[3]; Tdst[0][1] = q_src[4]; Tdst[0][2] = q_src[5];  Tdst[0][3] = 0.0;
    Tdst[1][0] = q_src[6]; Tdst[1][1] = q_src[7]; Tdst[1][2] = q_src[8]; Tdst[1][3] = 0.0;
    Tdst[2][0] = q_src[9]; Tdst[2][1] = q_src[10]; Tdst[2][2] = q_src[11]; Tdst[2][3] = 0.0;
    Tdst[3][0] = q_src[0]; Tdst[3][1] = q_src[1]; Tdst[3][2] = q_src[2];  Tdst[3][3] = 1.0;
}

mat3 mat3PolarDecomposition(const mat3& src){
    Eigen::JacobiSVD<mat3> svd(src, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    mat3 U = svd.matrixU();
    mat3 V = svd.matrixV();

    return U * V.transpose();
}

void getJ(const vec3& src, matJ& dst){
    dst.block<3, 3>(0, 0).setIdentity();
    dst.block<3, 3>(0, 3) = mat3::Identity() * src[0];
    dst.block<3, 3>(0, 6) = mat3::Identity() * src[1];
    dst.block<3, 3>(0, 9) = mat3::Identity() * src[2];
}

void getJR(const vec3& src, const mat3& R, matJ& dst){
    dst.block<3, 3>(0, 0) = R;
    dst.block<3, 3>(0, 3) = R * src[0];
    dst.block<3, 3>(0, 6) = R * src[1];
    dst.block<3, 3>(0, 9) = R * src[2];
}

scalar R2angle(const mat3& Raxis, const vec3& axis){
    scalar cosR = (Raxis.trace() - 1.0) / 2.0;
    cosR = std::min(cosR, 1.0);
    cosR = std::max(cosR, -1.0);
    
    scalar angle = std::acos(cosR);
    
    mat3 Raxis_temp = angle2R(axis, angle);
    mat3 m = Raxis - Raxis_temp;
    
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            if(std::fabs(m(i, j)) > 1.0e-3){
                return -angle;
            }
        }
    }
    return angle;
}

mat3 angle2R(const vec3& axis, const scalar angle){
    mat3 skewAxis = vec2skewSym(axis);
    mat3 R_axis = mat3::Identity() + sin(angle) * skewAxis + (1 - cos(angle)) * skewAxis * skewAxis;
    return R_axis;
}

mat3 vec2skewSym(const vec3& src){
    mat3 result;
    result.setZero();
    result(0, 1) = -src[2];
    result(0, 2) = src[1];
    result(1, 0) = src[2];
    result(1, 2) = -src[0];
    result(2, 0) = -src[1];
    result(2, 1) = src[0];

    return result;
}

vec3 skewSym2vec(const mat3& src){
    return vec3(src(2, 1), src(0, 2), src(1, 0));
}

bool AABBIntersectWithTolerance(const Eigen::AlignedBox3d& box0, const Eigen::AlignedBox3d& box1, scalar tolerance){
    Eigen::AlignedBox3d box0Expanded = box0;
    box0Expanded.min() -= vec3::Constant(tolerance);
    box0Expanded.max() += vec3::Constant(tolerance);

    return box0Expanded.intersects(box1);
}

}
}