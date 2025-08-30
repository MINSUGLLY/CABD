#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <map>
#include <stack>
#include <algorithm>
#include <unordered_set>
#include <memory>
#include <chrono>
#include <thread>
#include <mutex>

#include <CLI/CLI.hpp>

#include <spdlog/spdlog.h>

#include <mkl.h>
#define EIGEN_USE_MKL_ALL
#define EIGEN_DONT_PARALLELIZE

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/AABB.h>

#include <polyscope/polyscope.h>
#include <polyscope/messages.h>
#include <polyscope/surface_mesh.h>

#include <cnpy.h>

#include <nlohmann/json.hpp>

#include <ipc/ipc.hpp>

namespace CABD{

using scalar = double;

using vec2 = Eigen::Vector<scalar, 2>;
using vec3 = Eigen::Vector<scalar, 3>;
using vec4 = Eigen::Vector<scalar, 4>;
using vec6 = Eigen::Vector<scalar, 6>;
using vec9 = Eigen::Vector<scalar, 9>;
using vec12 = Eigen::Vector<scalar, 12>;
using vecX = Eigen::VectorX<scalar>;

using vec2i = Eigen::Vector<int, 2>;
using vec3i = Eigen::Vector<int, 3>;
using vec4i = Eigen::Vector<int, 4>;
using vec6i = Eigen::Vector<int, 6>;
using vec9i = Eigen::Vector<int, 9>;
using vec12i = Eigen::Vector<int, 12>;
using vecXi = Eigen::VectorXi;

using mat2 = Eigen::Matrix<scalar, 2, 2>;
using mat3 = Eigen::Matrix<scalar, 3, 3>;
using mat4 = Eigen::Matrix<scalar, 4, 4>;
using mat6 = Eigen::Matrix<scalar, 6, 6>;
using mat9 = Eigen::Matrix<scalar, 9, 9>;
using mat12 = Eigen::Matrix<scalar, 12, 12>;
using matJ = Eigen::Matrix<scalar, 3, 12>;
using matX = Eigen::MatrixX<scalar>;
using matXi = Eigen::MatrixXi;

const scalar pi = 3.141592653589793;

enum class simType{
    forward = 0,
    backward = 1,
};

namespace utils{

template<typename T> bool InputScalar(const char* label, T* value){
    if constexpr (std::is_same_v<T, float>) {
        return ImGui::InputFloat(label, value);
    } else if constexpr (std::is_same_v<T, scalar>) {
        return ImGui::InputDouble(label, value);
    } else {
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported scalar type");
        return false;
    }
}

template<typename T> bool InputScalar2(const char* label, T* value){
    if constexpr (std::is_same_v<T, float>) {
        return ImGui::InputFloat2(label, value);
    } else if constexpr (std::is_same_v<T, scalar>) {
        return ImGui::InputDouble2(label, value);
    } else {
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported scalar type");
        return false;
    }
}

template<typename T> bool InputScalar3(const char* label, T* value){
    if constexpr (std::is_same_v<T, float>) {
        return ImGui::InputFloat3(label, value);
    } else if constexpr (std::is_same_v<T, scalar>) {
        return ImGui::InputDouble3(label, value);
    } else {
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported scalar type");
        return false;
    }
}

void loadNpy(matX& Mat, const std::string& path);
void loadNpy(matXi& Mat, const std::string& path);
void saveNpy(const matX& Mat, const std::string& path);

vec3 jsonV3(const nlohmann::json& src);
mat3 jsonM3(const nlohmann::json& src);

mat3 rpy2R(const vec3& src);

void q2Rp(const Eigen::Ref<const vec12> q_src, mat3& R_dst, vec3& p_dst);
void Rp2q(Eigen::Ref<vec12> q_dst, const mat3& R_src, const vec3& p_src);

vec3 str2vec3(const std::string& src);
std::string padWithZeros(const std::string& numStr, int totalLength);
std::string path2filename(const std::string& path);

void q2glm(const Eigen::Ref<vec12>& q_src, glm::mat4& Tdst);

mat3 mat3PolarDecomposition(const mat3& src);

void getJ(const vec3& src, matJ& dst);
void getJR(const vec3& src, const mat3& R, matJ& dst);

mat3 angle2R(const vec3& axis, const scalar angle);
scalar R2angle(const mat3& Raxis, const vec3& axis);

mat3 vec2skewSym(const vec3& src);
vec3 skewSym2vec(const mat3& src);

bool AABBIntersectWithTolerance(const Eigen::AlignedBox3d& box0, const Eigen::AlignedBox3d& box1, scalar tolerance);
}
}