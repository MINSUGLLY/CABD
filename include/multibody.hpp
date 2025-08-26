#pragma once

#include "object.hpp"

namespace labd{

class Joint{
public:
    Joint(int idx0Input, int idx1Input, const vec3& pos0Input, const vec3& pos1Input);

    int jointId;
    void setIdx(int idx){jointId = idx;}

    int idx0;
    int idx1;

    vec3 pos0;
    vec3 pos1;
};

class Multibody{
public:
    Multibody(const std::string& csvPath);
    
    int idx;
    void setIdx(int idxInput){idx = idxInput;}

    std::vector<std::shared_ptr<Object>> linksList;
    std::vector<std::shared_ptr<Joint>> jointsList;
    void appendLink(std::shared_ptr<Object> object, const vec3 initp = vec3::Zero(), const mat3 initR = mat3::Identity());
    void appendJoint(std::shared_ptr<Joint> joint);

    std::vector<int> EEID;
    std::vector<vec12> qt;
};

}