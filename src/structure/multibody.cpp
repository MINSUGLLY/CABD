#include "multibody.hpp"

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "tinygltf/tiny_gltf.h"

namespace CABD{

Joint::Joint(int idx0Input, int idx1Input, const vec3& pos0Input, const vec3& pos1Input){
    idx0 = idx0Input;
    idx1 = idx1Input;
    pos0 = pos0Input;
    pos1 = pos1Input;
}

Multibody::Multibody(const std::string& configPath, std::map<std::string, std::shared_ptr<Object>>& geoCache){
    std::ifstream file(configPath);
    nlohmann::json j;
    file >> j;

    int linkCount = 0;
    for(const auto& linkConfig : j["links"]){
        spdlog::info("load link config {} / {}", ++linkCount, j["links"].size());
        mat3 initR = utils::rpy2R(utils::jsonV3(linkConfig["rpy"]));
        vec3 initp = utils::jsonV3(linkConfig["xyz"]);
        if(geoCache.find(linkConfig["configPath"]) != geoCache.end()){
            spdlog::info("use cached geometry {}", linkConfig["configPath"]);
            appendLink(std::make_shared<Object>(geoCache[linkConfig["configPath"]]), initp, initR);
        }
        else{
            appendLink(std::make_shared<Object>(linkConfig["configPath"]), initp, initR);
            geoCache[linkConfig["configPath"]] = linksList.back();
        }
        linksList.back()->fixed = linkConfig["static"];
        if(linksList.back()->fixed){
            utils::q2Rp(linksList.back()->q, linksList.back()->Rt, linksList.back()->pt);
            utils::q2Rp(linksList.back()->q, linksList.back()->Rd, linksList.back()->pd);
        }
    }

    for(const auto& jointConfig : j["joints"]){
        int idx0 = jointConfig["obj0"];
        int idx1 = jointConfig["obj1"];
        vec3 pos0 = utils::jsonV3(jointConfig["pos0"]);
        vec3 pos1 = utils::jsonV3(jointConfig["pos1"]);
        appendJoint(std::make_shared<Joint>(idx0, idx1, pos0, pos1));
    }
}

void Multibody::appendLink(std::shared_ptr<Object> object, const vec3 initp, const mat3 initR){
    object->setIdx(linksList.size());
    utils::Rp2q(object->q, initR, initp + initR * object->massCenter);
    linksList.push_back(object);
}

void Multibody::appendJoint(std::shared_ptr<Joint> joint){
    joint->setIdx(jointsList.size());
    joint->pos0 += -linksList[joint->idx0]->massCenter;
    joint->pos1 += -linksList[joint->idx1]->massCenter;
    jointsList.push_back(joint);
    linksList[joint->idx0]->begList.push_back(3 * joint->jointId);
    linksList[joint->idx1]->begList.push_back(3 * joint->jointId);
}

}

