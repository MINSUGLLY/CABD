#include "simulator.hpp"

namespace labd{

void Simulator::initViewer(){
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::XFront);

    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::state::lengthScale = 10.0;
    polyscope::state::boundingBox = std::tuple<glm::vec3, glm::vec3>{ {-1.0, -1.0, 0.0}, {1.0, 1.0, 1.0} };

    polyscope::options::alwaysRedraw = true;

    polyscope::options::maxFPS = 60;
    polyscope::options::enableVSync = true;

    polyscope::state::userCallback = [&](){
        {   
            {
                std::lock_guard<std::mutex> lock(mtx);
                updateMeshTransform();
            }

            // UI
            if(ImGui::Button("FORWARD")){
                type = simType::forward;
                spdlog::info("FROWARD MODE");
            }
            ImGui::SameLine();
            if(ImGui::Button("BACKWARD")){
                type = simType::backward;
                spdlog::info("BACKWARD MODE");
            }

            if(type == simType::forward){
                if(ImGui::InputInt("frames", &framesTotal)){
                    framesCount = framesTotal;
                };
                if(ImGui::Button("MOVE")){
                    framesCount = 0;
                    
                    // update move information
                    for(const auto& multibody : multibodiesList){
                        for(const auto& link : multibody->linksList){
                            if(link->fixed){
                                vec3 p;
                                mat3 R;
                                utils::q2Rp(link->q, R, p);
                                link->dp = (link->pd - p) / framesTotal;
                            }
                        }
                    }
                }

                utils::InputScalar3("gravity", gravity.data());

                if(ImGui::TreeNode("external force")){
                    for(const auto& object : objsList){
                        if(ImGui::TreeNode((std::string("force obj") + utils::padWithZeros(std::to_string(object->idx), 4)).c_str())){
                            utils::InputScalar3((std::string("obj") + utils::padWithZeros(std::to_string(object->idx), 4) + " force").c_str(), object->forceE.data());
                            utils::InputScalar3((std::string("obj") + utils::padWithZeros(std::to_string(object->idx), 4) + " torque").c_str(), object->torqueE.data());
                            ImGui::TreePop();
                        }
                    }
                    ImGui::TreePop();
                }

                if(ImGui::TreeNode("static link")){
                    for(const auto& multibody : multibodiesList){
                        for(const auto& link : multibody->linksList){
                            if(link->fixed){
                                utils::InputScalar3((std::string("link") + utils::padWithZeros(std::to_string(link->idx), 4)).c_str(), link->pd.data());
                            }
                        }
                    }
                    ImGui::TreePop();
                }
            }
        }
    };
}

void Simulator::registerMesh(){
    for(int objIdx = 0; objIdx < objsList.size(); ++objIdx){
        objsList[objIdx]->setMeshIdx(meshesList.size());
        meshesList.push_back(polyscope::registerSurfaceMesh(std::string("obj") + utils::padWithZeros(std::to_string(objIdx), 4), objsList[objIdx]->meshTV, objsList[objIdx]->meshTF));
    }
    for(int multibodyIdx = 0; multibodyIdx < multibodiesList.size(); ++multibodyIdx){
        glm::vec3 armColor(0.9, 0.9, 0.9);
        for(int objIdx = 0; objIdx < multibodiesList[multibodyIdx]->linksList.size(); ++objIdx){
            multibodiesList[multibodyIdx]->linksList[objIdx]->setMeshIdx(meshesList.size());
            meshesList.push_back(polyscope::registerSurfaceMesh(std::string("arm") + utils::padWithZeros(std::to_string(multibodyIdx), 4) + "link" + utils::padWithZeros(std::to_string(objIdx), 4), multibodiesList[multibodyIdx]->linksList[objIdx]->meshTV, multibodiesList[multibodyIdx]->linksList[objIdx]->meshTF));
            meshesList[meshesList.size() - 1]->setSurfaceColor(armColor);
        }
    }
}

void Simulator::updateMeshTransform(){
    glm::mat4 T;
    for(int objIdx = 0; objIdx < objsList.size(); ++objIdx){
        utils::q2glm(objsList[objIdx]->q, T);
        meshesList[objsList[objIdx]->meshIdx]->setTransform(T);
    }
    for(int multibodyIdx = 0; multibodyIdx < multibodiesList.size(); ++multibodyIdx){
        for(int linkIdx = 0; linkIdx < multibodiesList[multibodyIdx]->linksList.size(); linkIdx++){
            utils::q2glm(multibodiesList[multibodyIdx]->linksList[linkIdx]->q, T);
            meshesList[multibodiesList[multibodyIdx]->linksList[linkIdx]->meshIdx]->setTransform(T);
        }
    }
}

}