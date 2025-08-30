#include "simulator.hpp"

namespace CABD{

void Simulator::loadConfig(const std::string& configPath){
    spdlog::info("load scene config ...");
    std::ifstream file(configPath);
    nlohmann::json j;
    file >> j;

    gravity = utils::jsonV3(j["gravity"]);
    dt = j["dt"];
    convergenceEpsilon = j["convergenceEpsilon"];
    alphaGD = j["alphaGD"];

    framesTotal = j["framesTotal"];
    framesCount = framesTotal;

    // object
    int objCount = 0;
    spdlog::info("load object ...");
    for(const auto& objConfig : j["object"]){
        spdlog::info("load object config {} / {}", ++objCount, j["object"].size());
        mat3 initR = utils::jsonM3(objConfig["initR"]);
        vec3 initp = utils::jsonV3(objConfig["initp"]);
        if(geoCache.find(objConfig["configPath"]) != geoCache.end()){
            spdlog::info("use cached geometry {}", objConfig["configPath"]);
            appendObject(std::make_shared<Object>(geoCache[objConfig["configPath"]]), initp, initR);
        }
        else{
            appendObject(std::make_shared<Object>(objConfig["configPath"]), initp, initR);
            geoCache[objConfig["configPath"]] = objsList.back();
        }
    }
    totalObjNum = objsList.size();

    // multibody
    int multibodyCount = 0;
    spdlog::info("load multibody ...");
    for(const auto& multibodyConfig : j["multibody"]){
        spdlog::info("load multibody config {} / {}", ++objCount, j["multibody"].size());
        appendMultibody(std::make_shared<Multibody>(multibodyConfig["configPath"], geoCache));
        totalObjNum += multibodiesList.back()->linksList.size();
    }
}

void Simulator::launch(){
    polyscope::init();

    initViewer();
    registerMesh();
    updateMeshTransform();

    polyscope::show();
    polyscope::shutdown();

    simRunning = false;
}

void Simulator::appendObject(std::shared_ptr<Object> object, const vec3 initp, const mat3 initR){
    object->setIdx(objsList.size());
    utils::Rp2q(object->q, initR, initp + initR * object->massCenter);
    objsList.push_back(object);
}

void Simulator::appendMultibody(std::shared_ptr<Multibody> multibody){
    multibody->setIdx(multibodiesList.size());
    multibodiesList.push_back(multibody);
}

void Simulator::simLoop(){
    while(simRunning){
        if(framesCount < framesTotal){
            if(framesCount == 0){
                spdlog::trace("move start");
                start_time = std::chrono::high_resolution_clock::now();
            }
            advance();
            framesCount++;
            spdlog::trace("frame {} / {}", framesCount, framesTotal);
            if(framesCount == framesTotal){
                spdlog::trace("move end");
                end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end_time - start_time;
                spdlog::info("total: {}s, average: {}s, FPS: {}", duration.count(), duration.count() / framesTotal, framesTotal / duration.count());
            }
        }    
        else{
            std::this_thread::sleep_for(std::chrono::milliseconds(1)); 
        }
    }
}

void Simulator::saveState(const std::string& path){
    matX result(totalObjNum, 12);
    int objCount = 0;
    for(const auto& obj : objsList){
        result.row(objCount) = obj->q;
        objCount++;
    }
    for(const auto& multibody : multibodiesList){
        for(const auto& link : multibody->linksList){
            result.row(objCount) = link->q;
            objCount++;
        }
    }

    utils::saveNpy(result, path);
}

}