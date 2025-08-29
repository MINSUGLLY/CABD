#pragma once

#include "multibody.hpp"

namespace CABD{

class Simulator{
public:
    Simulator() = default;
    void loadConfig(const std::string& configPath);
    void launch();
    void simLoop();

    std::vector<std::shared_ptr<Object>> objsList;
    std::vector<std::shared_ptr<Multibody>> multibodiesList;

    std::map<std::string, std::shared_ptr<Object>> geoCache;

private:
    simType type = simType::forward;

    vec3 gravity;
    scalar dt;
    scalar convergenceEpsilon;
    scalar alphaGD;
    
    int framesTotal;
    int framesCount;

    std::chrono::_V2::system_clock::time_point start_time;
    std::chrono::_V2::system_clock::time_point end_time;

    std::mutex mtx;
    std::condition_variable cv;
    bool renderFinish = false;
    bool simulationFinish = true;

    bool simRunning = true;

    void initViewer();
    void registerMesh();
    void updateMeshTransform();
    std::vector<polyscope::SurfaceMesh*> meshesList;

    void appendObject(std::shared_ptr<Object> object, const vec3 initp = vec3::Zero(), const mat3 initR = mat3::Identity());
    void appendMultibody(std::shared_ptr<Multibody> multibody);

    void advance();
    void objectAdvance(std::shared_ptr<Object> object);
    void multibodyAdvance(std::shared_ptr<Multibody> multibody);
};
    
}