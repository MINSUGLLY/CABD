#include "simulator.hpp"

using namespace CABD;

int main(int argc, char** argv){
    spdlog::set_level(spdlog::level::trace);

    CLI::App app{"CABD Simulator"};
    
    std::string config;
    std::string savePath;
    app.add_option("-c,--config", config);
    app.add_option("-s,--save", savePath);

    CLI11_PARSE(app, argc, argv);

    if(config.empty()){
        spdlog::error("please provide config path");
        return 0;
    }

    Simulator sim;
    sim.loadConfig("../config/scene/" + config + ".json");

    if(!savePath.empty()){
        sim.savePath = std::string("../results/") + savePath + std::string("/");
        sim.saveInterval = static_cast<int>(1.0 / 60.0 / sim.dt);
    }

    std::thread simThread(std::bind(&Simulator::simLoop, &sim));
    sim.launch();
    simThread.join();
}

