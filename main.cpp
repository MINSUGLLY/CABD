#include "simulator.hpp"

using namespace CABD;

int main(int argc, char** argv){
    spdlog::set_level(spdlog::level::trace);

    CLI::App app{"CABD Simulator"};
    
    std::string config;
    app.add_option("-c,--config", config);
    
    CLI11_PARSE(app, argc, argv);

    if(config.empty()){
        spdlog::error("please provide config path");
        return 0;
    }

    Simulator sim;
    sim.loadConfig("../config/scene/" + config + ".json");

    std::thread simThread(std::bind(&Simulator::simLoop, &sim));
    sim.launch();
    simThread.join();
}

