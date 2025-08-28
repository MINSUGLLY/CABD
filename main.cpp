#include "simulator.hpp"

using namespace labd;

int main(){
    // spdlog::set_level(spdlog::level::info);

    Simulator sim;

    /*---- config ----*/

    // sim.loadConfig("../config/scene/ball.json");
    // sim.loadConfig("../config/scene/cube.json");
    // sim.loadConfig("../config/scene/rod.json");
    // sim.loadConfig("../config/scene/link.json");
    // sim.loadConfig("../config/scene/link_circle.json");

    // sim.loadConfig("../config/scene/chain.json");
    // sim.loadConfig("../config/scene/ring.json");
    sim.loadConfig("../config/scene/grid.json");

    // sim.loadConfig("../config/scene/deform_test.json");
    // sim.loadConfig("../config/scene/rod_static.json");

    /*---- end ----*/

    std::thread simThread(std::bind(&Simulator::simLoop, &sim));
    sim.launch();
    simThread.join();
}

