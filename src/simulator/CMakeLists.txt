file(GLOB SIMULATORSOURCES ./*.cpp)
add_library(Sim STATIC ${SIMULATORSOURCES})
target_link_libraries(Sim PUBLIC Structure)