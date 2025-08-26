file(GLOB STRUCTURESOURCES ./*.cpp)
add_library(Structure STATIC ${STRUCTURESOURCES})
target_link_libraries(Structure PUBLIC Utils)