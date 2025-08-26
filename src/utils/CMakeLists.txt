find_package(MKL CONFIG REQUIRED)

file(GLOB UTILSSOURCES ./*.cpp)

add_library(Utils STATIC ${UTILSSOURCES})
target_link_libraries(Utils PUBLIC 
                      Eigen3::Eigen 
                      igl::core
                      polyscope 
                      cnpy
                      ipc_toolkit)

target_compile_options(Utils PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(Utils PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(Utils PUBLIC $<LINK_ONLY:MKL::MKL>)