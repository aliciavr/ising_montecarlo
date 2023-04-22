cmake_minimum_required(VERSION 3.24)
project(ising_montecarlo)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED COMPONENTS RIO Net MathCore Hist Physics Core ROOTDataFrame)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

add_executable(ising_montecarlo main.cpp)
target_link_libraries(ising_montecarlo PUBLIC ROOT::RIO ROOT::Net ROOT::MathCore ROOT::Hist ROOT::Physics ROOT::Core ROOT::ROOTDataFrame)

