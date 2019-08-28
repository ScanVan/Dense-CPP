#Usage => cmake -DCMAKE_BUILD_TYPE=Debug .

cmake_minimum_required(VERSION 2.8.0)
project( Dense-CPP LANGUAGES CUDA CXX)
set(CMAKE_CXX_FLAGS "-std=c++11") 

find_package( OpenCV REQUIRED )
find_package( Eigen3 REQUIRED )

set(ROOT ./)
include_directories(${ROOT}/src)
include_directories(${ROOT}/third)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
include_directories(${ROOT}/third/nvof/Includes)
include_directories(${ROOT}/third/nvof/NvOF)
include_directories(${ROOT}/third/nvof/Utils)




set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


SET(GCC_COVERAGE_COMPILE_FLAGS "")
SET(GCC_COVERAGE_LINK_FLAGS "")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )


file(GLOB Processor_SRC
    "src/*.hpp"
    "src/*.cpp"
    "third/*.cpp"
#    "third/nvof/NvOF/*.cpp"
#    "third/nvof/Utils/*.cpp"
)
add_executable( Dense-CPP ${Processor_SRC} )
target_link_libraries( Dense-CPP ${OpenCV_LIBS} dl cuda)


add_custom_target(run
    COMMAND bin/Dense-CPP
    DEPENDS Dense-CPP
)
