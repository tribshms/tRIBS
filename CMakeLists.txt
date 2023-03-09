cmake_minimum_required(VERSION 3.24)
project(tRIBS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_executable(tRIBS main.cpp)
target_sources(tRIBS PRIVATE "main.cpp")
target_include_directories(tRIBS PRIVATE "${PROJECT_SOURCE_DIR}")

add_subdirectory("tArray")
add_subdirectory("tGraph")
add_subdirectory("tHydro")
add_subdirectory("tInOut")
add_subdirectory("tMesh")
add_subdirectory("tMeshElements")
add_subdirectory("tMeshList")
add_subdirectory("tRasTin")
add_subdirectory("tSimulator")
add_subdirectory("tStorm")
add_subdirectory("tCNode")
add_subdirectory("Mathutil")
add_subdirectory("tFlowNet")
add_subdirectory("tList")
add_subdirectory("tListInputData")
add_subdirectory("tParallel")
add_subdirectory("Utilities")
#

#
#target_link_libraries(tRIBS PRIVATE
#        tArray
#        tGraph
#        tHydro
#        tInOut
#        tMesh
#        tMeshElements
#        tMeshList
#        tRasTin
#        tSimulator
#        tStorm
#        tCNode
#        Mathutil
#        tFlowNet
#        tList
#        tListInputData
#        tParallel
#        Utilities)
#
#target_include_directories(tRIBS PUBLIC
#        "${PROJECT_BINARY_DIR}"
#        "${PROJECT_SOURCE_DIR}/tArray"
#        "${PROJECT_SOURCE_DIR}/tGraph"
#        "${PROJECT_SOURCE_DIR}/tHydro"
#        "${PROJECT_SOURCE_DIR}/tInOut"
#        "${PROJECT_SOURCE_DIR}/tMesh"
#        "${PROJECT_SOURCE_DIR}/tMeshElements"
#        "${PROJECT_SOURCE_DIR}/tMeshList"
#        "${PROJECT_SOURCE_DIR}/tRasTin"
#        "${PROJECT_SOURCE_DIR}/tSimulator"
#        "${PROJECT_SOURCE_DIR}/tStorm"
#        "${PROJECT_SOURCE_DIR}/tCNode"
#        "${PROJECT_SOURCE_DIR}/Mathutil"
#        "${PROJECT_SOURCE_DIR}/tFlowNet"
#        "${PROJECT_SOURCE_DIR}/tList"
#        "${PROJECT_SOURCE_DIR}/tListInputData"
#        "${PROJECT_SOURCE_DIR}/tParallel"
#        "${PROJECT_SOURCE_DIR}/Utilities")
#

