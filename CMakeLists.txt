# Required CMake arguments
cmake_minimum_required(VERSION 3.20)
#set standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY") #needed for testing cross compilation see: https://stackoverflow.com/questions/53633705/cmake-the-c-compiler-is-not-able-to-compile-a-simple-test-program

#set(parallel ON) # currently set through -Dparallel, but can be manually set by uncommenting
set(cxx_flags "-O2") #release mode overwrites with "-O3"
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${cxx_flags}") # -O2 appears to be the fastest option
set(CMAKE_INCLUDE_CURRENT_DIR ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

#set compiler if necessary, this has to be done before project()
#set (CMAKE_CXX_COMPILER <specify compiler here>)
#if(CMAKE_CXX_COMPILER_LOADED)
#    message(STATUS "Configured for C++, with compiler: ")
#    message(STATUS "${CMAKE_CXX_COMPILER_VERSION}}")
#endif()

#################### BELOW SPECIFIES COMPILATION FOR SERIAL OR PARALLEL VERSION OF TRIBS ###############################
if(parallel)
    set(exe "tRIBSpar")
    project(${exe}
            DESCRIPTION "Parallel Distribution of the TIN-based Real-time Integrated Basin Simulator Distributed Hydrological Model"
            VERSION 5.2.0
    )
    add_executable(${exe} src/main.cpp)

    message(WARNING "Compiling parallel version of tRIBS")

    target_compile_definitions(${exe}
            PRIVATE
            PARALLEL_TRIBS
            PRIVATE
            GRAPH_TRIBS
            PRIVATE
            LINUX_32
    )

    find_package(MPI REQUIRED) #Note: changes to MPI path etc may require reset of CMakeCache.txt
    if (MPI_FOUND) #Note:this should find whichever implementation of MPI is installed (i.e. Open MPI, Intel, mpich)
        include_directories(${MPI_INCLUDE_PATH})
        target_link_libraries(${exe} PUBLIC MPI::MPI_CXX)
        set(CMAKE_CXX_COMPILER mpicxx)
    else (MPI_FOUND)
        message(SEND_ERROR "This application cannot compile without MPI")
    endif (MPI_FOUND)

    target_sources(
            ${exe} PRIVATE src/main.cpp
            src/Headers/Classes.h
            src/Headers/Definitions.h
            src/Headers/Inclusions.h
            src/Headers/TemplDefinitions.h
            src/Headers/globalFns.cpp
            src/Headers/globalFns.h
            src/Headers/globalIO.h
            src/Mathutil/geometry.h
            src/Mathutil/mathutil.cpp
            src/Mathutil/mathutil.h
            src/Mathutil/predicates.cpp
            src/Mathutil/predicates.h
            src/tArray/tArray.h
            src/tArray/tMatrix.cpp
            src/tArray/tMatrix.h
            src/tCNode/tCNode.cpp
            src/tCNode/tCNode.h
            src/tFlowNet/tFlowNet.cpp
            src/tFlowNet/tFlowNet.h
            src/tFlowNet/tFlowResults.cpp
            src/tFlowNet/tFlowResults.h
            src/tFlowNet/tKinemat.cpp
            src/tFlowNet/tKinemat.h
            src/tFlowNet/tResData.cpp
            src/tFlowNet/tResData.h
            src/tFlowNet/tReservoir.cpp
            src/tFlowNet/tReservoir.h
            src/tGraph/tGraph.cpp
            src/tGraph/tGraph.h
            src/tGraph/tGraphNode.cpp
            src/tGraph/tGraphNode.h
            src/tHydro/tEvapoTrans.cpp
            src/tHydro/tEvapoTrans.h
            src/tHydro/tHydroMet.cpp
            src/tHydro/tHydroMet.h
            src/tHydro/tHydroMetConvert.cpp
            src/tHydro/tHydroMetConvert.h
            src/tHydro/tHydroMetStoch.cpp
            src/tHydro/tHydroMetStoch.h
            src/tHydro/tHydroModel.cpp
            src/tHydro/tHydroModel.h
            src/tHydro/tIntercept.cpp
            src/tHydro/tIntercept.h
            src/tHydro/tSnowPack.cpp
            src/tHydro/tSnowPack.h
            src/tHydro/tWaterBalance.cpp
            src/tHydro/tWaterBalance.h
            src/tInOut/tInputFile.cpp
            src/tInOut/tInputFile.h
            src/tInOut/tOstream.cpp
            src/tInOut/tOstream.h
            src/tInOut/tOutput.cpp
            src/tInOut/tOutput.h
            src/tList/tList.cpp
            src/tList/tList.h
            src/tListInputData/tListInputData.cpp
            src/tListInputData/tListInputData.h
            src/tMesh/heapsort.h
            src/tMesh/tMesh.cpp
            src/tMesh/tMesh.h
            src/tMesh/tTriangulator.cpp
            src/tMesh/tTriangulator.h
            src/tMeshElements/meshElements.cpp
            src/tMeshElements/meshElements.h
            src/tMeshList/tMeshList.h
            src/tParallel/tParallel.cpp
            src/tParallel/tParallel.h
            src/tParallel/tTimer.cpp
            src/tParallel/tTimer.h
            src/tParallel/tTimings.cpp
            src/tParallel/tTimings.h
            src/tPtrList/tPtrList.cpp
            src/tPtrList/tPtrList.h
            src/tRasTin/tInvariant.cpp
            src/tRasTin/tInvariant.h
            src/tRasTin/tRainGauge.cpp
            src/tRasTin/tRainGauge.h
            src/tRasTin/tRainfall.cpp
            src/tRasTin/tRainfall.h
            src/tRasTin/tResample.cpp
            src/tRasTin/tResample.h
            src/tRasTin/tShelter.cpp
            src/tRasTin/tShelter.h
            src/tRasTin/tVariant.cpp
            src/tRasTin/tVariant.h
            src/tSimulator/tControl.cpp
            src/tSimulator/tControl.h
            src/tSimulator/tPreProcess.cpp
            src/tSimulator/tPreProcess.h
            src/tSimulator/tRestart.cpp
            src/tSimulator/tRestart.h
            src/tSimulator/tRunTimer.cpp
            src/tSimulator/tRunTimer.h
            src/tSimulator/tSimul.cpp
            src/tSimulator/tSimul.h
            src/tStorm/tStorm.cpp
            src/tStorm/tStorm.h
            ${exe} PUBLIC tests/tHydroTests/ThydroModelTests.cpp
    )

else()
    set(exe "tRIBS")
    project(${exe}
            DESCRIPTION "Serial Distribution of the TIN-based Real-time Integrated Basin Simulator Distributed Hydrological Model"
            VERSION 5.2.0)
    add_executable(${exe} src/main.cpp)
    message(WARNING "Compiling serial version of tRIBS")

    if(APPLE) # not sure though if loading libraries are any different for linux--windows don't know
        target_compile_definitions(${exe}
                PRIVATE
                MAC
        )
    else()
        target_compile_definitions(${exe}
                PRIVATE
                LINUX_32
        )
    endif()

    target_sources(
            ${exe} PRIVATE src/main.cpp
            src/Headers/Classes.h
            src/Headers/Definitions.h
            src/Headers/Inclusions.h
            src/Headers/TemplDefinitions.h
            src/Headers/globalFns.cpp
            src/Headers/globalFns.h
            src/Headers/globalIO.h
            src/Mathutil/geometry.h
            src/Mathutil/mathutil.cpp
            src/Mathutil/mathutil.h
            src/Mathutil/predicates.cpp
            src/Mathutil/predicates.h
            src/tArray/tArray.h
            src/tArray/tMatrix.cpp
            src/tArray/tMatrix.h
            src/tCNode/tCNode.cpp
            src/tCNode/tCNode.h
            src/tFlowNet/tFlowNet.cpp
            src/tFlowNet/tFlowNet.h
            src/tFlowNet/tFlowResults.cpp
            src/tFlowNet/tFlowResults.h
            src/tFlowNet/tKinemat.cpp
            src/tFlowNet/tKinemat.h
            src/tFlowNet/tResData.cpp
            src/tFlowNet/tResData.h
            src/tFlowNet/tReservoir.cpp
            src/tFlowNet/tReservoir.h
            src/tGraph/tGraph.cpp
            src/tGraph/tGraph.h
            src/tGraph/tGraphNode.cpp
            src/tGraph/tGraphNode.h
            src/tHydro/tEvapoTrans.cpp
            src/tHydro/tEvapoTrans.h
            src/tHydro/tHydroMet.cpp
            src/tHydro/tHydroMet.h
            src/tHydro/tHydroMetConvert.cpp
            src/tHydro/tHydroMetConvert.h
            src/tHydro/tHydroMetStoch.cpp
            src/tHydro/tHydroMetStoch.h
            src/tHydro/tHydroModel.cpp
            src/tHydro/tHydroModel.h
            src/tHydro/tIntercept.cpp
            src/tHydro/tIntercept.h
            src/tHydro/tSnowPack.cpp
            src/tHydro/tSnowPack.h
            src/tHydro/tWaterBalance.cpp
            src/tHydro/tWaterBalance.h
            src/tInOut/tInputFile.cpp
            src/tInOut/tInputFile.h
            src/tInOut/tOstream.cpp
            src/tInOut/tOstream.h
            src/tInOut/tOutput.cpp
            src/tInOut/tOutput.h
            src/tList/tList.cpp
            src/tList/tList.h
            src/tListInputData/tListInputData.cpp
            src/tListInputData/tListInputData.h
            src/tMesh/heapsort.h
            src/tMesh/tMesh.cpp
            src/tMesh/tMesh.h
            src/tMesh/tTriangulator.cpp
            src/tMesh/tTriangulator.h
            src/tMeshElements/meshElements.cpp
            src/tMeshElements/meshElements.h
            src/tMeshList/tMeshList.h
            src/tPtrList/tPtrList.cpp
            src/tPtrList/tPtrList.h
            src/tRasTin/tInvariant.cpp
            src/tRasTin/tInvariant.h
            src/tRasTin/tRainGauge.cpp
            src/tRasTin/tRainGauge.h
            src/tRasTin/tRainfall.cpp
            src/tRasTin/tRainfall.h
            src/tRasTin/tResample.cpp
            src/tRasTin/tResample.h
            src/tRasTin/tShelter.cpp
            src/tRasTin/tShelter.h
            src/tRasTin/tVariant.cpp
            src/tRasTin/tVariant.h
            src/tSimulator/tControl.cpp
            src/tSimulator/tControl.h
            src/tSimulator/tPreProcess.cpp
            src/tSimulator/tPreProcess.h
            src/tSimulator/tRestart.cpp
            src/tSimulator/tRestart.h
            src/tSimulator/tRunTimer.cpp
            src/tSimulator/tRunTimer.h
            src/tSimulator/tSimul.cpp
            src/tSimulator/tSimul.h
            src/tStorm/tStorm.cpp
            src/tStorm/tStorm.h
            tests/tHydroTests/ThydroModelTests.cpp
    )

endif()

# can be passed via the command line with -Dpackaging=ON
#set(packaging ON)

if (packaging)

    if (NOT DEFINED os)
        set(os "macOS/Silicon") # options: macOS/Silicon, macOS/Intel, or Linux
    endif ()

    install(
        TARGETS ${PROJECT_NAME}
    )

    set(CPACK_VERBOSE ON)
    set(CPACK_PACKAGE_NAME "TIN-based Real-time Integrated Basin Simulator")
    set(CPACK_PACKAGE_FILE_NAME "TIN-based Real-time Integrated Basin Simulator-5.2.0")
    set(CPACK_PACKAGE_VENDOR "tRIBS developers")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE.txt")
    # Set the custom script for the STGZ generator
    #list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/packaging/utilities/CPack.STGZ_Header.sh.in")
    #set(CPACK_STGZ_HEADER_FILE "${CMAKE_SOURCE_DIR}/packaging/utilities/CPack.STGZ_Header.sh.in")

    set(CPACK_GENERATOR "STGZ")
    set(CPACK_PACKAGE_DIRECTORY "${CMAKE_SOURCE_DIR}/packaging/${os}")
    set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_SOURCE_DIR}/packaging/${os}")

    include(CPack)
    cpack_add_component(CPACK_COMPONENTS_ALL)

endif ()

enable_testing ()

include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(
        tHydroTests
        ${PROJECT_SOURCE_DIR}/tests/tHydrotests/THydroModelTests.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroModel.h
        ${PROJECT_SOURCE_DIR}/src/Headers/Classes.h
        ${PROJECT_SOURCE_DIR}/src/Headers/Inclusions.h
        ${PROJECT_SOURCE_DIR}/src/Headers/Definitions.h
        ${PROJECT_SOURCE_DIR}/src/Headers/TemplDefinitions.h
        ${PROJECT_SOURCE_DIR}/src/Headers/Classes.h
        ${PROJECT_SOURCE_DIR}/src/Headers/Definitions.h
        ${PROJECT_SOURCE_DIR}/src/Headers/Inclusions.h
        ${PROJECT_SOURCE_DIR}/src/Headers/TemplDefinitions.h
        ${PROJECT_SOURCE_DIR}/src/Headers/globalFns.cpp
        ${PROJECT_SOURCE_DIR}/src/Headers/globalFns.h
        ${PROJECT_SOURCE_DIR}/src/Headers/globalIO.h
        ${PROJECT_SOURCE_DIR}/src/Mathutil/geometry.h
        ${PROJECT_SOURCE_DIR}/src/Mathutil/mathutil.cpp
        ${PROJECT_SOURCE_DIR}/src/Mathutil/mathutil.h
        ${PROJECT_SOURCE_DIR}/src/Mathutil/predicates.cpp
        ${PROJECT_SOURCE_DIR}/src/Mathutil/predicates.h
        ${PROJECT_SOURCE_DIR}/src/tArray/tArray.h
        ${PROJECT_SOURCE_DIR}/src/tArray/tMatrix.cpp
        ${PROJECT_SOURCE_DIR}/src/tArray/tMatrix.h
        ${PROJECT_SOURCE_DIR}/src/tCNode/tCNode.cpp
        ${PROJECT_SOURCE_DIR}/src/tCNode/tCNode.h
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tFlowNet.cpp
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tFlowNet.h
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tFlowResults.cpp
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tFlowResults.h
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tKinemat.cpp
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tKinemat.h
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tResData.cpp
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tResData.h
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tReservoir.cpp
        ${PROJECT_SOURCE_DIR}/src/tFlowNet/tReservoir.h
        ${PROJECT_SOURCE_DIR}/src/tGraph/tGraph.cpp
        ${PROJECT_SOURCE_DIR}/src/tGraph/tGraph.h
        ${PROJECT_SOURCE_DIR}/src/tGraph/tGraphNode.cpp
        ${PROJECT_SOURCE_DIR}/src/tGraph/tGraphNode.h
        ${PROJECT_SOURCE_DIR}/src/tHydro/tEvapoTrans.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tEvapoTrans.h
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroMet.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroMet.h
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroMetConvert.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroMetConvert.h
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroMetStoch.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroMetStoch.h
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroModel.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tHydroModel.h
        ${PROJECT_SOURCE_DIR}/src/tHydro/tIntercept.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tIntercept.h
        ${PROJECT_SOURCE_DIR}/src/tHydro/tSnowPack.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tSnowPack.h
        ${PROJECT_SOURCE_DIR}/src/tHydro/tWaterBalance.cpp
        ${PROJECT_SOURCE_DIR}/src/tHydro/tWaterBalance.h
        ${PROJECT_SOURCE_DIR}/src/tInOut/tInputFile.cpp
        ${PROJECT_SOURCE_DIR}/src/tInOut/tInputFile.h
        ${PROJECT_SOURCE_DIR}/src/tInOut/tOstream.cpp
        ${PROJECT_SOURCE_DIR}/src/tInOut/tOstream.h
        ${PROJECT_SOURCE_DIR}/src/tInOut/tOutput.cpp
        ${PROJECT_SOURCE_DIR}/src/tInOut/tOutput.h
        ${PROJECT_SOURCE_DIR}/src/tList/tList.cpp
        ${PROJECT_SOURCE_DIR}/src/tList/tList.h
        ${PROJECT_SOURCE_DIR}/src/tListInputData/tListInputData.cpp
        ${PROJECT_SOURCE_DIR}/src/tListInputData/tListInputData.h
        ${PROJECT_SOURCE_DIR}/src/tMesh/heapsort.h
        ${PROJECT_SOURCE_DIR}/src/tMesh/tMesh.cpp
        ${PROJECT_SOURCE_DIR}/src/tMesh/tMesh.h
        ${PROJECT_SOURCE_DIR}/src/tMesh/tTriangulator.cpp
        ${PROJECT_SOURCE_DIR}/src/tMesh/tTriangulator.h
        ${PROJECT_SOURCE_DIR}/src/tMeshElements/meshElements.cpp
        ${PROJECT_SOURCE_DIR}/src/tMeshElements/meshElements.h
        ${PROJECT_SOURCE_DIR}/src/tMeshList/tMeshList.h
        ${PROJECT_SOURCE_DIR}/src/tPtrList/tPtrList.cpp
        ${PROJECT_SOURCE_DIR}/src/tPtrList/tPtrList.h
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tInvariant.cpp
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tInvariant.h
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tRainGauge.cpp
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tRainGauge.h
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tRainfall.cpp
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tRainfall.h
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tResample.cpp
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tResample.h
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tShelter.cpp
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tShelter.h
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tVariant.cpp
        ${PROJECT_SOURCE_DIR}/src/tRasTin/tVariant.h
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tControl.cpp
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tControl.h
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tPreProcess.cpp
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tPreProcess.h
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tRestart.cpp
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tRestart.h
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tRunTimer.cpp
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tRunTimer.h
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tSimul.cpp
        ${PROJECT_SOURCE_DIR}/src/tSimulator/tSimul.h
        ${PROJECT_SOURCE_DIR}/src/tStorm/tStorm.cpp
        ${PROJECT_SOURCE_DIR}/src/tStorm/tStorm.h
)
target_link_libraries(
  tHydroTests
  GTest::gtest_main
)

include(GoogleTest)
gtest_add_tests(TARGET tHydroTests)

add_subdirectory(tests)
