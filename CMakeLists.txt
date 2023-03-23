cmake_minimum_required(VERSION 3.24)
project(tRIBSgrid)
add_executable(tRIBSgrid main.cpp)
set(CMAKE_CXX_STANDARD 17)

# setup for compiling with parallelization feature requires MPI
#set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/open-mpi/4.1.5/lib")
find_package(MPI REQUIRED) #Note: changes to MPI path etc may require reset of CMakeCache.txt
if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(tRIBSgrid PUBLIC MPI::MPI_CXX)
    set(CMAKE_CXX_COMPILER mpicxx)
else (MPI_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)

# used while testing and developing
set(CMAKE_CXX_EXTENSIONS OFF)
set(DCMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -g") # different options available for using Sanitizer

# currently
if(APPLE) # not sure though if loading libraries are any different for linux--windows don't know
    target_compile_definitions(tRIBSgrid
            PRIVATE
            MAC
            )
endif()

# project files
set(CMAKE_INCLUDE_CURRENT_DIR ON)
target_sources(
        tRIBSgrid PRIVATE main.cpp
        Headers/Classes.h
        Headers/Definitions.h
        Headers/Inclusions.h
        Headers/TemplDefinitions.h
        Headers/globalFns.cpp
        Headers/globalFns.h
        Headers/globalIO.h
        Mathutil/geometry.h
        Mathutil/mathutil.cpp
        Mathutil/mathutil.h
        Mathutil/predicates.cpp
        Mathutil/predicates.h
        #Utilities/InitialGW.cpp
        #Utilities/RainInputCheck.cpp
        #Utilities/RunsTracker.cpp
        tArray/tArray.cpp
        tArray/tArray.h
        tArray/tMatrix.cpp
        tArray/tMatrix.h
        tCNode/tCNode.cpp
        tCNode/tCNode.h
        tFlowNet/tFlowNet.cpp
        tFlowNet/tFlowNet.h
        tFlowNet/tFlowResults.cpp
        tFlowNet/tFlowResults.h
        tFlowNet/tKinemat.cpp
        tFlowNet/tKinemat.h
        tFlowNet/tResData.cpp
        tFlowNet/tResData.h
        tFlowNet/tReservoir.cpp
        tFlowNet/tReservoir.h
        tGraph/tGraph.cpp
        tGraph/tGraph.h
        tGraph/tGraphNode.cpp
        tGraph/tGraphNode.h
        tHydro/tEvapoTrans.cpp
        tHydro/tEvapoTrans.h
        tHydro/tHydroMet.cpp
        tHydro/tHydroMet.h
        tHydro/tHydroMetConvert.cpp
        tHydro/tHydroMetConvert.h
        tHydro/tHydroMetStoch.cpp
        tHydro/tHydroMetStoch.h
        tHydro/tHydroModel.cpp
        tHydro/tHydroModel.h
        tHydro/tIntercept.cpp
        tHydro/tIntercept.h
        tHydro/tSnowIntercept.cpp
        tHydro/tSnowIntercept.h
        tHydro/tSnowPack.cpp
        tHydro/tSnowPack.h
        tHydro/tWaterBalance.cpp
        tHydro/tWaterBalance.h
        tInOut/tInputFile.cpp
        tInOut/tInputFile.h
        tInOut/tOstream.cpp
        tInOut/tOstream.h
        tInOut/tOutput.cpp
        tInOut/tOutput.h
        tList/tList.cpp
        tList/tList.h
        tListInputData/tListInputData.cpp
        tListInputData/tListInputData.h
        tMesh/heapsort.h
        tMesh/tMesh.cpp
        tMesh/tMesh.h
        tMesh/tTriangulator.cpp
        tMesh/tTriangulator.h
        tMeshElements/meshElements.cpp
        tMeshElements/meshElements.h
        tMeshList/tMeshList.cpp
        tMeshList/tMeshList.h
         tParallel/tParallel.cpp
         tParallel/tParallel.h
         tParallel/tTimer.cpp
         tParallel/tTimer.h
         tParallel/tTimings.cpp
         tParallel/tTimings.h
        tPtrList/tPtrList.cpp
        tPtrList/tPtrList.h
        tRasTin/tInvariant.cpp
        tRasTin/tInvariant.h
        tRasTin/tRainGauge.cpp
        tRasTin/tRainGauge.h
        tRasTin/tRainfall.cpp
        tRasTin/tRainfall.h
        tRasTin/tResample.cpp
        tRasTin/tResample.h
        tRasTin/tShelter.cpp
        tRasTin/tShelter.h
        tRasTin/tVariant.cpp
        tRasTin/tVariant.h
        tSimulator/tControl.cpp
        tSimulator/tControl.h
        tSimulator/tPreProcess.cpp
        tSimulator/tPreProcess.h
        tSimulator/tRestart.cpp
        tSimulator/tRestart.h
        tSimulator/tRunTimer.cpp
        tSimulator/tRunTimer.h
        tSimulator/tSimul.cpp
        tSimulator/tSimul.h
        tStorm/tStorm.cpp
        tStorm/tStorm.h
)