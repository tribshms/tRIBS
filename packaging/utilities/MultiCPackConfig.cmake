include("cmake-build-parallel/CPackConfig.cmake")
include("cmake-build-serial/CPackConfig.cmake")

set(CPACK_INSTALL_CMAKE_PROJECTS
        "cmake-build-serial;tRIBS;ALL;/"
        "cmake-build-parallel;tRIBSpar;ALL;/"
)

