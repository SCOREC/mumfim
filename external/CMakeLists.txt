set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL YES)

include(FetchContent)
SET(CMAKE_CXX_FLAGS "")

FetchContent_Declare(catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG 037ddbc75cc5e58b93cf5a010a94b32333ad824d)

FetchContent_MakeAvailable(catch2)

FetchContent_Declare(ModelTraits
        GIT_REPOSITORY https://github.com/jacobmerson/model-traits
        GIT_TAG 8c5e9296435f25e552ee6d20108072239024585c
        )
FetchContent_MakeAvailable(ModelTraits)

FetchContent_Declare(LAS
        GIT_REPOSITORY https://github.com/wrtobin/las
        GIT_TAG 09c0c79dd40407271be3388a761b7f2b6ee2521f
        )
FetchContent_MakeAvailable(LAS)

set(YAML_CPP_INSTALL ON CACHE INTERNAL "")
FetchContent_Declare(yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 9a3624205e8774953ef18f57067b3426c1c5ada6)
FetchContent_MakeAvailable(yaml-cpp)

FetchContent_Declare(AMSI
        GIT_REPOSITORY https://github.com/jacobmerson/amsi
        GIT_TAG 75f238f944426e9bb05344fca41091bbd77a55d1
        )
FetchContent_MakeAvailable(AMSI)

