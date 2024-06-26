cmake_minimum_required (VERSION 3.6) # Not sure what the minimum is required
project (a5_base CXX C)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

# Use modern C++
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

include_directories("${PROJECT_BINARY_DIR}" libs src)

set(DEFINITIONS _USE_MATH_DEFINES=1 GLM_FORCE_CXX14=1
    IMGUI_IMPL_OPENGL_LOADER_CUSTOM="glad/glad.h")

if(UNIX)
    # setup warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")

    # add extra checks for developer or debug builds
    if(CMAKE_BUILD_TYPE STREQUAL "" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_ASSERTIONS")
    endif()
endif()

find_package(OpenGL REQUIRED)
set(LIBRARIES ${LIBRARIES} ${OPENGL_gl_LIBRARY})

# GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(libs/glfw)
set(LIBRARIES ${LIBRARIES} glfw)

file(GLOB sources src/*.cpp src/*.h src/*.hpp src/*.tpp libs/*.h libs/*.hpp libs/*.cpp libs/*.c libs/imgui/*.h libs/imgui/*.cpp)

file(GLOB_RECURSE models RELATIVE ${CMAKE_SOURCE_DIR} models/*)
foreach(file ${models})
    configure_file(${file} ${file} COPYONLY)
endforeach(file)

add_executable(${PROJECT_NAME} ${sources} ${example_source})
target_link_libraries(${PROJECT_NAME} ${LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDES})
target_compile_definitions(${PROJECT_NAME} PRIVATE ${DEFINITIONS})
