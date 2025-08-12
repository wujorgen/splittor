cmake_minimum_required(VERSION 3.15...4.0)

project(
	splittor
	VERSION 0.0
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_library(
	libsplittor
    src/splittor.cpp
    src/InputReader.cpp
    src/InputProcessor.cpp
	src/SolverPoisson.cpp
)

target_include_directories(libsplittor PUBLIC src)

add_executable(splittor src/main.cpp)
target_link_libraries(splittor libsplittor)

file(GLOB_RECURSE TEST_SOURCES "unit_tests/*.cpp")
message("${TEST_SOURCES}")
add_executable(testsplittor ${TEST_SOURCES})
target_link_libraries(testsplittor libsplittor gtest_main)

# shenanigans to link Eigen
find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
	target_link_libraries(libsplittor Eigen3::Eigen)
	#target_link_libraries(testsplittor Eigen3::Eigen)
else()
	if (APPLE)
		target_include_directories(libsplittor PUBLIC /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/)
	elseif(UNIX)
		target_include_directories(libsplittor PUBLIC /usr/local/include/eigen-3.4.0)
	endif()
endif()

target_compile_options(splittor PRIVATE -march=native)


enable_testing()
include(GoogleTest)
gtest_discover_tests(testsplittor)