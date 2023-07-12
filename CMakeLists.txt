cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project("cppBWQT")

# Set required C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O0 -lfftw3 -g")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
file( GLOB SRC_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file( GLOB LIB_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} lib/*.cpp)
file( GLOB LIB_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} lib/*.h)
add_library(bwlib ${LIB_SOURCES} ${LIB_HEADERS})

#########################################
#Link all library files to sources codes
foreach(srcSource ${SRC_SOURCES})
	get_filename_component(srcName ${srcSource} NAME_WE)
	add_executable(${srcName} ${srcSource})
	target_link_libraries(${srcName} bwlib)
	set_target_properties(${srcName} PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endforeach()

############################################
#Enable testing and compile Boost test files
enable_testing()
set(CTEST_OUTPUT_ON_FAILURE 1)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories (${Boost_INCLUDE_DIRS})

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)

foreach(testSrc ${TEST_SRCS})
	get_filename_component(testName ${testSrc} NAME_WE)
	add_executable(${testName} ${testSrc})
	target_link_libraries(${testName} Boost::unit_test_framework)
	target_link_libraries(${testName} bwlib)
	set_target_properties(${testName} PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin)
	add_test(NAME ${testName}
		 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin
	 	 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testName} )
endforeach(testSrc)	 
