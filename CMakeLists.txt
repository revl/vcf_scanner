cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

project(vcf_scanner VERSION 0.1.0 LANGUAGES CXX)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)

add_library(${PROJECT_NAME}
	src/vcf_scanner.cc
)

target_include_directories(${PROJECT_NAME} PUBLIC
	${PROJECT_SOURCE_DIR}/include)

if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
	include(CTest)

	if(BUILD_TESTING)
		add_subdirectory(tests)
	endif()
endif()
