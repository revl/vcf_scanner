cmake_minimum_required(VERSION 3.1)

project(vcf_scanner VERSION 0.1.0 LANGUAGES CXX)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)

add_library(${PROJECT_NAME}
	src/vcf_scanner.cc
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/Modules)

if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
	include(CTest)

	if(BUILD_TESTING)
		option(BUILD_COVERAGE "Enable code coverage report" OFF)

		if(BUILD_COVERAGE)
			include(CodeCoverage)

			append_coverage_compiler_flags()
		endif()

		add_subdirectory(tests)
	endif()
endif()
