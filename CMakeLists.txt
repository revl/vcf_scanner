cmake_minimum_required(VERSION 3.6)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

project(vcf_scanner VERSION 0.1.0 LANGUAGES CXX)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)

add_library(${PROJECT_NAME}
	src/vcf_scanner.cc
)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake/Modules)

if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
	set_target_properties(${PROJECT_NAME}
		PROPERTIES EXPORT_NAME VCFScanner)

	include(GNUInstallDirs)

	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
		DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	)

	install(TARGETS ${PROJECT_NAME}
		EXPORT VCFScannerTargets
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
	)

	install(EXPORT VCFScannerTargets
		FILE VCFScannerTargets.cmake
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VCFScanner
	)

	include(CMakePackageConfigHelpers)
	write_basic_package_version_file(
		${CMAKE_CURRENT_BINARY_DIR}/VCFScannerConfigVersion.cmake
		VERSION ${PROJECT_VERSION}
		COMPATIBILITY SameMajorVersion
	)

	configure_package_config_file(
		${PROJECT_SOURCE_DIR}/cmake/VCFScannerConfig.cmake.in
		${CMAKE_CURRENT_BINARY_DIR}/VCFScannerConfig.cmake
		INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VCFScanner
	)

	install(FILES ${CMAKE_CURRENT_BINARY_DIR}/VCFScannerConfig.cmake
		${CMAKE_CURRENT_BINARY_DIR}/VCFScannerConfigVersion.cmake
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VCFScanner
	)

	export(EXPORT VCFScannerTargets
		FILE ${CMAKE_CURRENT_BINARY_DIR}/VCFScannerTargets.cmake
	)

	export(PACKAGE VCFScanner)

	include(CTest)

	if(BUILD_TESTING)
		option(BUILD_COVERAGE "Enable code coverage report" OFF)

		if(BUILD_COVERAGE)
			include(CodeCoverage)

			append_coverage_compiler_flags()
		endif()

		add_subdirectory(tests)
		add_subdirectory(examples)
	endif()

	include(InstallRequiredSystemLibraries)
	set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
	set(CPACK_PACKAGE_VERSION_MAJOR "${b_VERSION_MAJOR}")
	set(CPACK_PACKAGE_VERSION_MINOR "${b_VERSION_MINOR}")
	include(CPack)
endif()
