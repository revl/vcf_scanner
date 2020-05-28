add_library(catch2 catch_main.cc)

foreach(TEST_NAME
	vcf_tokenizer_test
	vcf_scanner_test
)
	add_executable(${TEST_NAME} ${TEST_NAME}.cc)
	target_link_libraries(${TEST_NAME} ${PROJECT_NAME} catch2)
	add_test(${TEST_NAME} ${TEST_NAME})

	list(APPEND UNIT_TESTS ${TEST_NAME})
endforeach(TEST_NAME)

add_executable(dump_vcf dump_vcf.cc)
target_link_libraries(dump_vcf ${PROJECT_NAME})

setup_target_for_coverage_lcov(
	NAME coverage
	EXECUTABLE ctest -j ${PROCESSOR_COUNT}
	DEPENDENCIES ${UNIT_TESTS}
	EXCLUDE "/usr/*"
		"${PROJECT_SOURCE_DIR}/tests/catch.hh"
		"${PROJECT_SOURCE_DIR}/tests/dump_vcf.cc"
)
