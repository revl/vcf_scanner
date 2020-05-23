foreach(TEST_NAME
	vcf_tokenizer_test
	vcf_scanner_test
)
	add_executable(${TEST_NAME} ${TEST_NAME}.cc catch_main.cc)
	target_link_libraries(${TEST_NAME} ${PROJECT_NAME})
	add_test(${TEST_NAME} ${TEST_NAME})
endforeach(TEST_NAME)

add_executable(dump_vcf dump_vcf.cc)
target_link_libraries(dump_vcf ${PROJECT_NAME})
