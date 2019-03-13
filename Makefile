CXXFLAGS = -O0 -ggdb -Wall -Wextra -Werror -pedantic --coverage

test_vcfscanner: test_vcfscanner.cc vcfscanner.cc vcfscanner.hh tokenizer.hh
	g++ -o $@ ${CXXFLAGS} test_vcfscanner.cc vcfscanner.cc

dump_vcf: dump_vcf.cc vcfscanner.cc vcfscanner.hh tokenizer.hh
	g++ -o $@ ${CXXFLAGS} dump_vcf.cc vcfscanner.cc

test_tokenizer: test_tokenizer.cc tokenizer.hh
	g++ -o $@ ${CXXFLAGS} $<
