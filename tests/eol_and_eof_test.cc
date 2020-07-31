#include "test_plan.hh"

struct Test_case {
    std::string vcf;
    std::vector<Test_check> test_plan;
};

static std::vector<Test_case> test_cases_sensitive_to_newline_at_eof = {
        // Unexpected EOF in the header
        {
                R"(##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">)",
                {{"^",
                        "E:Unexpected end of file while parsing "
                        "VCF file header"}}},

        // Test line counting when there is a newline
        // after the header line.
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
)",
                {{"^", ""}, {"@", "@3"}, {".", ""}}},

        // Test line counting when there is no newline
        // after the header line.
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO)",
                {{"^", ""}, {"@", "@2"}, {".", ""}}},
};

static std::vector<Test_case> test_cases_insensitive_to_newline_at_eof = {
        // Not a VCF file
        {R"(text
file)",
                {{"^", "E:VCF files must start with '##fileformat'"},
                        {".", "!EOF"}}},

        // Invalid meta-information line
        {"##fileformat=VCFv4.0\nKEY",
                {{"^", "E:Malformed meta-information line"}}},

        // Invalid meta-information line (no double-dash prefix)
        {"##fileformat=VCFv4.0\nKEY=VALUE",
                {{"^", "E:Malformed meta-information line"}}},

        // Missing header line
        {R"(##fileformat=VCFv4.0
1	100000	.	C	G	.	.	.)",
                {{"^", "E:Malformed meta-information line"}}},

        // Incomplete header line
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER)",
                {{"^", "E:Malformed VCF header line"}}},

        // Incorrect column name
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFORM)",
                {{"^", "E:Malformed VCF header line"}}},

        // File with no data lines
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO)",
                {{"^", ""}, {"HM*", ""}, {"HG", "no genotypes"}, {"HS*", "."},
                        {".", ""}}},

        // FORMAT in the header line, but no samples
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT)",
                {{"^", ""}, {"HM*", ""}, {"HG", "with genotypes"}, {"HS*", "."},
                        {".", ""}}},

        // The simplest of headers and a few samples
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3)",
                {{"^", ""}, {"HF", "[VCFv4.0]"}, {"HM*", ""},
                        {"HG", "with genotypes"}, {"HS*", "[S1,S2,S3]"},
                        {".", ""}}},

        // clear_line is OK at EOF
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO)",
                {{"^", ""}, {".", ""}, {";", ";"}, {".", ""}}},

        // Test many things at once.
        {R"(##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
1	100000	rs123;rs456	C	G	10	.	.	GT	0|1	1/.	1/0
2	200000	.	C	G,T	.	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT	0|0	0|1	1|2)",
                {
                        {"^", ""},
                        {"HF", "[VCFv4.0]"},
                        {"HM*",
                                "FORMAT=<ID=GT,Number=1,Type=String,"
                                "Description=\"Genotype\">\n"},
                        {"HS#", "S#=3"},
                        {"@", "@4"},
                        {"L", "L:1@100000"},
                        {"#", "ID:[rs123,rs456]"},
                        {"A", "R:C;A:G"},
                        {"Q", "Q:10"},
                        {"GF", "GF:OK"},
                        {"GC", "GT:OK"},
                        {"GT", "GT:[0,1]"},
                        {"GA", "GT:AVAIL"},
                        {"GT", "GT:[1,-1]"},
                        {"GA", "GT:AVAIL"},
                        {"GT", "GT:[1,0]"},
                        {"GA", "GT:NO MORE"},
                        {";", ";"},
                        {"@", "@5"},
                        {"L", "L:2@200000"},
                        {"A", "R:C;A:[G,T]"},
                        {"Q", "Q:"},
                        {"F", "F:PASS"},
                        {"I", "I:[NS=3,DP=14,AF=0.5,DB,H2]"},
                        {";", ";"},
                        {".", ""},

                }},

        // Missing a mandatory field
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100000	.	C
1	100000	.	C	G	.	.	.
1	100000	.	C	G)",
                {
                        {"^", ""},
                        {"@", "@3"},
                        {"A", "E:Missing mandatory VCF field \"ALT\""},
                        {"@", "@4"},
                        {"F", "F:."},
                        {";", ";"},
                        {"@", "@5"},
                        {"F", "E:Missing mandatory VCF field \"QUAL\""},
                }}};

TEST_CASE("With or without newline at EOF")
{
    for (const auto& tc : test_cases_sensitive_to_newline_at_eof) {
        run_test_case_with_and_without_cr(tc.vcf, tc.test_plan);
    }
}

TEST_CASE("With and without newline at EOF")
{
    for (const auto& tc : test_cases_insensitive_to_newline_at_eof) {
        run_test_case_with_and_without_cr(tc.vcf, tc.test_plan);
        run_test_case_with_and_without_cr(tc.vcf + '\n', tc.test_plan);
    }
}
