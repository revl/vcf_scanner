#include "test_plan.hh"

TEST_CASE("Unexpected EOF in the header")
{
    run_test_case_with_and_without_cr(R"(##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">)",
            {{"^",
                    "E:Unexpected end of file while parsing "
                    "VCF file header"}});
}

TEST_CASE("Line counting when there is a newline after the header line")
{
    run_test_case_with_and_without_cr(R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
)",
            {{"^", ""}, {"@", "@3"}, {".", ""}});
}

static void run_test_cases_insensitive_to_newline_at_eof(
        const std::string& vcf, const std::vector<Test_check>& test_plan)
{
    run_test_case_with_and_without_cr(vcf, test_plan);
    run_test_case_with_and_without_cr(vcf + '\n', test_plan);
}

TEST_CASE("Not a VCF file")
{
    run_test_cases_insensitive_to_newline_at_eof("text\nfile",
            {{"^", "E:VCF files must start with '##fileformat'"},
                    {".", "!EOF"}});
}

TEST_CASE("Invalid meta-information line")
{
    run_test_cases_insensitive_to_newline_at_eof("##fileformat=VCFv4.0\nKEY",
            {{"^", "E:Malformed meta-information line"}});
}

TEST_CASE("Invalid meta-information line (no double-dash prefix)")
{
    run_test_cases_insensitive_to_newline_at_eof(
            "##fileformat=VCFv4.0\nKEY=VALUE",
            {{"^", "E:Malformed meta-information line"}});
}

TEST_CASE("Missing header line")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
1	100000	.	C	G	.	.	.)",
            {{"^", "E:Malformed meta-information line"}});
}

TEST_CASE("Incomplete header line")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER)",
            {{"^", "E:Malformed VCF header line"}});
}

TEST_CASE("Incorrect column name")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFORM)",
            {{"^", "E:Malformed VCF header line"}});
}

TEST_CASE("File with no data lines")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO)",
            {{"^", ""}, {"HM*", ""}, {"HG", "no genotypes"}, {"HS*", "."},
                    {".", ""}});
}

TEST_CASE("FORMAT in the header line, but no samples")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT)",
            {{"^", ""}, {"HM*", ""}, {"HG", "with genotypes"}, {"HS*", "."},
                    {".", ""}});
}

TEST_CASE("The simplest of headers and a few samples")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3)",
            {{"^", ""}, {"HF", "[VCFv4.0]"}, {"HM*", ""},
                    {"HG", "with genotypes"}, {"HS*", "[S1,S2,S3]"},
                    {".", ""}});
}

TEST_CASE("clear_line is OK at EOF")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO)",
            {{"^", ""}, {".", ""}, {";", ";"}, {".", ""}});
}

TEST_CASE("Test many things at once")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
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

            });
}

TEST_CASE("Missing a mandatory field")
{
    run_test_cases_insensitive_to_newline_at_eof(R"(##fileformat=VCFv4.0
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
            });
}
