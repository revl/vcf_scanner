#include "test_plan.hh"

static std::vector<Test_case> test_cases_sensitive_to_newline_at_eof = {
        // Unexpected EOF in the header
        {
                R"(##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">)",
                "", "E:Unexpected end of file while parsing VCF file header\n"},
        // Test line counting when there is a newline
        // after the header line.
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
)",
                "@ .", "@3\n"},
        // Test line counting when there is no newline
        // after the header line.
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO)",
                "@ .", "@2\n"},
};

static std::vector<Test_case> test_cases_insensitive_to_newline_at_eof = {
        // Not a VCF file
        {R"(text
file)",
                ".", R"(E:VCF files must start with '##fileformat'
)"},
        // Invalid meta-information line
        {R"(##fileformat=VCFv4.0
KEY)",
                ".", R"(E:Malformed meta-information line
)"},
        // Invalid meta-information line (no double-dash prefix)
        {R"(##fileformat=VCFv4.0
KEY=VALUE)",
                ".", R"(E:Malformed meta-information line
)"},
        // Missing header line
        {R"(##fileformat=VCFv4.0
1	100000	.	C	G	.	.	.)",
                ".", R"(E:Malformed meta-information line
)"},
        // Incomplete header line
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER)",
                ".", R"(E:Malformed VCF header line
)"},
        // Incorrect column name
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFORM)",
                ".", R"(E:Malformed VCF header line
)"},
        // File with no data lines
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO)",
                "HM* HG HS* .", R"(no genotypes
)"},
        // FORMAT in the header line, but no samples
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT)",
                "HM* HG HS* .", R"(with genotypes
)"},
        // The simplest of headers and a few samples
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3)",
                "HF HM* HG HS* .", R"([VCFv4.0]
with genotypes
S1
S2
S3
)"},
        // clear_line is OK at EOF
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO)",
                ". ; .", R"(;
)"},
        // Test many things at once.
        {R"(##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
1	100000	rs123;rs456	C	G	10	.	.	GT	0|1	1/.	1/0
2	200000	.	C	G,T	.	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT	0|0	0|1	1|2)",
                "HF HM* HS# @ L # A Q GF GC GT GA GT GA GT GA ; "
                "@ L A Q F I ; .",
                R"([VCFv4.0]
FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
S#=3
@4
L:1@100000
ID:[rs123,rs456]
R:C;A:G
Q:10
GF:OK
GT:OK
GT:[0,1]
GT:AVAIL
GT:[1,-1]
GT:AVAIL
GT:[1,0]
GT:NO MORE
;
@5
L:2@200000
R:C;A:[G,T]
Q:
F:PASS
I:[NS=3,DP=14,AF=0.5,DB,H2]
;
)"},
        // Missing a mandatory field
        {R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100000	.	C
1	100000	.	C	G	.	.	.
1	100000	.	C	G)",
                "@ A @ F ; @ F", R"(@3
E:Missing mandatory VCF field "ALT"
@4
F:.
;
@5
E:Missing mandatory VCF field "QUAL"
)"}};

TEST_CASE("With or without newline at EOF")
{
    for (const auto& tc : test_cases_sensitive_to_newline_at_eof) {
        run_test_case_with_and_without_cr(tc, tc.vcf);
    }
}

TEST_CASE("With and without newline at EOF")
{
    for (const auto& tc : test_cases_insensitive_to_newline_at_eof) {
        run_test_case_with_and_without_cr(tc, tc.vcf);
        run_test_case_with_and_without_cr(tc, tc.vcf + '\n');
    }
}
