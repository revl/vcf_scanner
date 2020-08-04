#include "test_plan.hh"

TEST_CASE("List parsing")
{
    run_test_case_with_and_without_cr(R"(##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	1000	.	C	.	.	.	.
1	1000	ID1;.	C	G,.	.	F1;.	.
1	1000	.;ID1	C	.,G	.	.;F1	.)",
            {
                    {"^", ""},
                    {"#", "ID:[]"},
                    {"A", "R:C;A:[]"},
                    {"F", "F:[]"},
                    {";", ";"},
                    {"#", "ID:[ID1]"},
                    {"A", "R:C;A:[G]"},
                    {"F", "F:[F1]"},
                    {";", ";"},
                    {"#", "ID:[ID1]"},
                    {"A", "R:C;A:[G]"},
                    {"F", "F:[F1]"},
            });
}
