#include <vcf_scanner/vcf_scanner.hh>

#include <sstream>

#include "catch.hh"

using Dump = std::stringstream;

struct Test_case {
    std::string vcf;
    std::string test_plan;
    std::string expected_result;
};

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

class Test_reader
{
public:
    Test_reader(const std::string& vcf, size_t cs) :
        vcf_data(vcf),
        current_ptr(vcf_data.data()),
        eof_ptr(current_ptr + vcf_data.length()),
        chunk_size(cs)
    {}

    VCF_parsing_event read_and_feed(VCF_scanner& vcf_scanner)
    {
        for (;;) {
            size_t buf_size = eof_ptr - current_ptr;
            if (buf_size > chunk_size) {
                buf_size = chunk_size;
            }
            VCF_parsing_event pe = vcf_scanner.feed(current_ptr, buf_size);
            current_ptr += buf_size;
            if (pe != VCF_parsing_event::need_more_data) {
                return pe;
            }
        }
    }

private:
    const std::string vcf_data;
    const char* current_ptr;
    const char* eof_ptr;
    const size_t chunk_size;
};

bool update_dump(Dump& dump, VCF_scanner& vcf_scanner, Test_reader& test_reader,
        VCF_parsing_event pe)
{
    if (pe == VCF_parsing_event::need_more_data) {
        pe = test_reader.read_and_feed(vcf_scanner);
    }
    if (pe == VCF_parsing_event::error) {
        dump << "E:" << vcf_scanner.get_error() << std::endl;
        return false;
    }
    if (pe == VCF_parsing_event::ok_with_warnings) {
        for (const auto& warning : vcf_scanner.get_warnings()) {
            dump << "W:" << warning.line_number << warning.warning_message
                 << std::endl;
        }
    }
    return true;
}

bool dump_issues_and_clear_line(Dump& dump, VCF_scanner& vcf_scanner,
        Test_reader& test_reader, VCF_parsing_event pe)
{
    if (!update_dump(dump, vcf_scanner, test_reader, pe)) {
        update_dump(dump, vcf_scanner, test_reader, vcf_scanner.clear_line());
        return false;
    }
    return true;
}

const char* dump_meta_info(Dump& dump, const VCF_header::Meta_info& meta_info,
        const char* test_plan)
{
    switch (*test_plan) {
    case '*':
        ++test_plan;
        for (const auto& kv : meta_info) {
            for (const auto& v : kv.second) {
                dump << kv.first << '=' << v << std::endl;
            }
        }
        return test_plan;
    case '{':
        do {
            size_t key_len = strcspn(++test_plan, ",}");
            if (key_len > 0) {
                std::string key(test_plan, key_len);
                auto key_iter = meta_info.find(key);
                if (key_iter == meta_info.end()) {
                    dump << key << ": NOT FOUND" << std::endl;
                } else {
                    for (const auto& v : key_iter->second) {
                        dump << key << '=' << v << std::endl;
                    }
                }
                test_plan += key_len;
            }
        } while (*test_plan == ',');
        ++test_plan;
    }
    return test_plan;
}

const char* dump_header(
        Dump& dump, const VCF_header& vcf_header, const char* test_plan)
{
    switch (*test_plan) {
    case 'F':
        ++test_plan;
        dump << "[" << vcf_header.get_file_format_version() << ']' << std::endl;
        break;
    case 'M':
        test_plan =
                dump_meta_info(dump, vcf_header.get_meta_info(), test_plan + 1);
        break;
    case 'S':
        switch (*++test_plan) {
        case '*':
            ++test_plan;
            for (const auto& s : vcf_header.get_sample_ids()) {
                dump << s << std::endl;
            }
            break;
        case '#':
            ++test_plan;
            dump << "S#=" << vcf_header.get_sample_ids().size() << std::endl;
            break;
        }
        break;
    case 'G':
        ++test_plan;
        dump << (vcf_header.has_genotype_info() ? "with genotypes" :
                                                  "no genotypes")
             << std::endl;
    }
    return test_plan;
}

template <typename T>
void dump_list(Dump& dump, const T& l)
{
    if (l.empty()) {
        dump << ".";
    } else if (l.size() == 1) {
        dump << l.front();
    } else {
        auto iter = l.begin();
        dump << '[' << *iter;
        while (++iter != l.end()) {
            dump << ',' << *iter;
        }
        dump << ']';
    }
}

const char* dump_genotype(Dump& dump, VCF_scanner& vcf_scanner,
        Test_reader& test_reader, const char* test_plan)
{
    switch (*test_plan) {
    case 'F':
        ++test_plan;
        if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                    vcf_scanner.parse_genotype_format())) {
            dump << "GF:OK" << std::endl;
        }
        break;
    case 'C':
        ++test_plan;
        dump << (vcf_scanner.capture_gt() ? "GT:OK" : "GT:NOT FOUND")
             << std::endl;
        break;
    case 'T':
        ++test_plan;
        if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                    vcf_scanner.parse_genotype())) {
            dump << "GT:";
            dump_list(dump, vcf_scanner.get_gt());
            dump << std::endl;
        }
        break;
    case 'A':
        ++test_plan;
        dump << (vcf_scanner.genotype_available() ? "GT:AVAIL" : "GT:NO MORE")
             << std::endl;
    }
    return test_plan;
}

void ineterpret_test_plan(const char* test_plan, Dump& dump,
        VCF_scanner& vcf_scanner, Test_reader& test_reader)
{
    for (;;) {
        switch (*test_plan) {
        case '.':
            ++test_plan;
            if (!vcf_scanner.at_eof()) {
                dump << "!EOF" << std::endl;
            }
            break;
        case ' ':
            ++test_plan;
            break;
        case '@':
            ++test_plan;
            dump << '@' << vcf_scanner.get_line_number() << std::endl;
            break;
        case 'H':
            test_plan =
                    dump_header(dump, vcf_scanner.get_header(), test_plan + 1);
            break;
        case 'L':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_loc())) {
                dump << "L:" << vcf_scanner.get_chrom() << '@'
                     << vcf_scanner.get_pos() << std::endl;
            }
            break;
        case '#':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_ids())) {
                dump << "ID:";
                dump_list(dump, vcf_scanner.get_ids());
                dump << std::endl;
            }
            break;
        case 'A':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_alleles())) {
                dump << "R:" << vcf_scanner.get_ref() << ";A:";
                dump_list(dump, vcf_scanner.get_alts());
                dump << std::endl;
            }
            break;
        case 'Q':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_quality())) {
                dump << "Q:" << vcf_scanner.get_quality_as_string()
                     << std::endl;
            }
            break;
        case 'F':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_filters())) {
                dump << "F:";
                dump_list(dump, vcf_scanner.get_filters());
                dump << std::endl;
            }
            break;
        case 'I':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_info())) {
                dump << "I:";
                dump_list(dump, vcf_scanner.get_info());
                dump << std::endl;
            }
            break;
        case 'G':
            test_plan = dump_genotype(
                    dump, vcf_scanner, test_reader, test_plan + 1);
            break;
        case ';':
            ++test_plan;
            update_dump(
                    dump, vcf_scanner, test_reader, vcf_scanner.clear_line());
            dump << ';' << std::endl;
            break;
        default:
            return;
        }
    }
}

static void run_test_case_with_all_buffer_sizes(
        const Test_case& tc, const std::string& vcf)
{
    for (size_t buf_size = 1; buf_size <= vcf.length(); ++buf_size) {
        Test_reader test_reader(vcf, buf_size);

        Dump dump;

        VCF_scanner vcf_scanner;

        if (update_dump(dump, vcf_scanner, test_reader,
                    VCF_parsing_event::need_more_data)) {
            ineterpret_test_plan(
                    tc.test_plan.c_str(), dump, vcf_scanner, test_reader);
        }

        CHECK(tc.expected_result == dump.str());
    }
}

static void run_test_case_with_and_without_cr(
        const Test_case& tc, const std::string& vcf)
{
    run_test_case_with_all_buffer_sizes(tc, vcf);

    std::string vcf_with_crs = vcf;

    size_t start_pos = 0;
    while ((start_pos = vcf_with_crs.find('\n', start_pos)) !=
            std::string::npos) {
        vcf_with_crs.insert(start_pos, 1, '\r');
        start_pos += 2;
    }

    run_test_case_with_all_buffer_sizes(tc, vcf_with_crs);
}

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
