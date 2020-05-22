#include <vcf_scanner/vcf_scanner.hh>

#include "test_case.h"

#include <sstream>

using TDump = std::stringstream;

struct STestCase {
    std::string vcf;
    std::string test_plan;
    std::string expected_result;
};

static std::vector<STestCase> s_TestCasesSensitiveToNewlineAtEOF = {
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

static std::vector<STestCase> s_TestCasesInsensitiveToNewlineAtEOF = {
        // Not a VCF file
        {R"(text
file)",
                ".", R"(E:VCF file must start with '##fileformat'
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
        // ClearLine is OK at EOF
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

class CTestReader
{
public:
    CTestReader(const std::string& vcf, size_t buf_size) :
        m_VCFData(vcf),
        m_CurrentPtr(m_VCFData.data()),
        m_EOFPtr(m_CurrentPtr + m_VCFData.length()),
        m_BufSize(buf_size)
    {}

    vcf::CVCFScanner::EParsingEvent ReadAndFeed(vcf::CVCFScanner& vcf_scanner)
    {
        for (;;) {
            size_t buf_size = m_EOFPtr - m_CurrentPtr;
            if (buf_size > m_BufSize)
                buf_size = m_BufSize;
            vcf::CVCFScanner::EParsingEvent pe =
                    vcf_scanner.Feed(m_CurrentPtr, buf_size);
            m_CurrentPtr += buf_size;
            if (pe != vcf::CVCFScanner::eNeedMoreData)
                return pe;
        }
    }

private:
    std::string m_VCFData;
    const char* m_CurrentPtr;
    const char* m_EOFPtr;
    size_t m_BufSize;
};

bool s_UpdateDump(TDump& dump, vcf::CVCFScanner& vcf_scanner,
        CTestReader& test_reader, vcf::CVCFScanner::EParsingEvent pe)
{
    if (pe == vcf::CVCFScanner::eNeedMoreData) {
        pe = test_reader.ReadAndFeed(vcf_scanner);
    }
    if (pe == vcf::CVCFScanner::eError) {
        dump << "E:" << vcf_scanner.GetError() << std::endl;
        return false;
    }
    if (pe == vcf::CVCFScanner::eOKWithWarnings) {
        for (const auto& warning : vcf_scanner.GetWarnings()) {
            dump << "W:" << warning.line_number << warning.warning_message
                 << std::endl;
        }
    }
    return true;
}

bool s_DumpIssuesAndClearLine(TDump& dump, vcf::CVCFScanner& vcf_scanner,
        CTestReader& test_reader, vcf::CVCFScanner::EParsingEvent pe)
{
    if (!s_UpdateDump(dump, vcf_scanner, test_reader, pe)) {
        s_UpdateDump(dump, vcf_scanner, test_reader, vcf_scanner.ClearLine());
        return false;
    }
    return true;
}

const char* s_DumpMetaInfo(TDump& dump,
        const vcf::CVCFHeader::TMetaInfo& meta_info, const char* test_plan)
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

const char* s_DumpHeader(
        TDump& dump, const vcf::CVCFHeader& vcf_header, const char* test_plan)
{
    switch (*test_plan) {
    case 'F':
        ++test_plan;
        dump << "[" << vcf_header.GetFileFormat() << ']' << std::endl;
        break;
    case 'M':
        test_plan =
                s_DumpMetaInfo(dump, vcf_header.GetMetaInfo(), test_plan + 1);
        break;
    case 'S':
        switch (*++test_plan) {
        case '*':
            ++test_plan;
            for (const auto& s : vcf_header.GetSampleIDs()) {
                dump << s << std::endl;
            }
            break;
        case '#':
            ++test_plan;
            dump << "S#=" << vcf_header.GetSampleIDs().size() << std::endl;
            break;
        }
        break;
    case 'G':
        ++test_plan;
        dump << (vcf_header.HasGenotypeInfo() ? "with genotypes" :
                                                "no genotypes")
             << std::endl;
    }
    return test_plan;
}

template <typename T>
void s_DumpList(TDump& dump, const T& l)
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

const char* s_DumpGenotype(TDump& dump, vcf::CVCFScanner& vcf_scanner,
        CTestReader& test_reader, const char* test_plan)
{
    switch (*test_plan) {
    case 'F':
        ++test_plan;
        if (s_DumpIssuesAndClearLine(dump, vcf_scanner, test_reader,
                    vcf_scanner.ParseGenotypeFormat())) {
            dump << "GF:OK" << std::endl;
        }
        break;
    case 'C':
        ++test_plan;
        dump << (vcf_scanner.CaptureGT() ? "GT:OK" : "GT:NOT FOUND")
             << std::endl;
        break;
    case 'T':
        ++test_plan;
        if (s_DumpIssuesAndClearLine(dump, vcf_scanner, test_reader,
                    vcf_scanner.ParseGenotype())) {
            dump << "GT:";
            s_DumpList(dump, vcf_scanner.GetGT());
            dump << std::endl;
        }
        break;
    case 'A':
        ++test_plan;
        dump << (vcf_scanner.GenotypeAvailable() ? "GT:AVAIL" : "GT:NO MORE")
             << std::endl;
    }
    return test_plan;
}

void s_IneterpretTestPlan(const char* test_plan, TDump& dump,
        vcf::CVCFScanner& vcf_scanner, CTestReader& test_reader)
{
    for (;;) {
        switch (*test_plan) {
        case '.':
            ++test_plan;
            if (!vcf_scanner.AtEOF()) {
                dump << "!EOF" << std::endl;
            }
            break;
        case ' ':
            ++test_plan;
            break;
        case '@':
            ++test_plan;
            dump << '@' << vcf_scanner.GetLineNumber() << std::endl;
            break;
        case 'H':
            test_plan =
                    s_DumpHeader(dump, vcf_scanner.GetHeader(), test_plan + 1);
            break;
        case 'L':
            ++test_plan;
            if (s_DumpIssuesAndClearLine(dump, vcf_scanner, test_reader,
                        vcf_scanner.ParseLoc())) {
                dump << "L:" << vcf_scanner.GetChrom() << '@'
                     << vcf_scanner.GetPos() << std::endl;
            }
            break;
        case '#':
            ++test_plan;
            if (s_DumpIssuesAndClearLine(dump, vcf_scanner, test_reader,
                        vcf_scanner.ParseIDs())) {
                dump << "ID:";
                s_DumpList(dump, vcf_scanner.GetIDs());
                dump << std::endl;
            }
            break;
        case 'A':
            ++test_plan;
            if (s_DumpIssuesAndClearLine(dump, vcf_scanner, test_reader,
                        vcf_scanner.ParseAlleles())) {
                dump << "R:" << vcf_scanner.GetRef() << ";A:";
                s_DumpList(dump, vcf_scanner.GetAlts());
                dump << std::endl;
            }
            break;
        case 'Q':
            ++test_plan;
            if (s_DumpIssuesAndClearLine(dump, vcf_scanner, test_reader,
                        vcf_scanner.ParseQuality())) {
                dump << "Q:" << vcf_scanner.GetQuality() << std::endl;
            }
            break;
        case 'F':
            ++test_plan;
            if (s_DumpIssuesAndClearLine(dump, vcf_scanner, test_reader,
                        vcf_scanner.ParseFilters())) {
                dump << "F:";
                s_DumpList(dump, vcf_scanner.GetFilters());
                dump << std::endl;
            }
            break;
        case 'I':
            ++test_plan;
            if (s_DumpIssuesAndClearLine(dump, vcf_scanner, test_reader,
                        vcf_scanner.ParseInfo())) {
                dump << "I:";
                s_DumpList(dump, vcf_scanner.GetInfo());
                dump << std::endl;
            }
            break;
        case 'G':
            test_plan = s_DumpGenotype(
                    dump, vcf_scanner, test_reader, test_plan + 1);
            break;
        case ';':
            ++test_plan;
            s_UpdateDump(
                    dump, vcf_scanner, test_reader, vcf_scanner.ClearLine());
            dump << ';' << std::endl;
            break;
        default:
            return;
        }
    }
}

/* TEST_CASE(debug)
{
    STestCase tc = {R"(##fileformat=VCFv4.0
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
)"};

    size_t buf_size = tc.vcf.length();

    CTestReader test_reader(tc.vcf, buf_size);

    TDump dump;

    vcf::CVCFScanner vcf_scanner;

    if (s_UpdateDump(
                dump, vcf_scanner, test_reader,
vcf::CVCFScanner::eNeedMoreData)) { s_IneterpretTestPlan( tc.test_plan.c_str(),
dump, vcf_scanner, test_reader);
    }

    std::string result = dump.str();
    if (tc.expected_result != result) {
        FILE* expected = fopen("EXPECTED", "wt");
        fwrite(tc.expected_result.data(), 1, tc.expected_result.length(),
                expected);
        fclose(expected);
        FILE* actual = fopen("ACTUAL", "wt");
        fwrite(result.data(), 1, result.length(), actual);
        fclose(actual);
    }
} */

static bool s_RunTestCaseWithAllBufferSizes(
        const STestCase& tc, const std::string& vcf)
{
    for (size_t buf_size = 1; buf_size <= vcf.length(); ++buf_size) {
        CTestReader test_reader(vcf, buf_size);

        TDump dump;

        vcf::CVCFScanner vcf_scanner;

        if (s_UpdateDump(dump, vcf_scanner, test_reader,
                    vcf::CVCFScanner::eNeedMoreData)) {
            s_IneterpretTestPlan(
                    tc.test_plan.c_str(), dump, vcf_scanner, test_reader);
        }

        std::string result = dump.str();
        if (tc.expected_result != result) {
            CHECK(tc.expected_result == result);
            FILE* expected = fopen("EXPECTED", "wt");
            fwrite(tc.expected_result.data(), 1, tc.expected_result.length(),
                    expected);
            fclose(expected);
            FILE* actual = fopen("ACTUAL", "wt");
            fwrite(result.data(), 1, result.length(), actual);
            fclose(actual);
            return false;
        }
    }
    return true;
}

static bool s_RunTestCaseWithAndWithoutCR(
        const STestCase& tc, const std::string& vcf)
{
    if (!s_RunTestCaseWithAllBufferSizes(tc, vcf)) {
        return false;
    }

    std::string vcf_with_crs = vcf;

    size_t start_pos = 0;
    while ((start_pos = vcf_with_crs.find('\n', start_pos)) !=
            std::string::npos) {
        vcf_with_crs.insert(start_pos, 1, '\r');
        start_pos += 2;
    }

    return s_RunTestCaseWithAllBufferSizes(tc, vcf_with_crs);
}

TEST_CASE(with_or_without_newline_at_eof)
{
    for (const auto& tc : s_TestCasesSensitiveToNewlineAtEOF) {
        if (!s_RunTestCaseWithAndWithoutCR(tc, tc.vcf)) {
            break;
        }
    }
}

TEST_CASE(with_and_without_newline_at_eof)
{
    for (const auto& tc : s_TestCasesInsensitiveToNewlineAtEOF) {
        if (!s_RunTestCaseWithAndWithoutCR(tc, tc.vcf) ||
                !s_RunTestCaseWithAndWithoutCR(tc, tc.vcf + '\n')) {
            break;
        }
    }
}
