#ifndef TEST_PLAN__HH
#define TEST_PLAN__HH

#include <vcf_scanner/vcf_scanner.hh>

#include "catch.hh"

#include <sstream>

namespace {

using Dump = std::stringstream;

struct Test_case {
    std::string vcf;
    std::string test_plan;
    std::string expected_result;
};

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
    std::string chrom;
    unsigned pos;
    std::vector<std::string> ids;
    std::string ref;
    std::vector<std::string> alts;
    std::string quality_str;
    bool quality_is_missing;
    std::vector<std::string> filters;

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
                        vcf_scanner.parse_loc(&chrom, &pos))) {
                dump << "L:" << chrom << '@' << pos << std::endl;
            }
            break;
        case '#':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_ids(&ids))) {
                dump << "ID:";
                dump_list(dump, ids);
                dump << std::endl;
            }
            break;
        case 'A':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_alleles(&ref, &alts))) {
                dump << "R:" << ref << ";A:";
                dump_list(dump, alts);
                dump << std::endl;
            }
            break;
        case 'Q':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_quality(
                                &quality_str, &quality_is_missing))) {
                dump << "Q:" << quality_str << std::endl;
            }
            break;
        case 'F':
            ++test_plan;
            if (dump_issues_and_clear_line(dump, vcf_scanner, test_reader,
                        vcf_scanner.parse_filters(&filters))) {
                dump << "F:";
                dump_list(dump, filters);
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

void run_test_case_with_all_buffer_sizes(
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

void run_test_case_with_and_without_cr(
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

}

#endif /* !defined(TEST_PLAN__HH) */
