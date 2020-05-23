#include <vcf_scanner/vcf_scanner.hh>

static char buffer[1];
static size_t buffer_size;

static void read_buffer(FILE* input)
{
    buffer_size = fread(buffer, 1, sizeof(buffer), input);
}

static bool parse_to_completion(
        VCF_scanner::Parsing_event pe, VCF_scanner& vcf_scanner, FILE* input)
{
    while (pe == VCF_scanner::need_more_data) {
        read_buffer(input);
        pe = vcf_scanner.feed(buffer, buffer_size);
    }

    if (pe == VCF_scanner::error) {
        return false;
    }

    if (pe == VCF_scanner::ok_with_warnings) {
        for (const auto& warning : vcf_scanner.get_warnings())
            std::cerr << "Warning: " << warning.warning_message << std::endl;
    }

    return true;
}

static bool parse_data_line(VCF_scanner& vcf_scanner, FILE* input)
{
    std::cout << vcf_scanner.get_line_number() << ':';

    const char* sep;

    if (!parse_to_completion(vcf_scanner.parse_loc(), vcf_scanner, input)) {
        return false;
    }
    std::cout << vcf_scanner.get_chrom() << '\t' << vcf_scanner.get_pos();

    std::vector<std::string> ids;
    if (!parse_to_completion(vcf_scanner.parse_ids(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.get_ids().empty()) {
        sep = "\t";
        for (const auto& id : vcf_scanner.get_ids()) {
            std::cout << sep << id;
            sep = ",";
        }
        std::cout << '\t';
    } else
        std::cout << "\t.\t";

    if (!parse_to_completion(vcf_scanner.parse_alleles(), vcf_scanner, input)) {
        return false;
    }
    std::cout << vcf_scanner.get_ref();
    if (!vcf_scanner.get_alts().empty()) {
        sep = "\t";
        for (const auto& alt : vcf_scanner.get_alts()) {
            std::cout << sep << alt;
            sep = ",";
        }
        std::cout << '\t';
    } else
        std::cout << "\t.\t";

    if (!parse_to_completion(vcf_scanner.parse_quality(), vcf_scanner, input)) {
        return false;
    }
    std::string quality = vcf_scanner.get_quality();
    if (!quality.empty())
        std::cout << quality;
    else
        std::cout << ".";

    if (!parse_to_completion(vcf_scanner.parse_filters(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.get_filters().empty()) {
        sep = "\t";
        for (const auto& filter : vcf_scanner.get_filters()) {
            std::cout << sep << filter;
            sep = ";";
        }
    } else
        std::cout << "\t.";

    if (!parse_to_completion(vcf_scanner.parse_info(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.get_info().empty()) {
        sep = "\t";
        for (const auto& info_item : vcf_scanner.get_info()) {
            std::cout << sep << info_item;
            sep = ";";
        }
        std::cout << '\t';
    } else
        std::cout << "\t.\t";

    if (vcf_scanner.get_header().has_genotype_info()) {
        if (!parse_to_completion(
                    vcf_scanner.parse_genotype_format(), vcf_scanner, input)) {
            return false;
        }

        if (!vcf_scanner.capture_gt()) {
            std::cout << std::endl;
            std::cout << "\tERR: no GT key" << std::endl;
            return true;
        }

        std::cout << "GT";

        while (vcf_scanner.genotype_available()) {
            if (!parse_to_completion(
                        vcf_scanner.parse_genotype(), vcf_scanner, input)) {
                return false;
            }
            sep = "\t";
            for (auto allele : vcf_scanner.get_gt()) {
                std::cout << sep << allele;
                sep = vcf_scanner.is_phased_gt() ? "|" : "/";
            }
        }
    }

    std::cout << std::endl;

    return true;
}

int main(int argc, const char* argv[])
{
    if (argc != 2) {
        fprintf(stderr, "Usage %s VCF_FILE\n", *argv);
        return 2;
    }

    FILE* input = fopen(argv[1], "rb");
    if (input == nullptr) {
        perror(argv[1]);
        return 1;
    }

    VCF_scanner vcf_scanner;

    VCF_scanner::Parsing_event pe;

    // Read the header
    do
        read_buffer(input);
    while ((pe = vcf_scanner.feed(buffer, buffer_size)) ==
            VCF_scanner::need_more_data);

    if (pe != VCF_scanner::ok) {
        std::cerr << vcf_scanner.get_error() << std::endl;
        return 1;
    }

    if (pe == VCF_scanner::ok_with_warnings) {
        for (const auto& warning : vcf_scanner.get_warnings())
            std::cerr << "Warning: " << warning.warning_message << std::endl;
        pe = VCF_scanner::ok;
    }

    std::cout << "##fileformat="
              << vcf_scanner.get_header().get_file_format_version()
              << std::endl;

    for (const auto& kv : vcf_scanner.get_header().get_meta_info()) {
        for (const auto& v : kv.second) {
            std::cout << "##" << kv.first << '=' << v << std::endl;
        }
    }

    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& v : vcf_scanner.get_header().get_sample_ids()) {
        std::cout << '\t' << v;
    }
    std::cout << std::endl;

    while (!vcf_scanner.at_eof()) {
        if (!parse_data_line(vcf_scanner, input)) {
            std::cout << std::endl;
            std::cerr << "<-ERR@" << vcf_scanner.get_line_number() << ": "
                      << vcf_scanner.get_error() << std::endl;
        }

        parse_to_completion(vcf_scanner.clear_line(), vcf_scanner, input);
    }

    return 0;
}
