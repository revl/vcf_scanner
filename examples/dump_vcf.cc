// This example parses the specified VCF file and prints the extracted data
// to the standard output stream.

#include <vcf_scanner/vcf_scanner.hh>

#include <iostream>

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

    static std::array<char, 1024 * 1024> buffer;

    VCF_scanner vcf_scanner;

    auto parse_to_completion = [&](VCF_parsing_event pe) {
        while (pe == VCF_parsing_event::need_more_data) {
            pe = vcf_scanner.feed(buffer.data(),
                    fread(buffer.data(), 1, buffer.size(), input));
        }

        if (pe == VCF_parsing_event::error) {
            return false;
        }

        if (pe == VCF_parsing_event::ok_with_warnings) {
            for (const auto& warning : vcf_scanner.get_warnings()) {
                std::cerr << "Warning: " << warning.warning_message
                          << std::endl;
            }
        }

        return true;
    };

    VCF_header header;

    if (!parse_to_completion(vcf_scanner.parse_header(&header))) {
        return 1;
    }

    std::cout << "##fileformat=" << header.get_file_format_version()
              << std::endl;

    for (const auto& kv : header.get_meta_info()) {
        for (const auto& v : kv.second) {
            std::cout << "##" << kv.first << '=' << v << std::endl;
        }
    }

    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& v : header.get_sample_ids()) {
        std::cout << '\t' << v;
    }
    std::cout << std::endl;

    std::string chrom;
    unsigned pos;
    std::vector<std::string> ids;
    std::string ref;
    std::vector<std::string> alts;
    std::string quality_str;
    std::vector<std::string> filters;

    auto parse_data_line = [&] {
        const char* sep;

        if (!parse_to_completion(vcf_scanner.parse_loc(&chrom, &pos))) {
            return false;
        }
        std::cout << chrom << '\t' << pos;

        std::vector<std::string> ids;
        if (!parse_to_completion(vcf_scanner.parse_ids(&ids))) {
            return false;
        }
        if (!ids.empty()) {
            sep = "\t";
            for (const auto& id : ids) {
                std::cout << sep << id;
                sep = ",";
            }
            std::cout << '\t';
        } else {
            std::cout << "\t.\t";
        }

        if (!parse_to_completion(vcf_scanner.parse_alleles(&ref, &alts))) {
            return false;
        }
        std::cout << ref;
        if (!alts.empty()) {
            sep = "\t";
            for (const auto& alt : alts) {
                std::cout << sep << alt;
                sep = ",";
            }
            std::cout << '\t';
        } else {
            std::cout << "\t.\t";
        }

        if (!parse_to_completion(vcf_scanner.parse_quality(&quality_str))) {
            return false;
        }
        if (!quality_str.empty()) {
            std::cout << quality_str;
        } else {
            std::cout << ".";
        }

        if (!parse_to_completion(vcf_scanner.parse_filters(&filters))) {
            return false;
        }
        if (!filters.empty()) {
            sep = "\t";
            for (const auto& filter : filters) {
                std::cout << sep << filter;
                sep = ";";
            }
        } else {
            std::cout << "\t.";
        }

        if (!parse_to_completion(vcf_scanner.parse_info())) {
            return false;
        }
        if (!vcf_scanner.get_info().empty()) {
            sep = "\t";
            for (const auto& info_item : vcf_scanner.get_info()) {
                std::cout << sep << info_item;
                sep = ";";
            }
            std::cout << '\t';
        } else {
            std::cout << "\t.\t";
        }

        if (header.has_genotype_info()) {
            if (!parse_to_completion(vcf_scanner.parse_genotype_format())) {
                return false;
            }

            if (!vcf_scanner.capture_gt()) {
                std::cout << std::endl;
                std::cerr << "\tERR: no GT key" << std::endl;
                return true;
            }

            std::cout << "GT";

            while (vcf_scanner.genotype_available()) {
                if (!parse_to_completion(vcf_scanner.parse_genotype())) {
                    return false;
                }
                sep = "\t";
                for (auto allele : vcf_scanner.get_gt()) {
                    std::cout << sep;
                    if (allele < 0) {
                        std::cout << '.';
                    } else {
                        std::cout << allele;
                    }
                    sep = vcf_scanner.is_phased_gt() ? "|" : "/";
                }
            }
        }

        std::cout << std::endl;

        return true;
    };

    while (!vcf_scanner.at_eof()) {
        if (!parse_data_line()) {
            std::cout << std::endl;
            std::cerr << "<-ERR@" << vcf_scanner.get_line_number() << ": "
                      << vcf_scanner.get_error() << std::endl;
        }

        parse_to_completion(vcf_scanner.clear_line());
    }
}
