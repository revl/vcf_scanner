#include <vcf_scanner/vcf_scanner.hh>

static char buffer[1];
static size_t buffer_size;

static void s_ReadBuffer(FILE* input)
{
    buffer_size = fread(buffer, 1, sizeof(buffer), input);
}

static bool s_ParseToCompletion(vcf::CVCFScanner::EParsingEvent pe,
        vcf::CVCFScanner& vcf_scanner, FILE* input)
{
    while (pe == vcf::CVCFScanner::eNeedMoreData) {
        s_ReadBuffer(input);
        pe = vcf_scanner.Feed(buffer, buffer_size);
    }

    if (pe == vcf::CVCFScanner::eError) {
        return false;
    }

    if (pe == vcf::CVCFScanner::eOKWithWarnings) {
        for (const auto& warning : vcf_scanner.GetWarnings())
            std::cerr << "Warning: " << warning.warning_message << std::endl;
    }

    return true;
}

static bool ParseDataLine(vcf::CVCFScanner& vcf_scanner, FILE* input)
{
    std::cout << vcf_scanner.GetLineNumber() << ':';

    const char* sep;

    if (!s_ParseToCompletion(vcf_scanner.ParseLoc(), vcf_scanner, input)) {
        return false;
    }
    std::cout << vcf_scanner.GetChrom() << '\t' << vcf_scanner.GetPos();

    std::vector<std::string> ids;
    if (!s_ParseToCompletion(vcf_scanner.ParseIDs(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.GetIDs().empty()) {
        sep = "\t";
        for (const auto& id : vcf_scanner.GetIDs()) {
            std::cout << sep << id;
            sep = ",";
        }
        std::cout << '\t';
    } else
        std::cout << "\t.\t";

    if (!s_ParseToCompletion(vcf_scanner.ParseAlleles(), vcf_scanner, input)) {
        return false;
    }
    std::cout << vcf_scanner.GetRef();
    if (!vcf_scanner.GetAlts().empty()) {
        sep = "\t";
        for (const auto& alt : vcf_scanner.GetAlts()) {
            std::cout << sep << alt;
            sep = ",";
        }
        std::cout << '\t';
    } else
        std::cout << "\t.\t";

    if (!s_ParseToCompletion(vcf_scanner.ParseQuality(), vcf_scanner, input)) {
        return false;
    }
    std::string quality = vcf_scanner.GetQuality();
    if (!quality.empty())
        std::cout << quality;
    else
        std::cout << ".";

    if (!s_ParseToCompletion(vcf_scanner.ParseFilters(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.GetFilters().empty()) {
        sep = "\t";
        for (const auto& filter : vcf_scanner.GetFilters()) {
            std::cout << sep << filter;
            sep = ";";
        }
    } else
        std::cout << "\t.";

    if (!s_ParseToCompletion(vcf_scanner.ParseInfo(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.GetInfo().empty()) {
        sep = "\t";
        for (const auto& info_item : vcf_scanner.GetInfo()) {
            std::cout << sep << info_item;
            sep = ";";
        }
        std::cout << '\t';
    } else
        std::cout << "\t.\t";

    if (vcf_scanner.GetHeader().HasGenotypeInfo()) {
        if (!s_ParseToCompletion(
                    vcf_scanner.ParseGenotypeFormat(), vcf_scanner, input)) {
            return false;
        }

        if (!vcf_scanner.CaptureGT()) {
            std::cout << std::endl;
            std::cout << "\tERR: no GT key" << std::endl;
            return true;
        }

        std::cout << "GT";

        while (vcf_scanner.GenotypeAvailable()) {
            if (!s_ParseToCompletion(
                        vcf_scanner.ParseGenotype(), vcf_scanner, input)) {
                return false;
            }
            sep = "\t";
            for (auto allele : vcf_scanner.GetGT()) {
                std::cout << sep << allele;
                sep = vcf_scanner.IsPhasedGT() ? "|" : "/";
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

    vcf::CVCFScanner vcf_scanner;

    vcf::CVCFScanner::EParsingEvent pe;

    // Read the header
    do
        s_ReadBuffer(input);
    while ((pe = vcf_scanner.Feed(buffer, buffer_size)) ==
            vcf::CVCFScanner::eNeedMoreData);

    if (pe != vcf::CVCFScanner::eOK) {
        std::cerr << vcf_scanner.GetError() << std::endl;
        return 1;
    }

    if (pe == vcf::CVCFScanner::eOKWithWarnings) {
        for (const auto& warning : vcf_scanner.GetWarnings())
            std::cerr << "Warning: " << warning.warning_message << std::endl;
        pe = vcf::CVCFScanner::eOK;
    }

    std::cout << "##fileformat=" << vcf_scanner.GetHeader().GetFileFormat()
              << std::endl;

    for (const auto& kv : vcf_scanner.GetHeader().GetMetaInfo()) {
        for (const auto& v : kv.second) {
            std::cout << "##" << kv.first << '=' << v << std::endl;
        }
    }

    std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& v : vcf_scanner.GetHeader().GetSampleIDs()) {
        std::cout << '\t' << v;
    }
    std::cout << std::endl;

    while (!vcf_scanner.AtEOF()) {
        if (!ParseDataLine(vcf_scanner, input)) {
            std::cout << std::endl;
            std::cerr << "<-ERR@" << vcf_scanner.GetLineNumber() << ": "
                      << vcf_scanner.GetError() << std::endl;
        }

        s_ParseToCompletion(vcf_scanner.ClearLine(), vcf_scanner, input);
    }

    return 0;
}
