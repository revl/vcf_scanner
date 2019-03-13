#include "vcfscanner.hh"

static char buffer[1];
static size_t buffer_size;

static void s_ReadBuffer(FILE* input)
{
    buffer_size = fread(buffer, 1, sizeof(buffer), input);
}

static bool s_ParseToCompletion(
        CVCFScanner::EParsingEvent pe, CVCFScanner& vcf_scanner, FILE* input)
{
    while (pe == CVCFScanner::eNeedMoreData) {
        s_ReadBuffer(input);
        pe = vcf_scanner.Feed(buffer, buffer_size);
    }

    if (pe == CVCFScanner::eError) {
        return false;
    }

    if (pe == CVCFScanner::eOKWithWarnings) {
        for (const auto& warning : vcf_scanner.GetWarnings())
            cerr << "Warning: " << warning.warning_message << endl;
    }

    return true;
}

static bool ParseDataLine(CVCFScanner& vcf_scanner, FILE* input)
{
    cout << vcf_scanner.GetLineNumber() << ':';

    const char* sep;

    if (!s_ParseToCompletion(vcf_scanner.ParseLoc(), vcf_scanner, input)) {
        return false;
    }
    cout << vcf_scanner.GetChrom() << '\t' << vcf_scanner.GetPos();

    vector<string> ids;
    if (!s_ParseToCompletion(vcf_scanner.ParseIDs(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.GetIDs().empty()) {
        sep = "\t";
        for (const auto& id : vcf_scanner.GetIDs()) {
            cout << sep << id;
            sep = ",";
        }
        cout << '\t';
    } else
        cout << "\t.\t";

    if (!s_ParseToCompletion(vcf_scanner.ParseAlleles(), vcf_scanner, input)) {
        return false;
    }
    cout << vcf_scanner.GetRef();
    if (!vcf_scanner.GetAlts().empty()) {
        sep = "\t";
        for (const auto& alt : vcf_scanner.GetAlts()) {
            cout << sep << alt;
            sep = ",";
        }
        cout << '\t';
    } else
        cout << "\t.\t";

    if (!s_ParseToCompletion(vcf_scanner.ParseQuality(), vcf_scanner, input)) {
        return false;
    }
    string quality = vcf_scanner.GetQuality();
    if (!quality.empty())
        cout << quality;
    else
        cout << ".";

    if (!s_ParseToCompletion(vcf_scanner.ParseFilters(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.GetFilters().empty()) {
        sep = "\t";
        for (const auto& filter : vcf_scanner.GetFilters()) {
            cout << sep << filter;
            sep = ";";
        }
    } else
        cout << "\t.";

    if (!s_ParseToCompletion(vcf_scanner.ParseInfo(), vcf_scanner, input)) {
        return false;
    }
    if (!vcf_scanner.GetInfo().empty()) {
        sep = "\t";
        for (const auto& info_item : vcf_scanner.GetInfo()) {
            cout << sep << info_item;
            sep = ";";
        }
        cout << '\t';
    } else
        cout << "\t.\t";

    if (vcf_scanner.GetHeader().HasGenotypeInfo()) {
        if (!s_ParseToCompletion(
                    vcf_scanner.ParseGenotypeFormat(), vcf_scanner, input)) {
            return false;
        }

        if (!vcf_scanner.CaptureGT()) {
            cout << endl;
            cout << "\tERR: no GT key" << endl;
            return true;
        }

        cout << "GT";

        while (vcf_scanner.GenotypeAvailable()) {
            if (!s_ParseToCompletion(
                        vcf_scanner.ParseGenotype(), vcf_scanner, input)) {
                return false;
            }
            sep = "\t";
            for (auto allele : vcf_scanner.GetGT()) {
                cout << sep << allele;
                sep = vcf_scanner.IsPhasedGT() ? "|" : "/";
            }
        }
    }

    cout << endl;

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

    CVCFScanner vcf_scanner;

    CVCFScanner::EParsingEvent pe;

    // Read the header
    do
        s_ReadBuffer(input);
    while ((pe = vcf_scanner.Feed(buffer, buffer_size)) ==
            CVCFScanner::eNeedMoreData);

    if (pe != CVCFScanner::eOK) {
        cerr << vcf_scanner.GetError().m_ErrorMessage << endl;
        return 1;
    }

    if (pe == CVCFScanner::eOKWithWarnings) {
        for (const auto& warning : vcf_scanner.GetWarnings())
            cerr << "Warning: " << warning.warning_message << endl;
        pe = CVCFScanner::eOK;
    }

    cout << "##fileformat=" << vcf_scanner.GetHeader().GetFileFormat() << endl;

    for (const auto& kv : vcf_scanner.GetHeader().GetMetaInfo()) {
        for (const auto& v : kv.second) {
            cout << "##" << kv.first << '=' << v << endl;
        }
    }

    cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const auto& v : vcf_scanner.GetHeader().GetSampleIDs()) {
        cout << '\t' << v;
    }
    cout << endl;

    while (!vcf_scanner.AtEOF()) {
        if (!ParseDataLine(vcf_scanner, input)) {
            cout << endl;
            cerr << "<-ERR@" << vcf_scanner.GetLineNumber() << ": "
                 << vcf_scanner.GetError().m_ErrorMessage << endl;
        }

        s_ParseToCompletion(vcf_scanner.ClearLine(), vcf_scanner, input);
    }

    return 0;
}
