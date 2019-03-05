/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 */

#include "vcfscanner.hh"

#define NUMBER_OF_MANDATORY_COLUMNS 8

void CVCFScanner::SetNewInputBuffer(const char* buffer, ssize_t buffer_size)
{
    m_Tokenizer.SetNewBuffer(buffer, buffer_size);
}

static const char* const s_HeaderLineColumns[NUMBER_OF_MANDATORY_COLUMNS + 2] =
        {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                "FORMAT", "GENOTYPE"};

CVCFScanner::EParsingEvent CVCFScanner::Feed(
        const char* buffer, ssize_t buffer_size)
{
    m_Tokenizer.SetNewBuffer(buffer, buffer_size);

    switch (m_ParsingState) {
    case ePeekAfterEOL:
        x_ResetDataLine();
        return eOK;

    default:
        if (m_ParsingState < eChrom) {
            return x_ParseHeader();
        }
    }
    return eNeedMoreData;
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParseHeader()
{
    switch (m_ParsingState) {
    case eFileFormatVersion:
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindNewline())) {
            return eNeedMoreData;
        }

        {
            CTempString key, value;

            if (!m_Tokenizer.GetKeyValue(&key, &value) ||
                    key != "##fileformat") {
                return x_HeaderError("VCF file must start with '##fileformat'");
            }

            m_Header.m_FileFormat = value;
        }

    ParseMetaInfoKey:
        m_ParsingState = eMetaInfoKey;
        /* FALL THROUGH */

    case eMetaInfoKey:
        if (!m_Tokenizer.PrepareTokenOrAccumulate(
                    m_Tokenizer.FindNewlineOrTabOrEquals())) {
            return eNeedMoreData;
        }

        switch (m_Tokenizer.GetTokenTerm()) {
        case '\t':
            if (m_Tokenizer.GetToken() != s_HeaderLineColumns[0]) {
                return x_InvalidMetaInfoLineError();
            }
            m_HeaderLineColumnOK = 1;
            goto ParseHeaderColumns;
        case '\n':
            return x_InvalidMetaInfoLineError();
        case EOF:
            return x_UnexpectedEOFInHeader();
        }

        // Found an equals sign - save the key and proceed
        // to parsing the value.
        {
            const CTempString& key = m_Tokenizer.GetToken();
            if (key.length() < 3 || key[0] != '#' || key[1] != '#') {
                return x_InvalidMetaInfoLineError();
            }
            m_CurrentMetaInfoKey = key.substr(2);
        }

        m_ParsingState = eMetaInfoValue;
        /* FALL THROUGH */

    case eMetaInfoValue:
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindNewline())) {
            return eNeedMoreData;
        }

        if (m_Tokenizer.GetTokenTerm() == EOF) {
            return x_UnexpectedEOFInHeader();
        }

        m_Header.m_MetaInfo[m_CurrentMetaInfoKey].push_back(
                m_Tokenizer.GetToken());

        // Go back to parsing the next key.
        goto ParseMetaInfoKey;

    ParseHeaderColumns:
        m_ParsingState = eHeaderLineColumns;
        /* FALL THROUGH */

    case eHeaderLineColumns:
        do {
            if (!m_Tokenizer.PrepareTokenOrAccumulate(
                        m_Tokenizer.FindNewlineOrTab())) {
                return eNeedMoreData;
            }

            if (m_Tokenizer.GetToken() !=
                    s_HeaderLineColumns[m_HeaderLineColumnOK]) {
                return x_InvalidHeaderLineError();
            }

            ++m_HeaderLineColumnOK;

            if (m_Tokenizer.TokenIsLast())
                switch (m_HeaderLineColumnOK) {
                case NUMBER_OF_MANDATORY_COLUMNS + 1:
                    // The FORMAT field is present,
                    // but there are no samples.
                    m_Header.m_GenotypeInfoPresent = true;
                    /* FALL THROUGH */
                case NUMBER_OF_MANDATORY_COLUMNS:
                    return eOK;
                default:
                    return x_InvalidHeaderLineError();
                }

            // The current token ends with a tab.
            // Parse the next header line column.
        } while (m_HeaderLineColumnOK <= NUMBER_OF_MANDATORY_COLUMNS);

        m_Header.m_GenotypeInfoPresent = true;
        m_ParsingState = eSampleIDs;
        /* FALL THROUGH */

    case eSampleIDs:
        do {
            if (!m_Tokenizer.PrepareTokenOrAccumulate(
                        m_Tokenizer.FindNewlineOrTab())) {
                return eNeedMoreData;
            }

            m_Header.m_SampleIDs.push_back(m_Tokenizer.GetToken());
        } while (m_Tokenizer.GetTokenTerm() == '\t');

        if (m_Tokenizer.BufferIsEmpty() && !m_Tokenizer.AtEOF()) {
            m_ParsingState = ePeekAfterEOL;
            return eNeedMoreData;
        }

        x_ResetDataLine();

        return eOK;
    }

    return x_HeaderNotParsedError();
}

#define PARSE_STRING_FIELD(target_state)                                       \
    if (!m_Tokenizer.PrepareTokenOrAccumulate(                                 \
                m_Tokenizer.FindNewlineOrTab())) {                             \
        return eNeedMoreData;                                                  \
    }                                                                          \
    if (m_Tokenizer.TokenIsLast()) {                                           \
        return x_MissingMandatoryFieldError(                                   \
                s_HeaderLineColumns[target_state - eChrom]);                   \
    }

CVCFScanner::EParsingEvent CVCFScanner::ParseLoc()
{
    if (m_ParsingState < eChrom) {
        assert(false && "Header parsing failed");
        return x_HeaderNotParsedError();
    }

    if (m_ParsingState > ePos) {
        assert(false && "Must call ClearLine() before ParseLoc()");
        return eError;
    }

    if (m_ParsingState == eChrom) {
        PARSE_STRING_FIELD(ePos);
        m_Chrom = m_Tokenizer.GetToken();

        m_ParsingState = ePos;
        m_Pos = 0;
        m_NumberLen = 0;
    }

    switch (m_Tokenizer.ParseUnsignedInt(&m_Pos, &m_NumberLen)) {
    case CVCFTokenizer::eEndOfBuffer:
        return eNeedMoreData;
    case CVCFTokenizer::eIntegerOverflow:
        return x_DataLineError("Integer overflow in the POS column");
    default /* case CVCFTokenizer::eEndOfNumber */:
        break;
    }

    if (m_NumberLen == 0) {
        return x_DataLineError("Missing an integer in the POS column");
    }

    if (m_Tokenizer.GetTokenTerm() != '\t') {
        return x_DataLineError("Invalid data line format");
    }

    m_ParsingState = eID;
    m_IDs.clear();
    return eOK;
}

#define PARSE_LIST_FIELD(target_state, container, character_set)               \
    do {                                                                       \
        if (!m_Tokenizer.PrepareTokenOrAccumulate(                             \
                    m_Tokenizer.FindCharFromSet(character_set))) {             \
            return eNeedMoreData;                                              \
        }                                                                      \
        if (m_Tokenizer.TokenIsLast()) {                                       \
            return x_MissingMandatoryFieldError(                               \
                    s_HeaderLineColumns[target_state - eChrom]);               \
        }                                                                      \
        if (!m_Tokenizer.TokenIsDot()) {                                       \
            container.push_back(m_Tokenizer.GetToken());                       \
        }                                                                      \
    } while (m_Tokenizer.GetTokenTerm() != '\t');

#define SKIP_TO_FIELD(field)                                                   \
    if (m_ParsingState < eChrom) {                                             \
        assert(false && "Header parsing failed");                              \
        return x_HeaderNotParsedError();                                       \
    }                                                                          \
    do {                                                                       \
        if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTab())) {          \
            return eNeedMoreData;                                              \
        }                                                                      \
        if (m_Tokenizer.TokenIsLast()) {                                       \
            return x_MissingMandatoryFieldError(                               \
                    s_HeaderLineColumns[m_ParsingState - eChrom + 1]);         \
        }                                                                      \
    } while (++m_ParsingState < field)

CVCFScanner::EParsingEvent CVCFScanner::ParseIDs()
{
    if (m_ParsingState > eID) {
        assert(false && "Must call ClearLine() before ParseIDs()");
        return eError;
    }

    if (m_ParsingState < eID) {
        SKIP_TO_FIELD(eID);
        m_IDs.clear();
    }

    PARSE_LIST_FIELD(eRef, m_IDs, m_Tokenizer.m_NewlineOrTabOrSemicolon);

    m_ParsingState = eRef;
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseAlleles()
{
    if (m_ParsingState > eAlt) {
        assert(false && "Must call ClearLine() before ParseAlleles()");
        return eError;
    }

    if (m_ParsingState < eRef) {
        SKIP_TO_FIELD(eRef);
    }

    if (m_ParsingState == eRef) {
        PARSE_STRING_FIELD(eAlt);
        m_Ref = m_Tokenizer.GetToken();

        m_ParsingState = eAlt;
        m_Alts.clear();
    }

    PARSE_LIST_FIELD(eQuality, m_Alts, m_Tokenizer.m_NewlineOrTabOrComma);

    m_AllelesParsed = true;

    m_ParsingState = eQuality;
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseQuality()
{
    if (m_ParsingState > eQuality) {
        assert(false && "Must call ClearLine() before ParseQuality()");
        return eError;
    }

    if (m_ParsingState < eQuality) {
        SKIP_TO_FIELD(eQuality);
    }

    PARSE_STRING_FIELD(eFilter);
    if (!m_Tokenizer.TokenIsDot())
        m_Quality = m_Tokenizer.GetToken();
    else
        m_Quality.clear();

    m_ParsingState = eFilter;
    m_Filters.clear();

    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseFilters()
{
    if (m_ParsingState > eFilter) {
        assert(false && "Must call ClearLine() before ParseFilters()");
        return eError;
    }

    if (m_ParsingState < eFilter) {
        SKIP_TO_FIELD(eFilter);
        m_Filters.clear();
    }

    PARSE_LIST_FIELD(
            eInfoField, m_Filters, m_Tokenizer.m_NewlineOrTabOrSemicolon);

    m_ParsingState = eInfoField;
    m_Info.clear();
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseInfo()
{
    if (m_ParsingState > eInfoField) {
        assert(false && "Must call ClearLine() before ParseInfo()");
        return eError;
    }

    if (m_ParsingState < eInfoField) {
        SKIP_TO_FIELD(eInfoField);
        m_Info.clear();
    }

    do {
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindCharFromSet(
                    m_Tokenizer.m_NewlineOrTabOrSemicolon))) {
            return eNeedMoreData;
        }
        if (m_Tokenizer.TokenIsLast()) {
            m_ParsingState = eEndOfDataLine;
            return eOK;
        }
        if (!m_Tokenizer.TokenIsDot()) {
            m_Info.push_back(m_Tokenizer.GetToken());
        }
    } while (m_Tokenizer.GetTokenTerm() != '\t');

    m_ParsingState = eGenotypeFormat;
    m_GenotypeKeyPositions.Clear();
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseGenotypeFormat()
{
    if (m_ParsingState > eGenotypeFormat) {
        assert(false && "Must call ClearLine() before ParseGenotypeFormat()");
        return eError;
    }

    if (m_ParsingState < eGenotypeFormat) {
        SKIP_TO_FIELD(eGenotypeFormat);
        m_GenotypeKeyPositions.Clear();
    }

    do {
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindCharFromSet(
                    m_Tokenizer.m_NewlineOrTabOrColon))) {
            return eNeedMoreData;
        }
        if (m_Tokenizer.TokenIsLast()) {
            m_ParsingState = eEndOfDataLine;
            if (m_Header.m_SampleIDs.empty()) {
                return eOK;
            }
            return x_DataLineError("No genotype information present");
        }
        string key = m_Tokenizer.GetToken();
        auto key_iter = m_FormatKeys.insert(key).first;
        if (key == "GT") {
            if (m_GenotypeKeyPositions.number_of_positions != 0) {
                // TODO Generate a warning: GT must be the first key.
            }
            m_GenotypeKeyPositions.GT =
                    ++m_GenotypeKeyPositions.number_of_positions;
        } else
            m_GenotypeKeyPositions.other_keys[key_iter->c_str()] =
                    ++m_GenotypeKeyPositions.number_of_positions;
    } while (m_Tokenizer.GetTokenTerm() != '\t');

    m_ParsingState = eGenotypes;
    x_ClearGenotypeValues();
    return eOK;
}

bool CVCFScanner::CaptureGT(vector<int>* alleles)
{
    unsigned gt_index = m_GenotypeKeyPositions.GT;
    if (gt_index == 0) {
        return false;
    }
    --gt_index;
    auto gt_value = x_AllocGenotypeValue(gt_index);
    gt_value->data_type = eGT;
    gt_value->int_vector = alleles;
    alleles->clear();
    return true;
}

const char* CVCFScanner::x_ParseGT(vector<int>* int_vector)
{
    int_vector->clear();

    const string& token = m_Tokenizer.GetToken();

    size_t len = token.length();

    if (len == 0) {
        return "Empty GT value";
    }

    const char* ptr = token.data();
    unsigned digit, allele;

    for (;; ++ptr, --len) {
        if (*ptr == '.') {
            int_vector->push_back(-1);
            ++ptr;
            --len;
        } else {
            if ((allele = (unsigned) *ptr - '0') > 9) {
                break;
            }

            while (--len > 0 && (digit = (unsigned) *++ptr - '0') <= 9) {
                if (allele > (UINT_MAX / 10) ||
                        (allele == (UINT_MAX / 10) && digit > UINT_MAX % 10)) {
                    // TODO ERR+="in genotype info for the sample #" +
                    //           str(m_CurrentGenotypeValueIndex)
                    return "Integer overflow in allele index";
                }

                allele = allele * 10 + digit;
            }

            int_vector->push_back((int) allele);

            if (m_AllelesParsed && allele > m_Alts.size()) {
                return "Allele index exceeds the number of alleles";
            }
        }
        if (len == 0) {
            return nullptr;
        }
        switch (*ptr) {
        case '/':
            m_PhasedGT = false;
            continue;
        case '|':
            m_PhasedGT = true;
            continue;
        }
        break;
    }
    return "Invalid character in GT value";
}

CVCFScanner::EParsingEvent CVCFScanner::ParseGenotype()
{
    if (m_ParsingState != eGenotypes) {
        assert(false &&
                "Must call ParseGenotypeFormat() before ParseGenotype()");
        return eError;
    }

    if (m_CurrentGenotypeValueIndex >=
            m_GenotypeKeyPositions.number_of_positions) {
        return x_DataLineError(
                "The number of genotype fields exceeds the number of samples");
    }

    SGenotypeValue* value =
            m_GenotypeValues.data() + m_CurrentGenotypeValueIndex;

    const char* error_message;

    do {
        if (value->flag == nullptr) {
            if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTabOrColon())) {
                return eNeedMoreData;
            }
            if (m_Tokenizer.TokenIsLast()) {
                m_ParsingState = eEndOfDataLine;
                return eOK;
            }
        } else {
            if (!m_Tokenizer.PrepareTokenOrAccumulate(
                        m_Tokenizer.FindCharFromSet(
                                m_Tokenizer.m_NewlineOrTabOrColon))) {
                return eNeedMoreData;
            }

            if (m_Tokenizer.TokenIsLast()) {
                m_ParsingState = eEndOfDataLine;
            }

            switch (value->data_type) {
            // TODO case eInteger:
            // TODO case eFloat:
            // TODO case eCharacter:
            // TODO case eString:
            case eGT:
                error_message = x_ParseGT(value->int_vector);
                if (error_message != nullptr) {
                    return x_DataLineError(error_message);
                }
                break;
            default /* eFlag - impossible type for genotype info */:
                break;
            }

            if (m_Tokenizer.TokenIsLast()) {
                return eOK;
            }
        }
        if (m_Tokenizer.GetTokenTerm() == '\t') {
            return eOK;
        }

        ++value;
    } while (++m_CurrentGenotypeValueIndex <
            m_GenotypeKeyPositions.number_of_positions);

    return x_DataLineError("Too many genotype info fields");
}

CVCFScanner::EParsingEvent CVCFScanner::ClearLine()
{
    if (m_Tokenizer.AtEOF()) {
        return eOK;
    }

    if (m_ParsingState != ePeekAfterEOL) {
        if (m_ParsingState != eEndOfDataLine &&
                !m_Tokenizer.SkipToken(m_Tokenizer.FindNewline())) {
            return eNeedMoreData;
        }

        if (m_Tokenizer.BufferIsEmpty()) {
            m_ParsingState = ePeekAfterEOL;
            return eNeedMoreData;
        }
    }

    x_ResetDataLine();

    return eOK;
}

static char buffer[1];
static size_t buffer_size;

static void s_ReadBuffer(FILE* input)
{
    buffer_size = fread(buffer, 1, sizeof(buffer), input);
}

static void s_ReadAndSetNewInputBuffer(CVCFScanner& vcf_scanner, FILE* input)
{
    vcf_scanner.SetNewInputBuffer(
            buffer, fread(buffer, 1, sizeof(buffer), input));
}

#define RETRY_UNTIL_OK_OR_ERROR(method)                                        \
    while ((pe = vcf_scanner.method()) == CVCFScanner::eNeedMoreData) {        \
        s_ReadAndSetNewInputBuffer(vcf_scanner, input);                        \
    }                                                                          \
    if (pe == CVCFScanner::eOKWithWarnings) {                                  \
        for (const auto& warning : vcf_scanner.GetWarnings())                  \
            cerr << "Warning: " << warning.warning_message << endl;            \
        pe = CVCFScanner::eOK;                                                 \
    }

#define RETRY_UNTIL_OK_OR_SKIP_LINE(method)                                    \
    RETRY_UNTIL_OK_OR_ERROR(method);                                           \
    if (pe == CVCFScanner::eError) {                                           \
        cerr << "<-ERR@" << vcf_scanner.GetLineNumber() << ": "                \
             << vcf_scanner.GetError().m_ErrorMessage << endl;                 \
        return;                                                                \
    }

static void ParseDataLine(CVCFScanner& vcf_scanner, FILE* input)
{
    CVCFScanner::EParsingEvent pe;

    cout << vcf_scanner.GetLineNumber() << ':';

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseLoc);
    cout << vcf_scanner.GetChrom() << '\t' << vcf_scanner.GetPos();

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseIDs);
    const char* sep;
    if (!vcf_scanner.GetIDs().empty()) {
        sep = "\t";
        for (const auto& id : vcf_scanner.GetIDs()) {
            cout << sep << id;
            sep = ",";
        }
        cout << '\t';
    } else
        cout << "\t.\t";

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseAlleles);
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

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseQuality);
    string quality = vcf_scanner.GetQuality();
    if (!quality.empty())
        cout << quality;
    else
        cout << ".";

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseFilters);
    if (!vcf_scanner.GetFilters().empty()) {
        sep = "\t";
        for (const auto& filter : vcf_scanner.GetFilters()) {
            cout << sep << filter;
            sep = ";";
        }
    } else
        cout << "\t.";

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseInfo);
    if (!vcf_scanner.GetInfo().empty()) {
        sep = "\t";
        for (const auto& info_item : vcf_scanner.GetInfo()) {
            cout << sep << info_item;
            sep = ";";
        }
        cout << '\t';
    } else
        cout << "\t.\t";

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseGenotypeFormat);
    cout << "GT";

    vector<int> alleles;

    if (!vcf_scanner.CaptureGT(&alleles)) {
        cout << "\tERR: no GT key" << endl;
        return;
    }

    while (vcf_scanner.GenotypeAvailable()) {
        RETRY_UNTIL_OK_OR_SKIP_LINE(ParseGenotype);
        sep = "\t";
        for (auto allele : alleles) {
            cout << sep << allele;
            sep = vcf_scanner.IsPhasedGT() ? "|" : "/";
        }
    }

    cout << endl;
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
        ParseDataLine(vcf_scanner, input);

        while (vcf_scanner.ClearLine() == CVCFScanner::eNeedMoreData)
            s_ReadAndSetNewInputBuffer(vcf_scanner, input);
    }

    return 0;
}
