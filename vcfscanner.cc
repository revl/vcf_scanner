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


void CVCFScanner::SetNewInputBuffer(const char* buffer, ssize_t buffer_size)
{
    m_Buffer = buffer;
    m_BufferSize = buffer_size;

    m_Tokenizer.SetNewBuffer(buffer, buffer_size);
}

const char* const CVCFScanner::m_HeaderLineColumns[] = {"#CHROM", "POS", "ID",
        "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GENOTYPE"};

CVCFScanner::EParsingEvent CVCFScanner::ParseHeader()
{
    switch (m_HeaderParsingState) {
    case eFileFormatVersion:
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindNewline()))
            return eNeedMoreData;

        {
            CTempString key, value;

            if (!m_Tokenizer.GetKeyValue(&key, &value) ||
                    key != "##fileformat") {
                return x_HeaderError("VCF file must start with '##fileformat'");
            }

            m_Header.m_FileFormat = value;
        }

    ParseMetaInfoKey:
        m_HeaderParsingState = eMetaInfoKey;
        /* FALL THROUGH */

    case eMetaInfoKey:
        if (!m_Tokenizer.PrepareTokenOrAccumulate(
                    m_Tokenizer.FindNewlineOrTabOrEquals()))
            return eNeedMoreData;

        switch (m_Tokenizer.GetTokenTerm()) {
        case '\t':
            if (m_Tokenizer.GetToken() != m_HeaderLineColumns[0])
                goto InvalidMetaInfoLine;
            m_HeaderLineColumnOK = 1;
            goto ParseHeaderColumns;
        case '\n':
            goto InvalidMetaInfoLine;
        case EOF:
            goto UnexpectedEOF;
        }

        // Found an equals sign - save the key and proceed
        // to parsing the value.
        {
            const CTempString& key = m_Tokenizer.GetToken();
            if (key.length() < 3 || key[0] != '#' || key[1] != '#')
                goto InvalidMetaInfoLine;
            m_CurrentMetaInfoKey = key.substr(2);
        }

        m_HeaderParsingState = eMetaInfoValue;
        /* FALL THROUGH */

    case eMetaInfoValue:
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindNewline()))
            return eNeedMoreData;

        if (m_Tokenizer.GetTokenTerm() == EOF)
            goto UnexpectedEOF;

        m_Header.m_MetaInfo[m_CurrentMetaInfoKey].push_back(
                m_Tokenizer.GetToken());

        // Go back to parsing the next key.
        goto ParseMetaInfoKey;

    ParseHeaderColumns:
        m_HeaderParsingState = eHeaderLineColumns;
        /* FALL THROUGH */

    case eHeaderLineColumns:
        do {
            if (!m_Tokenizer.PrepareTokenOrAccumulate(
                        m_Tokenizer.FindNewlineOrTab()))
                return eNeedMoreData;

            if (m_Tokenizer.GetToken() !=
                    m_HeaderLineColumns[m_HeaderLineColumnOK])
                goto InvalidHeaderLine;

            ++m_HeaderLineColumnOK;

            switch (m_Tokenizer.GetTokenTerm()) {
            case '\n':
            case EOF:
                switch (m_HeaderLineColumnOK) {
                case eNumberOfMandatoryColumns + 1:
                    m_Header.m_GenotypeInfoPresent = true;
                    /* FALL THROUGH */
                case eNumberOfMandatoryColumns:
                    return eOK;
                default:
                    goto InvalidHeaderLine;
                }
            }

            // The current token ends with a tab.
            // Parse the next header line column.
        } while (m_HeaderLineColumnOK <= eNumberOfMandatoryColumns);

        m_Header.m_GenotypeInfoPresent = true;
        m_HeaderParsingState = eSampleIDs;
        /* FALL THROUGH */

    case eSampleIDs:
        do {
            if (!m_Tokenizer.PrepareTokenOrAccumulate(
                        m_Tokenizer.FindNewlineOrTab()))
                return eNeedMoreData;

            m_Header.m_SampleIDs.push_back(m_Tokenizer.GetToken());
        } while (m_Tokenizer.GetTokenTerm() == '\t');

        return eOK;
    }

UnexpectedEOF:
    return x_HeaderError(
            "Unexpected end of file while parsing VCF file header");

InvalidMetaInfoLine:
    return x_HeaderError("Malformed meta-information line");

InvalidHeaderLine:
    return x_HeaderError("Malformed VCF header line");
}

bool CVCFScanner::Rewind()
{
    if (m_Tokenizer.AtEOF())
        return true;

    if (m_Tokenizer.BufferIsEmpty())
        return false;

    if (m_DataLineParsingState != eEndOfDataLine &&
            !m_Tokenizer.SkipToken(m_Tokenizer.FindNewline()))
        return false;

    m_DataLineParsingState = eChrom;
    m_AllelesParsed = false;
    return true;
}

#define PARSE_STRING_FIELD(target_state)                                       \
    if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindNewlineOrTab())) \
        return eNeedMoreData;                                                  \
    if (m_Tokenizer.TokenIsLast()) {                                           \
        return x_MissingMandatoryFieldError(                                   \
                m_HeaderLineColumns[target_state]);                            \
    }

CVCFScanner::EParsingEvent CVCFScanner::ParseLoc()
{
    if (m_DataLineParsingState > ePos) {
        assert(false && "Must call Rewind() before ParseLoc()");
        return eError;
    }

    if (m_DataLineParsingState == eChrom) {
        PARSE_STRING_FIELD(eChrom);
        m_Chrom = m_Tokenizer.GetToken();

        m_DataLineParsingState = ePos;
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

    if (m_NumberLen == 0)
        return x_DataLineError("Missing an integer in the POS column");

    if (m_Tokenizer.GetTokenTerm() != '\t')
        return x_DataLineError("Invalid data line format");

    m_DataLineParsingState = eID;
    m_IDs.clear();
    return eOK;
}

#define PARSE_LIST_FIELD(target_state, container, character_set)               \
    do {                                                                       \
        if (!m_Tokenizer.PrepareTokenOrAccumulate(                             \
                    m_Tokenizer.FindCharFromSet(character_set)))               \
            return eNeedMoreData;                                              \
        if (m_Tokenizer.TokenIsLast())                                         \
            return x_MissingMandatoryFieldError(                               \
                    m_HeaderLineColumns[target_state]);                        \
        if (!m_Tokenizer.TokenIsDot())                                         \
            container.push_back(m_Tokenizer.GetToken());                       \
    } while (m_Tokenizer.GetTokenTerm() != '\t');

#define SKIP_TO_FIELD(field)                                                   \
    do {                                                                       \
        if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTab()))            \
            return eNeedMoreData;                                              \
        if (m_Tokenizer.TokenIsLast())                                         \
            return x_MissingMandatoryFieldError(                               \
                    m_HeaderLineColumns[m_DataLineParsingState + 1]);          \
    } while (++m_DataLineParsingState < field)

CVCFScanner::EParsingEvent CVCFScanner::ParseIDs()
{
    if (m_DataLineParsingState > eID) {
        assert(false && "Must call Rewind() before ParseIDs()");
        return eError;
    }

    if (m_DataLineParsingState < eID) {
        SKIP_TO_FIELD(eID);
        m_IDs.clear();
    }

    PARSE_LIST_FIELD(eRef, m_IDs, m_Tokenizer.m_NewlineOrTabOrSemicolon);

    m_DataLineParsingState = eRef;
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseAlleles()
{
    if (m_DataLineParsingState > eAlt) {
        assert(false && "Must call Rewind() before ParseAlleles()");
        return eError;
    }

    if (m_DataLineParsingState < eRef) {
        SKIP_TO_FIELD(eRef);
    }

    if (m_DataLineParsingState == eRef) {
        PARSE_STRING_FIELD(eRef);
        m_Ref = m_Tokenizer.GetToken();

        m_DataLineParsingState = eAlt;
        m_Alts.clear();
    }

    PARSE_LIST_FIELD(eQuality, m_Alts, m_Tokenizer.m_NewlineOrTabOrComma);

    m_AllelesParsed = true;

    m_DataLineParsingState = eQuality;
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseQuality()
{
    if (m_DataLineParsingState > eQuality) {
        assert(false && "Must call Rewind() before ParseQuality()");
        return eError;
    }

    if (m_DataLineParsingState < eQuality) {
        SKIP_TO_FIELD(eQuality);
    }

    PARSE_STRING_FIELD(eQuality);
    if (!m_Tokenizer.TokenIsDot())
        m_Quality = m_Tokenizer.GetToken();
    else
        m_Quality.clear();

    m_DataLineParsingState = eFilter;
    m_Filters.clear();

    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseFilters()
{
    if (m_DataLineParsingState > eFilter) {
        assert(false && "Must call Rewind() before ParseFilters()");
        return eError;
    }

    if (m_DataLineParsingState < eFilter) {
        SKIP_TO_FIELD(eFilter);
        m_Filters.clear();
    }

    PARSE_LIST_FIELD(
            eInfoField, m_Filters, m_Tokenizer.m_NewlineOrTabOrSemicolon);

    m_DataLineParsingState = eInfoField;
    m_Info.clear();
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseInfo()
{
    if (m_DataLineParsingState > eInfoField) {
        assert(false && "Must call Rewind() before ParseInfo()");
        return eError;
    }

    if (m_DataLineParsingState < eInfoField) {
        SKIP_TO_FIELD(eInfoField);
        m_Info.clear();
    }

    do {
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindCharFromSet(
                    m_Tokenizer.m_NewlineOrTabOrSemicolon)))
            return eNeedMoreData;
        if (m_Tokenizer.TokenIsLast()) {
            m_DataLineParsingState = eEndOfDataLine;
            return eOK;
        }
        if (!m_Tokenizer.TokenIsDot())
            m_Info.push_back(m_Tokenizer.GetToken());
    } while (m_Tokenizer.GetTokenTerm() != '\t');

    m_DataLineParsingState = eGenotypeFormat;
    m_GenotypeKeyPositions.Clear();
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseGenotypeFormat()
{
    if (m_DataLineParsingState > eGenotypeFormat) {
        assert(false && "Must call Rewind() before ParseGenotypeFormat()");
        return eError;
    }

    if (m_DataLineParsingState < eGenotypeFormat) {
        SKIP_TO_FIELD(eGenotypeFormat);
        m_GenotypeKeyPositions.Clear();
    }

    do {
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindCharFromSet(
                    m_Tokenizer.m_NewlineOrTabOrColon)))
            return eNeedMoreData;
        if (m_Tokenizer.TokenIsLast()) {
            m_DataLineParsingState = eEndOfDataLine;
            if (m_Header.m_SampleIDs.empty())
                return eOK;
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

    m_DataLineParsingState = eGenotypes;
    x_ClearGenotypeValues();
    return eOK;
}

bool CVCFScanner::CaptureGT(vector<int>* alleles)
{
    unsigned gt_index = m_GenotypeKeyPositions.GT;
    if (gt_index == 0)
        return false;
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

    if (len == 0)
        return "Empty GT value";

    const char* ptr = token.data();
    unsigned digit, allele;

    for (;; ++ptr, --len) {
        if (*ptr == '.') {
            int_vector->push_back(-1);
            ++ptr;
            --len;
        } else {
            if ((allele = (unsigned) *ptr - '0') > 9)
                break;

            while (--len > 0 && (digit = (unsigned) *++ptr - '0') <= 9) {
                if (allele > (UINT_MAX / 10) ||
                        (allele == (UINT_MAX / 10) && digit > UINT_MAX % 10))
                    return "Integer overflow in allele index";

                allele = allele * 10 + digit;
            }

            int_vector->push_back((int) allele);

            if (m_AllelesParsed && allele > m_Alts.size())
                return "Allele index exceeds the number of alleles";
        }
        if (len == 0)
            return nullptr;
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
    if (m_DataLineParsingState != eGenotypes) {
        assert(false &&
                "Must call ParseGenotypeFormat() before ParseGenotype()");
        return eError;
    }

    if (m_CurrentGenotypeValueIndex >=
            m_GenotypeKeyPositions.number_of_positions)
        return x_DataLineError(
                "The number of genotype fields exceeds the number of samples");

    SGenotypeValue* value =
            m_GenotypeValues.data() + m_CurrentGenotypeValueIndex;

    const char* error_message;

    do {
        if (value->flag == nullptr) {
            if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTabOrColon()))
                return eNeedMoreData;
            if (m_Tokenizer.TokenIsLast()) {
                m_DataLineParsingState = eEndOfDataLine;
                return eOK;
            }
        } else {
            if (!m_Tokenizer.PrepareTokenOrAccumulate(
                        m_Tokenizer.FindCharFromSet(
                                m_Tokenizer.m_NewlineOrTabOrColon)))
                return eNeedMoreData;

            if (m_Tokenizer.TokenIsLast())
                m_DataLineParsingState = eEndOfDataLine;

            switch (value->data_type) {
            // TODO case eInteger:
            // TODO case eFloat:
            // TODO case eCharacter:
            // TODO case eString:
            case eGT:
                error_message = x_ParseGT(value->int_vector);
                if (error_message != nullptr)
                    return x_DataLineError(error_message);
                break;
            default /* eFlag - impossible type for genotype info */:
                break;
            }

            if (m_Tokenizer.TokenIsLast())
                return eOK;
        }
        if (m_Tokenizer.GetTokenTerm() == '\t')
            return eOK;

        ++value;
    } while (++m_CurrentGenotypeValueIndex <
            m_GenotypeKeyPositions.number_of_positions);

    return x_DataLineError("Too many genotype info fields");

    /* do {
        if (!m_Tokenizer.PrepareTokenOrAccumulate(
                    m_Tokenizer.FindCharFromSet(
                            m_Tokenizer.m_NewlineTabColonSlashBar)))
            return eNeedMoreData;
        if (m_Tokenizer.TokenIsLast()) {
            m_DataLineParsingState = eEndOfDataLine;
        }

        unsigned allele_index;
        if (!m_Tokenizer.GetTokenAsUInt(&allele_index)) {
            if (m_Tokenizer.GetToken() != ".")
                return x_DataLineError(
                        "GT values must be either numbers or dots");
            value.int_vector->push_back(-1);
        } else {
            if (m_AllelesParsed && m_UIntAcc > m_Alts.size())
                return x_DataLineError(
                        "Allele index exceeds the number of alleles");
            value.int_vector->push_back((int) m_UIntAcc);
        }
    } while (m_Tokenizer.GetTokenTerm() != '\t');

    switch (m_Tokenizer.ParseUnsignedInt(&m_UIntAcc, &m_NumberLen)) {
    case CVCFTokenizer::eEndOfBuffer:
        return eNeedMoreData;
    case CVCFTokenizer::eIntegerOverflow:
        // TODO ERR+=in genotype info for sample m_CurrentGenotypeValueIndex
        return x_DataLineError("Integer overflow in allele index");
    default / * case CVCFTokenizer::eEndOfNumber * /:
        break;
    }
    if (m_NumberLen == 0) {
        if (m_Tokenizer.GetTokenTerm() != '.')
            value.int_vector->push_back(-1);
    } else {
        if (m_AllelesParsed && m_UIntAcc > m_Alts.size())
            return x_DataLineError(
                    "Allele index exceeds the number of alleles");
        value.int_vector->push_back((int) m_UIntAcc);
    } */

    /* if (!m_Tokenizer.PrepareTokenOrAccumulate(
                m_Tokenizer.FindCharFromSet(m_Tokenizer.m_NewlineOrTab)))
        return eNeedMoreData;

    switch (m_Tokenizer.GetTokenTerm()) {
    case EOF:
    case '\n':
        m_DataLineParsingState = eEndOfDataLine;
    }

    return eOK; */

    /* EOLHandler:

        m_DataLineParsingState = eEndOfDataLine; */
}

static const char vcf_file[] = R"(##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S-1	S-2	S-3
1	100000	rs123;rs333	C	G	10	.	.	GT	0|1	1/.	1/0
2	200000	.	C	G,T	.	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT	0|0	0|1	1|2)";

static const char* vcf_eof_ptr = vcf_file + sizeof(vcf_file) - 1;
static const char* current_ptr = vcf_file;

static ssize_t s_ReadVCFFile(void* target_buffer, size_t buffer_size)
{
    size_t remaining_size = vcf_eof_ptr - current_ptr;
    if (buffer_size > remaining_size)
        buffer_size = remaining_size;
    memcpy(target_buffer, current_ptr, buffer_size);
    current_ptr += buffer_size;
    return buffer_size;
}

static char buffer[1];

#define SET_NEW_BUFFER()                                                       \
    vcf_scanner.SetNewInputBuffer(buffer, s_ReadVCFFile(buffer, sizeof(buffer)))

#define RETRY_UNTIL_OK_OR_ERROR(method)                                        \
    while ((pe = vcf_scanner.method()) == CVCFScanner::eNeedMoreData)          \
        SET_NEW_BUFFER();                                                      \
    if (pe == CVCFScanner::eOKWithWarnings) {                                  \
        for (const auto& warning : vcf_scanner.GetWarnings())                  \
            cerr << "Warning: " << warning.warning_message << endl;            \
        pe = CVCFScanner::eOK;                                                 \
    }

#define RETRY_UNTIL_OK_OR_SKIP_LINE(method)                                    \
    RETRY_UNTIL_OK_OR_ERROR(method);                                           \
    if (pe == CVCFScanner::eError) {                                           \
        cerr << "file.vcf:" << vcf_scanner.GetLineNumber() << ": "             \
             << vcf_scanner.GetError().m_ErrorMessage << endl                  \
             << endl;                                                          \
        return;                                                                \
    }

static void ParseDataLine(CVCFScanner& vcf_scanner)
{
    CVCFScanner::EParsingEvent pe;

    cout << '#' << vcf_scanner.GetLineNumber() << endl;

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseLoc);
    cout << "CHROM:[" << vcf_scanner.GetChrom() << ']' << endl
         << "POS:[" << vcf_scanner.GetPos() << ']' << endl;

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseIDs);
    for (const auto& id : vcf_scanner.GetIDs())
        cout << "ID:[" << id << ']' << endl;

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseAlleles);
    cout << "REF:[" << vcf_scanner.GetRef() << ']' << endl;
    for (const auto& alt : vcf_scanner.GetAlts())
        cout << "ALT:[" << alt << ']' << endl;

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseQuality);
    string quality = vcf_scanner.GetQuality();
    if (!quality.empty())
        cout << "QUAL:[" << quality << ']' << endl;

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseFilters);
    for (const auto& filter : vcf_scanner.GetFilters())
        cout << "FILTER:[" << filter << ']' << endl;

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseInfo);
    for (const auto& info_item : vcf_scanner.GetInfo())
        cout << "INFO:[" << info_item << ']' << endl;

    RETRY_UNTIL_OK_OR_SKIP_LINE(ParseGenotypeFormat);

    vector<int> alleles;

    if (!vcf_scanner.CaptureGT(&alleles)) {
        cout << "ERR: no GT key" << endl;
        return;
    }

    while (vcf_scanner.GenotypeAvailable()) {
        RETRY_UNTIL_OK_OR_SKIP_LINE(ParseGenotype);
        cout << "---" << endl;
        for (auto allele : alleles)
            cout << "GT:[" << allele << ']' << endl;
    }

    cout << endl;
}

int main()
{
    CVCFScanner vcf_scanner;

    SET_NEW_BUFFER();

    CVCFScanner::EParsingEvent pe;

    // Read the header
    RETRY_UNTIL_OK_OR_ERROR(ParseHeader);

    if (pe != CVCFScanner::eOK) {
        cerr << vcf_scanner.GetError().m_ErrorMessage << endl;
        return 1;
    }

    for (const auto& kv : vcf_scanner.GetHeader().GetMetaInfo()) {
        cout << '[' << kv.first << ']' << endl;
        for (const auto& v : kv.second) {
            cout << v << endl;
        }
    }

    cout << endl << "[Sample IDs]" << endl;
    for (const auto& v : vcf_scanner.GetHeader().GetSampleIDs()) {
        cout << v << endl;
    }

    cout << endl << "[Data Lines]" << endl;
    for (;;) {
        while (!vcf_scanner.Rewind())
            SET_NEW_BUFFER();
        if (vcf_scanner.AtEOF())
            break;

        ParseDataLine(vcf_scanner);
    }

    return 0;
}
