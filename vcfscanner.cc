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

static const char* const s_HeaderLineColumns[NUMBER_OF_MANDATORY_COLUMNS + 2] =
        {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                "GENOTYPE"};

#define PARSE_STRING(target_state)                                             \
    if (!m_Tokenizer.PrepareTokenOrAccumulate(                                 \
                m_Tokenizer.FindNewlineOrTab())) {                             \
        return eNeedMoreData;                                                  \
    }                                                                          \
    if (m_Tokenizer.TokenIsLast()) {                                           \
        return x_MissingMandatoryFieldError(                                   \
                s_HeaderLineColumns[target_state - eChrom]);                   \
    }                                                                          \
    m_ParsingState = target_state;

#define PARSE_STRING_LIST(target_state, container, character_set)              \
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
            if (m_NextListIndex < container.size()) {                          \
                container.at(m_NextListIndex) = m_Tokenizer.GetToken();        \
            } else {                                                           \
                container.push_back(m_Tokenizer.GetToken());                   \
            }                                                                  \
            ++m_NextListIndex;                                                 \
        }                                                                      \
    } while (m_Tokenizer.GetTokenTerm() != '\t');                              \
    container.resize(m_NextListIndex);                                         \
    m_ParsingState = target_state;

#define SKIP_TO_STATE(target_state)                                            \
    if (m_ParsingState < eChrom) {                                             \
        assert(false && "VCF header must be parsed first");                    \
        return eError;                                                         \
    }                                                                          \
    if (m_ParsingState > target_state) {                                       \
        assert(false && "ClearLine must call be called first");                \
        return eError;                                                         \
    }                                                                          \
    while (m_ParsingState < target_state) {                                    \
        if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTab())) {          \
            m_FieldsToSkip = target_state - m_ParsingState;                    \
            m_ParsingState = target_state;                                     \
            return eNeedMoreData;                                              \
        }                                                                      \
        if (m_Tokenizer.TokenIsLast()) {                                       \
            return x_MissingMandatoryFieldError(                               \
                    s_HeaderLineColumns[m_ParsingState - eChrom + 1]);         \
        }                                                                      \
        ++m_ParsingState;                                                      \
    }

CVCFScanner::EParsingEvent CVCFScanner::x_SkipToState(
        CVCFScanner::EParsingState target_state)
{
    // LCOV_EXCL_START
    if (m_ParsingState < eChrom) {
        assert(false && "VCF header must be parsed first");
        return eError;
    }
    if (m_ParsingState > target_state) {
        assert(false && "ClearLine must call be called first");
        return eError;
    }
    // LCOV_EXCL_STOP

    while (m_ParsingState < target_state) {
        if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTab())) {
            m_FieldsToSkip = target_state - m_ParsingState;
            m_ParsingState = target_state;
            return eNeedMoreData;
        }
        if (m_Tokenizer.TokenIsLast()) {
            return x_MissingMandatoryFieldError(
                    s_HeaderLineColumns[m_ParsingState - eChrom + 1]);
        }
        ++m_ParsingState;
    }
    return eOK;
}

#define PARSE_CHROM()                                                          \
    PARSE_STRING(ePos);                                                        \
    m_Chrom = m_Tokenizer.GetToken();

#define PARSE_REF()                                                            \
    PARSE_STRING(eAlt);                                                        \
    m_Ref = m_Tokenizer.GetToken();

CVCFScanner::EParsingEvent CVCFScanner::Feed(
        const char* buffer, ssize_t buffer_size)
{
    m_Tokenizer.SetNewBuffer(buffer, buffer_size);

    if (m_ParsingState == eGenotypes) {
        return x_ParseGenotype();
    }

    if (m_ParsingState <= ePos) {
        if (m_ParsingState < eChrom) {
            return x_ParseHeader();
        }
        if (m_ParsingState == eChrom) {
            PARSE_CHROM();
        }
        return x_ParsePos();
    } else {
        for (; m_FieldsToSkip > 0; --m_FieldsToSkip) {
            if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTab())) {
                return eNeedMoreData;
            }
            if (m_Tokenizer.TokenIsLast()) {
                unsigned missing_field =
                        m_ParsingState - eChrom + 1 - m_FieldsToSkip;
                m_FieldsToSkip = 0;
                return x_MissingMandatoryFieldError(
                        s_HeaderLineColumns[missing_field]);
            }
        }
    }

    switch (m_ParsingState) {
    case eID:
        return x_ParseIDs();
    case eRef:
        PARSE_REF();
        /* FALL THROUGH */
    case eAlt:
        return x_ParseAlts();
    case eQuality:
        return x_ParseQuality();
    case eFilter:
        return x_ParseFilters();
    case eInfoField:
        return x_ParseInfo();
    case eGenotypeFormat:
        return x_ParseGenotypeFormat();
    case eClearLine:
        if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewline())) {
            return eNeedMoreData;
        }
        if (m_Tokenizer.BufferIsEmpty() && !m_Tokenizer.AtEOF()) {
            m_ParsingState = ePeekAfterEOL;
            return eNeedMoreData;
        }
        /* FALL THROUGH */
    case ePeekAfterEOL:
        x_ResetDataLine();
        return eOK;
    }

    return eError; // LCOV_EXCL_LINE
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

        if (m_Tokenizer.TokenIsLast()) {
            return x_InvalidMetaInfoLineError();
        }

        if (m_Tokenizer.GetTokenTerm() == '\t') {
            if (m_Tokenizer.GetToken().substr(1) != s_HeaderLineColumns[0]) {
                return x_InvalidMetaInfoLineError();
            }
            m_HeaderLineColumnOK = 1;
            goto ParseHeaderLine;
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
            return x_HeaderError(
                    "Unexpected end of file while parsing VCF file header");
        }

        m_Header.m_MetaInfo[m_CurrentMetaInfoKey].push_back(
                m_Tokenizer.GetToken());

        // Go back to parsing the next key.
        goto ParseMetaInfoKey;

    ParseHeaderLine:
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

            if (m_Tokenizer.TokenIsLast()) {
                if (m_HeaderLineColumnOK < NUMBER_OF_MANDATORY_COLUMNS) {
                    return x_InvalidHeaderLineError();
                }
                if (m_HeaderLineColumnOK > NUMBER_OF_MANDATORY_COLUMNS) {
                    // The FORMAT field is present,
                    // but there are no samples.
                    m_Header.m_GenotypeInfoPresent = true;
                }
                goto EndOfHeaderLine;
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
    }

EndOfHeaderLine:
    if (m_Tokenizer.BufferIsEmpty() && !m_Tokenizer.AtEOF()) {
        m_ParsingState = ePeekAfterEOL;
        return eNeedMoreData;
    }

    x_ResetDataLine();

    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseLoc()
{
    // LCOV_EXCL_START
    if (m_ParsingState != eChrom) {
        if (m_ParsingState < eChrom) {
            assert(false && "VCF header must be parsed first");
            return eError;
        }

        assert(false && "Must call ClearLine() before ParseLoc()");
        return eError;
    }
    // LCOV_EXCL_STOP

    m_Pos = 0;
    m_NumberLen = 0;

    PARSE_CHROM();

    return x_ParsePos();
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParsePos()
{
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
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseIDs()
{
    m_NextListIndex = 0;

    SKIP_TO_STATE(eID);

    return x_ParseIDs();
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParseIDs()
{
    PARSE_STRING_LIST(eRef, m_IDs, m_Tokenizer.m_NewlineOrTabOrSemicolon);
    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseAlleles()
{
    m_NextListIndex = 0;

    SKIP_TO_STATE(eRef);

    PARSE_REF();

    return x_ParseAlts();
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParseAlts()
{
    PARSE_STRING_LIST(eQuality, m_Alts, m_Tokenizer.m_NewlineOrTabOrComma);

    m_AllelesParsed = true;

    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseQuality()
{
    SKIP_TO_STATE(eQuality);

    return x_ParseQuality();
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParseQuality()
{
    PARSE_STRING(eFilter);
    if (!m_Tokenizer.TokenIsDot())
        m_Quality = m_Tokenizer.GetToken();
    else
        m_Quality.clear();

    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseFilters()
{
    m_NextListIndex = 0;

    EParsingEvent pe = x_SkipToState(eFilter);
    if (pe != eOK) {
        return pe;
    }

    return x_ParseFilters();
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParseFilters()
{
    PARSE_STRING_LIST(
            eInfoField, m_Filters, m_Tokenizer.m_NewlineOrTabOrSemicolon);

    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseInfo()
{
    m_Info.clear();

    SKIP_TO_STATE(eInfoField);

    return x_ParseInfo();
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParseInfo()
{
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

    return eOK;
}

CVCFScanner::EParsingEvent CVCFScanner::ParseGenotypeFormat()
{
    m_GenotypeKeyPositions.Clear();

    SKIP_TO_STATE(eGenotypeFormat);

    return x_ParseGenotypeFormat();
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParseGenotypeFormat()
{
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

    x_ClearGenotypeValues();
    m_ParsingState = eGenotypes;
    return eOK;
}

bool CVCFScanner::CaptureGT()
{
    unsigned gt_index = m_GenotypeKeyPositions.GT;
    if (gt_index == 0) {
        return false;
    }
    --gt_index;
    auto gt_value = x_AllocGenotypeValue(gt_index);
    gt_value->data_type = eGT;
    gt_value->int_vector = &m_GT;
    return true;
}

const char* CVCFScanner::x_ParseGT()
{
    m_GT.clear();

    const string& token = m_Tokenizer.GetToken();

    size_t len = token.length();

    if (len == 0) {
        return "Empty GT value";
    }

    const char* ptr = token.data();
    unsigned digit, allele;

    for (;; ++ptr, --len) {
        if (*ptr == '.') {
            m_GT.push_back(-1);
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

            m_GT.push_back((int) allele);

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
    // LCOV_EXCL_START
    if (m_ParsingState != eGenotypes) {
        assert(false &&
                "ParseGenotypeFormat must be called before ParseGenotype");
        return eError;
    }
    // LCOV_EXCL_STOP

    if (m_CurrentGenotypeFieldIndex >= m_Header.m_SampleIDs.size()) {
        return x_DataLineError(
                "The number of genotype fields exceeds the number of samples");
    }

    m_CurrentGenotypeValueIndex = 0;

    return x_ParseGenotype();
}

CVCFScanner::EParsingEvent CVCFScanner::x_ParseGenotype()
{
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
                error_message = x_ParseGT();
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
            ++m_CurrentGenotypeFieldIndex;
            return eOK;
        }

        ++value;
    } while (++m_CurrentGenotypeValueIndex <
            m_GenotypeKeyPositions.number_of_positions);

    return x_DataLineError("Too many genotype info fields");
}

CVCFScanner::EParsingEvent CVCFScanner::ClearLine()
{
    if (!m_Tokenizer.AtEOF()) {
        if (m_ParsingState != ePeekAfterEOL) {
            if (m_ParsingState != eEndOfDataLine &&
                    !m_Tokenizer.SkipToken(m_Tokenizer.FindNewline())) {
                m_ParsingState = eClearLine;
                return eNeedMoreData;
            }

            if (m_Tokenizer.BufferIsEmpty()) {
                m_ParsingState = ePeekAfterEOL;
                return eNeedMoreData;
            }
        }
    }

    x_ResetDataLine();

    return eOK;
}
