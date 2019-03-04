#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <cassert>
#include <stdio.h>

#include "tokenizer.hh"

using namespace std;

class CErrorReport
{
public:
    string m_ErrorMessage;
};

class CVCFScanner;

// CVCFHeader contains information extracted from the VCF header.
class CVCFHeader
{
public:
    typedef vector<string> TMetaInfoLines;
    typedef map<string, TMetaInfoLines> TMetaInfo;
    typedef vector<string> TSampleIDs;

    // GetFileFormat returns the file format version.
    const string& GetFileFormat() const
    {
        return m_FileFormat;
    }

    const TMetaInfo& GetMetaInfo() const
    {
        return m_MetaInfo;
    }

    const TSampleIDs& GetSampleIDs() const
    {
        return m_SampleIDs;
    }

private:
    string m_FileFormat;
    TMetaInfo m_MetaInfo;
    bool m_GenotypeInfoPresent = false;
    TSampleIDs m_SampleIDs;

    friend class CVCFScanner;
};

struct CVCFWarning {
    unsigned line_number;
    string warning_message;
};

// CVCFScanner parses VCF (Variant Call Format) files.
//
// First, CVCFScanner parses the header in its entirety,
// and then it parses data lines one by one.
//
// All header information is kept by the parser in its member variables.
// That includes sample IDs from the header line.
//
// The fields of the data lines, however, are never stored internally.
// CVCFScanner is merely a tokenizer. Each data line field is discarded
// as soon as the caller proceeds to parsing the next field.
//
// CVCFScanner does not have a stream reading loop inside. It relies on
// the client code to provide the input data. As a result, it never blocks
// on I/O operations. This allows for using a separate thread to read data
// into a new buffer while the main thread is parsing a previously read
// buffer. Alternatively to reading the input file into memory, the whole
// file can be memory-mapped.
//
// Before parsing begins, or when a parsing function returns eNeedMoreData,
// a new buffer with input data must be supplied to the parser by calling
// SetNewInputBuffer(). The buffer must not be freed or overwritten until
// eNeedMoreData is received again or the client code chooses not to continue
// parsing.
class CVCFScanner
{
public:
    CVCFScanner() {}

    // SetNewInputBuffer sets the next buffer to parse. The parser stores
    // the buffer pointer internally. Freeing the buffer will cause a
    // segmentation fault. A buffer of zero size is treated as an EOF
    // condition.
    void SetNewInputBuffer(const char* buffer, ssize_t buffer_size);

    enum EParsingEvent {
        eNeedMoreData, // The parser needs a new input buffer to
                       // continue parsing.

        eOK, // The current token (the VCF header or a data field) has
             // been successfully parsed. Header meta-information or
             // the data field value is now available for retrieval.

        eOKWithWarnings, // The token has been successfully parsed, but
                         // parser encountered issues during parsing. Use
                         // GetWarnings() to retrieve the warning messages.

        eError // A parsing error has occurred. Use GetError() to get
               // the error report.
               //
               // If the error happened while parsing the VCF header,
               // this parser instance can no longer be used.
               //
               // If the error happened while parsing a data line, there
               // is an option to ignore it and skip to the next line
               // by calling Rewind().
    };

    // GetLineNumber returns the current line number in the
    // input VCF file. The returned value is one-based.
    unsigned GetLineNumber() const
    {
        return m_Tokenizer.GetLineNumber();
    }

    vector<CVCFWarning> GetWarnings() const
    {
        return m_Warnings;
    }

    CErrorReport GetError() const
    {
        return m_ErrorReport;
    }

    // ParseHeader begins or continues parsing the header.
    // It returns eOK when the entire header has been parsed.
    EParsingEvent ParseHeader();

    // GetHeader returns the VCF header structure parsed by
    // ParseHeader().
    const CVCFHeader& GetHeader() const
    {
        return m_Header;
    }

    // Rewind determines whether end-of-file has been reached.
    // It also positions the "reading head" of the parser at
    // the beginning of the next data line. If the previous
    // line was skipped due to an error or due to the lack of
    // interest in it from the client code, this operation may
    // require more data to be read.
    // Rewind returns false if the parser needs more input data
    // to detect the end-of-file condition. Supply more data by
    // calling SetNewInputBuffer until Rewind returns true.
    bool Rewind();

    // AtEOF returns true if the entire input stream has been
    // successfully parsed. The method returns false if the
    // VCF file has at least one more data line to parse.
    bool AtEOF() const
    {
        return m_Tokenizer.AtEOF();
    }

    // ParseLoc parses the CHROM and the POS fields.
    EParsingEvent ParseLoc();
    // GetChrom returns the CHROM field parsed by ParseLoc.
    const string& GetChrom() const
    {
        return m_Chrom;
    }
    // GetPos returns the POS field parsed by ParseLoc.
    unsigned GetPos() const
    {
        return m_Pos;
    }

    // ParseIDs parses the ID field.
    EParsingEvent ParseIDs();
    // GetID returns the IDs parsed by ParsePos.
    const vector<string>& GetIDs() const
    {
        return m_IDs;
    }

    // ParseRef parses the REF and the ALT fields.
    EParsingEvent ParseAlleles();
    // GetRef returns the REF field parsed by ParseAlleles.
    const string& GetRef() const
    {
        return m_Ref;
    }
    // GetAlts returns the ALT field parsed by ParseAlleles.
    const vector<string>& GetAlts() const
    {
        return m_Alts;
    }

    // ParseQuality parses the QUAL field.
    EParsingEvent ParseQuality();
    // GetQuality returns the QUAL field parsed by ParseQuality.
    string GetQuality() const
    {
        return m_Quality;
    }

    // ParseFilters parses the FILTER field.
    EParsingEvent ParseFilters();
    // GetFilter returns the FILTER field parsed by ParseFilters.
    // The word "PASS" is returned when the current record passed
    // all filters.
    const vector<string>& GetFilters() const
    {
        return m_Filters;
    }

    // ParseInfo parses the INFO key-value pairs.
    EParsingEvent ParseInfo();
    // GetInfo returns the INFO field parsed by ParseInfo.
    const vector<string>& GetInfo() const
    {
        return m_Info;
    }

    // ParseGenotypeFormat parses the genotype format keys.
    EParsingEvent ParseGenotypeFormat();

    // GetGenotypeFormat returns the FORMAT field parsed by
    // ParseGenotypeFormat.
    // TODO const TGenotypeFormat& GetGenotypeFormat() const;

    // CaptureGT binds or re-binds a client-side variable to receive
    // GT values as they are parsed by the ParseGenotype method below.
    // Use the IsPhasedGT method to check whether the sample was phased.
    // CaptureGT returns false and does nothing if the GT key was not
    // specified in the FORMAT field.
    // Set alleles to NULL to stop receiving the GT values.
    bool CaptureGT(vector<int>* alleles);

    // TODO bool CaptureString(const char* key, string* value);
    // TODO bool CaptureStrings(const char* key, vector<string>* values);
    // TODO bool CaptureInt(const char* key, int* value);
    // TODO bool CaptureInts(const char* key, vector<int>* values);

    // ParseGenotype sequentially parses genotype fields.
    EParsingEvent ParseGenotype();

    // IsPhasedGT returns true if the parsed sample was phased.
    bool IsPhasedGT() const
    {
        return m_PhasedGT;
    }

    // GenotypeAvailable returns true if at least one more genotype
    // field is available. The caller has an option to either use
    // this method or count the retrieved genotypes to determine
    // when to stop parsing the current data line.
    bool GenotypeAvailable() const
    {
        return m_Tokenizer.GetTokenTerm() == '\t';
    }

private:
    enum EHeaderParsingState {
        eFileFormatVersion,
        eMetaInfoKey,
        eMetaInfoValue,
        eHeaderLineColumns,
        eSampleIDs
    } m_HeaderParsingState = eFileFormatVersion;

    string m_CurrentMetaInfoKey;

    enum EColumn {
        eChrom,
        ePos,
        eID,
        eRef,
        eAlt,
        eQuality,
        eFilter,
        eInfoField,
        eGenotypeFormat,
        eGenotypes,
        eEndOfDataLine,
        // FORMAT is the first optional column
        eNumberOfMandatoryColumns = eGenotypeFormat
    };
    int m_DataLineParsingState = eEndOfDataLine;

    static const char* const m_HeaderLineColumns[];

    unsigned m_HeaderLineColumnOK;

    EParsingEvent x_HeaderError(const char* error_message)
    {
        m_ErrorReport.m_ErrorMessage = error_message;
        return eError;
    }

    EParsingEvent x_DataLineError(const string& msg)
    {
        m_ErrorReport.m_ErrorMessage = msg;
        return eError;
    }

    EParsingEvent x_MissingMandatoryFieldError(const char* field_name)
    {
        m_DataLineParsingState = eEndOfDataLine;

        string msg = "Missing mandatory VCF field \"";
        msg += field_name;
        msg += '"';
        return x_DataLineError(msg);
    }

    const char* m_Buffer;
    size_t m_BufferSize = 0;

    vector<CVCFWarning> m_Warnings;
    CErrorReport m_ErrorReport;

    CVCFTokenizer m_Tokenizer;

    CVCFHeader m_Header;

    unsigned m_NumberLen;

    string m_Chrom;
    unsigned m_Pos;
    vector<string> m_IDs;
    string m_Ref;
    vector<string> m_Alts;
    bool m_AllelesParsed;
    string m_Quality;
    vector<string> m_Filters;
    vector<string> m_Info;

    // TODO Implement the INFO and FORMAT type definitions in the header.
    set<string> m_FormatKeys;

    struct CLessCStr {
        bool operator()(const char* left, const char* right) const
        {
            return strcmp(left, right) < 0;
        }
    };

    // Positions of the reserved genotype keys in the FORMAT field.
    struct SGenotypeKeyPositions {
        unsigned number_of_positions;
        unsigned GT;
        map<const char*, unsigned, CLessCStr> other_keys;

        void Clear()
        {
            GT = number_of_positions = 0;
            other_keys.clear();
        }
    } m_GenotypeKeyPositions;

    enum EDataType { eInteger, eFloat, eFlag, eCharacter, eString, eGT };

    enum ENumberOfValues {
        eScalar,
        eOnePerAlt,
        eOnePerAllele,
        eOnePerGenotype,
        eUnbound,
        eExactNumber
    };

    struct SGenotypeValue {
        EDataType data_type;
        unsigned number_of_values;
        union {
            bool* flag;
            int* int_scalar;
            vector<int>* int_vector;
            string* string_scalar;
            vector<string>* string_vector;
            char* char_scalar;
            vector<char>* char_vector;
        };
    };

    vector<SGenotypeValue> m_GenotypeValues;
    unsigned m_CurrentGenotypeValueIndex;
    union {
        int m_IntAcc;
        unsigned m_UIntAcc;
    };
    unsigned m_Ploidy;
    bool m_PhasedGT;

    void x_ClearGenotypeValues()
    {
        memset(m_GenotypeValues.data(), 0,
                (char*) &*m_GenotypeValues.end() -
                        (char*) m_GenotypeValues.data());
        m_CurrentGenotypeValueIndex = 0;
        m_Ploidy = 0;
        m_NumberLen = 0;
    }

    SGenotypeValue* x_AllocGenotypeValue(unsigned index)
    {
        if (m_GenotypeValues.size() <= index) {
            size_t old_size = m_GenotypeValues.size();
            m_GenotypeValues.resize(index + 1);
            char* new_struct_ptr = (char*) (m_GenotypeValues.data() + old_size);
            memset(new_struct_ptr, 0,
                    (char*) &*m_GenotypeValues.end() - (char*) new_struct_ptr);
        }
        return m_GenotypeValues.data() + index;
    }

    const char* x_ParseGT(vector<int>* int_vector);
};

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
            string key, value;

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
        m_CurrentMetaInfoKey = m_Tokenizer.GetToken();

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
