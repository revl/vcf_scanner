#include <string>
#include <vector>
#include <map>
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
        eOK, // The current token (the VCF header or a data field) has
             // been successfully parsed. Header meta-information or
             // the data field value is now available for retrieval.

        eNeedMoreData, // The parser needs a new input buffer to
                       // continue parsing.

        eWarning, // The parser has generated a warning. Use
                  // GetWarning() to retrieve the warning message.

        eError // A parsing error has occurred. Use GetError() to get
               // the error report.
               //
               // If the error happened while parsing the VCF header,
               // this parser instance can no longer be used.
               //
               // If the error happened while parsing a data line, there
               // is an option to ignore it and skip to the next line
               // by calling ClearLine().
    };

    // GetLineNumber returns the current line number in the
    // input VCF file. The returned value is one-based.
    unsigned GetLineNumber() const
    {
        return m_Tokenizer.GetLineNumber();
    }

    CErrorReport GetWarning() const
    {
        return CErrorReport{"hello"};
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
    EParsingEvent Rewind();

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

    // ParseGenotypeFormats parses the genotype format keys.
    EParsingEvent ParseGenotypeFormats();
    // GetGenotypeFormat returns the FORMAT field parsed by
    // ParseGenotypeFormat.
    const vector<string>& GetGenotypeFormats() const
    {
        return m_GenotypeFormats;
    }

    // ParseGenotype sequentially parses genotype fields.
    EParsingEvent ParseGenotype();
    // GetGenotype returns the genotype field parsed by ParseGenotype.
    string GetGenotype() const
    {
        return m_Tokenizer.GetToken();
    }
    // HasAnotherGenotype returns true if another genotype value is available.
    // ParseGenotype followed by GetGenotype must be called to retrieve it.
    bool HasAnotherGenotype() const
    {
        return m_Tokenizer.GetTokenTerm() == '\t';
    }

    // Finishes parsing of the current line. This method must
    // be called repeatedly until it returns eOK thus confirming
    // that the end of the current line has been reached.
    EParsingEvent ClearLine();

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

    CErrorReport m_ErrorReport;

    CVCFTokenizer m_Tokenizer;

    CVCFHeader m_Header;

    unsigned m_NumberLen;

    string m_Chrom;
    unsigned m_Pos;
    vector<string> m_IDs;
    string m_Ref;
    vector<string> m_Alts;
    string m_Quality;
    vector<string> m_Filters;
    vector<string> m_Info;
    vector<string> m_GenotypeFormats;
    vector<string> m_GenotypeInfo;
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

CVCFScanner::EParsingEvent CVCFScanner::Rewind()
{
    if (m_Tokenizer.AtEOF())
        return eOK;

    if (m_Tokenizer.BufferIsEmpty())
        return eNeedMoreData;

    if (m_DataLineParsingState != eEndOfDataLine &&
            !m_Tokenizer.SkipToken(m_Tokenizer.FindNewline()))
        return eNeedMoreData;

    m_DataLineParsingState = eChrom;
    return eOK;
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

#define REALLY_DO_SKIP_TO_FIELD(field)                                         \
    if (m_DataLineParsingState < field) {                                      \
        do {                                                                   \
            if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTab()))        \
                return eNeedMoreData;                                          \
            if (m_Tokenizer.TokenIsLast())                                     \
                return x_MissingMandatoryFieldError(                           \
                        m_HeaderLineColumns[m_DataLineParsingState + 1]);      \
        } while (++m_DataLineParsingState < field);                            \
    } else {                                                                   \
        assert(false && "Must call Rewind() first");                           \
        return eError;                                                         \
    }

CVCFScanner::EParsingEvent CVCFScanner::ParseIDs()
{
    if (m_DataLineParsingState > eID) {
        assert(false && "Must call Rewind() before ParseIDs()");
        return eError;
    }

    if (m_DataLineParsingState < eID) {
        REALLY_DO_SKIP_TO_FIELD(eID);
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
        REALLY_DO_SKIP_TO_FIELD(eRef);
    }

    if (m_DataLineParsingState == eRef) {
        PARSE_STRING_FIELD(eRef);
        m_Ref = m_Tokenizer.GetToken();

        m_DataLineParsingState = eAlt;
        m_Alts.clear();
    }

    PARSE_LIST_FIELD(eQuality, m_Alts, m_Tokenizer.m_NewlineOrTabOrComma);

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
        REALLY_DO_SKIP_TO_FIELD(eQuality);
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
        REALLY_DO_SKIP_TO_FIELD(eFilter);
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
    if (m_DataLineParsingState != eInfoField) {
        REALLY_DO_SKIP_TO_FIELD(eInfoField);
        m_Info.clear();
    }

    for (;;) {
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindCharFromSet(
                    m_Tokenizer.m_NewlineOrTabOrSemicolon)))
            return eNeedMoreData;
        switch (m_Tokenizer.GetTokenTerm()) {
        case EOF:
        case '\n':
            m_DataLineParsingState = eEndOfDataLine;
            return eOK;
        case '\t':
            m_DataLineParsingState = eGenotypeFormat;
            m_GenotypeFormats.clear();
            return eOK;
        }
        if (!m_Tokenizer.TokenIsDot())
            m_Info.push_back(m_Tokenizer.GetToken());
    }
}

CVCFScanner::EParsingEvent CVCFScanner::ParseGenotypeFormats()
{
    if (m_DataLineParsingState != eGenotypeFormat) {
        REALLY_DO_SKIP_TO_FIELD(eGenotypeFormat);
        m_GenotypeFormats.clear();
    }

    for (;;) {
        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindCharFromSet(
                    m_Tokenizer.m_NewlineOrTabOrColon)))
            return eNeedMoreData;
        switch (m_Tokenizer.GetTokenTerm()) {
        case EOF:
        case '\n':
            m_DataLineParsingState = eEndOfDataLine;
            return eOK;
        case '\t':
            m_DataLineParsingState = eGenotypes;
            m_GenotypeInfo.clear();
            return eOK;
        }
        if (!m_Tokenizer.TokenIsDot())
            m_GenotypeFormats.push_back(m_Tokenizer.GetToken());
    }
}

#define SKIP_TO_FIELD(new_state)                                               \
    for (; m_DataLineParsingState < new_state; ++m_DataLineParsingState) {     \
        if (!m_Tokenizer.SkipToken(m_Tokenizer.FindNewlineOrTab()))            \
            return eNeedMoreData;                                              \
        switch (m_Tokenizer.GetTokenTerm()) {                                  \
        case EOF:                                                              \
        case '\n':                                                             \
            return x_MissingMandatoryFieldError(                               \
                    m_HeaderLineColumns[new_state]);                           \
        }                                                                      \
    }

CVCFScanner::EParsingEvent CVCFScanner::ParseGenotype()
{
    assert(m_DataLineParsingState <= eGenotypes &&
            "Must Rewind() before ParseGenotype()");

    SKIP_TO_FIELD(eGenotypes);

    if (!m_Tokenizer.PrepareTokenOrAccumulate(
                m_Tokenizer.FindCharFromSet(m_Tokenizer.m_NewlineOrTab)))
        return eNeedMoreData;

    switch (m_Tokenizer.GetTokenTerm()) {
    case EOF:
    case '\n':
        m_DataLineParsingState = eEndOfDataLine;
    }

    return eOK;
}

static const char vcf_file[] =
        "##fileformat=VCFv4.0\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\r\n"
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	"
        "INFO	"
        "FORMAT	S-1	S-2	S-3\n"
        "1	100000	rs123;rs333	C	G	10	.	"
        ".	"
        "GT	0|M	1/.	1/0\n"
        "2	200000	.	C	G,T	.	PASS	"
        "NS=3;DP=14;AF=0.5;DB;H2	"
        "GT	0|0	0|1	1|E";

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

#define RETRY_UNTIL_OK_OR_ERROR(vcf_scanner_call)                              \
    for (;;) {                                                                 \
        switch (pe = (vcf_scanner_call)) {                                     \
        case CVCFScanner::eWarning:                                            \
            cerr << vcf_scanner.GetWarning().m_ErrorMessage << endl;           \
            continue;                                                          \
        case CVCFScanner::eNeedMoreData:                                       \
            vcf_scanner.SetNewInputBuffer(                                     \
                    buffer, s_ReadVCFFile(buffer, sizeof(buffer)));            \
            continue;                                                          \
        default:                                                               \
            break;                                                             \
        }                                                                      \
        break;                                                                 \
    }

#define RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner_call)                \
    RETRY_UNTIL_OK_OR_ERROR(vcf_scanner_call);                                 \
    if (pe == CVCFScanner::eError) {                                           \
        cerr << "file.vcf:" << vcf_scanner.GetLineNumber() << ": "             \
             << vcf_scanner.GetError().m_ErrorMessage << endl                  \
             << endl;                                                          \
        continue;                                                              \
    }

int main()
{
    char buffer[1];

    CVCFScanner vcf_scanner;

    vcf_scanner.SetNewInputBuffer(
            buffer, s_ReadVCFFile(buffer, sizeof(buffer)));

    CVCFScanner::EParsingEvent pe;

    // Read the header
    RETRY_UNTIL_OK_OR_ERROR(vcf_scanner.ParseHeader());

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
        RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.Rewind());
        if (vcf_scanner.AtEOF())
            break;

        cout << '#' << vcf_scanner.GetLineNumber() << endl;

        RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseLoc());
        cout << "CHROM:[" << vcf_scanner.GetChrom() << ']' << endl
             << "POS:[" << vcf_scanner.GetPos() << ']' << endl;

        RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseIDs());
        for (const auto& id : vcf_scanner.GetIDs())
            cout << "ID:[" << id << ']' << endl;

        RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseAlleles());
        cout << "REF:[" << vcf_scanner.GetRef() << ']' << endl;
        for (const auto& alt : vcf_scanner.GetAlts())
            cout << "ALT:[" << alt << ']' << endl;

        /* RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseQuality());
        cout << "QUAL:["
             << vcf_scanner.GetQuality() << ']' << endl;

        RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseFilter());
        cout << "FILTER:["
             << vcf_scanner.GetFilter();
        while (vcf_scanner.HasAnotherFilter()) {
            RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseFilter());
            cout << " and " << vcf_scanner.GetFilter();
        }
        cout << ']' << endl;

        RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseInfo());
        cout << "INFO:["
             << vcf_scanner.GetInfo();
        while (vcf_scanner.HasMoreInfo()) {
            RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseInfo());
            cout << " and " << vcf_scanner.GetInfo();
        }
        cout << ']' << endl;

        RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(
                vcf_scanner.ParseGenotypeFormat());
        cout << "FORMAT:["
             << vcf_scanner.GetGenotypeFormat();
        while (vcf_scanner.HasAnotherGenotypeFormat()) {
            RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(
                    vcf_scanner.ParseGenotypeFormat());
            cout << " and " << vcf_scanner.GetGenotypeFormat();
        }
        cout << ']' << endl; */

        RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseGenotype());
        cout << "GENOTYPE:[" << vcf_scanner.GetGenotype();
        while (vcf_scanner.HasAnotherGenotype()) {
            RETRY_UNTIL_OK_AND_SKIP_LINE_ON_ERROR(vcf_scanner.ParseGenotype());
            cout << " and " << vcf_scanner.GetGenotype();
        }
        cout << ']' << endl;

        cout << endl;
    }

    /* // Read data lines
        case CVCFScanner::eInfoField:
            break;
        case CVCFScanner::eGenotypeField:
            break;
        case CVCFScanner::eEndOfDataLine:
            break;
        case CVCFScanner::eEndOfStream:
            return 0;
        }
    }

    pe = vcf_scanner.NextEvent(); */

    return 0;
}
