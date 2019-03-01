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
    // segmentation fault. An buffer of size zero is treated as an EOF
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

    // HasAnotherDataLine returns true if the VCF file has more
    // data lines to parse and false if the entire stream has been
    // successfully parsed.
    bool HasAnotherDataLine() const;

    // GetLineNumber returns the current one-based line number
    // in the input VCF file.
    unsigned GetLineNumber() const
    {
        return m_Tokenizer.GetLineNumber();
    }

    // ParseChrom parses the CHROM field.
    EParsingEvent ParseChrom();
    // GetChrom returns the CHROM field parsed by ParseChrom.
    const string& GetChrom() const;

    // ParsePos parses the POS field.
    EParsingEvent ParsePos();
    // GetPos returns the POS field parsed by ParsePos.
    int GetPos() const;

    // ParseID parses the ID field.
    EParsingEvent ParseID();
    // GetID returns the POS field parsed by ParsePos.
    // A dot (".") is returned when no ID is available.
    string GetID() const;
    // HasAnotherID returns true if another ID is available for
    // this record. In which case, ParseID() must be called again,
    // whereupon GetID() will return the next value.
    bool HasAnotherID() const;

    // ParseRef parses the REF field.
    EParsingEvent ParseRef();
    // GetRef returns the REF field parsed by ParseRef.
    string GetRef() const;

    // ParseAlt parses the ALT field.
    EParsingEvent ParseAlt();
    // GetAlt returns the ALT field parsed by ParseAlt.
    // A dot (".") is returned when there are no alternative alleles.
    string GetAlt() const;
    // HasAnotherAlt returns true if another ALT is available.
    // ParseAlt followed by GetAlt must be called to retrieve it.
    bool HasAnotherAlt() const;

    // ParseQuality parses the QUAL field.
    EParsingEvent ParseQuality();
    // GetQuality returns the QUAL field parsed by ParseQuality.
    string GetQuality() const;

    // ParseFilter parses the FILTER field.
    EParsingEvent ParseFilter();
    // GetFilter returns the FILTER field parsed by ParseFilter.
    // The word "PASS" is returned when the current record passed
    // all filters. If filters have not been applied, GetFilter
    // returns a dot (".").
    string GetFilter() const;
    // HasAnotherFilter returns true if another FILTER is available.
    // ParseFilter followed by GetFilter must be called to retrieve it.
    bool HasAnotherFilter() const;

    // ParseInfo parses the INFO key-value pairs.
    EParsingEvent ParseInfo();
    // GetInfo returns the INFO field parsed by ParseInfo.
    // If the INFO field is missing, GetInfo returns a dot (".").
    string GetInfo() const;
    // HasMoreInfo returns true if another INFO field is available.
    // ParseInfo followed by GetInfo must be called to retrieve it.
    bool HasMoreInfo() const;

    // ParseFormat parses the genotype format keys.
    EParsingEvent ParseFormat();
    // GetFormat returns the FORMAT field parsed by ParseFormat.
    string GetFormat() const;
    // HasAnotherFormat returns true if another FORMAT is available.
    // ParseFormat followed by GetFormat must be called to retrieve it.
    bool HasAnotherFormat() const;

    // ParseGenotype sequentially parses genotype fields.
    EParsingEvent ParseGenotype();
    // GetGenotype returns the genotype field parsed by ParseGenotype.
    string GetGenotype() const;
    // HasAnotherGenotype returns true if another genotype value is available.
    // ParseGenotype followed by GetGenotype must be called to retrieve it.
    bool HasAnotherGenotype() const;

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
        eSampleIDs,
    } m_HeaderParsingState = eFileFormatVersion;

    /* enum {
        eChrom,
        ePos,
        eID,
        eRef,
        eAlt,
        eQuality,
        ePassedAllFormats,
        eFailedFormat,
        eInfoField,
        eGenotypeField,
        eEndOfDataLine,
        eEndOfStream,
    } */

    string m_CurrentMetaInfoKey;

    unsigned m_HeaderLineColumnOK;

    EParsingEvent x_HeaderError(const char* error_message)
    {
        m_ErrorReport.m_ErrorMessage = error_message;
        return eError;
    }

    const char* m_Buffer;
    size_t m_BufferSize = 0;

    CErrorReport m_ErrorReport;

    CVCFTokenizer m_Tokenizer;

    CVCFHeader m_Header;
};

void CVCFScanner::SetNewInputBuffer(const char* buffer, ssize_t buffer_size)
{
    m_Buffer = buffer;
    m_BufferSize = buffer_size;

    m_Tokenizer.SetNewBuffer(buffer, buffer_size);
}

static constexpr unsigned kNumberOfMandatoryColumns = 8;

static const char* s_HeaderLineColumns[kNumberOfMandatoryColumns + 1] = {
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT"};

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
            if (m_Tokenizer.GetToken() != s_HeaderLineColumns[0])
                goto InvalidMetaInfoLine;
            m_HeaderLineColumnOK = 1;
            goto ParseHeaderColumns;
        case '\n':
            goto InvalidMetaInfoLine;
        case EOF:
            goto UnexpectedEOF;
        }

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
                    s_HeaderLineColumns[m_HeaderLineColumnOK])
                goto InvalidHeaderLine;

            ++m_HeaderLineColumnOK;

            switch (m_Tokenizer.GetTokenTerm()) {
            case '\n':
            case EOF:
                switch (m_HeaderLineColumnOK) {
                case kNumberOfMandatoryColumns + 1:
                    m_Header.m_GenotypeInfoPresent = true;
                    /* FALL THROUGH */
                case kNumberOfMandatoryColumns:
                    return eOK;
                default:
                    goto InvalidHeaderLine;
                }
            }

            // The current token ends with a tab.
        } while (m_HeaderLineColumnOK <= kNumberOfMandatoryColumns);

        m_Header.m_GenotypeInfoPresent = true;
        m_HeaderParsingState = eSampleIDs;
        /* FALL THROUGH */

    case eSampleIDs:
        for (;;) {
            if (!m_Tokenizer.PrepareTokenOrAccumulate(
                        m_Tokenizer.FindNewlineOrTab()))
                return eNeedMoreData;

            m_Header.m_SampleIDs.push_back(m_Tokenizer.GetToken());

            switch (m_Tokenizer.GetTokenTerm()) {
            case '\n':
            case EOF:
                return eOK;
            }
        }
    }

UnexpectedEOF:
    return x_HeaderError(
            "Unexpected end of file while parsing VCF file header");

InvalidMetaInfoLine:
    return x_HeaderError("Malformed meta-information line");

InvalidHeaderLine:
    return x_HeaderError("Malformed VCF header line");
}

static const char vcf_file[] =
        "##fileformat=VCFv4.0\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\r\n"
        "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	"
        "INFO	"
        "FORMAT	S-1	S-2	S-3\n"
        "1	100000	.	C	G	.	.	.	"
        "GT	0|M	1/.	1/0\n"
        "2	200000	.	C	G	.	.	.	"
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
    do {                                                                       \
        vcf_scanner.SetNewInputBuffer(                                         \
                buffer, s_ReadVCFFile(buffer, sizeof(buffer)));                \
        while ((pe = vcf_scanner_call) == CVCFScanner::eWarning)               \
            cerr << vcf_scanner.GetWarning().m_ErrorMessage << endl;           \
    } while (pe == CVCFScanner::eNeedMoreData)

int main()
{
    char buffer[1];

    CVCFScanner vcf_scanner;

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

    cout << "[Sample IDs]" << endl;
    for (const auto& v : vcf_scanner.GetHeader().GetSampleIDs()) {
        cout << v << endl;
    }

    /* // Read data lines
    for (;;) {
        switch (vcf_scanner.NextEvent()) {
        case CVCFScanner::eNeedMoreData:
            vcf_scanner.SetNewInputBuffer(
                buffer, s_ReadVCFFile(buffer, sizeof(buffer)));
            break;
        case CVCFScanner::eWarning:
            cerr << vcf_scanner.GetWarning().m_ErrorMessage << endl;
            break;
        case CVCFScanner::eError:
            cerr << vcf_scanner.GetWarning().m_ErrorMessage << endl;
            vcf_scanner.SkipToNextDataLine();
            break;
        case CVCFScanner::eChrom:
            break;
        case CVCFScanner::ePos:
            break;
        case CVCFScanner::eID:
            break;
        case CVCFScanner::eRef:
            break;
        case CVCFScanner::eAlt:
            break;
        case CVCFScanner::eQuality:
            break;
        case CVCFScanner::ePassedAllFormats:
            break;
        case CVCFScanner::eFailedFormat:
            break;
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
