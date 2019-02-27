#include <string>
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
    // GetFileFormat returns the file format version.
    const string& GetFileFormat() const
    {
        return m_FileFormat;
    }

private:
    string m_FileFormat;
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
// on I/O operations.
//
// Before parsing begins, or when a parsing function returns eNeedMoreData,
// a new buffer with input data must be supplied to the parser by calling
// SetNewInputBuffer(). The buffer must not be freed or overwritten until
// eNeedMoreData is received again or the client code chooses not to continue
// parsing.
class CVCFScanner
{
public:
    CVCFScanner();

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
        return m_LineNumber;
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
    enum ECurrentState {
        eParseVersion,
        eMetaInfoLineParsing,
        eHeaderLineParsing,
        eHeaderError,
    } m_CurrentState;

    bool x_IsParsingHeader() const
    {
        return m_CurrentState <= eHeaderError;
    }

    bool x_IsAtEOF() const
    {
        return m_BufferSize == 0;
    }

    EParsingEvent x_HeaderError(const char* error_message)
    {
        m_ErrorReport.m_ErrorMessage = error_message;
        m_CurrentState = eHeaderError;
        return m_LastParsingEvent = eError;
    }

    EParsingEvent m_LastParsingEvent;

    const char* m_Buffer;
    size_t m_BufferSize;

    unsigned m_LineNumber;
    CErrorReport m_ErrorReport;

    CVCFTokenizer m_Tokenizer;

    CVCFHeader m_Header;
};

CVCFScanner::CVCFScanner() :
    m_CurrentState(eParseVersion),
    m_LastParsingEvent(eNeedMoreData),
    m_BufferSize(0),
    m_LineNumber(0)
{}

void CVCFScanner::SetNewInputBuffer(const char* buffer, ssize_t buffer_size)
{
    assert(m_LastParsingEvent == eNeedMoreData);

    /* if (!m_SeamBuffer.empty() && buffer_size <= 0) {
        m_ErrorReport.m_ErrorMessage = x_IsParsingHeader() ?
                "Unexpected EOF while parsing VCF header" :
                "Unexpected EOF while parsing a data line";
        m_LastParsingEvent = eError;
    } */

    m_Buffer = buffer;
    m_BufferSize = buffer_size;

    m_Tokenizer.SetNewBuffer(buffer, buffer_size);
}

CVCFScanner::EParsingEvent CVCFScanner::ParseHeader()
{
    assert(x_IsParsingHeader());

    // Extract \n (header, skip line) ?+ (\t (all data fields) ?+ ; (id, filter,
    // info) or , (alt) or : (format, genotype))

    char meta_info_prefix[2];

    switch (m_CurrentState) {
    case eParseVersion:
        if (x_IsAtEOF())
            goto UnexpectedEOF;

        if (!m_Tokenizer.PrepareTokenOrAccumulate(m_Tokenizer.FindNewline()))
            return eNeedMoreData;

        {
            string key, value;

            if (!m_Tokenizer.GetKeyValue(&key, &value) ||
                    key != "##fileformat") {
                return x_HeaderError("VCF file must start with '##fileformat'");
            }

            m_Header.m_FileFormat = value;

            m_CurrentState = eMetaInfoLineParsing;
        }
        break;

    case eMetaInfoLineParsing:
        if (x_IsAtEOF())
            goto UnexpectedEOF;

        for (;;) {
            if (!m_Tokenizer.Peek(meta_info_prefix, 2))
                return eNeedMoreData;

            if (meta_info_prefix[0] != '#')
                goto InvalidMetaInfoLinePrefix;

            if (meta_info_prefix[1] == '#') {
                if (m_Tokenizer.PrepareTokenOrAccumulate(
                            m_Tokenizer.FindNewline())) {
                    cout << m_Tokenizer.GetToken() << endl;
                }

                return eNeedMoreData;
            } else if (meta_info_prefix[1] == 'C') {
                m_CurrentState = eHeaderLineParsing;
                break;
            } else
                goto InvalidMetaInfoLinePrefix;
        }

    case eHeaderLineParsing:

    // Cannot recover from errors in the header.
    case eHeaderError:
        break;
    }
    return eError;

UnexpectedEOF:
    return x_HeaderError(
            "Unexpected end of file while parsing VCF file header");

InvalidMetaInfoLinePrefix:
    return x_HeaderError("Invalid meta-information line prefix");
}

static const char vcf_file[] = R"(##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FORMAT	INFO	FORMAT	S-1	S-2	S-3
1	100000	.	C	G	.	.	.	GT	0|M	1/.	1/0
2	200000	.	C	G	.	.	.	GT	0|0	0|1	1|E)";

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

#define ERR(what)                                                              \
    {                                                                          \
        cerr << what << endl;                                                  \
        return 1;                                                              \
    }

int main()
{
    char buffer[1];

    CVCFScanner vcf_scanner;

    CVCFScanner::EParsingEvent pe;

    // Read the header
    do {
        vcf_scanner.SetNewInputBuffer(
                buffer, s_ReadVCFFile(buffer, sizeof(buffer)));

        while ((pe = vcf_scanner.ParseHeader()) == CVCFScanner::eWarning)
            cerr << vcf_scanner.GetWarning().m_ErrorMessage << endl;
    } while (pe == CVCFScanner::eNeedMoreData);

    if (pe != CVCFScanner::eOK) {
        cerr << vcf_scanner.GetError().m_ErrorMessage << endl;
        return 1;
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
