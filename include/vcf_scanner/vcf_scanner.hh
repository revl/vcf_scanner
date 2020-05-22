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

#ifndef VCF_SCANNER__HH
#define VCF_SCANNER__HH

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <cassert>
#include <stdio.h>

#include "impl/vcf_tokenizer.hh"

namespace vcf {

class CVCFScanner;

// CVCFHeader contains information extracted from the VCF header.
class CVCFHeader
{
public:
    typedef std::vector<std::string> TMetaInfoLines;
    typedef std::map<std::string, TMetaInfoLines> TMetaInfo;
    typedef std::vector<std::string> TSampleIDs;

    // GetFileFormat returns the file format version.
    const std::string& GetFileFormat() const
    {
        return m_FileFormat;
    }

    const TMetaInfo& GetMetaInfo() const
    {
        return m_MetaInfo;
    }

    bool HasGenotypeInfo() const
    {
        return m_GenotypeInfoPresent;
    }

    const TSampleIDs& GetSampleIDs() const
    {
        return m_SampleIDs;
    }

private:
    std::string m_FileFormat;
    TMetaInfo m_MetaInfo;
    bool m_GenotypeInfoPresent = false;
    TSampleIDs m_SampleIDs;

    friend class CVCFScanner;
};

struct CVCFWarning {
    unsigned line_number;
    std::string warning_message;
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
// Feed(). The buffer must not be freed or overwritten until eNeedMoreData
// is received again or the client code chooses not to continue parsing.
class CVCFScanner
{
public:
    CVCFScanner() {}

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
               // by calling ClearLine().
    };

    // Feed supplies a chunk of input data to this parser either
    // when the parser has just been created and is in the process
    // of parsing the VCF header or when a previously called method
    // returned eNeedMoreData.
    //
    // The parser stores the buffer pointer internally. Freeing the
    // buffer before the parser returns eNeedMoreData from any of its
    // methods will cause a segmentation fault. A buffer of zero size
    // is treated as an EOF condition.
    //
    // Feed resumes parsing of the previously requested token and
    // returns eOK when the entire token has been parsed.
    EParsingEvent Feed(const char* buffer, ssize_t buffer_size);

    // GetLineNumber returns the current line number in the
    // input VCF file before parsing the next token. The line
    // number will increase after the last token on the current
    // line has been parsed. The returned value is one-based.
    unsigned GetLineNumber() const
    {
        return m_Tokenizer.GetLineNumber();
    }

    std::vector<CVCFWarning> GetWarnings() const
    {
        return m_Warnings;
    }

    std::string GetError() const
    {
        return m_ErrorReport;
    }

    // GetHeader returns the VCF header, which becomes available
    // once the last of the initial series of calls to Feed
    // returns eOK.
    const CVCFHeader& GetHeader() const
    {
        return m_Header;
    }

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
    const std::string& GetChrom() const
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
    // GetID returns the IDs parsed by ParseIDs.
    const std::vector<std::string>& GetIDs() const
    {
        return m_IDs;
    }

    // ParseRef parses the REF and the ALT fields.
    EParsingEvent ParseAlleles();
    // GetRef returns the REF field parsed by ParseAlleles.
    const std::string& GetRef() const
    {
        return m_Ref;
    }
    // GetAlts returns the ALT field parsed by ParseAlleles.
    const std::vector<std::string>& GetAlts() const
    {
        return m_Alts;
    }

    // ParseQuality parses the QUAL field.
    EParsingEvent ParseQuality();
    // GetQuality returns the QUAL field parsed by ParseQuality.
    std::string GetQuality() const
    {
        return m_Quality;
    }

    // ParseFilters parses the FILTER field.
    EParsingEvent ParseFilters();
    // GetFilter returns the FILTER field parsed by ParseFilters.
    // The word "PASS" is returned when the current record passed
    // all filters.
    const std::vector<std::string>& GetFilters() const
    {
        return m_Filters;
    }

    // ParseInfo parses the INFO key-value pairs.
    EParsingEvent ParseInfo();
    // GetInfo returns the INFO field parsed by ParseInfo.
    const std::vector<std::string>& GetInfo() const
    {
        return m_Info;
    }

    // ParseGenotypeFormat parses the genotype format keys.
    EParsingEvent ParseGenotypeFormat();

    // GetGenotypeFormat returns the FORMAT field parsed by
    // ParseGenotypeFormat.
    // TODO const TGenotypeFormat& GetGenotypeFormat() const;

    // CaptureGT enables parsing of GT values in the ParseGenotype method.
    // CaptureGT returns false and does nothing if the GT key was not
    // specified in the FORMAT field.
    bool CaptureGT();

    // TODO bool CaptureString(const char* key, std::string* value);
    // TODO bool CaptureStrings(const char* key, std::vector<std::string>*
    // values);
    // TODO bool CaptureInt(const char* key, int* value);
    // TODO bool CaptureInts(const char* key, std::vector<int>* values);

    // ParseGenotype sequentially parses genotype fields.
    EParsingEvent ParseGenotype();

    // GetGT returns the GT values parsed by ParseGenotype in
    // its previous iteration.
    const std::vector<int>& GetGT() const
    {
        return m_GT;
    }
    // IsPhasedGT returns true if the parsed sample was phased.
    bool IsPhasedGT() const
    {
        return m_PhasedGT;
    }

    // GenotypeAvailable returns true if at least one more genotype
    // field is available. The caller has an option to either use
    // this method or count the retrieved genotypes to determine
    // when the last genotype on the current data line has been parsed.
    bool GenotypeAvailable() const
    {
        return m_Tokenizer.GetTokenTerm() == '\t';
    }

    // ClearLine skips the remaining part of the current data line.
    // This operation may require more data to be read. The client
    // code must call this method after parsing each line even if
    // the line has been parsed to the end, because this method also
    // determines whether end-of-file has been reached.
    EParsingEvent ClearLine();

private:
    enum EParsingState {
        eFileFormatVersion,
        eMetaInfoKey,
        eMetaInfoValue,
        eHeaderLineColumns,
        eSampleIDs,
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
        eClearLine,
        ePeekAfterEOL
    };
    int m_ParsingState = eFileFormatVersion;

    int m_FieldsToSkip = 0;

    std::string m_CurrentMetaInfoKey;

    unsigned m_HeaderLineColumnOK;

    EParsingEvent x_HeaderError(const char* error_message)
    {
        m_ErrorReport = error_message;
        return eError;
    }

    EParsingEvent x_InvalidMetaInfoLineError()
    {
        return x_HeaderError("Malformed meta-information line");
    }

    EParsingEvent x_InvalidHeaderLineError()
    {
        return x_HeaderError("Malformed VCF header line");
    }

    EParsingEvent x_DataLineError(const std::string& msg)
    {
        m_ErrorReport = msg;
        return eError;
    }

    EParsingEvent x_MissingMandatoryFieldError(const char* field_name)
    {
        m_ParsingState = eEndOfDataLine;

        std::string msg = "Missing mandatory VCF field \"";
        msg += field_name;
        msg += '"';
        return x_DataLineError(msg);
    }

    std::vector<CVCFWarning> m_Warnings;
    std::string m_ErrorReport;

    CVCFTokenizer m_Tokenizer;

    CVCFHeader m_Header;

    unsigned m_NumberLen;
    size_t m_NextListIndex;

    std::string m_Chrom;
    unsigned m_Pos;
    std::vector<std::string> m_IDs;
    std::string m_Ref;
    std::vector<std::string> m_Alts;
    bool m_AllelesParsed;
    std::string m_Quality;
    std::vector<std::string> m_Filters;
    std::vector<std::string> m_Info;

    void x_ResetDataLine()
    {
        m_ParsingState = eChrom;
        m_AllelesParsed = false;
    }

    // TODO Implement the INFO and FORMAT type definitions in the header.
    std::set<std::string> m_FormatKeys;

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
        std::map<const char*, unsigned, CLessCStr> other_keys;

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
            std::vector<int>* int_vector;
            std::string* string_scalar;
            std::vector<std::string>* string_vector;
            char* char_scalar;
            std::vector<char>* char_vector;
        };
    };

    unsigned m_CurrentGenotypeFieldIndex;

    std::vector<SGenotypeValue> m_GenotypeValues;
    unsigned m_CurrentGenotypeValueIndex;

    std::vector<int> m_GT;
    bool m_PhasedGT;

    void x_ClearGenotypeValues()
    {
        memset(m_GenotypeValues.data(), 0,
                (char*) &*m_GenotypeValues.end() -
                        (char*) m_GenotypeValues.data());
        m_CurrentGenotypeFieldIndex = 0;
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

    EParsingEvent x_SkipToState(EParsingState target_state);

    EParsingEvent x_ParseHeader();
    EParsingEvent x_ParsePos();
    EParsingEvent x_ParseIDs();
    EParsingEvent x_ParseAlts();
    EParsingEvent x_ParseQuality();
    EParsingEvent x_ParseFilters();
    EParsingEvent x_ParseInfo();
    EParsingEvent x_ParseGenotypeFormat();
    EParsingEvent x_ParseGenotype();

    const char* x_ParseGT();
};

} /* namespace vcf */

#endif /* !defined(VCF_SCANNER__HH) */
