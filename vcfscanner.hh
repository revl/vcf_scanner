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

#pragma once

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
