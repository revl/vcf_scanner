/*
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
 */

#ifndef VCF_SCANNER__HH
#define VCF_SCANNER__HH

// These constants returned by most VCF_scanner methods to indicate the result
// of a parsing operation.
enum class VCF_parsing_event {
    need_more_data, // The parser needs a new input buffer to
                    // continue parsing. See VCF_scanner::feed().

    ok, // The current token (the VCF header or a data field) has
        // been successfully parsed. Header meta-information or
        // the data field value is now available for retrieval.

    ok_with_warnings, // The token has been successfully parsed, but
                      // parser encountered issues during parsing. Use
                      // VCF_scanner::get_warnings() to retrieve the warning
                      // messages.

    error // A parsing error has occurred. Use VCF_scanner::get_error()
          // to get the error message.
          //
          // If the error happened while parsing the VCF header,
          // this parser instance can no longer be used.
          //
          // If the error happened while parsing a data line, there
          // is an option to ignore it and skip to the next line
          // by calling VCF_scanner::clear_line().
};

#include "impl/header.hh"

// Metadata extracted from the VCF header.
class VCF_header final : public VCF_header_impl
{
public:
    typedef std::map<std::string, std::vector<std::string>> Meta_info;

    // Returns the VCF version of the current input file.
    const std::string& get_file_format_version() const
    {
        return file_format_version;
    }

    const Meta_info& get_meta_info() const
    {
        return meta_info;
    }

    bool has_genotype_info() const
    {
        return genotype_info_present;
    }

    const std::vector<std::string>& get_sample_ids() const
    {
        return sample_ids;
    }
};

// TODO FIXME Not used yet.
struct VCF_warning {
    unsigned line_number;
    std::string warning_message;
};

#include "impl/scanner.hh"

// Parser of VCF (Variant Call Format) files.
//
// This class parses and returns the header first, and then it parses
// data lines one by one.
//
// All header information is kept by the parser in its member variables.
// That includes sample IDs from the header line.
//
// The fields of the data lines, however, are never stored internally.
// This class does not store the parsed data. Each data line field is
// discarded as soon as the caller proceeds to parsing the next field.
//
// The parser does not have a stream reading loop inside. It relies on
// the client code to provide the input data. As a result, it never blocks
// on I/O operations. This allows for using a separate thread to read data
// into a new buffer while the main thread is parsing a previously read
// buffer. Alternatively to reading the input file into memory, the whole
// file can be memory-mapped.
//
// Before parsing begins, or when a parsing function returns 'need_more_data',
// a new buffer with input data must be supplied to the parser by calling
// feed(). The buffer must not be freed or overwritten until 'need_more_data'
// is received again or the client code chooses not to continue parsing.
class VCF_scanner final : public VCF_scanner_impl
{
public:
    VCF_scanner() = default;

    // Supplies a chunk of input data to this parser either
    // when the parser has just been created and is in the process
    // of parsing the VCF header or when a previously called method
    // returned 'need_more_data'.
    //
    // The parser stores the buffer pointer internally. Freeing the
    // buffer before the parser returns 'need_more_data' from any of its
    // methods will cause a segmentation fault. A buffer of zero size
    // is treated as an EOF condition.
    //
    // The method resumes parsing of the previously requested token and
    // returns 'ok' when the entire token has been parsed.
    VCF_parsing_event feed(const char* buffer, ssize_t buffer_size)
    {
        return feed_impl(buffer, buffer_size);
    }

    // Returns the current line number in the input VCF file before parsing the
    // next token. The line number will increase after the last token on the
    // current line has been parsed. The returned value is one-based.
    unsigned get_line_number() const
    {
        return tokenizer.get_line_number();
    }

    // TODO FIXME Not used yet.
    std::vector<VCF_warning> get_warnings() const
    {
        return warnings;
    }

    // Returns the description of the error that caused parsing to fail. Use
    // get_line_number() to get the line number where the error has occurred.
    std::string get_error() const
    {
        return error_message;
    }

    // Returns the VCF header, which becomes available once the last of the
    // initial series of calls to Feed returns 'ok'.
    const VCF_header& get_header() const
    {
        return header;
    }

    // Returns true if the entire input stream has been successfully parsed.
    // The method returns false if the VCF file has at least one more data line
    // to parse.
    bool at_eof() const
    {
        return tokenizer.at_eof();
    }

    // Parses the CHROM and the POS fields and stores the parsed values into
    // the variables pointed to by 'chrom' and 'pos'.  The lifespan of those
    // variables must exceed this 'parse_loc()' call as well as all 'feed()'
    // calls that may be required to finish parsing the CHROM and POS fields.
    VCF_parsing_event parse_loc(std::string* chrom, unsigned* pos)
    {
        return parse_loc_impl(chrom, pos);
    }

    // Parses the ID field into the 'ids' array.  The lifespan of the array
    // must exceed this 'parse_ids()' call as well as all 'feed()' calls that
    // may be required to finish parsing the ID field.
    VCF_parsing_event parse_ids(std::vector<std::string>* ids)
    {
        return parse_ids_impl(ids);
    }

    // Parses the REF and the ALT fields.  The lifespan of 'ref' and 'alts'
    // must exceed this 'parse_alleles()' call as well as all 'feed()' calls
    // that may be required to finish parsing the REF and ALT fields.
    VCF_parsing_event parse_alleles(
            std::string* ref, std::vector<std::string>* alts)
    {
        return parse_alleles_impl(ref, alts);
    }

    // Parses the QUAL field.
    VCF_parsing_event parse_quality()
    {
        return parse_quality_impl();
    }
    // Returns whether the QUAL field contains the MISSING value ('.').
    bool quality_is_missing() const
    {
        return quality.empty();
    }
    // Returns the QUAL value. This method must not be called if the value
    // is missing.
    float get_quality() const
    {
        assert(!quality_is_missing());

        return std::stof(quality);
    }
    // Returns the original string representation of the QUAL value as it
    // appears in the VCF file or an empty string if the value is missing.
    const std::string& get_quality_as_string() const
    {
        return quality;
    }

    // Parses the FILTER field.
    VCF_parsing_event parse_filters()
    {
        return parse_filters_impl();
    }

    // Returns the FILTER field parsed by parse_filters().
    // The word "PASS" is returned when the current record passed
    // all filters.
    const std::vector<std::string>& get_filters() const
    {
        return filters;
    }

    // Parses the INFO key-value pairs.
    VCF_parsing_event parse_info()
    {
        return parse_info_impl();
    }

    // Returns the INFO field parsed by parse_info().
    const std::vector<std::string>& get_info() const
    {
        return info;
    }

    // Parses the genotype format keys.
    VCF_parsing_event parse_genotype_format()
    {
        return parse_genotype_format_impl();
    }

    // Returns the FORMAT field parsed by parse_genotype_format().
    // TODO const Genotype_format& get_genotype_format() const;

    // Enables parsing of GT values in the parse_genotype method.
    // Returns false and does nothing if the GT key was not
    // specified in the FORMAT field.
    bool capture_gt()
    {
        return capture_gt_impl();
    }

    // TODO bool capture_string(const char* key, std::string* value);
    // TODO bool capture_strings(const char* key, std::vector<std::string>*
    // values);
    // TODO bool capture_int(const char* key, int* value);
    // TODO bool capture_ints(const char* key, std::vector<int>* values);

    // Parses genotype fields one by one.
    VCF_parsing_event parse_genotype()
    {
        return parse_genotype_impl();
    }

    // Returns the GT values parsed by parse_genotype in
    // its previous iteration.
    const std::vector<int>& get_gt() const
    {
        return gt;
    }
    // Returns true if the parsed sample was phased.
    bool is_phased_gt() const
    {
        return phased_gt;
    }

    // Returns true if at least one more genotype field is available.
    // The caller has an option to either use this method or count the
    // retrieved genotypes to determine when the last genotype on the
    // current data line has been parsed.
    bool genotype_available() const
    {
        return tokenizer.get_terminator() == '\t';
    }

    // Skips the remaining part of the current data line.
    // This operation may require more data to be read. The client
    // code must call this method after parsing each line even if
    // the line has been parsed to the end, because this method also
    // determines whether end-of-file has been reached.
    VCF_parsing_event clear_line()
    {
        return clear_line_impl();
    }
};

#endif /* !defined(VCF_SCANNER__HH) */
