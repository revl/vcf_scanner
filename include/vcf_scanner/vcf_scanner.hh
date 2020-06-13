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

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <cassert>
#include <stdio.h>

#include "impl/vcf_tokenizer.hh"

class VCF_scanner;

// Metadata extracted from the VCF header.
class VCF_header
{
public:
    typedef std::vector<std::string> Meta_info_lines;
    typedef std::map<std::string, Meta_info_lines> Meta_info;
    typedef std::vector<std::string> Sample_IDs;

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

    const Sample_IDs& get_sample_ids() const
    {
        return sample_ids;
    }

private:
    std::string file_format_version;
    Meta_info meta_info;
    bool genotype_info_present = false;
    Sample_IDs sample_ids;

    friend class VCF_scanner;
};

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
class VCF_scanner
{
public:
    VCF_scanner() = default;

    enum Parsing_event {
        need_more_data, // The parser needs a new input buffer to
                        // continue parsing. See feed().

        ok, // The current token (the VCF header or a data field) has
            // been successfully parsed. Header meta-information or
            // the data field value is now available for retrieval.

        ok_with_warnings, // The token has been successfully parsed, but
                          // parser encountered issues during parsing. Use
                          // GetWarnings() to retrieve the warning messages.

        error // A parsing error has occurred. Use get_error() to get
              // the error message.
              //
              // If the error happened while parsing the VCF header,
              // this parser instance can no longer be used.
              //
              // If the error happened while parsing a data line, there
              // is an option to ignore it and skip to the next line
              // by calling clear_line().
    };

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
    // returns eOK when the entire token has been parsed.
    Parsing_event feed(const char* buffer, ssize_t buffer_size);

    // Returns the current line number in the input VCF file before
    // parsing the next token. The line number will increase after the
    // last token on the current line has been parsed. The returned value
    // is one-based.
    unsigned get_line_number() const
    {
        return tokenizer.get_line_number();
    }

    // TODO FIXME Not used yet.
    struct Warning {
        unsigned line_number;
        std::string warning_message;
    };

    // TODO FIXME Not used yet.
    std::vector<Warning> get_warnings() const
    {
        return warnings;
    }

    // Returns the description of the error that caused parsing to fail. Use
    // get_line_number() to get the line number where the error has occurred.
    std::string get_error() const
    {
        return error_message;
    }

    // Returns the VCF header, which becomes available
    // once the last of the initial series of calls to Feed
    // returns eOK.
    const VCF_header& get_header() const
    {
        return header;
    }

    // Returns true if the entire input stream has been
    // successfully parsed. The method returns false if the
    // VCF file has at least one more data line to parse.
    bool at_eof() const
    {
        return tokenizer.at_eof();
    }

    // Parses the CHROM and the POS fields.
    Parsing_event parse_loc();
    // Returns the CHROM field parsed by parse_loc.
    const std::string& get_chrom() const
    {
        return chrom;
    }
    // Returns the POS field parsed by parse_loc.
    unsigned get_pos() const
    {
        return pos;
    }

    // Parses the ID field.
    Parsing_event parse_ids();
    // Returns the IDs parsed by parse_ids().
    const std::vector<std::string>& get_ids() const
    {
        return ids;
    }

    // Parses the REF and the ALT fields.
    Parsing_event parse_alleles();
    // Returns the REF field parsed by parse_alleles.
    const std::string& get_ref() const
    {
        return ref;
    }
    // Returns the ALT field parsed by parse_alleles().
    const std::vector<std::string>& get_alts() const
    {
        return alts;
    }

    // Parses the QUAL field.
    Parsing_event parse_quality();
    // Returns the QUAL field parsed by parse_quality().
    std::string get_quality() const
    {
        return quality;
    }

    // Parses the FILTER field.
    Parsing_event parse_filters();
    // Returns the FILTER field parsed by parse_filters().
    // The word "PASS" is returned when the current record passed
    // all filters.
    const std::vector<std::string>& get_filters() const
    {
        return filters;
    }

    // Parses the INFO key-value pairs.
    Parsing_event parse_info();
    // Returns the INFO field parsed by parse_info().
    const std::vector<std::string>& get_info() const
    {
        return info;
    }

    // Parses the genotype format keys.
    Parsing_event parse_genotype_format();

    // Returns the FORMAT field parsed by parse_genotype_format().
    // TODO const TGenotypeFormat& get_genotype_format() const;

    // Enables parsing of GT values in the parse_genotype method.
    // Returns false and does nothing if the GT key was not
    // specified in the FORMAT field.
    bool capture_gt();

    // TODO bool capture_string(const char* key, std::string* value);
    // TODO bool capture_strings(const char* key, std::vector<std::string>*
    // values);
    // TODO bool capture_int(const char* key, int* value);
    // TODO bool capture_ints(const char* key, std::vector<int>* values);

    // Parses genotype fields one by one.
    Parsing_event parse_genotype();

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
    Parsing_event clear_line();

private:
    enum State {
        parsing_fileformat,
        parsing_metainfo_key,
        parsing_metainfo_value,
        parsing_header_line_columns,
        parsing_sample_ids,
        parsing_chrom,
        parsing_pos,
        parsing_id,
        parsing_ref,
        parsing_alt,
        parsing_quality,
        parsing_filter,
        parsing_info_field,
        parsing_genotype_format,
        parsing_genotypes,
        end_of_data_line,
        skipping_to_next_line,
        peeking_beyond_newline
    };
    int state = parsing_fileformat;

    int fields_to_skip = 0;

    std::string current_meta_info_key;

    unsigned header_line_column_ok;

    Parsing_event header_error(const char* err_msg)
    {
        error_message = err_msg;
        return error;
    }

    Parsing_event invalid_meta_info_line_error()
    {
        return header_error("Malformed meta-information line");
    }

    Parsing_event invalid_header_line_error()
    {
        return header_error("Malformed VCF header line");
    }

    Parsing_event data_line_error(const std::string& err_msg)
    {
        error_message = err_msg;
        return error;
    }

    Parsing_event missing_mandatory_field_error(const char* field_name)
    {
        state = end_of_data_line;

        std::string err_msg = "Missing mandatory VCF field \"";
        err_msg += field_name;
        err_msg += '"';
        return data_line_error(err_msg);
    }

    std::vector<Warning> warnings;
    std::string error_message;

    VCF_tokenizer tokenizer;

    VCF_header header;

    size_t next_list_index;
    unsigned number_len;

    std::string chrom;
    unsigned pos;
    std::vector<std::string> ids;
    std::string ref;
    std::vector<std::string> alts;
    bool alleles_parsed;
    std::string quality;
    std::vector<std::string> filters;
    std::vector<std::string> info;

    void reset_data_line()
    {
        state = parsing_chrom;
        alleles_parsed = false;
    }

    // TODO Implement the INFO and FORMAT type definitions in the header.
    std::set<std::string> format_keys;

    struct Strcmp {
        bool operator()(const char* left, const char* right) const
        {
            return strcmp(left, right) < 0;
        }
    };

    // Positions of the reserved genotype keys in the FORMAT field.
    struct Genotype_key_positions {
        unsigned number_of_positions;
        unsigned gt;
        std::map<const char*, unsigned, Strcmp> other_keys;

        void clear()
        {
            gt = number_of_positions = 0;
            other_keys.clear();
        }
    } genotype_key_positions;

    enum Data_type {
        vcf_integer,
        vcf_float,
        vcf_flag,
        vcf_character,
        vcf_string,
        vcf_gt
    };

    enum Number_of_values {
        scalar,
        one_per_alt,
        one_per_allele,
        one_per_genotype,
        unbound,
        exact_number
    };

    struct Genotype_value {
        Data_type data_type;
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

    unsigned current_genotype_field_index;

    unsigned current_genotype_value_index;
    std::vector<Genotype_value> genotype_values;

    std::vector<int> gt;
    bool phased_gt;

    void reset_genotype_values()
    {
        memset(genotype_values.data(), 0,
                (char*) &*genotype_values.end() -
                        (char*) genotype_values.data());
        current_genotype_field_index = 0;
        number_len = 0;
    }

    Genotype_value* alloc_genotype_value(unsigned index)
    {
        if (genotype_values.size() <= index) {
            size_t old_size = genotype_values.size();
            genotype_values.resize(index + 1);
            char* new_struct_ptr = (char*) (genotype_values.data() + old_size);
            memset(new_struct_ptr, 0,
                    (char*) &*genotype_values.end() - (char*) new_struct_ptr);
        }
        return genotype_values.data() + index;
    }

    Parsing_event parse_string(State target_state);
    Parsing_event parse_string_list(State target_state,
            std::vector<std::string>& container, const bool* character_set);
    Parsing_event skip_to_state(State target_state);

    Parsing_event continue_parsing_header();
    Parsing_event continue_parsing_pos();
    Parsing_event continue_parsing_ids();
    Parsing_event continue_parsing_alts();
    Parsing_event continue_parsing_quality();
    Parsing_event continue_parsing_filters();
    Parsing_event continue_parsing_info();
    Parsing_event continue_parsing_genotype_format();
    Parsing_event continue_parsing_genotype();

    const char* parse_gt();
};

#endif /* !defined(VCF_SCANNER__HH) */
