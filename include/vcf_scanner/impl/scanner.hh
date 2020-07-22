// This header contains implementation details.
// It is not meant to be included directly.

#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <stdio.h>

#include "tokenizer.hh"

class VCF_scanner_impl
{
protected:
    VCF_parsing_event feed_impl(const char* buffer, ssize_t buffer_size)
    {
        tokenizer.set_new_buffer(buffer, buffer_size);

        if (state == parsing_genotypes) {
            return continue_parsing_genotype();
        }

        if (state <= parsing_pos) {
            if (state < parsing_chrom) {
                return continue_parsing_header();
            }
            if (state == parsing_chrom) {
                const VCF_parsing_event pe = parse_string(parsing_pos);
                if (pe != VCF_parsing_event::ok) {
                    return pe;
                }
                *output.loc.chrom = tokenizer.get_token();
            }
            return continue_parsing_pos();
        }

        for (; fields_to_skip > 0; --fields_to_skip) {
            if (!tokenizer.skip_token(tokenizer.find_newline_or_tab())) {
                return VCF_parsing_event::need_more_data;
            }
            if (tokenizer.at_eol()) {
                unsigned missing_field =
                        state - parsing_chrom + 1 - fields_to_skip;
                fields_to_skip = 0;
                return missing_mandatory_field_error(missing_field);
            }
        }

        switch (state) {
        case parsing_id:
            return continue_parsing_ids();
        case parsing_ref:
            // Parsing of 'ref' and 'alts' is requested by a single
            // method parse_alleles(). Once 'ref' has been parsed,
            // proceed to parsing 'alts'.
            {
                const VCF_parsing_event pe = parse_string(parsing_alt);
                if (pe != VCF_parsing_event::ok) {
                    return pe;
                }
                ref = tokenizer.get_token();
            }
            /* FALL THROUGH */
        case parsing_alt:
            return continue_parsing_alts();
        case parsing_quality:
            return continue_parsing_quality();
        case parsing_filter:
            return continue_parsing_filters();
        case parsing_info_field:
            return continue_parsing_info();
        case parsing_genotype_format:
            return continue_parsing_genotype_format();
        case skipping_to_next_line:
            if (!tokenizer.skip_token(tokenizer.find_newline())) {
                return VCF_parsing_event::need_more_data;
            }
            if (tokenizer.buffer_is_empty() && !tokenizer.at_eof()) {
                state = peeking_beyond_newline;
                return VCF_parsing_event::need_more_data;
            }
            /* FALL THROUGH */
        case peeking_beyond_newline:
            reset_state_for_next_data_line();
            return VCF_parsing_event::ok;
        }

        return VCF_parsing_event::error; // LCOV_EXCL_LINE
    }

    // Parses the CHROM and the POS fields.
    VCF_parsing_event parse_loc_impl(std::string* chrom, unsigned* pos)
    {
        // LCOV_EXCL_START
        if (state != parsing_chrom) {
            if (state < parsing_chrom) {
                assert(false && "VCF header must be parsed first");
                return VCF_parsing_event::error;
            }

            assert(false && "Must call clear_line() before parse_loc()");
            return VCF_parsing_event::error;
        }
        // LCOV_EXCL_STOP

        output.loc.chrom = chrom;
        *(output.loc.pos = pos) = 0;
        number_len = 0;

        const VCF_parsing_event pe = parse_string(parsing_pos);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }
        *output.loc.chrom = tokenizer.get_token();

        return continue_parsing_pos();
    }

    VCF_parsing_event parse_ids_impl()
    {
        next_list_index = 0;

        const VCF_parsing_event pe = skip_to_state(parsing_id);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }

        return continue_parsing_ids();
    }

    VCF_parsing_event parse_alleles_impl()
    {
        next_list_index = 0;

        VCF_parsing_event pe = skip_to_state(parsing_ref);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }

        // Parsing of 'ref' and 'alts' is requested by a single
        // method parse_alleles(). Once 'ref' has been parsed,
        // proceed to parsing 'alts'.
        pe = parse_string(parsing_alt);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }
        ref = tokenizer.get_token();

        return continue_parsing_alts();
    }

    VCF_parsing_event parse_quality_impl()
    {
        const VCF_parsing_event pe = skip_to_state(parsing_quality);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }

        return continue_parsing_quality();
    }

    VCF_parsing_event parse_filters_impl()
    {
        next_list_index = 0;

        VCF_parsing_event pe = skip_to_state(parsing_filter);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }

        return continue_parsing_filters();
    }

    VCF_parsing_event parse_info_impl()
    {
        info.clear();

        const VCF_parsing_event pe = skip_to_state(parsing_info_field);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }

        return continue_parsing_info();
    }

    VCF_parsing_event parse_genotype_format_impl()
    {
        genotype_key_positions.clear();

        const VCF_parsing_event pe = skip_to_state(parsing_genotype_format);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }

        return continue_parsing_genotype_format();
    }

    bool capture_gt_impl()
    {
        unsigned gt_index = genotype_key_positions.gt;
        if (gt_index == 0) {
            return false;
        }
        --gt_index;
        auto* gt_value = alloc_genotype_value(gt_index);
        gt_value->data_type = vcf_gt;
        gt_value->int_vector = &gt;
        return true;
    }

    VCF_parsing_event parse_genotype_impl()
    {
        // LCOV_EXCL_START
        if (state != parsing_genotypes) {
            assert(false &&
                    "parse_genotype_format must be called before "
                    "parse_genotype");
            return VCF_parsing_event::error;
        }
        // LCOV_EXCL_STOP

        if (current_genotype_field_index >= header.sample_ids.size()) {
            return data_line_error(
                    "The number of genotype fields exceeds "
                    "the number of samples");
        }

        current_genotype_value_index = 0;

        return continue_parsing_genotype();
    }

    VCF_parsing_event clear_line_impl()
    {
        if (!tokenizer.at_eof()) {
            if (state != peeking_beyond_newline) {
                if (state != end_of_data_line &&
                        !tokenizer.skip_token(tokenizer.find_newline())) {
                    state = skipping_to_next_line;
                    return VCF_parsing_event::need_more_data;
                }

                if (tokenizer.buffer_is_empty()) {
                    state = peeking_beyond_newline;
                    return VCF_parsing_event::need_more_data;
                }
            }
        }

        reset_state_for_next_data_line();

        return VCF_parsing_event::ok;
    }

protected:
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

    VCF_parsing_event header_error(const char* err_msg)
    {
        error_message = err_msg;
        return VCF_parsing_event::error;
    }

    VCF_parsing_event invalid_meta_info_line_error()
    {
        return header_error("Malformed meta-information line");
    }

    VCF_parsing_event invalid_header_line_error()
    {
        return header_error("Malformed VCF header line");
    }

    VCF_parsing_event data_line_error(const std::string& err_msg)
    {
        error_message = err_msg;
        return VCF_parsing_event::error;
    }

    static constexpr unsigned number_of_mandatory_columns = 8;

    static const char* get_header_line_column(unsigned field_index)
    {
        static const char* const columns[number_of_mandatory_columns + 2] = {
                "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                "FORMAT", "GENOTYPE"};

        return columns[field_index];
    }

    VCF_parsing_event missing_mandatory_field_error(unsigned field_index)
    {
        state = end_of_data_line;

        std::string err_msg = "Missing mandatory VCF field \"";
        err_msg += get_header_line_column(field_index);
        err_msg += '"';
        return data_line_error(err_msg);
    }

    std::vector<VCF_warning> warnings;
    std::string error_message;

    VCF_tokenizer tokenizer;

    VCF_header header;

    size_t next_list_index;
    unsigned number_len;

    union {
        struct {
            std::string* chrom;
            unsigned* pos;
        } loc;
    } output;
    std::vector<std::string> ids;
    std::string ref;
    std::vector<std::string> alts;
    bool alleles_parsed;
    std::string quality;
    std::vector<std::string> filters;
    std::vector<std::string> info;

    void reset_state_for_next_data_line()
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

    VCF_parsing_event parse_string(State target_state)
    {
        if (!tokenizer.prepare_token_or_accumulate(
                    tokenizer.find_newline_or_tab())) {
            return VCF_parsing_event::need_more_data;
        }
        if (tokenizer.at_eol()) {
            return missing_mandatory_field_error(target_state - parsing_chrom);
        }
        state = target_state;
        return VCF_parsing_event::ok;
    }

    VCF_parsing_event parse_string_list(State target_state,
            std::vector<std::string>& container,
            const std::array<bool, 256>& character_set)
    {
        do {
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_char_from_set(character_set))) {
                return VCF_parsing_event::need_more_data;
            }
            if (tokenizer.at_eol()) {
                return missing_mandatory_field_error(
                        target_state - parsing_chrom);
            }
            if (!tokenizer.token_is_dot()) {
                if (next_list_index < container.size()) {
                    container.at(next_list_index) = tokenizer.get_token();
                } else {
                    container.push_back(tokenizer.get_token());
                }
                ++next_list_index;
            }
        } while (tokenizer.get_terminator() != '\t');
        container.resize(next_list_index);
        state = target_state;
        return VCF_parsing_event::ok;
    }

    VCF_parsing_event skip_to_state(State target_state)
    {
        // LCOV_EXCL_START
        if (state < parsing_chrom) {
            assert(false && "VCF header must be parsed first");
            return VCF_parsing_event::error;
        }
        if (state > target_state) {
            assert(false && "clear_line() must call be called first");
            return VCF_parsing_event::error;
        }
        // LCOV_EXCL_STOP

        while (state < target_state) {
            if (!tokenizer.skip_token(tokenizer.find_newline_or_tab())) {
                fields_to_skip = target_state - state;
                state = target_state;
                return VCF_parsing_event::need_more_data;
            }
            if (tokenizer.at_eol()) {
                return missing_mandatory_field_error(state - parsing_chrom + 1);
            }
            ++state;
        }
        return VCF_parsing_event::ok;
    }

    VCF_parsing_event continue_parsing_header()
    {
        switch (state) {
        case parsing_fileformat:
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_newline())) {
                return VCF_parsing_event::need_more_data;
            }

            {
                VCF_string_view key, value;

                if (!tokenizer.get_key_value(&key, &value) ||
                        key != "##fileformat") {
                    return header_error(
                            "VCF files must start with '##fileformat'");
                }

                header.file_format_version = value;
            }

        parse_meta_info_key:
            state = parsing_metainfo_key;
            /* FALL THROUGH */

        case parsing_metainfo_key:
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_newline_or_tab_or_equals())) {
                return VCF_parsing_event::need_more_data;
            }

            if (tokenizer.at_eol()) {
                return invalid_meta_info_line_error();
            }

            if (tokenizer.get_terminator() == '\t') {
                if (tokenizer.get_token().substr(1) !=
                        get_header_line_column(0)) {
                    return invalid_meta_info_line_error();
                }
                header_line_column_ok = 1;
                goto parse_header_line;
            }

            // Found an equals sign - save the key and proceed
            // to parsing the value.
            {
                const VCF_string_view& key = tokenizer.get_token();
                if (key.length() < 3 || key[0] != '#' || key[1] != '#') {
                    return invalid_meta_info_line_error();
                }
                current_meta_info_key = key.substr(2);
            }

            state = parsing_metainfo_value;
            /* FALL THROUGH */

        case parsing_metainfo_value:
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_newline())) {
                return VCF_parsing_event::need_more_data;
            }

            if (tokenizer.get_terminator() == EOF) {
                return header_error(
                        "Unexpected end of file while parsing VCF file header");
            }

            header.add_meta_info(
                    std::move(current_meta_info_key), tokenizer.get_token());

            // Go back to parsing the next key.
            goto parse_meta_info_key;

        parse_header_line:
            state = parsing_header_line_columns;
            /* FALL THROUGH */

        case parsing_header_line_columns:
            do {
                if (!tokenizer.prepare_token_or_accumulate(
                            tokenizer.find_newline_or_tab())) {
                    return VCF_parsing_event::need_more_data;
                }

                if (tokenizer.get_token() !=
                        get_header_line_column(header_line_column_ok)) {
                    return invalid_header_line_error();
                }

                ++header_line_column_ok;

                if (tokenizer.at_eol()) {
                    if (header_line_column_ok < number_of_mandatory_columns) {
                        return invalid_header_line_error();
                    }
                    if (header_line_column_ok > number_of_mandatory_columns) {
                        // The FORMAT field is present,
                        // but there are no samples.
                        header.genotype_info_present = true;
                    }
                    goto end_of_header_line;
                }

                // The current token ends with a tab.
                // Parse the next header line column.
            } while (header_line_column_ok <= number_of_mandatory_columns);

            header.genotype_info_present = true;
            state = parsing_sample_ids;
            /* FALL THROUGH */

        case parsing_sample_ids:
            do {
                if (!tokenizer.prepare_token_or_accumulate(
                            tokenizer.find_newline_or_tab())) {
                    return VCF_parsing_event::need_more_data;
                }

                header.sample_ids.push_back(tokenizer.get_token());
            } while (tokenizer.get_terminator() == '\t');
        }

    end_of_header_line:
        if (tokenizer.buffer_is_empty() && !tokenizer.at_eof()) {
            state = peeking_beyond_newline;
            return VCF_parsing_event::need_more_data;
        }

        reset_state_for_next_data_line();

        return VCF_parsing_event::ok;
    }

    VCF_parsing_event continue_parsing_pos()
    {
        switch (tokenizer.parse_uint(output.loc.pos, &number_len)) {
        case VCF_tokenizer::end_of_buffer:
            return VCF_parsing_event::need_more_data;
        case VCF_tokenizer::integer_overflow:
            return data_line_error("Integer overflow in the POS column");
        default /* case VCF_tokenizer::end_of_number */:
            break;
        }

        if (number_len == 0) {
            return data_line_error("Missing an integer in the POS column");
        }

        if (tokenizer.get_terminator() != '\t') {
            return data_line_error("Invalid data line format");
        }

        state = parsing_id;
        return VCF_parsing_event::ok;
    }

    VCF_parsing_event continue_parsing_ids()
    {
        return parse_string_list(
                parsing_ref, ids, tokenizer.newline_or_tab_or_semicolon);
    }

    VCF_parsing_event continue_parsing_alts()
    {
        const VCF_parsing_event pe = parse_string_list(
                parsing_quality, alts, tokenizer.newline_or_tab_or_comma);

        if (pe == VCF_parsing_event::ok) {
            alleles_parsed = true;
        }

        return pe;
    }

    VCF_parsing_event continue_parsing_quality()
    {
        const VCF_parsing_event pe = parse_string(parsing_filter);
        if (pe != VCF_parsing_event::ok) {
            return pe;
        }
        if (!tokenizer.token_is_dot()) {
            quality = tokenizer.get_token();
        } else {
            quality.clear();
        }

        return VCF_parsing_event::ok;
    }

    VCF_parsing_event continue_parsing_filters()
    {
        return parse_string_list(parsing_info_field, filters,
                tokenizer.newline_or_tab_or_semicolon);
    }

    VCF_parsing_event continue_parsing_info()
    {
        do {
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_char_from_set(
                                tokenizer.newline_or_tab_or_semicolon))) {
                return VCF_parsing_event::need_more_data;
            }
            if (tokenizer.at_eol()) {
                state = end_of_data_line;
                return VCF_parsing_event::ok;
            }
            if (!tokenizer.token_is_dot()) {
                info.push_back(tokenizer.get_token());
            }
        } while (tokenizer.get_terminator() != '\t');

        state = parsing_genotype_format;

        return VCF_parsing_event::ok;
    }

    VCF_parsing_event continue_parsing_genotype_format()
    {
        do {
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_char_from_set(
                                tokenizer.newline_or_tab_or_colon))) {
                return VCF_parsing_event::need_more_data;
            }
            if (tokenizer.at_eol()) {
                state = end_of_data_line;
                if (header.sample_ids.empty()) {
                    return VCF_parsing_event::ok;
                }
                return data_line_error("No genotype information present");
            }
            std::string key = tokenizer.get_token();
            auto key_iter = format_keys.insert(key).first;
            if (key == "GT") {
                if (genotype_key_positions.number_of_positions != 0) {
                    // TODO Generate a warning: GT must be the first key.
                }
                genotype_key_positions.gt =
                        ++genotype_key_positions.number_of_positions;
            } else {
                genotype_key_positions.other_keys[key_iter->c_str()] =
                        ++genotype_key_positions.number_of_positions;
            }
        } while (tokenizer.get_terminator() != '\t');

        reset_genotype_values();
        state = parsing_genotypes;
        return VCF_parsing_event::ok;
    }

    VCF_parsing_event continue_parsing_genotype()
    {
        Genotype_value* value =
                genotype_values.data() + current_genotype_value_index;

        do {
            if (value->flag == nullptr) {
                if (!tokenizer.skip_token(
                            tokenizer.find_newline_or_tab_or_colon())) {
                    return VCF_parsing_event::need_more_data;
                }
                if (tokenizer.at_eol()) {
                    state = end_of_data_line;
                    return VCF_parsing_event::ok;
                }
            } else {
                if (!tokenizer.prepare_token_or_accumulate(
                            tokenizer.find_char_from_set(
                                    tokenizer.newline_or_tab_or_colon))) {
                    return VCF_parsing_event::need_more_data;
                }

                if (tokenizer.at_eol()) {
                    state = end_of_data_line;
                }

                switch (value->data_type) {
                // TODO case vcf_integer:
                // TODO case vcf_float:
                // TODO case vcf_character:
                // TODO case vcf_string:
                case vcf_gt:
                    // Hi
                    {
                        const char* err_msg = parse_gt();
                        if (err_msg != nullptr) {
                            return data_line_error(err_msg);
                        }
                    }
                    break;
                default /* vcf_flag - impossible type for genotype info */:
                    break;
                }

                if (tokenizer.at_eol()) {
                    return VCF_parsing_event::ok;
                }
            }

            if (tokenizer.get_terminator() == '\t') {
                ++current_genotype_field_index;
                return VCF_parsing_event::ok;
            }

            ++value;
        } while (++current_genotype_value_index <
                genotype_key_positions.number_of_positions);

        return data_line_error("Too many genotype info fields");
    }

    const char* parse_gt()
    {
        gt.clear();

        const std::string& token = tokenizer.get_token();

        size_t len = token.length();

        if (len == 0) {
            return "Empty GT value";
        }

        const char* ptr = token.data();
        unsigned digit, allele;

        for (;; ++ptr, --len) {
            if (*ptr == '.') {
                gt.push_back(-1);
                ++ptr;
                --len;
            } else {
                if ((allele = (unsigned) *ptr - '0') > 9) {
                    break;
                }

                while (--len > 0 && (digit = (unsigned) *++ptr - '0') <= 9) {
                    if (allele > (UINT_MAX / 10) ||
                            (allele == (UINT_MAX / 10) &&
                                    digit > UINT_MAX % 10)) {
                        // TODO ERR+="in genotype info for the sample #" +
                        //           str(current_genotype_value_index)
                        return "Integer overflow in allele index";
                    }

                    allele = allele * 10 + digit;
                }

                gt.push_back((int) allele);

                if (alleles_parsed && allele > alts.size()) {
                    return "Allele index exceeds the number of alleles";
                }
            }
            if (len == 0) {
                return nullptr;
            }
            switch (*ptr) {
            case '/':
                phased_gt = false;
                continue;
            case '|':
                phased_gt = true;
                continue;
            }
            break;
        }
        return "Invalid character in GT value";
    }
};
