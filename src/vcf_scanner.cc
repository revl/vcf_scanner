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

#include <vcf_scanner/vcf_scanner.hh>

#define NUMBER_OF_MANDATORY_COLUMNS 8

static const char* const header_line_columns[NUMBER_OF_MANDATORY_COLUMNS + 2] =
        {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                "GENOTYPE"};

#define PARSE_STRING(target_state)                                             \
    if (!tokenizer.prepare_token_or_accumulate(                                \
                tokenizer.find_newline_or_tab())) {                            \
        return need_more_data;                                                 \
    }                                                                          \
    if (tokenizer.token_is_last()) {                                           \
        return missing_mandatory_field_error(                                  \
                header_line_columns[target_state - parsing_chrom]);            \
    }                                                                          \
    state = target_state;

#define PARSE_STRING_LIST(target_state, container, character_set)              \
    do {                                                                       \
        if (!tokenizer.prepare_token_or_accumulate(                            \
                    tokenizer.find_char_from_set(character_set))) {            \
            return need_more_data;                                             \
        }                                                                      \
        if (tokenizer.token_is_last()) {                                       \
            return missing_mandatory_field_error(                              \
                    header_line_columns[target_state - parsing_chrom]);        \
        }                                                                      \
        if (!tokenizer.token_is_dot()) {                                       \
            if (next_list_index < container.size()) {                          \
                container.at(next_list_index) = tokenizer.get_token();         \
            } else {                                                           \
                container.push_back(tokenizer.get_token());                    \
            }                                                                  \
            ++next_list_index;                                                 \
        }                                                                      \
    } while (tokenizer.get_terminator() != '\t');                              \
    container.resize(next_list_index);                                         \
    state = target_state;

VCF_scanner::Parsing_event VCF_scanner::skip_to_state(
        VCF_scanner::State target_state)
{
    // LCOV_EXCL_START
    if (state < parsing_chrom) {
        assert(false && "VCF header must be parsed first");
        return error;
    }
    if (state > target_state) {
        assert(false && "clear_line() must call be called first");
        return error;
    }
    // LCOV_EXCL_STOP

    while (state < target_state) {
        if (!tokenizer.skip_token(tokenizer.find_newline_or_tab())) {
            fields_to_skip = target_state - state;
            state = target_state;
            return need_more_data;
        }
        if (tokenizer.token_is_last()) {
            return missing_mandatory_field_error(
                    header_line_columns[state - parsing_chrom + 1]);
        }
        ++state;
    }
    return ok;
}

#define PARSE_CHROM()                                                          \
    PARSE_STRING(parsing_pos);                                                 \
    chrom = tokenizer.get_token();

#define PARSE_REF()                                                            \
    PARSE_STRING(parsing_alt);                                                 \
    ref = tokenizer.get_token();

VCF_scanner::Parsing_event VCF_scanner::feed(
        const char* buffer, ssize_t buffer_size)
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
            PARSE_CHROM();
        }
        return continue_parsing_pos();
    }

    for (; fields_to_skip > 0; --fields_to_skip) {
        if (!tokenizer.skip_token(tokenizer.find_newline_or_tab())) {
            return need_more_data;
        }
        if (tokenizer.token_is_last()) {
            unsigned missing_field = state - parsing_chrom + 1 - fields_to_skip;
            fields_to_skip = 0;
            return missing_mandatory_field_error(
                    header_line_columns[missing_field]);
        }
    }

    switch (state) {
    case parsing_id:
        return continue_parsing_ids();
    case parsing_ref:
        PARSE_REF();
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
            return need_more_data;
        }
        if (tokenizer.buffer_is_empty() && !tokenizer.at_eof()) {
            state = peeking_beyond_newline;
            return need_more_data;
        }
        /* FALL THROUGH */
    case peeking_beyond_newline:
        reset_data_line();
        return ok;
    }

    return error; // LCOV_EXCL_LINE
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_header()
{
    switch (state) {
    case parsing_fileformat:
        if (!tokenizer.prepare_token_or_accumulate(tokenizer.find_newline())) {
            return need_more_data;
        }

        {
            VCF_string_view key, value;

            if (!tokenizer.get_key_value(&key, &value) ||
                    key != "##fileformat") {
                return header_error("VCF files must start with '##fileformat'");
            }

            header.file_format_version = value;
        }

    parse_meta_info_key:
        state = parsing_metainfo_key;
        /* FALL THROUGH */

    case parsing_metainfo_key:
        if (!tokenizer.prepare_token_or_accumulate(
                    tokenizer.find_newline_or_tab_or_equals())) {
            return need_more_data;
        }

        if (tokenizer.token_is_last()) {
            return invalid_meta_info_line_error();
        }

        if (tokenizer.get_terminator() == '\t') {
            if (tokenizer.get_token().substr(1) != header_line_columns[0]) {
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
        if (!tokenizer.prepare_token_or_accumulate(tokenizer.find_newline())) {
            return need_more_data;
        }

        if (tokenizer.get_terminator() == EOF) {
            return header_error(
                    "Unexpected end of file while parsing VCF file header");
        }

        header.meta_info[current_meta_info_key].push_back(
                tokenizer.get_token());

        // Go back to parsing the next key.
        goto parse_meta_info_key;

    parse_header_line:
        state = parsing_header_line_columns;
        /* FALL THROUGH */

    case parsing_header_line_columns:
        do {
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_newline_or_tab())) {
                return need_more_data;
            }

            if (tokenizer.get_token() !=
                    header_line_columns[header_line_column_ok]) {
                return invalid_header_line_error();
            }

            ++header_line_column_ok;

            if (tokenizer.token_is_last()) {
                if (header_line_column_ok < NUMBER_OF_MANDATORY_COLUMNS) {
                    return invalid_header_line_error();
                }
                if (header_line_column_ok > NUMBER_OF_MANDATORY_COLUMNS) {
                    // The FORMAT field is present,
                    // but there are no samples.
                    header.genotype_info_present = true;
                }
                goto end_of_header_line;
            }

            // The current token ends with a tab.
            // Parse the next header line column.
        } while (header_line_column_ok <= NUMBER_OF_MANDATORY_COLUMNS);

        header.genotype_info_present = true;
        state = parsing_sample_ids;
        /* FALL THROUGH */

    case parsing_sample_ids:
        do {
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_newline_or_tab())) {
                return need_more_data;
            }

            header.sample_ids.push_back(tokenizer.get_token());
        } while (tokenizer.get_terminator() == '\t');
    }

end_of_header_line:
    if (tokenizer.buffer_is_empty() && !tokenizer.at_eof()) {
        state = peeking_beyond_newline;
        return need_more_data;
    }

    reset_data_line();

    return ok;
}

VCF_scanner::Parsing_event VCF_scanner::parse_loc()
{
    // LCOV_EXCL_START
    if (state != parsing_chrom) {
        if (state < parsing_chrom) {
            assert(false && "VCF header must be parsed first");
            return error;
        }

        assert(false && "Must call clear_line() before parse_loc()");
        return error;
    }
    // LCOV_EXCL_STOP

    pos = 0;
    number_len = 0;

    PARSE_CHROM();

    return continue_parsing_pos();
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_pos()
{
    switch (tokenizer.parse_uint(&pos, &number_len)) {
    case VCF_tokenizer::end_of_buffer:
        return need_more_data;
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
    return ok;
}

VCF_scanner::Parsing_event VCF_scanner::parse_ids()
{
    next_list_index = 0;

    const Parsing_event pe = skip_to_state(parsing_id);
    if (pe != ok) {
        return pe;
    }

    return continue_parsing_ids();
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_ids()
{
    PARSE_STRING_LIST(parsing_ref, ids, tokenizer.newline_or_tab_or_semicolon);
    return ok;
}

VCF_scanner::Parsing_event VCF_scanner::parse_alleles()
{
    next_list_index = 0;

    const Parsing_event pe = skip_to_state(parsing_ref);
    if (pe != ok) {
        return pe;
    }

    PARSE_REF();

    return continue_parsing_alts();
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_alts()
{
    PARSE_STRING_LIST(parsing_quality, alts, tokenizer.newline_or_tab_or_comma);

    alleles_parsed = true;

    return ok;
}

VCF_scanner::Parsing_event VCF_scanner::parse_quality()
{
    const Parsing_event pe = skip_to_state(parsing_quality);
    if (pe != ok) {
        return pe;
    }

    return continue_parsing_quality();
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_quality()
{
    PARSE_STRING(parsing_filter);
    if (!tokenizer.token_is_dot()) {
        quality = tokenizer.get_token();
    } else {
        quality.clear();
    }

    return ok;
}

VCF_scanner::Parsing_event VCF_scanner::parse_filters()
{
    next_list_index = 0;

    Parsing_event pe = skip_to_state(parsing_filter);
    if (pe != ok) {
        return pe;
    }

    return continue_parsing_filters();
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_filters()
{
    PARSE_STRING_LIST(
            parsing_info_field, filters, tokenizer.newline_or_tab_or_semicolon);

    return ok;
}

VCF_scanner::Parsing_event VCF_scanner::parse_info()
{
    info.clear();

    const Parsing_event pe = skip_to_state(parsing_info_field);
    if (pe != ok) {
        return pe;
    }

    return continue_parsing_info();
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_info()
{
    do {
        if (!tokenizer.prepare_token_or_accumulate(tokenizer.find_char_from_set(
                    tokenizer.newline_or_tab_or_semicolon))) {
            return need_more_data;
        }
        if (tokenizer.token_is_last()) {
            state = end_of_data_line;
            return ok;
        }
        if (!tokenizer.token_is_dot()) {
            info.push_back(tokenizer.get_token());
        }
    } while (tokenizer.get_terminator() != '\t');

    state = parsing_genotype_format;

    return ok;
}

VCF_scanner::Parsing_event VCF_scanner::parse_genotype_format()
{
    genotype_key_positions.clear();

    const Parsing_event pe = skip_to_state(parsing_genotype_format);
    if (pe != ok) {
        return pe;
    }

    return continue_parsing_genotype_format();
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_genotype_format()
{
    do {
        if (!tokenizer.prepare_token_or_accumulate(tokenizer.find_char_from_set(
                    tokenizer.newline_or_tab_or_colon))) {
            return need_more_data;
        }
        if (tokenizer.token_is_last()) {
            state = end_of_data_line;
            if (header.sample_ids.empty()) {
                return ok;
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
    return ok;
}

bool VCF_scanner::capture_gt()
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

const char* VCF_scanner::parse_gt()
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
                        (allele == (UINT_MAX / 10) && digit > UINT_MAX % 10)) {
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

VCF_scanner::Parsing_event VCF_scanner::parse_genotype()
{
    // LCOV_EXCL_START
    if (state != parsing_genotypes) {
        assert(false &&
                "parse_genotype_format must be called before parse_genotype");
        return error;
    }
    // LCOV_EXCL_STOP

    if (current_genotype_field_index >= header.sample_ids.size()) {
        return data_line_error(
                "The number of genotype fields exceeds the number of samples");
    }

    current_genotype_value_index = 0;

    return continue_parsing_genotype();
}

VCF_scanner::Parsing_event VCF_scanner::continue_parsing_genotype()
{
    Genotype_value* value =
            genotype_values.data() + current_genotype_value_index;

    do {
        if (value->flag == nullptr) {
            if (!tokenizer.skip_token(
                        tokenizer.find_newline_or_tab_or_colon())) {
                return need_more_data;
            }
            if (tokenizer.token_is_last()) {
                state = end_of_data_line;
                return ok;
            }
        } else {
            if (!tokenizer.prepare_token_or_accumulate(
                        tokenizer.find_char_from_set(
                                tokenizer.newline_or_tab_or_colon))) {
                return need_more_data;
            }

            if (tokenizer.token_is_last()) {
                state = end_of_data_line;
            }

            switch (value->data_type) {
            // TODO case vcf_integer:
            // TODO case vcf_float:
            // TODO case vcf_character:
            // TODO case vcf_string:
            case vcf_gt: {
                const char* err_msg = parse_gt();
                if (err_msg != nullptr) {
                    return data_line_error(err_msg);
                }
            } break;
            default /* vcf_flag - impossible type for genotype info */:
                break;
            }

            if (tokenizer.token_is_last()) {
                return ok;
            }
        }

        if (tokenizer.get_terminator() == '\t') {
            ++current_genotype_field_index;
            return ok;
        }

        ++value;
    } while (++current_genotype_value_index <
            genotype_key_positions.number_of_positions);

    return data_line_error("Too many genotype info fields");
}

VCF_scanner::Parsing_event VCF_scanner::clear_line()
{
    if (!tokenizer.at_eof()) {
        if (state != peeking_beyond_newline) {
            if (state != end_of_data_line &&
                    !tokenizer.skip_token(tokenizer.find_newline())) {
                state = skipping_to_next_line;
                return need_more_data;
            }

            if (tokenizer.buffer_is_empty()) {
                state = peeking_beyond_newline;
                return need_more_data;
            }
        }
    }

    reset_data_line();

    return ok;
}
