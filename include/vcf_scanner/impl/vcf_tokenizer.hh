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

#ifndef VCF_TOKENIZER__HH
#define VCF_TOKENIZER__HH

#include <string>
#include <iostream>
#include <climits>
#include <cstring>

typedef std::string VCF_string_view;

// Tokenizer for VCF streams. This class is not meant to be used directly.
class VCF_tokenizer
{
public:
    void set_new_buffer(const char* buffer, size_t buffer_size)
    {
        current_ptr = buffer;

        eof_reached = (remaining_size = buffer_size) == 0;
    }

    bool buffer_is_empty() const
    {
        return remaining_size == 0;
    }

    bool at_eof() const
    {
        return eof_reached;
    }

    const char* find_newline() const
    {
        return (const char*) memchr(current_ptr, '\n', remaining_size);
    }

private:
    static const char* find_char_from_set(
            const char* buffer, size_t buffer_size, const bool* character_set)
    {
        for (; buffer_size > 0; ++buffer, --buffer_size) {
            if (character_set[(unsigned char) *buffer]) {
                return buffer;
            }
        }

        return nullptr;
    }

public:
    const char* find_char_from_set(const bool* character_set) const
    {
        return find_char_from_set(current_ptr, remaining_size, character_set);
    }

    const char* find_newline_or_tab() const
    {
        return find_char_from_set(current_ptr, remaining_size, newline_or_tab);
    }

    const char* find_newline_or_tab_or_equals() const
    {
        return find_char_from_set(
                current_ptr, remaining_size, newline_or_tab_or_equals);
    }

    const char* find_newline_or_tab_or_semicolon() const
    {
        return find_char_from_set(
                current_ptr, remaining_size, newline_or_tab_or_semicolon);
    }

    const char* find_newline_or_tab_or_comma() const
    {
        return find_char_from_set(
                current_ptr, remaining_size, newline_or_tab_or_comma);
    }

    const char* find_newline_or_tab_or_colon() const
    {
        return find_char_from_set(
                current_ptr, remaining_size, newline_or_tab_or_colon);
    }

private:
    void set_terminator(int term)
    {
        terminator = term;
    }

    void set_terminator_and_inc_line_num_if_newline(int term)
    {
        set_terminator(term);

        if (term == '\n') {
            ++line_number;
        }
    }

    void advance_by(size_t number_of_bytes)
    {
        current_ptr += number_of_bytes;
        remaining_size -= number_of_bytes;
    }

public:
    enum Int_parsing_result { end_of_number, integer_overflow, end_of_buffer };

    Int_parsing_result parse_uint(unsigned* number, unsigned* number_len)
    {
        if (remaining_size == 0) {
            if (eof_reached) {
                set_terminator(EOF);
                return end_of_number;
            }
            return end_of_buffer;
        }

        unsigned digit;

        do {
            if ((digit = (unsigned) *current_ptr - '0') > 9) {
                set_terminator_and_inc_line_num_if_newline(
                        (unsigned char) *current_ptr);
                ++current_ptr;
                --remaining_size;
                return end_of_number;
            }

            if (*number > (UINT_MAX / 10) ||
                    (*number == (UINT_MAX / 10) && digit > UINT_MAX % 10)) {
                return integer_overflow;
            }

            *number = *number * 10 + digit;
            ++*number_len;

            ++current_ptr;
        } while (--remaining_size > 0);

        return end_of_buffer;
    }

    bool prepare_token_or_accumulate(const char* const end_of_token)
    {
        if (end_of_token == nullptr) {
            if (!eof_reached) {
                if (accumulating) {
                    accumulator.append(current_ptr, remaining_size);
                } else {
                    accumulating = true;
                    accumulator.assign(current_ptr, remaining_size);
                }

                return false;
            }

            // End of file has been reached. Return the accumulated
            // bytes as the last token.

            set_terminator(EOF);
            if (!accumulating) {
                token.clear();
            } else {
                accumulating = false;
                token = accumulator;
            }
            return true;
        }

        set_terminator_and_inc_line_num_if_newline(
                (unsigned char) *end_of_token);

        const size_t token_len = end_of_token - current_ptr;

        if (!accumulating) {
            if (token_len > 0) {
                token.assign(current_ptr,
                        *end_of_token == '\n' && end_of_token[-1] == '\r' ?
                                token_len - 1 :
                                token_len);
            } else {
                token.clear();
            }
        } else {
            accumulating = false;
            if (token_len > 0) {
                accumulator.append(current_ptr,
                        *end_of_token == '\n' && end_of_token[-1] == '\r' ?
                                token_len - 1 :
                                token_len);
            } else if (*end_of_token == '\n' && accumulator.length() > 0 &&
                    accumulator.back() == '\r') {
                accumulator.pop_back();
            }

            token = accumulator;
        }

        advance_by(token_len + 1);

        return true;
    }

    bool skip_token(const char* const end_of_token)
    {
        accumulating = false;

        if (end_of_token == nullptr) {
            if (!eof_reached) {
                return false;
            }

            // End of file has been reached.
            set_terminator(EOF);
            return true;
        }

        set_terminator_and_inc_line_num_if_newline(
                (unsigned char) *end_of_token);

        const size_t skipped_len = end_of_token - current_ptr;

        advance_by(skipped_len + 1);

        return true;
    }

    const VCF_string_view& get_token() const
    {
        return token;
    }

    bool get_token_as_uint(unsigned* number) const
    {
        unsigned len = (unsigned int) token.size();

        if (len == 0) {
            return false;
        }

        *number = 0;

        const char* ptr = token.data();
        unsigned digit;

        do {
            if ((digit = (unsigned) *ptr - '0') > 9) {
                return false;
            }

            if (*number > (UINT_MAX / 10) ||
                    (*number == (UINT_MAX / 10) && digit > UINT_MAX % 10)) {
                return false;
            }

            *number = *number * 10 + digit;
            ++ptr;
        } while (--len > 0);

        return true;
    }

    bool token_is_dot() const
    {
        return token.empty() || (token.length() == 1 && token.front() == '.');
    }

    bool token_is_last() const
    {
        return terminator == '\n' || terminator == EOF;
    }

    bool get_key_value(
            VCF_string_view* key, VCF_string_view* value, char delim = '=')
    {
        size_t delim_pos = token.find(delim);

        if (delim_pos == std::string::npos) {
            return false;
        }

        if (key == nullptr) {
            if (value != nullptr) {
                ++delim_pos;

                value->assign(
                        token.data() + delim_pos, token.length() - delim_pos);
            }

            return true;
        }

        if (value == nullptr) {
            key->assign(token.data(), delim_pos);

            return true;
        }

        key->assign(token.data(), delim_pos);

        ++delim_pos;

        value->assign(token.data() + delim_pos, token.length() - delim_pos);

        return true;
    }

    unsigned get_line_number() const
    {
        return line_number;
    }

    int get_terminator() const
    {
        return terminator;
    }

public:
    VCF_tokenizer()
    {
        memset(newline_or_tab, 0, sizeof(newline_or_tab));
        newline_or_tab['\n'] = true;
        newline_or_tab['\t'] = true;

        memset(newline_or_tab_or_equals, 0, sizeof(newline_or_tab_or_equals));
        newline_or_tab_or_equals['\n'] = true;
        newline_or_tab_or_equals['\t'] = true;
        newline_or_tab_or_equals['='] = true;

        memset(newline_or_tab_or_semicolon, 0,
                sizeof(newline_or_tab_or_semicolon));
        newline_or_tab_or_semicolon['\n'] = true;
        newline_or_tab_or_semicolon['\t'] = true;
        newline_or_tab_or_semicolon[';'] = true;

        memset(newline_or_tab_or_comma, 0, sizeof(newline_or_tab_or_comma));
        newline_or_tab_or_comma['\n'] = true;
        newline_or_tab_or_comma['\t'] = true;
        newline_or_tab_or_comma[','] = true;

        memset(newline_or_tab_or_colon, 0, sizeof(newline_or_tab_or_colon));
        newline_or_tab_or_colon['\n'] = true;
        newline_or_tab_or_colon['\t'] = true;
        newline_or_tab_or_colon[':'] = true;

        memset(newline_tab_colon_slash_bar, 0,
                sizeof(newline_tab_colon_slash_bar));
        newline_tab_colon_slash_bar['\n'] = true;
        newline_tab_colon_slash_bar['\t'] = true;
        newline_tab_colon_slash_bar[':'] = true;
        newline_tab_colon_slash_bar['/'] = true;
        newline_tab_colon_slash_bar['|'] = true;
    }

private:
    unsigned line_number = 1;
    int terminator;

    const char* current_ptr;
    size_t remaining_size;
    bool eof_reached;

    bool accumulating = false;
    std::string accumulator;

    VCF_string_view token;

public:
    // For parsing the meta-information lines
    // as well as the first token of the header line
    bool newline_or_tab_or_equals[256];
    // For extracting CHROM, POS, REF, or QUAL fields,
    // or skipping any other field
    bool newline_or_tab[256];
    // For extracting ID, FILTER, or INFO
    bool newline_or_tab_or_semicolon[256];
    // For extracting ALT
    bool newline_or_tab_or_comma[256];
    // For extracting FORMAT or GENOTYPE
    bool newline_or_tab_or_colon[256];
    // For extracting the GT values
    bool newline_tab_colon_slash_bar[256];
};

#endif /* !defined(VCF_TOKENIZER__HH) */
