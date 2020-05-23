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

#ifndef VCF_TOKENIZER__HH
#define VCF_TOKENIZER__HH

#include <string>
#include <iostream>
#include <climits>
#include <cstring>

typedef std::string CTempString;

namespace vcf {

// CVCFTokenizer finds tokens in the input VCF stream.
class CVCFTokenizer
{
public:
    void SetNewBuffer(const char* buffer, size_t buffer_size)
    {
        m_CurrentPtr = (m_RemainingSize = buffer_size) > 0 ? buffer : nullptr;
    }

    bool BufferIsEmpty() const
    {
        return m_RemainingSize == 0;
    }

    bool AtEOF() const
    {
        return m_CurrentPtr == nullptr;
    }

    const char* FindNewline() const
    {
        return (const char*) memchr(m_CurrentPtr, '\n', m_RemainingSize);
    }

private:
    static const char* x_FindCharFromSet(
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
    const char* FindCharFromSet(const bool* character_set) const
    {
        return x_FindCharFromSet(m_CurrentPtr, m_RemainingSize, character_set);
    }

    const char* FindNewlineOrTab() const
    {
        return x_FindCharFromSet(m_CurrentPtr, m_RemainingSize, m_NewlineOrTab);
    }

    const char* FindNewlineOrTabOrEquals() const
    {
        return x_FindCharFromSet(
                m_CurrentPtr, m_RemainingSize, m_NewlineOrTabOrEquals);
    }

    const char* FindNewlineOrTabOrSemicolon() const
    {
        return x_FindCharFromSet(
                m_CurrentPtr, m_RemainingSize, m_NewlineOrTabOrSemicolon);
    }

    const char* FindNewlineOrTabOrComma() const
    {
        return x_FindCharFromSet(
                m_CurrentPtr, m_RemainingSize, m_NewlineOrTabOrComma);
    }

    const char* FindNewlineOrTabOrColon() const
    {
        return x_FindCharFromSet(
                m_CurrentPtr, m_RemainingSize, m_NewlineOrTabOrColon);
    }

private:
    void x_SetTokenTermAndPossiblyIncrementLineNumber(int token_term)
    {
        if ((m_TokenTerm = token_term) == '\n') {
            ++m_LineNumber;
        }
    }

    void x_AdvanceBy(size_t number_of_bytes)
    {
        m_CurrentPtr += number_of_bytes;
        m_RemainingSize -= number_of_bytes;
    }

public:
    enum EIntParsingResult { eEndOfNumber, eIntegerOverflow, eEndOfBuffer };

    EIntParsingResult ParseUnsignedInt(unsigned* number, unsigned* number_len)
    {
        if (m_RemainingSize == 0) {
            if (m_CurrentPtr == nullptr) {
                x_SetTokenTermAndPossiblyIncrementLineNumber(EOF);
                return eEndOfNumber;
            }
            return eEndOfBuffer;
        }

        unsigned digit;

        do {
            if ((digit = (unsigned) *m_CurrentPtr - '0') > 9) {
                x_SetTokenTermAndPossiblyIncrementLineNumber(
                        (unsigned char) *m_CurrentPtr);
                ++m_CurrentPtr;
                --m_RemainingSize;
                return eEndOfNumber;
            }

            if (*number > (UINT_MAX / 10) ||
                    (*number == (UINT_MAX / 10) && digit > UINT_MAX % 10)) {
                return eIntegerOverflow;
            }

            *number = *number * 10 + digit;
            ++*number_len;

            ++m_CurrentPtr;
        } while (--m_RemainingSize > 0);

        return eEndOfBuffer;
    }

    bool PrepareTokenOrAccumulate(const char* const end_of_token)
    {
        if (end_of_token == nullptr) {
            if (m_CurrentPtr != nullptr) { // Check for EOF
                if (m_Accumulating) {
                    m_Accumulator.append(m_CurrentPtr, m_RemainingSize);
                } else {
                    m_Accumulating = true;
                    m_Accumulator.assign(m_CurrentPtr, m_RemainingSize);
                }

                return false;
            }

            // EOF has been reached. Return whatever has been
            // accumulated as the last token.

            x_SetTokenTermAndPossiblyIncrementLineNumber(EOF);
            if (!m_Accumulating) {
                m_Token.clear();
            } else {
                m_Accumulating = false;
                m_Token = m_Accumulator;
            }
            return true;
        }

        x_SetTokenTermAndPossiblyIncrementLineNumber(
                (unsigned char) *end_of_token);

        const size_t token_len = end_of_token - m_CurrentPtr;

        if (!m_Accumulating) {
            if (token_len > 0) {
                m_Token.assign(m_CurrentPtr,
                        *end_of_token == '\n' && end_of_token[-1] == '\r' ?
                                token_len - 1 :
                                token_len);
            } else {
                m_Token.clear();
            }
        } else {
            m_Accumulating = false;
            if (token_len > 0) {
                m_Accumulator.append(m_CurrentPtr,
                        *end_of_token == '\n' && end_of_token[-1] == '\r' ?
                                token_len - 1 :
                                token_len);
            } else if (*end_of_token == '\n' && m_Accumulator.length() > 0 &&
                    m_Accumulator.back() == '\r') {
                m_Accumulator.pop_back();
            }

            m_Token = m_Accumulator;
        }

        x_AdvanceBy(token_len + 1);

        return true;
    }

    bool SkipToken(const char* const end_of_token)
    {
        m_Accumulating = false;

        if (end_of_token == nullptr) {
            if (m_CurrentPtr != nullptr) { // Check for EOF
                return false;
            }

            // EOF has been reached.
            x_SetTokenTermAndPossiblyIncrementLineNumber(EOF);
            return true;
        }

        x_SetTokenTermAndPossiblyIncrementLineNumber(
                (unsigned char) *end_of_token);

        const size_t skipped_len = end_of_token - m_CurrentPtr;

        x_AdvanceBy(skipped_len + 1);

        return true;
    }

    const CTempString& GetToken() const
    {
        return m_Token;
    }

    bool GetTokenAsUInt(unsigned* number) const
    {
        unsigned len = (unsigned int) m_Token.size();

        if (len == 0) {
            return false;
        }

        *number = 0;

        const char* ptr = m_Token.data();
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

    bool TokenIsDot() const
    {
        return m_Token.empty() || m_Token == ".";
    }

    bool TokenIsLast() const
    {
        return m_TokenTerm == '\n' || m_TokenTerm == EOF;
    }

    bool GetKeyValue(CTempString* key, CTempString* value, char delim = '=')
    {
        size_t delim_pos = m_Token.find(delim);

        if (delim_pos == std::string::npos) {
            return false;
        }

        if (key == nullptr) {
            if (value != nullptr) {
                ++delim_pos;

                value->assign(m_Token.data() + delim_pos,
                        m_Token.length() - delim_pos);
            }

            return true;
        }

        if (value == nullptr) {
            key->assign(m_Token.data(), delim_pos);

            return true;
        }

        key->assign(m_Token.data(), delim_pos);

        ++delim_pos;

        value->assign(m_Token.data() + delim_pos, m_Token.length() - delim_pos);

        return true;
    }

    unsigned GetLineNumber() const
    {
        return m_LineNumber;
    }

    int GetTokenTerm() const
    {
        return m_TokenTerm;
    }

public:
    CVCFTokenizer()
    {
        memset(m_NewlineOrTab, 0, sizeof(m_NewlineOrTab));
        m_NewlineOrTab['\n'] = true;
        m_NewlineOrTab['\t'] = true;

        memset(m_NewlineOrTabOrEquals, 0, sizeof(m_NewlineOrTabOrEquals));
        m_NewlineOrTabOrEquals['\n'] = true;
        m_NewlineOrTabOrEquals['\t'] = true;
        m_NewlineOrTabOrEquals['='] = true;

        memset(m_NewlineOrTabOrSemicolon, 0, sizeof(m_NewlineOrTabOrSemicolon));
        m_NewlineOrTabOrSemicolon['\n'] = true;
        m_NewlineOrTabOrSemicolon['\t'] = true;
        m_NewlineOrTabOrSemicolon[';'] = true;

        memset(m_NewlineOrTabOrComma, 0, sizeof(m_NewlineOrTabOrComma));
        m_NewlineOrTabOrComma['\n'] = true;
        m_NewlineOrTabOrComma['\t'] = true;
        m_NewlineOrTabOrComma[','] = true;

        memset(m_NewlineOrTabOrColon, 0, sizeof(m_NewlineOrTabOrColon));
        m_NewlineOrTabOrColon['\n'] = true;
        m_NewlineOrTabOrColon['\t'] = true;
        m_NewlineOrTabOrColon[':'] = true;

        memset(m_NewlineTabColonSlashBar, 0, sizeof(m_NewlineTabColonSlashBar));
        m_NewlineTabColonSlashBar['\n'] = true;
        m_NewlineTabColonSlashBar['\t'] = true;
        m_NewlineTabColonSlashBar[':'] = true;
        m_NewlineTabColonSlashBar['/'] = true;
        m_NewlineTabColonSlashBar['|'] = true;
    }

private:
    unsigned m_LineNumber = 1;
    int m_TokenTerm = '\0';

    const char* m_CurrentPtr = nullptr;
    size_t m_RemainingSize = 0;

    std::string m_Accumulator;
    bool m_Accumulating = false;

    CTempString m_Token;

public:
    // For parsing the meta-information lines
    // as well as the first token of the header line
    bool m_NewlineOrTabOrEquals[256];
    // For extracting CHROM, POS, REF, or QUAL fields,
    // or skipping any other field
    bool m_NewlineOrTab[256];
    // For extracting ID, FILTER, or INFO
    bool m_NewlineOrTabOrSemicolon[256];
    // For extracting ALT
    bool m_NewlineOrTabOrComma[256];
    // For extracting FORMAT or GENOTYPE
    bool m_NewlineOrTabOrColon[256];
    // For extracting the GT values
    bool m_NewlineTabColonSlashBar[256];
};

} /* namespace vcf */

#endif /* !defined(VCF_TOKENIZER__HH) */