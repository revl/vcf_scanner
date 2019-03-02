#include <string>
#include <string.h>
#include <iostream>
#include <limits.h>

// For the EOF definition
#include <stdio.h>

using namespace std;

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

    const char* FindNewline()
    {
        return (const char*) memchr(m_CurrentPtr, '\n', m_RemainingSize);
    }

private:
    static const char* x_FindCharFromSet(
            const char* buffer, size_t buffer_size, const bool* character_set)
    {
        for (; buffer_size > 0; ++buffer, --buffer_size)
            if (character_set[(unsigned char) *buffer])
                return buffer;

        return nullptr;
    }

public:
    const char* FindCharFromSet(const bool* character_set)
    {
        return x_FindCharFromSet(m_CurrentPtr, m_RemainingSize, character_set);
    }

    const char* FindNewlineOrTab()
    {
        return x_FindCharFromSet(m_CurrentPtr, m_RemainingSize, m_NewlineOrTab);
    }

    const char* FindNewlineOrTabOrEquals()
    {
        return x_FindCharFromSet(
                m_CurrentPtr, m_RemainingSize, m_NewlineOrTabOrEquals);
    }

    const char* FindNewlineOrTabOrSemicolon()
    {
        return x_FindCharFromSet(
                m_CurrentPtr, m_RemainingSize, m_NewlineOrTabOrSemicolon);
    }

    const char* FindNewlineOrTabOrComma()
    {
        return x_FindCharFromSet(
                m_CurrentPtr, m_RemainingSize, m_NewlineOrTabOrComma);
    }

    const char* FindNewlineOrTabOrColon()
    {
        return x_FindCharFromSet(
                m_CurrentPtr, m_RemainingSize, m_NewlineOrTabOrColon);
    }

private:
    void x_SetTokenTermAndPossiblyIncrementLineNumber(int token_term)
    {
        if (m_TokenTerm == '\n')
            ++m_LineNumber;
        m_TokenTerm = token_term;
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
            if (AtEOF()) {
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
                return eEndOfNumber;
            }

            if (*number > (UINT_MAX / 10) ||
                    (*number == (UINT_MAX / 10) && digit > UINT_MAX % 10))
                return eIntegerOverflow;

            *number = *number * 10 + digit;
            ++*number_len;

            ++m_CurrentPtr;
        } while (--m_RemainingSize > 0);

        return eEndOfBuffer;
    }

    bool PrepareTokenOrAccumulate(const char* const end_of_token)
    {
        if (end_of_token == nullptr) {
            if (!AtEOF()) {
                x_SetTokenTermAndPossiblyIncrementLineNumber('\0');

                if (m_Accumulating)
                    m_Accumulator.append(m_CurrentPtr, m_RemainingSize);
                else {
                    m_Accumulating = true;
                    m_Accumulator.assign(m_CurrentPtr, m_RemainingSize);
                }

                return false;
            }

            // EOF has been reached. Return whatever has been
            // accumulated as the last token.

            x_SetTokenTermAndPossiblyIncrementLineNumber(EOF);
            if (!m_Accumulating)
                m_Token.clear();
            else {
                m_Accumulating = false;
                m_Token = m_Accumulator;
            }
            return true;
        }

        x_SetTokenTermAndPossiblyIncrementLineNumber(
                (unsigned char) *end_of_token);

        const size_t token_len = end_of_token - m_CurrentPtr;

        if (!m_Accumulating)
            if (token_len > 0)
                m_Token.assign(m_CurrentPtr,
                        *end_of_token == '\n' && end_of_token[-1] == '\r' ?
                                token_len - 1 :
                                token_len);
            else
                m_Token.clear();
        else {
            m_Accumulating = false;
            if (token_len > 0) {
                m_Accumulator.append(m_CurrentPtr,
                        *end_of_token == '\n' && end_of_token[-1] == '\r' ?
                                token_len - 1 :
                                token_len);
            } else if (*end_of_token == '\n' && m_Accumulator.length() > 0 &&
                    m_Accumulator.back() == '\r')
                m_Accumulator.pop_back();

            m_Token = m_Accumulator;
        }

        x_AdvanceBy(token_len + 1);

        return true;
    }

    bool SkipToken(const char* const end_of_token)
    {
        m_Accumulating = false;

        if (end_of_token == nullptr) {
            if (!AtEOF()) {
                x_SetTokenTermAndPossiblyIncrementLineNumber('\0');
                return false;
            }
            // EOF has been reached.
            x_SetTokenTermAndPossiblyIncrementLineNumber(EOF);
            return true;
        }

        x_SetTokenTermAndPossiblyIncrementLineNumber(
                (unsigned char) *end_of_token);

        x_AdvanceBy(end_of_token + 1 - m_CurrentPtr);

        return true;
    }

    const string& GetToken() const
    {
        return m_Token;
    }

    bool GetKeyValue(string* key, string* value, char delim = '=')
    {
        size_t delim_pos = m_Token.find(delim);

        if (delim_pos == string::npos)
            return false;

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
    }

private:
    unsigned m_LineNumber = 1;
    int m_TokenTerm = '\0';

    const char* m_CurrentPtr = nullptr;
    size_t m_RemainingSize = 0;

    string m_Accumulator;
    bool m_Accumulating = false;

    string m_Token;

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
};
