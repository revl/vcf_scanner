#include <string>
#include <string.h>

using namespace std;

class CVCFTokenizer
{
public:
    CVCFTokenizer() : m_Accumulating(false)
    {
        memset(m_NewlineOrTab, 0, sizeof(m_NewlineOrTab));
        m_NewlineOrTab['\n'] = true;
        m_NewlineOrTab['\t'] = true;

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

    void SetNewBuffer(const char* buffer, size_t buffer_size)
    {
        m_CurrentPtr = buffer;
        m_RemainingSize = buffer_size;
    }

    bool Peek(char* buffer, size_t how_many_bytes)
    {
        if (m_Accumulating) {
            const size_t acc_len = m_Accumulator.length();
            if (acc_len >= how_many_bytes) {
                memcpy(buffer, m_Accumulator.data(), how_many_bytes);
                return true;
            }

            how_many_bytes -= acc_len;
            if (m_RemainingSize >= how_many_bytes) {
                memcpy(buffer, m_Accumulator.data(), acc_len);
                memcpy(buffer + acc_len, m_CurrentPtr, how_many_bytes);
                return true;
            }
        } else if (m_RemainingSize >= how_many_bytes) {
            memcpy(buffer, m_CurrentPtr, how_many_bytes);
            return true;
        }

        x_AccumulateRemainingBytes();
        return false;
    }

    const char* FindNewline()
    {
        return (const char*) memchr(m_CurrentPtr, '\n', m_RemainingSize);
    }

    const char* FindNewlineOrTab()
    {
        return x_FindCharFromSet(m_CurrentPtr, m_RemainingSize, m_NewlineOrTab);
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

    bool PrepareTokenOrAccumulate(const char* end_of_token)
    {
        if (end_of_token == nullptr) {
            x_AccumulateRemainingBytes();
            return false;
        }

        size_t token_len = end_of_token - m_CurrentPtr;

        if (!m_Accumulating)
            m_Token.assign(m_CurrentPtr, token_len);
        else {
            m_Accumulating = false;
            m_Accumulator.append(m_CurrentPtr, token_len);
            m_Token = m_Accumulator;
        }

        ++token_len;

        m_CurrentPtr += token_len;
        m_RemainingSize -= token_len;

        return true;
    }

    bool SkipToken(const char* end_of_token)
    {
        m_Accumulating = false;

        if (end_of_token == nullptr)
            return false;

        size_t skip_bytes = end_of_token + 1 - m_CurrentPtr;

        m_CurrentPtr += skip_bytes;
        m_RemainingSize -= skip_bytes;

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

private:
    static const char* x_FindCharFromSet(
            const char* buffer, size_t buffer_size, const bool* character_set)
    {
        while (buffer_size > 0) {
            if (character_set[(unsigned char) *buffer])
                return buffer;
            ++buffer;
            --buffer_size;
        }
        return nullptr;
    }

    void x_AccumulateRemainingBytes()
    {
        if (m_Accumulating)
            m_Accumulator.append(m_CurrentPtr, m_RemainingSize);
        else {
            m_Accumulating = true;
            m_Accumulator.assign(m_CurrentPtr, m_RemainingSize);
        }
    }

    const char* m_CurrentPtr;
    size_t m_RemainingSize;

    string m_Accumulator;
    bool m_Accumulating;

    string m_Token;

    // For extracting CHROM, POS, REF, or QUAL fields,
    // or skipping any other field
    bool m_NewlineOrTab[256];
    // For extracting ID, FORMAT, or INFO
    bool m_NewlineOrTabOrSemicolon[256];
    // FOR extracting ALT
    bool m_NewlineOrTabOrComma[256];
    // For extracting FORMAT or GENOTYPE
    bool m_NewlineOrTabOrColon[256];
};
