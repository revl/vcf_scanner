#include <stdlib.h>
#include <stdio.h>
#include "tokenizer.hh"

constexpr size_t kTestDataSize = 512 * 1024 * 1024;

int main()
{
    srand(1234);

    char* buf = new char[kTestDataSize];
    if (buf == nullptr)
        return 1;

    const char* end_of_buf = buf + kTestDataSize;
    for (char* b = buf; b < end_of_buf; ++b)
        *b = (char) rand();

    CVCFTokenizer tokenizer;

    size_t number_of_newlines = 0;
    size_t number_of_tabs = 0;
    size_t number_of_commas = 0;

    tokenizer.SetNewBuffer(buf, kTestDataSize);
    for (;;) {
        const char* end_of_token = tokenizer.FindNewlineOrTabOrComma();
        if (!tokenizer.SkipToken(end_of_token))
            break;
        if (*end_of_token == '\n')
            ++number_of_newlines;
        else if (*end_of_token == '\t')
            ++number_of_tabs;
        else if (*end_of_token == ',')
            ++number_of_commas;
        else
            puts("WTF");
    }

    printf("Number of newlines: %lu\n"
           "Number of tabs: %lu\n"
           "Number of commas: %lu\n",
            (long unsigned) number_of_newlines, (long unsigned) number_of_tabs,
            (long unsigned) number_of_commas);

    return 0;
}