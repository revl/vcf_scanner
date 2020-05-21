#include <vcf_scanner/impl/vcf_tokenizer.hh>

#include "test_case.h"

TEST_CASE(newline_no_newline)
{
    static const char test_data[] = "two\nlines";

    CVCFTokenizer tokenizer;

    // Start with a non-empty buffer
    tokenizer.SetNewBuffer(test_data, sizeof(test_data) - 1);
    CHECK(tokenizer.GetLineNumber() == 1);

    CHECK(!tokenizer.BufferIsEmpty());
    CHECK(!tokenizer.AtEOF());

    // Find the newline character
    const char* newline = tokenizer.FindNewline();
    REQUIRE(newline != nullptr);

    // Extract the token before the newline
    REQUIRE(tokenizer.PrepareTokenOrAccumulate(newline));
    CHECK(tokenizer.GetToken() == "two");
    CHECK(tokenizer.GetTokenTerm() == '\n');

    // The second line has started
    CHECK(tokenizer.GetLineNumber() == 2);

    // Confirm that there is no second newline
    newline = tokenizer.FindNewline();
    REQUIRE(newline == nullptr);
    REQUIRE(!tokenizer.PrepareTokenOrAccumulate(newline));

    // It is unknown whether EOF has been reached
    CHECK(!tokenizer.AtEOF());
    // The buffer is exhausted but the previous token
    // is not overwritten
    CHECK(tokenizer.GetToken() == "two");
    CHECK(tokenizer.GetTokenTerm() == '\n');

    // Simulate EOF condition
    tokenizer.SetNewBuffer("", 0);

    // The buffer is still empty
    CHECK(tokenizer.BufferIsEmpty());
    // And EOF condition is recognized
    CHECK(tokenizer.AtEOF());

    newline = tokenizer.FindNewline();
    REQUIRE(newline == nullptr);
    REQUIRE(tokenizer.PrepareTokenOrAccumulate(newline));

    CHECK(tokenizer.GetToken() == "lines");
    CHECK(tokenizer.GetTokenTerm() == EOF);
}

TEST_CASE(skipping)
{
    static const char test_data[] = "1\n2";

    CVCFTokenizer tokenizer;

    tokenizer.SetNewBuffer(test_data, sizeof(test_data) - 1);
    CHECK(tokenizer.GetLineNumber() == 1);

    CHECK(!tokenizer.BufferIsEmpty());
    CHECK(!tokenizer.AtEOF());

    // Find the first newline character
    const char* newline = tokenizer.FindNewline();
    REQUIRE(newline != nullptr);

    // Skip the first line
    REQUIRE(tokenizer.SkipToken(newline));
    CHECK(tokenizer.GetTokenTerm() == '\n');

    // The second line has started
    CHECK(tokenizer.GetLineNumber() == 2);

    // Confirm that there is no second newline
    newline = tokenizer.FindNewline();
    REQUIRE(newline == nullptr);
    REQUIRE(!tokenizer.SkipToken(newline));

    // It is unknown whether EOF has been reached
    CHECK(!tokenizer.AtEOF());
    CHECK(tokenizer.GetTokenTerm() == '\n');

    // Simulate EOF condition
    tokenizer.SetNewBuffer("", 0);
    newline = tokenizer.FindNewline();
    REQUIRE(newline == nullptr);
    REQUIRE(tokenizer.SkipToken(newline));

    CHECK(tokenizer.GetTokenTerm() == EOF);
}

TEST_CASE(empty_token)
{
    CVCFTokenizer tokenizer;

    tokenizer.SetNewBuffer("\t\n", 2);

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));
    CHECK(tokenizer.GetToken().empty());

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));
    CHECK(tokenizer.GetToken().empty());
}

static string s_Stitch3(CVCFTokenizer& tokenizer, const string& part1,
        const string& part2, const string& part3)
{
    tokenizer.SetNewBuffer(part1.data(), part1.length());
    REQUIRE(!tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));

    tokenizer.SetNewBuffer(part2.data(), part2.length());
    REQUIRE(!tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));

    tokenizer.SetNewBuffer(part3.data(), part3.length());
    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));

    return tokenizer.GetToken();
}

TEST_CASE(seams)
{
    CVCFTokenizer tokenizer;

    tokenizer.SetNewBuffer("", 0);
    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));
    CHECK(tokenizer.GetToken().empty());

    CHECK(s_Stitch3(tokenizer, "heads ", "and", " tails\n") ==
            "heads and tails");
    CHECK(s_Stitch3(tokenizer, "heads ", "and", " tails\r\n") ==
            "heads and tails");
    CHECK(s_Stitch3(tokenizer, "grid", "lock\r", "\n") == "gridlock");
    CHECK(s_Stitch3(tokenizer, "grid", "lock", "") == "gridlock");
}

TEST_CASE(key_value)
{
    CVCFTokenizer tokenizer;

    CTempString k, v;

    static const char kv[] = "key=value\nnokeyvalue\n";
    tokenizer.SetNewBuffer(kv, sizeof(kv) - 1);
    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));

    REQUIRE(tokenizer.GetKeyValue(&k, &v));
    CHECK(k == "key" && v == "value");

    k.clear();
    v.clear();
    REQUIRE(tokenizer.GetKeyValue(nullptr, &v));
    CHECK(v == "value");
    v.clear();
    REQUIRE(tokenizer.GetKeyValue(&k, nullptr));
    CHECK(k == "key");

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));
    REQUIRE(!tokenizer.GetKeyValue(&k, &v));
}

TEST_CASE(parse_unsigned_int)
{
    CVCFTokenizer tokenizer;

    static const char two_numbers[] = "\t12345-6789";
    tokenizer.SetNewBuffer(two_numbers, sizeof(two_numbers) - 1);

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));

    unsigned number = 0, number_len = 0;
    REQUIRE(tokenizer.ParseUnsignedInt(&number, &number_len) ==
            CVCFTokenizer::eEndOfNumber);
    CHECK(number == 12345 && number_len == 5);
    CHECK(tokenizer.GetTokenTerm() == '-');

    number = number_len = 0;
    REQUIRE(tokenizer.ParseUnsignedInt(&number, &number_len) ==
            CVCFTokenizer::eEndOfBuffer);
    CHECK(number == 6789 && number_len == 4);

    number = number_len = 0;
    REQUIRE(tokenizer.ParseUnsignedInt(&number, &number_len) ==
            CVCFTokenizer::eEndOfBuffer);
    CHECK(number == 0 && number_len == 0);

    static const char overflow[] = "4294967296";
    tokenizer.SetNewBuffer(overflow, sizeof(overflow) - 1);
    number = number_len = 0;
    REQUIRE(tokenizer.ParseUnsignedInt(&number, &number_len) ==
            CVCFTokenizer::eIntegerOverflow);

    tokenizer.SetNewBuffer("", 0);
    number = number_len = 0;
    REQUIRE(tokenizer.ParseUnsignedInt(&number, &number_len) ==
            CVCFTokenizer::eEndOfNumber);
    CHECK(number == 0 && number_len == 0);

    static const char test_data[] = "123456789\n4294967296\n\n100X\n";
    tokenizer.SetNewBuffer(test_data, sizeof(test_data) - 1);

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));
    REQUIRE(tokenizer.GetTokenAsUInt(&number));
    CHECK(number == 123456789);

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));
    CHECK(!tokenizer.GetTokenAsUInt(&number));

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));
    CHECK(!tokenizer.GetTokenAsUInt(&number));

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewlineOrTab()));
    CHECK(!tokenizer.GetTokenAsUInt(&number));
}

TEST_CASE(simple_checks)
{
    CVCFTokenizer tokenizer;

    static const char test_data[] = ".\n. \n";
    tokenizer.SetNewBuffer(test_data, sizeof(test_data) - 1);
    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewline()));

    CHECK(tokenizer.GetToken() == ".");
    CHECK(tokenizer.TokenIsDot());
    CHECK(tokenizer.TokenIsLast());

    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewline()));
    CHECK(!tokenizer.TokenIsDot());
    CHECK(tokenizer.TokenIsLast());

    tokenizer.SetNewBuffer(test_data, sizeof(test_data) - 1);
    REQUIRE(tokenizer.PrepareTokenOrAccumulate(tokenizer.FindNewline()));
    CHECK(tokenizer.TokenIsLast());
}
