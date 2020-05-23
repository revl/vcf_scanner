#include <vcf_scanner/impl/vcf_tokenizer.hh>

#include "catch.hh"

using Catch::Matchers::Equals;

TEST_CASE("Newline, no newline")
{
    static const char test_data[] = "two\nlines";

    VCF_tokenizer tokenizer;

    // Start with a non-empty buffer
    tokenizer.set_new_buffer(test_data, sizeof(test_data) - 1);
    CHECK(tokenizer.get_line_number() == 1);

    CHECK(!tokenizer.buffer_is_empty());
    CHECK(!tokenizer.at_eof());

    // Find the newline character
    const char* newline = tokenizer.find_newline();
    REQUIRE(newline != nullptr);

    // Extract the token before the newline
    REQUIRE(tokenizer.prepare_token_or_accumulate(newline));
    CHECK(tokenizer.get_token() == "two");
    CHECK(tokenizer.get_terminator() == '\n');

    // The second line has started
    CHECK(tokenizer.get_line_number() == 2);

    // Confirm that there is no second newline
    newline = tokenizer.find_newline();
    REQUIRE(newline == nullptr);
    REQUIRE(!tokenizer.prepare_token_or_accumulate(newline));

    // It is unknown whether EOF has been reached
    CHECK(!tokenizer.at_eof());
    // The buffer is exhausted but the previous token
    // is not overwritten
    CHECK(tokenizer.get_token() == "two");
    CHECK(tokenizer.get_terminator() == '\n');

    // Simulate EOF condition
    tokenizer.set_new_buffer("", 0);

    // The buffer is still empty
    CHECK(tokenizer.buffer_is_empty());
    // And EOF condition is recognized
    CHECK(tokenizer.at_eof());

    newline = tokenizer.find_newline();
    REQUIRE(newline == nullptr);
    REQUIRE(tokenizer.prepare_token_or_accumulate(newline));

    CHECK(tokenizer.get_token() == "lines");
    CHECK(tokenizer.get_terminator() == EOF);
}

TEST_CASE("Skipping")
{
    static const char test_data[] = "1\n2";

    VCF_tokenizer tokenizer;

    tokenizer.set_new_buffer(test_data, sizeof(test_data) - 1);
    CHECK(tokenizer.get_line_number() == 1);

    CHECK(!tokenizer.buffer_is_empty());
    CHECK(!tokenizer.at_eof());

    // Find the first newline character
    const char* newline = tokenizer.find_newline();
    REQUIRE(newline != nullptr);

    // Skip the first line
    REQUIRE(tokenizer.skip_token(newline));
    CHECK(tokenizer.get_terminator() == '\n');

    // The second line has started
    CHECK(tokenizer.get_line_number() == 2);

    // Confirm that there is no second newline
    newline = tokenizer.find_newline();
    REQUIRE(newline == nullptr);
    REQUIRE(!tokenizer.skip_token(newline));

    // It is unknown whether EOF has been reached
    CHECK(!tokenizer.at_eof());
    CHECK(tokenizer.get_terminator() == '\n');

    // Simulate EOF condition
    tokenizer.set_new_buffer("", 0);
    newline = tokenizer.find_newline();
    REQUIRE(newline == nullptr);
    REQUIRE(tokenizer.skip_token(newline));

    CHECK(tokenizer.get_terminator() == EOF);
}

TEST_CASE("Empty token")
{
    VCF_tokenizer tokenizer;

    tokenizer.set_new_buffer("\t\n", 2);

    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));
    CHECK(tokenizer.get_token().empty());

    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));
    CHECK(tokenizer.get_token().empty());
}

static std::string stitch3(VCF_tokenizer& tokenizer, const std::string& part1,
        const std::string& part2, const std::string& part3)
{
    tokenizer.set_new_buffer(part1.data(), part1.length());
    REQUIRE(!tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));

    tokenizer.set_new_buffer(part2.data(), part2.length());
    REQUIRE(!tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));

    tokenizer.set_new_buffer(part3.data(), part3.length());
    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));

    return tokenizer.get_token();
}

TEST_CASE("Seams")
{
    VCF_tokenizer tokenizer;

    tokenizer.set_new_buffer("", 0);
    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));
    CHECK(tokenizer.get_token().empty());

    CHECK(stitch3(tokenizer, "heads ", "and", " tails\n") == "heads and tails");
    CHECK(stitch3(tokenizer, "heads ", "and", " tails\r\n") ==
            "heads and tails");
    CHECK(stitch3(tokenizer, "grid", "lock\r", "\n") == "gridlock");
    CHECK(stitch3(tokenizer, "grid", "lock", "") == "gridlock");
}

TEST_CASE("Key-value")
{
    VCF_tokenizer tokenizer;

    VCF_string_view k, v;

    static const char kv[] = "key=value\nnokeyvalue\n";
    tokenizer.set_new_buffer(kv, sizeof(kv) - 1);
    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));

    REQUIRE(tokenizer.get_key_value(&k, &v));
    CHECK(k == "key");
    CHECK(v == "value");

    k.clear();
    v.clear();
    REQUIRE(tokenizer.get_key_value(nullptr, &v));
    CHECK(v == "value");
    v.clear();
    REQUIRE(tokenizer.get_key_value(&k, nullptr));
    CHECK(k == "key");

    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));
    REQUIRE(!tokenizer.get_key_value(&k, &v));
}

TEST_CASE("Parse unsigned int")
{
    VCF_tokenizer tokenizer;

    static const char two_numbers[] = "\t12345-6789";
    tokenizer.set_new_buffer(two_numbers, sizeof(two_numbers) - 1);

    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));

    unsigned number = 0, number_len = 0;
    REQUIRE(tokenizer.parse_uint(&number, &number_len) ==
            VCF_tokenizer::end_of_number);
    CHECK(number == 12345);
    CHECK(number_len == 5);
    CHECK(tokenizer.get_terminator() == '-');

    number = number_len = 0;
    REQUIRE(tokenizer.parse_uint(&number, &number_len) ==
            VCF_tokenizer::end_of_buffer);
    CHECK(number == 6789);
    CHECK(number_len == 4);

    number = number_len = 0;
    REQUIRE(tokenizer.parse_uint(&number, &number_len) ==
            VCF_tokenizer::end_of_buffer);
    CHECK(number == 0);
    CHECK(number_len == 0);

    static const char overflow[] = "4294967296";
    tokenizer.set_new_buffer(overflow, sizeof(overflow) - 1);
    number = number_len = 0;
    REQUIRE(tokenizer.parse_uint(&number, &number_len) ==
            VCF_tokenizer::integer_overflow);

    tokenizer.set_new_buffer("", 0);
    number = number_len = 0;
    REQUIRE(tokenizer.parse_uint(&number, &number_len) ==
            VCF_tokenizer::end_of_number);
    CHECK(number == 0);
    CHECK(number_len == 0);

    static const char test_data[] = "123456789\n4294967296\n\n100X\n";
    tokenizer.set_new_buffer(test_data, sizeof(test_data) - 1);

    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));
    REQUIRE(tokenizer.get_token_as_uint(&number));
    CHECK(number == 123456789);

    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));
    CHECK(!tokenizer.get_token_as_uint(&number));

    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));
    CHECK(!tokenizer.get_token_as_uint(&number));

    REQUIRE(tokenizer.prepare_token_or_accumulate(
            tokenizer.find_newline_or_tab()));
    CHECK(!tokenizer.get_token_as_uint(&number));
}

TEST_CASE("Simple checks")
{
    VCF_tokenizer tokenizer;

    static const char test_data[] = ".\n. \n";
    tokenizer.set_new_buffer(test_data, sizeof(test_data) - 1);
    REQUIRE(tokenizer.prepare_token_or_accumulate(tokenizer.find_newline()));

    CHECK(tokenizer.get_token() == ".");
    CHECK(tokenizer.token_is_dot());
    CHECK(tokenizer.token_is_last());

    REQUIRE(tokenizer.prepare_token_or_accumulate(tokenizer.find_newline()));
    CHECK(!tokenizer.token_is_dot());
    CHECK(tokenizer.token_is_last());

    tokenizer.set_new_buffer(test_data, sizeof(test_data) - 1);
    REQUIRE(tokenizer.prepare_token_or_accumulate(tokenizer.find_newline()));
    CHECK(tokenizer.token_is_last());
}
