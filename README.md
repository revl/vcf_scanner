# VCF_scanner - a blazingly fast Variant Call Format parser

## Features

*   The parser does not read the input VCF file stream and leaves it to the
    caller to load the input data into memory one buffer at a time.  As a
    result, the parser is non-blocking. The caller can choose to read input
    data in a separate thread so that reading and parsing happen in parallel.
*   The caller decides which VCF fields to parse. Fields that are not requested
    by the caller are skipped and not parsed.
*   Very few bytes in the input buffer are accessed more than once.
*   Memory is allocated frugally.
*   Exceptions are not used for error reporting.
*   The library is header-only with no dependencies outside the standard
    library.

## How to use

### Preparation

1.  Create an instance of `VCF_parser`.

2.  Implement a function that reads the next chunk of input data into a buffer.
    The buffer must be external to the function.

        bool parse_to_completion(
                VCF_parsing_event pe, VCF_scanner& vcf_scanner, FILE* input)
        {
            while (pe == VCF_parsing_event::need_more_data) {
                read_buffer(input);
                pe = vcf_scanner.feed(buffer, buffer_size);
            }

            if (pe == VCF_parsing_event::error) {
                return false;
            }

            if (pe == VCF_parsing_event::ok_with_warnings) {
                for (const auto& warning : vcf_scanner.get_warnings()) {
                    std::cerr << "Warning: " << warning.warning_message << std::endl;
                }
            }

            return true;
        }

### Parsing the header

2.  Supply input data to the parser using the `feed()` method until it returns
    `ok`, which means that the header has been successfully parsed and can be
    accessed using the `get_header()` method.

    Do not discard the current input buffer - it contains the beginning of the
    data lines.  Proceed to the data line parsing loop.

    If the `feed()` method returns `error` while parsing the header, the file
    cannot be parsed. Call `get_error()` to get the error message and
    `get_line_number()` to get the line number containing the error.

### Data line parsing loop

All data line fields are optional and can be omitted by not calling the
respective `parse_...` methods.

1.  Request parsing of the CHROM and POS fields by calling `parse_loc()`.
	1.  If the above method returns `need_more_data`, more input data must
	    be provided using the `feed()` method until it returns `ok`.
	2.  When `parse_loc()` or 

## To build a test coverage report

1.  Install `lcov`
2.  Configure, build, and run tests as follows:

        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_COVERAGE=TRUE ..
        make coverage

3.  Find the HTML report in the `coverage` directory:

        xdg-open coverage/index.html
