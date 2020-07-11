# VCF_scanner - a blazingly fast Variant Call Format parser

## Features

*   The parser itself does not read the input VCF file and leaves it to the
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

### Parsing the header

1.  Create an instance of `VCF_parser`.

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

## To build a test coverage report

1.  Install `lcov`
2.  Configure, build, and run tests as follows:

        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_COVERAGE=TRUE ..
        make coverage

3.  Find the HTML report in the `coverage` directory:

        xdg-open coverage/index.html
