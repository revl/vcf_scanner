# VCF_scanner - a blazingly fast Variant Call Format parser

## Features

*   The parser does not read the input VCF file stream and leaves it to the
    caller to load the input data into memory one buffer at a time.  As a
    result, the parser is non-blocking. The caller can choose to read input
    data in a separate thread so that reading and parsing happen in parallel.
*   The caller decides which VCF fields to parse. Fields that are not requested
    by the caller are skipped and not parsed.
*   Very few bytes in the input buffer are accessed more than once.
*   Memory is allocated frugally and reused whenever possible.
*   Exceptions are not used for error reporting.
*   The library is header-only with no dependencies outside the standard
    library.

## How to use

### Preparation and parsing the header

1.  Allocate a generous amount of memory for the input buffer. The buffer must
    not be changed or deleted after it's been passed to the parser using the
    `feed()` method (see below).

        char buffer[1024 * 1024];

2.  Create an instance of the parser.

        VCF_scanner vcf_scanner;

3.  Implement a function to read the input stream and feed the data into the
    parser.

        auto read_and_feed = [&]() -> VCF_parsing_event {
            size_t bytes_read;

            if (vcf_stream->at_eof())
                // The parser treats an empty input buffer as an indicator
                // that end of file has been reached.
                // Your stream reading function may already follow this
                // convention.
                bytes_read = 0;
            else
                bytes_read = vcf_stream->read(buffer, sizeof(buffer));

            // Continue parsing the token whose parsing was suspended
            // because the previous buffer was depleted.
            // Return the result of parsing with additional data.
            return vcf_scanner.feed(buffer, bytes_read);
        };

4.  Implement a function that repeatedly feeds the parser until the current
    token is completely parsed.

        auto parse_to_completion = [&](VCF_parsing_event pe) {
            while (pe == VCF_parsing_event::need_more_data)
                pe = read_and_feed();

            if (pe == VCF_parsing_event::error)
                terminate_parsing_with_error(vcf_scanner.get_error(),
                        vcf_scanner.get_line_number());

            if (pe == VCF_parsing_event::ok_with_warnings)
                show_warnings(vcf_scanner.get_warnings());
        };

5.  Read the initial input buffer and use the completion function to parse the
    header.

        parse_to_completion(read_and_feed());

        // The header has been successfully parsed and and can be
        // accessed using the get_header() method.

6.  Do not discard the current input buffer - it contains the beginning of the
    data lines.  Proceed to the data line parsing loop.

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
