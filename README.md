# VCF_scanner - a blazing fast library to parse large Variant Call Format files

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

This section is a step-by-step guide on how to use the library.  Additionally,
the `examples` directory contains compilable code that illustrates key
concepts, and the main header (`include/vcf_scanner/vcf_scanner.hh`) has plenty
of comments that provide further details.

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

4.  Implement a function to check the result of parsing a token and call the
    above function until the token is completely parsed.

        auto parse_to_completion = [&](VCF_parsing_event pe) {
            while (pe == VCF_parsing_event::need_more_data)
                pe = read_and_feed();

            if (pe == VCF_parsing_event::error)
                throw Parsing_error(vcf_scanner.get_error(),
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

All data line fields are optional and can be skipped by not calling the
respective `parse_...()` methods. However, the `clear_line()` method must be
called at the end of each line regardless of whether any fields were skipped.

The fields must be extracted in the exact order as defined by the VCF
specification.

1.  Define the variables for the field values outside the data reading loop to
    reduce unnecessary memory reallocation. Every variable in this list is
    valid only immedialtely after the respective `parse_...()` method returned
    `ok` or `ok_with_warnings`.

        std::string chrom;
        unsigned pos;
        std::vector<std::string> ids;
        std::string ref;
        std::vector<std::string> alts;
        std::string quality_str;
        bool quality_is_missing;
        std::vector<std::string> filters;

2.  Repeat until there are no more data lines left to read.

        while (!vcf_scanner.at_eof()) {

3.  Request parsing of the CHROM and POS fields by calling `parse_loc()`.

            parse_to_completion(vcf_scanner.parse_loc(&chrom, &pos));

4.  Parse the ID field.

            parse_to_completion(vcf_scanner.parse_ids(&ids));

5. The REF and ALT alleles are parsed by the `parse_alleles()` method.

            parse_to_completion(vcf_scanner.parse_alleles(&ref, &alts));

6.  Parse the QUAL field.

            parse_to_completion(vcf_scanner.parse_quality(
                &quality_str, &quality_is_missing));

            if (!quality_is_missing) {
                float quality = std::stof(quality_str);
            }

7.  Parse the FILTER field.

            parse_to_completion(vcf_scanner.parse_filters(&filters));

8.  Parse the INFO field.

            parse_to_completion(vcf_scanner.parse_info());

            std::vector<std::string> info = vcf_scanner.get_info();

9.  Parse genotype info.

            if (vcf_scanner.get_header().has_genotype_info()) {
                parse_to_completion(vcf_scanner.parse_genotype_format());

                if (vcf_scanner.capture_gt()) {
                    while (vcf_scanner.genotype_available()) {
                        parse_to_completion(vcf_scanner.parse_genotype());

                        for (auto allele_index : vcf_scanner.get_gt()) {
                            // 1. The MISSING value is represented by
                            //    allele_index < 0.
                            // 2. vcf_scanner.is_phased_gt() will return
                            //    true if the genotype is phased.
                        }
                    }
                }
            }

10. Skip to the next line by calling `clear_line()`.

            parse_to_completion(vcf_scanner.clear_line());
        }

## To build a test coverage report

1.  Install `lcov`
2.  Configure, build, and run tests as follows:

        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_COVERAGE=TRUE ..
        make coverage

3.  Find the HTML report in the `coverage` directory:

        xdg-open coverage/index.html
