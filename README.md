## Features

* This parser leaves it to the caller to load the input data into memory
  buffers. As a result, the parser is non-blocking. The caller can choose to
  read input data in a separate thread so that reading and parsing happens in
  parallel.
* The caller decides which VCF fields to parse. Fields that present no interest
  to the caller are skipped and not parsed.
* Very few bytes in the input buffer are accessed more than once.
* Memory is allocated frugally.
* Exceptions are not used for error reporting.
* The library is header-only with no external dependencies.

## To build a test coverage report

1.  Install `lcov`
2.  Configure, build, and run tests as follows:

        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_COVERAGE=TRUE ..
        make coverage

3.  Find the HTML report in the `coverage` directory:

        xdg-open coverage/index.html
