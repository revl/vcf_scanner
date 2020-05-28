## To build a test coverage report

1. Install `lcov`
2. Configure as follows:

        mkdir build
        cd build
        cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_COVERAGE=TRUE ..
        make coverage
        xdg-open coverage/index.html
