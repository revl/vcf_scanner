#ifndef TEST_CASE_H
#define TEST_CASE_H

#include <list>
#include <stdexcept>

class test_case
{
public:
    virtual ~test_case() {}

    test_case(const char* name) : test_name(name), failed_checks(0)
    {
        test_case::test_case_list.push_back(this);
    }

    virtual void run() const = 0;

    const char* const test_name;
    int failed_checks;

    static std::list<test_case*> test_case_list;
};

class abort_test_case : public exception
{
};

std::list<test_case*> test_case::test_case_list;

test_case* current_test_case = NULL;

#define TEST_CASE(class_name)                                                  \
    class test_##class_name : public test_case                                 \
    {                                                                          \
    public:                                                                    \
        test_##class_name(const char* name) : test_case(name) {}               \
        virtual void run() const;                                              \
    } static test_##class_name##_instance(#class_name);                        \
    void test_##class_name::run() const

#define CHECK(condition)                                                       \
    do                                                                         \
        if (!(condition)) {                                                    \
            fprintf(stderr,                                                    \
                    "%s:%d: error: "                                           \
                    "check \"%s\" failed\n",                                   \
                    __FILE__, __LINE__, #condition);                           \
            ++current_test_case->failed_checks;                                \
        }                                                                      \
    while (0)

#define REQUIRE(condition)                                                     \
    do                                                                         \
        if (!(condition)) {                                                    \
            fprintf(stderr,                                                    \
                    "%s:%d: error: required "                                  \
                    "check \"%s\" failed\n",                                   \
                    __FILE__, __LINE__, #condition);                           \
            throw abort_test_case();                                           \
        }                                                                      \
    while (0)

int main(int /*argc*/, char* /*argv*/ [])
{
    int failed_tests = 0;

    for (auto tc : test_case::test_case_list) {
        current_test_case = tc;
        try {
            current_test_case->run();

            if (current_test_case->failed_checks > 0)
                ++failed_tests; // LCOV_EXCL_LINE
        }
        // LCOV_EXCL_START
        catch (abort_test_case&) {
            ++failed_tests;
        }
        // LCOV_EXCL_STOP
    }

    return failed_tests;
}

#endif /* !defined(TEST_CASE_H) */
