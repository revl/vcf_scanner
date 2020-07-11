// This header contains implementation details.
// It is not meant to be included directly.

#include <string>
#include <map>
#include <vector>

#include "string_view.hh"

class VCF_scanner_impl;

class VCF_header_impl
{
protected:
    void add_meta_info(
            std::string&& meta_info_key, const VCF_string_view& meta_info_line)
    {
        meta_info[meta_info_key].push_back(meta_info_line);
    }

    std::string file_format_version;
    std::map<std::string, std::vector<std::string>> meta_info;
    bool genotype_info_present = false;
    std::vector<std::string> sample_ids;

    friend class VCF_scanner_impl;
};
