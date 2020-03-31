#pragma once
//std
#include <vector>
//seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/search/fm_index/all.hpp>
//types
using fm_index_t = seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection >;
using cursor_t = fm_index_t::cursor_type;
using text_t = std::vector<seqan3::dna5_vector>;
using id_type = seqan3::sequence_file_input<>::id_type;
using sequence_type = seqan3::sequence_file_input<>::sequence_type;
using size_type = fm_index_t::size_type;
//funcs
void print_cur(const cursor_t &cur, const text_t &texts);
//f is callable of the structure: f(const custom_struct_t &defaults, const cursor_t &cur, const std::vector<auto> &texts)
template <class Callable>
void preorder_cur(Callable f,auto &defaults, const cursor_t &cur, const text_t &texts, int level, const int max_depth =-1);
template <class Callable>
void apply_cur(Callable f,auto &defaults, const fm_index_t &index, const text_t &texts, const int max_depth = -1);
#define META_FILE_NAME "index.meta.json"
// void store_meta(const id_type &ids);
// void load_meta(id_type &ids);
#define INDEX_FILE_NAME "index.bin"
// void load_fm(fm_index_t &index, const std::filesystem::path path = INDEX_FILE_NAME);
// void store_fm(const fm_index_t&index,const std::filesystem::path path = INDEX_FILE_NAME,const bool is_verify=false);

// void store_fm(const fm_index_t &index,const bool is_verify);
// void load_fm(fm_index_t &index);