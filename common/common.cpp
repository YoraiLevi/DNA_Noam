#include <string>
#include <vector>
#include <seqan3/argument_parser/all.hpp> 
#include <seqan3/core/debug_stream.hpp>      // for debug_stream
#include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/search/fm_index/all.hpp>

#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include "common.h"


//exclude //include //iterator //ostream save to file
void print_cur(const cursor_t &cur, text_t &texts)
{
    // if(texts == nullptr) return;
    seqan3::debug_stream
        /* << " last rank: " << cur.last_rank()
        << " count: " << cur.count()
        << " len: " << cur.query_length()
        << " locate: " << cur.locate()
        << " str: "*/
        << cur.path_label(texts)
        << '\n';
}
template <class Callable>
void preorder_cur(Callable f,auto &defaults, const cursor_t &cur, text_t &texts, int level, const int max_depth)
{
    //the correct thing to do is probably writing an iterator and ++ing the iterator
    #ifndef NDEBUG
    print_cur(cur, texts);
    #endif
    //f is callable of the structure: f(const custom_struct_t &defaults, const cursor_t &cur, text_t &texts)
    f(defaults,cur,texts);
    if(level>max_depth) return;
    cursor_t temp_cur(cur); //copying
    if (!temp_cur.extend_right())
    { //if deepest then done.
        return;
    }
    do
    {
        preorder_cur(f,defaults,temp_cur, texts,level+1,max_depth);
    } while (temp_cur.cycle_back());
}
template <class Callable>
void apply_cur(Callable f,auto &defaults, const fm_index_t &index, text_t &texts, const int max_depth)
{
    auto cur = index.cursor(); //moves to level=1
    preorder_cur(f,defaults,cur,texts,1,max_depth);
}
void load_fm(fm_index_t &index, const std::filesystem::path path = INDEX_FILE_NAME)
{
  // std::filesystem::path path = INDEX_FILE_NAME;
      {
        std::ifstream is{path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
    }
}
void store_fm(const fm_index_t&index,const std::filesystem::path path = INDEX_FILE_NAME,const bool is_verify=false)
{   
    // std::filesystem::path path = INDEX_FILE_NAME;
    {
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
    }
    if(is_verify)
    {
        fm_index_t index2;
        load_fm(index2,path);
    if (index == index2)
        std::cout << "The indices are identical!\n";
    else
        std::cout << "The indices differ!\n";
    }
    
}


void store_meta(const std::vector<id_type> &ids, const std::filesystem::path path = META_FILE_NAME)
{
    {std::ofstream os{path};
    cereal::JSONOutputArchive oarchive{os};
    oarchive(CEREAL_NVP(ids));}
}

void load_meta(std::vector<id_type> &ids, const std::filesystem::path path = META_FILE_NAME)
{
        #ifndef NDEBUG
    std::cout << path<<std::endl;
    #endif
    {std::ifstream is{path};
    cereal::JSONInputArchive iarchive{is};
    iarchive(CEREAL_NVP(ids));
}
}

void load_records(const std::vector<std::filesystem::path> &file_paths, std::vector<sequence_type> &sequences, std::vector<id_type> &ids, const bool meta_only=false)
{
    for (auto &file_path : file_paths)
    {
        seqan3::sequence_file_input fin{file_path};
        for (auto &[seq, id, qual] : fin)
        {
            #ifndef NDEBUG
            seqan3::debug_stream <<"seq found:"<< id << std::endl;
            #endif
            if(!meta_only){
            sequences.push_back(std::move(seq));    
            }
            ids.push_back(std::move(id));
        }
    }
}
template<typename cmd_arguments>
cmd_arguments parse_args(auto name,int argc, char **argv)
{
    //throws -1 if failed
    seqan3::argument_parser myparser{name, argc, argv}; // initialise myparser
    cmd_arguments args{};
    initialise_argument_parser(myparser, args);

    //process argc argv into args
    try
    {
        myparser.parse(); // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const &ext) // catch user errors
    {
        seqan3::debug_stream << "[Error Message:] " << ext.what() << "\n"; // customise your error message
        throw ext;
    }
    return args;
}
class breadth_cur{
  public:
  breadth_cur(cursor_t root){
    queue.push(root);
  }
  cursor_t next(){
    cursor_t cur(this->queue.front());
    // temp_cur(cur);
    cursor_t& temp_cur = this->queue.front();
    if(temp_cur.extend_right()){
            //if has children
            do{
                //appending node children
                queue.push(temp_cur);
            } while(temp_cur.cycle_back());   
            }
    this->queue.pop();
    return cur;
  }
  bool empty(){
    return this->queue.empty();
  }
  private:  
    std::queue<cursor_t> queue;
};
class minmax_breadth_cur : public breadth_cur{
  public:
  minmax_breadth_cur(cursor_t root,fm_index_t::size_type min=-1,fm_index_t::size_type max=-1): breadth_cur(root), _min{min}, _max{max}{}
  cursor_t next(){
    cursor_t cur(this->queue.front());
    // temp_cur(cur);
    cursor_t& temp_cur = this->queue.front();
    if(temp_cur.extend_right()){
            //if has children
            do{
                //appending node children
                queue.push(temp_cur);
            } while(temp_cur.cycle_back());   
            }
    this->queue.pop();
    return cur;
  }
  bool empty(){
    return this->queue.empty();
  }
  private:  
    std::queue<cursor_t> queue;
    fm_index_t::size_type _min;
    fm_index_t::size_type _max;
};
// template <std::ranges::range text_t>
fm_index_t::size_type has_node(const auto &idx,auto &label){
    cursor_t cur{idx.cursor()};
    if(cur.extend_right(label))
        return cur.count();
    return 0;
}