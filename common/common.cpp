#include <string>
#include <vector>
#include <seqan3/argument_parser/all.hpp> 
#include <seqan3/core/debug_stream.hpp>      // for debug_stream
#include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/range/views/char_to.hpp>


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
    using size_type = fm_index_t::size_type;
    using cursor_type =  fm_index_t::cursor_type;
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

    size_type depth(){
    return this->queue.front().query_length();
  }
    bool fast_forward(size_type depth){
    //Moves internal head to the first node in {depth} depth
    size_type dept_atm(this->depth());
    if(dept_atm>depth) return false;
    while(!this->empty() && depth>this->depth()){
        this->next();
    }
    return true;
  }
  bool empty(){
    return this->queue.empty();
  }
  cursor_type front(){
    return this->queue.front();
  }
  private:  
    std::queue<cursor_t> queue;
};
// class depth_breadth_cur{
//   public:
//     using size_type = fm_index_t::size_type;
//     using cursor_type =  fm_index_t::cursor_type;
//     using value_type = std::tuple<size_type,cursor_type>;
//   depth_breadth_cur(cursor_t root){
//     this->queue.push(value_type(0,root));
//   }
//     depth_breadth_cur(value_type root){
//         this->queue.push(root);
//   }
//   size_type depth(){
//     return std::get<size_type>(this->queue.front());
//   }
//   bool fast_forward(size_type depth){
//     //Moves internal head to the first node in {depth} depth
//     size_type dept_atm(this->depth());
//     if(dept_atm>depth) return false;
//     while(!this->end() && depth>this->depth()){
//         this->next();
//     }
//     return true;
//   }
// value_type next(){
//     value_type cur_tuple{this->queue.front()}; //hopefully copy constructor;
//     value_type& temp_tuple = this->queue.front(); //work reference;
//     cursor_type& cur = std::get<cursor_type>(temp_tuple);
//     size_type depth = std::get<size_type>(temp_tuple);
//     if(cur.extend_right()){
//         do{
//             queue.push(value_type(depth,cur));
//         }while(cur.cycle_back());
//     }
//     this->queue.pop();
//     return cur_tuple;
// }

//   //   cursor_t cur(this->queue.front());
//   //   // temp_cur(cur);
//   //   cursor_t& temp_cur = this->queue.front();
//   //   if(temp_cur.extend_right()){
//   //           //if has children
//   //           do{
//   //               //appending node children
//   //               queue.push(temp_cur);
//   //           } while(temp_cur.cycle_back());   
//   //           }
//   //   this->queue.pop();
//   //   return cur;
//   // }
//   bool end(){
//     return this->queue.empty();
//   }
//   private:  
//     std::queue<value_type> queue;
// };
// template <std::ranges::range text_t>
fm_index_t::size_type has_node(const auto &idx,auto &label){
    cursor_t cur{idx.cursor()};
    if(cur.extend_right(label))
        return cur.count();
    return 0;
}

sequence_type str_to_dna5(const std::string &str){
    return str | seqan3::views::char_to<seqan3::dna5>;
}