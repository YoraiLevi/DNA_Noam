#include <string>
#include <fstream>
#include <iostream>
// #include <iostream>
#include <vector>
#include <algorithm>

#include <seqan3/search/configuration/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/argument_parser/all.hpp> 
#include <seqan3/search/fm_index/fm_index.hpp>

#include <common.cpp>
#define PROGRAM_NAME "BREATH-ITERATION-BENCHMARK"
#define VERSION "0.0.1"
#include <seqan3/core/detail/all.hpp>
class Timer{
public:
    Timer(){
        start = std::chrono::high_resolution_clock::now();
    }
    auto Stop()
    {
        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::high_resolution_clock::now();
        auto diffnans = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start);
        auto diffmics = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
        auto diffmils = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
        std::cout<<
        "milliseconds: "<<diffmils.count()<<
        " ;; microseconds: "<<diffmics.count()<<
        " ;; nanoseconds: "<< diffnans.count()<<
        std::endl;        
        return diffnans.count();
    }
    // template<typename T=std::chrono::microseconds>
     auto deltaMS()
    {
        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::high_resolution_clock::now();
        auto diffnans = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
        return diffnans.count();
    }
    virtual ~Timer(){
        Stop();
    }
private:
    std::chrono::time_point<std::chrono::system_clock> start;
};
struct cmd_arguments
{
    std::filesystem::path index_path{};
    std::vector<std::filesystem::path> fasta_path{};

    std::filesystem::path index2_path{};
    size_type max_depth{0};
    size_type min_depth{0};

};
void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
{
    parser.info.author = "Yorai Levi";
    parser.info.short_description = "Benchmarks different querying methods";
    parser.info.version = VERSION;
    parser.add_option(args.index_path, 'x', "index", "index of fasta input");
    parser.add_option(args.index2_path, 'y', "index2", "index2 of fasta input");
    
    parser.add_option(args.fasta_path, 'f', "file", "fasta1 input"); 
    parser.add_option(args.max_depth, '\0', "max", "max depth to search");
    parser.add_option(args.min_depth, '\0', "min", "min depth to search");

}
auto query_cursor(std::vector<sequence_type>& texts,fm_index_t &index){
    std::vector<size_type> v;
    for(auto &text : texts)
    {
        cursor_t cur = index.cursor();
        cur.extend_right(text);
        v.push_back(cur.count());
    }
return v;
}
auto query_search(std::vector<sequence_type>& texts,fm_index_t &index){
    std::vector<size_type> v(texts.size());      
    for (auto &&r : search(texts, index)){
        v[std::get<size_type>(r)]++;
    }
    return v;
}
auto query_search_parallel(std::vector<sequence_type>& texts,fm_index_t &index){
seqan3::configuration const search_config = 
seqan3::search_cfg::parallel{8} |
seqan3::search_cfg::max_error{seqan3::search_cfg::total{0},
seqan3::search_cfg::substitution{0},
seqan3::search_cfg::insertion{0},
seqan3::search_cfg::deletion{0}} |
seqan3::search_cfg::output{seqan3::search_cfg::text_position} |
seqan3::search_cfg::mode{seqan3::search_cfg::all};

    std::vector<size_type> v(texts.size());      
    for (auto &&r : seqan3::search(texts, index,search_config)){
        v[std::get<size_type>(r)]++;
    }
    return v;
}
template <typename T, std::size_t N = 5>
class Queue { 
    using block = std::array<T, N>;
    using iterator = typename block::iterator;

    public:
    Queue();
    void push(const T& value);
    void pop();
    T front() const;
    bool empty();

    private:
    std::list<block> storage;
    iterator first, last;
};

template <typename T, std::size_t N>
 bool Queue<T, N>::empty(){ 
    return first==last;
    // storage.push_back(block{});
    // first = std::begin(*std::begin(storage));
    // last  = first; // last == first <=> empty queue
}

template <typename T, std::size_t N>
Queue<T, N>::Queue() {
    storage.push_back(block{});
    first = std::begin(*std::begin(storage));
    last  = first;
}

template <typename T, std::size_t N>
void Queue<T, N>::push(const T& value) {
    if (last == std::end(*std::rbegin(storage))) {
        // std::cout << "requiring new storage" << '\n';
        storage.push_back(block{});
        last = std::begin(*std::rbegin(storage));
    }
    *last++ = value;
}

template <typename T, std::size_t N>
void Queue<T, N>::pop() {
    if (++first == std::end(*std::begin(storage))) {
        // std::cout << "freeing unused memory" << '\n';
        storage.pop_front();
        first = std::begin(*std::begin(storage));
    }
}

template <typename T, std::size_t N>
T Queue<T, N>::front() const {
    return *first;
}
class breadth_cur_q{
  public:
    using size_type = fm_index_t::size_type;
    using cursor_type =  fm_index_t::cursor_type;
  breadth_cur_q(cursor_t root){
    queue.push(root);
  }
  cursor_t next(){
    cursor_t cur(this->queue.front());
    // temp_cur(cur);
    cursor_t temp_cur = this->queue.front();
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
    Queue<cursor_t,100000> queue;
};


int main(int argc, char **argv)
{
    cmd_arguments args = parse_args<cmd_arguments>(PROGRAM_NAME,argc, argv);    
    std::cout << "Loading..."<<std::endl;
    fm_index_t idx;
    fm_index_t idx2;
    {
        load_fm(idx,args.index_path);
        load_fm(idx2,args.index2_path);
    }
    std::vector<sequence_type> texts{};
    {
    std::vector<id_type> ids{};
    load_records(args.fasta_path,texts,ids,false);
    }

    std::cout << "Done." << std::endl;
    std::cout << "Initing Memory..."<<std::endl;
    breadth_cur_q queue(idx);
    std::vector<sequence_type> v{};
    v.reserve(1000000);
    {
    Timer timer{};
    queue.fast_forward(args.min_depth);
    while(timer.deltaMS()<30*1000 && (args.max_depth==0 || queue.depth()<args.max_depth) && !queue.empty())
        v.push_back(queue.next().path_label(texts));
    }
    // for(auto &label:v)
    //     seqan3::debug_stream << label;
    std::cout << "Size: "<<v.size()<<";; depth: "<<queue.depth() << std::endl;

    std::cout << "Done." << std::endl;
    
    std::cout << "Processing...\n";
    {   

        std::cout<<"Starting extend_right:"<<std::endl;
        {Timer t{};
            std::vector<bool> counts(v.size());
            size_type i=0;
            for(auto &label:v){
                counts[i]=idx2.cursor().extend_right(label);
                i++;
            }
        }
        std::cout<<"Ended extend_right."<<std::endl;

        std::cout<<"Starting extend_right count:"<<std::endl;
        {Timer t{};
            std::vector<size_type> counts(v.size());
            size_type i=0;      
            for(auto &label:v){
                auto cur=idx2.cursor();
                cur.extend_right(label);
                counts[i]=cur.count();
                i++;
                }
        }
        std::cout<<"Ended extend_right count."<<std::endl;


        std::cout<<"Starting search single label:"<<std::endl;
        {Timer t{};
        std::vector<size_type> counts(v.size());      
        for(auto &label:v)
            for (auto &&r : seqan3::search(label, idx2)){
                counts[std::get<size_type>(r)]++;
            }
        }
        std::cout<<"Ended search single label"<<std::endl;
        

        std::cout<<"Starting search single label bool:"<<std::endl;
        {Timer t{};
        std::vector<bool> counts(v.size());      
        for(auto &label:v)
            for (auto &&r : seqan3::search(label, idx2)){
                counts[std::get<size_type>(r)]=true;
                break;
            }
        }
        std::cout<<"Ended search single label bool"<<std::endl;

        std::cout<<"Starting search all label:"<<std::endl;
        {Timer t{};
        std::vector<size_type> counts(v.size());      
        // for(auto &label:v)
            for (auto &&r : seqan3::search(v, idx2)){
                counts[std::get<size_type>(r)]++;
            }
        }
        std::cout<<"Ended search all label"<<std::endl;
        
        std::cout<<"Starting search all label bool:"<<std::endl;
        {Timer t{};
        std::vector<bool> counts(v.size());      
        // for(auto &label:v)
            for (auto &&r : seqan3::search(v, idx2)){
                counts[std::get<size_type>(r)]=true;
                // break;
            }
        }
        std::cout<<"Ended search all label bool"<<std::endl;

        seqan3::configuration const search_config = 
        seqan3::search_cfg::parallel{8} |
        seqan3::search_cfg::max_error{seqan3::search_cfg::total{0},
        seqan3::search_cfg::substitution{0},
        seqan3::search_cfg::insertion{0},
        seqan3::search_cfg::deletion{0}} |
        seqan3::search_cfg::output{seqan3::search_cfg::text_position} |
        seqan3::search_cfg::mode{seqan3::search_cfg::all};
        std::cout<<"Starting search parallel single label:"<<std::endl;
        {Timer t();
        std::vector<size_type> counts(v.size());      
        for(auto &label:v)
            for (auto &&r : seqan3::search(label, idx2,search_config)){
                counts[std::get<size_type>(r)]++;
            }
        }
        std::cout<<"Ended search parallel single label"<<std::endl;
        
        std::cout<<"Starting search parallel single label bool:"<<std::endl;
        {Timer t{};
        std::vector<bool> counts(v.size());      
        for(auto &label:v)
            for (auto &&r : seqan3::search(label, idx2,search_config)){
                counts[std::get<size_type>(r)]=true;
                break;
            }
        }
        std::cout<<"Ended search parallel single label bool"<<std::endl;

        std::cout<<"Starting search parallel all label:"<<std::endl;
        {Timer t{};
        std::vector<size_type> counts(v.size());      
        // for(auto &label:v)
            for (auto &&r : seqan3::search(v, idx2,search_config)){
                counts[std::get<size_type>(r)]++;
            }
        }
        std::cout<<"Ended search parallel all label"<<std::endl;
        
        std::cout<<"Starting search parallel all label bool:"<<std::endl;
        {Timer t{};
        std::vector<bool> counts(v.size());      
        // for(auto &label:v)
            for (auto &&r : seqan3::search(v, idx2,search_config)){
                counts[std::get<size_type>(r)]=true;
                // break;
            }
        }
        std::cout<<"Ended search parallel all label bool"<<std::endl;

    }
    std::cout << "Done." << std::endl;
    std::cout << "Closing stream...";
    std::cout << "Done." << std::endl;

    std::cout << "All Complete. GZ GZ" << std::endl;
}
        // {
        //     std::cout<<"Starting Cursor Test"<<std::endl;
        //     std::chrono::time_point<std::chrono::system_clock> past = std::chrono::high_resolution_clock::now();
        //     auto range = query_cursor(queries,head_idx);
        //     // seqan3::debug_stream<<range;
        //     std::chrono::time_point<std::chrono::system_clock> now = std::chrono::high_resolution_clock::now();
        //     auto diffnans = std::chrono::duration_cast<std::chrono::nanoseconds>(now-past);
        //     auto diffmics = std::chrono::duration_cast<std::chrono::microseconds>(now-past);
        //     auto diffmils = std::chrono::duration_cast<std::chrono::milliseconds>(now-past);
        //     std::cout<<"Counts: ";
        //     for(auto &count : range){
        //         std::cout<<count<<" ";
        //     }
        //     std::cout<<std::endl;
        //     std::cout<<
        //     "milliseconds: "<<diffmils.count()<<
        //     " ;; microseconds: "<<diffmics.count()<<
        //     " ;; nanoseconds: "<< diffnans.count()<<
        //     std::endl;
        // std::cout<<"Done Cursor Test"<<std::endl;
        // }
        // std::cout<<"\n\n"<<std::endl;
        // {
        //     std::cout<<"Starting Search Test"<<std::endl;

        //     std::chrono::time_point<std::chrono::system_clock> past = std::chrono::high_resolution_clock::now();
        //     auto range = query_search(queries,head_idx);
        //     std::chrono::time_point<std::chrono::system_clock> now = std::chrono::high_resolution_clock::now();
        //     auto diffnans = std::chrono::duration_cast<std::chrono::nanoseconds>(now-past);
        //     auto diffmics = std::chrono::duration_cast<std::chrono::microseconds>(now-past);
        //     auto diffmils = std::chrono::duration_cast<std::chrono::milliseconds>(now-past);
        //     std::cout<<"Counts: ";
        //     for(auto &count : range){
        //         std::cout<<count<<" ";
        //     }
        //     std::cout<<std::endl;
        //     std::cout<<
        //     "milliseconds: "<<diffmils.count()<<
        //     " ;; microseconds: "<<diffmics.count()<<
        //     " ;; nanoseconds: "<< diffnans.count()<<
        //     std::endl;
        //     std::cout<<"Done Search Test"<<std::endl;
        // }   
        //             std::cout<<"\n\n"<<std::endl;

        // {
        //     std::cout<<"Starting Search Parallel Test"<<std::endl;

        //     std::chrono::time_point<std::chrono::system_clock> past = std::chrono::high_resolution_clock::now();
        //     auto range = query_search_parallel(queries,head_idx);
        //     std::chrono::time_point<std::chrono::system_clock> now = std::chrono::high_resolution_clock::now();
        //     auto diffnans = std::chrono::duration_cast<std::chrono::nanoseconds>(now-past);
        //     auto diffmics = std::chrono::duration_cast<std::chrono::microseconds>(now-past);
        //     auto diffmils = std::chrono::duration_cast<std::chrono::milliseconds>(now-past);
        //     std::cout<<"Counts: ";
        //     for(auto &count : range){
        //         std::cout<<count<<" ";
        //     }
        //     std::cout<<std::endl;
        //     std::cout<<
        //     "milliseconds: "<<diffmils.count()<<
        //     " ;; microseconds: "<<diffmics.count()<<
        //     " ;; nanoseconds: "<< diffnans.count()<<
        //     std::endl;
        //     std::cout<<"Done Search Parallel Test"<<std::endl;

        // }
