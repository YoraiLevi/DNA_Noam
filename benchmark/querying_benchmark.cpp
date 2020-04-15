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
#define PROGRAM_NAME "QUERY-BENCHMARK"
#define VERSION "0.0.1"
#include <seqan3/core/detail/all.hpp>
struct cmd_arguments
{
    std::filesystem::path head_index_path{};
    std::vector<std::string> queries{};

};
void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
{
    parser.info.author = "Yorai Levi";
    parser.info.short_description = "Benchmarks peformance of cursor and search api of seqan3";
    parser.info.version = VERSION;
    parser.add_option(args.head_index_path, 'x', "index", "index of fasta input"); 
    parser.add_positional_option(args.queries, "List of queries to search");
}
auto query_cursor(std::vector<sequence_type>& texts,fm_index_t &index){
    std::vector<size_type> v;
    for(auto &text : texts)
{            cursor_t cur = index.cursor();
        cur.extend_right(text);
        v.push_back(cur.count());}
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
    for (auto &&r : search(texts, index,search_config)){
        v[std::get<size_type>(r)]++;
    }
    return v;
}
int main(int argc, char **argv)
{
    cmd_arguments args = parse_args<cmd_arguments>(PROGRAM_NAME,argc, argv);
    std::vector<sequence_type> queries{};
    queries.reserve(args.queries.size());
    for(auto &str: args.queries){
        queries.push_back(str_to_dna5(str));
        std::cout << str<<std::endl;
    }
    
    std::cout << "Loading..."<<std::endl;
    fm_index_t head_idx;
    {
        load_fm(head_idx,args.head_index_path);
    }
    std::cout << "Done." << std::endl;
    std::cout << "Processing...\n";
    {   
        {
            std::cout<<"Starting Cursor Test"<<std::endl;
            std::chrono::time_point<std::chrono::system_clock> past = std::chrono::high_resolution_clock::now();
            auto range = query_cursor(queries,head_idx);
            // seqan3::debug_stream<<range;
            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::high_resolution_clock::now();
            auto diffnans = std::chrono::duration_cast<std::chrono::nanoseconds>(now-past);
            auto diffmics = std::chrono::duration_cast<std::chrono::microseconds>(now-past);
            auto diffmils = std::chrono::duration_cast<std::chrono::milliseconds>(now-past);
            std::cout<<"Counts: ";
            for(auto &count : range){
                std::cout<<count<<" ";
            }
            std::cout<<std::endl;
            std::cout<<
            "milliseconds: "<<diffmils.count()<<
            " ;; microseconds: "<<diffmics.count()<<
            " ;; nanoseconds: "<< diffnans.count()<<
            std::endl;
        std::cout<<"Done Cursor Test"<<std::endl;
        }
        std::cout<<"\n\n"<<std::endl;
        {
            std::cout<<"Starting Search Test"<<std::endl;

            std::chrono::time_point<std::chrono::system_clock> past = std::chrono::high_resolution_clock::now();
            auto range = query_search(queries,head_idx);
            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::high_resolution_clock::now();
            auto diffnans = std::chrono::duration_cast<std::chrono::nanoseconds>(now-past);
            auto diffmics = std::chrono::duration_cast<std::chrono::microseconds>(now-past);
            auto diffmils = std::chrono::duration_cast<std::chrono::milliseconds>(now-past);
            std::cout<<"Counts: ";
            for(auto &count : range){
                std::cout<<count<<" ";
            }
            std::cout<<std::endl;
            std::cout<<
            "milliseconds: "<<diffmils.count()<<
            " ;; microseconds: "<<diffmics.count()<<
            " ;; nanoseconds: "<< diffnans.count()<<
            std::endl;
            std::cout<<"Done Search Test"<<std::endl;
        }   
                    std::cout<<"\n\n"<<std::endl;

        {
            std::cout<<"Starting Search Parallel Test"<<std::endl;

            std::chrono::time_point<std::chrono::system_clock> past = std::chrono::high_resolution_clock::now();
            auto range = query_search_parallel(queries,head_idx);
            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::high_resolution_clock::now();
            auto diffnans = std::chrono::duration_cast<std::chrono::nanoseconds>(now-past);
            auto diffmics = std::chrono::duration_cast<std::chrono::microseconds>(now-past);
            auto diffmils = std::chrono::duration_cast<std::chrono::milliseconds>(now-past);
            std::cout<<"Counts: ";
            for(auto &count : range){
                std::cout<<count<<" ";
            }
            std::cout<<std::endl;
            std::cout<<
            "milliseconds: "<<diffmils.count()<<
            " ;; microseconds: "<<diffmics.count()<<
            " ;; nanoseconds: "<< diffnans.count()<<
            std::endl;
            std::cout<<"Done Search Parallel Test"<<std::endl;

        }
    }
    std::cout << "Done." << std::endl;
    std::cout << "Closing stream...";
    std::cout << "Done." << std::endl;

    std::cout << "All Complete. GZ GZ" << std::endl;
}