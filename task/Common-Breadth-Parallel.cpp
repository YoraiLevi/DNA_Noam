#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <omp.h>
#include <array>
#include <stdio.h>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>  


#include <seqan3/argument_parser/all.hpp> 
#include <seqan3/search/fm_index/fm_index.hpp>
#include <common.cpp>
#define PROGRAM_NAME "Common-Breath-Serial"
#define VERSION "0.0.1"
#include <seqan3/core/detail/all.hpp>
struct cmd_arguments
{
    //input files
    std::filesystem::path dog{};
    std::filesystem::path wolf{};
    std::filesystem::path chimp{};
    std::filesystem::path human{};
    std::vector<std::filesystem::path> file_paths{};

    //search configuration
    size_type min_depth{0};
    size_type max_depth{0};

    //output configuration
    std::filesystem::path outTSV{"out.tsv"};

};
void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
{
    parser.info.author = "Yorai Levi";
    parser.info.short_description = "FM-INDEXES fasta formatted files";
    parser.info.version = VERSION;

    parser.add_option(args.max_depth, '\0', "max", "max depth to search");
    parser.add_option(args.min_depth, '\0', "min", "min depth to search");

    parser.add_option(args.dog, 'd', "dog", "dog index");
    parser.add_option(args.wolf, 'w', "wolf", "wolf index");
    parser.add_option(args.chimp, 'c', "chimp", "chimp index");
    parser.add_option(args.human, 'u', "human", "human index");

    parser.add_option(args.file_paths, 'f', "file", "original fasta input files that were indexed");
    // parser.add_option(args.tile_paths, 't', "tile", "original fasta input files that were indexed");
    parser.add_option(args.outTSV, '\0', "out-file-name", "name of file to output, as for this version please include .tsv as a file extension");
}

int main(int argc, char **argv)
{
    cmd_arguments args = parse_args<cmd_arguments>(PROGRAM_NAME,argc, argv);
    // fm_index_t idx;
    fm_index_t dog_idx,wolf_idx,chimp_idx,human_idx;
    // load_fm(idx);
    std::cout << "Loading...";
    load_fm(dog_idx,args.dog);
    load_fm(wolf_idx,args.wolf);
    load_fm(chimp_idx,args.chimp);
    load_fm(human_idx,args.human);

    std::vector<sequence_type> texts{};
    std::vector<id_type> ids{};
    load_records(args.file_paths,texts,ids,false);

    std::cout << "Done." << std::endl;

    std::cout << "Prepearing stream for writing...";
    
    std::ofstream os;
    seqan3::debug_stream_type<> deos{os};
    os.open(args.outTSV);
    deos << "Human Count" <<"\t"<< "Dog Count"<<"\t"<< "Path"<< std::endl;
    
    std::cout << "Done." << std::endl;


    std::cout << "Processing..."<<std::endl;
    {   
        fm_index_t::size_type total,count;
        breadth_cur queue(human_idx.cursor());
        //fast forward to desired min depth;
        if(args.min_depth>0)
            queue.fast_forward(args.min_depth);

        std::chrono::time_point<std::chrono::system_clock> past = std::chrono::system_clock::now();
        #pragma omp parallel
        while( (args.max_depth==0 || queue.depth()<args.max_depth) && !queue.empty()){
        cursor_t cur;
        #pragma omp critical
        {
            total++;
            count++;
           cur=queue.next();
        }

           // print_cur(cur,texts);
           auto label = cur.path_label(texts);
        #ifndef NDEBUG
        seqan3::debug_stream << \
        "queue depth:" << queue.depth() << "\t\t label:"<< queue.front().path_label(texts) << "\n" <<
        "depth:" << cur.query_length()/*queue.depth()*/ <<"\t\t label:" << label << "\n" << std::endl;
        #endif
        // #pragma omp single
        
        {
            std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
            auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now-past);
            // seqan3::debug_stream<<milliseconds.count();
            if(milliseconds.count()>30*1000){

                past = std::chrono::system_clock::now();
                auto givemetime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                seqan3::debug_stream << ctime(&givemetime);
                seqan3::debug_stream << "total: "<< total <<" count: "<<count<< " label: \n" << label;
                seqan3::debug_stream <<std::endl;
                count=0;
            }
            
        }
           {
            if(!chimp_idx.cursor().extend_right(label) && !wolf_idx.cursor().extend_right(label)){
                auto dog_cur = dog_idx.cursor();
                if(dog_cur.extend_right(label)){
                    auto dog_count{dog_cur.count()};
                    #pragma omp critical
                    deos << cur.count()<< "\t" << dog_count <<"\t" << label<< std::endl;   
                }    
            }
            
           }
           
        }
    }        
    std::cout << "Done." << std::endl;
    std::cout << "Closing stream...";
    os.close();
    std::cout << "Done." << std::endl;

    std::cout << "All Complete. GZ GZ" << std::endl;
}