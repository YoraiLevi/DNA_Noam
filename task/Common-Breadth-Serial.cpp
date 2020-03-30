#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <seqan3/argument_parser/all.hpp> 
#include <seqan3/search/fm_index/fm_index.hpp>
#include <common.cpp>
#define PROGRAM_NAME "Common-Breath-Serial"
#define VERSION "0.0.1"
#include <seqan3/core/detail/all.hpp>
struct cmd_arguments
{
    std::filesystem::path dog{};
    std::filesystem::path wolf{};
    std::filesystem::path chimp{};
    std::filesystem::path human{};
    std::vector<std::filesystem::path> file_paths{};

    std::vector<std::filesystem::path> tile_paths{};

    std::filesystem::path outTSV{"out.tsv"};

};
void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
{
    parser.info.author = "Yorai Levi";
    parser.info.short_description = "FM-INDEXES fasta formatted files";
    parser.info.version = VERSION;
    parser.add_option(args.dog, 'd', "dog", "dog index");
    parser.add_option(args.wolf, 'w', "wolf", "wolf index");
    parser.add_option(args.chimp, 'c', "chimp", "chimp index");
    parser.add_option(args.human, 'u', "human", "human index");

    parser.add_option(args.file_paths, 'f', "file", "original fasta input files that were indexed");
    parser.add_option(args.tile_paths, 't', "tile", "original fasta input files that were indexed");
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

    std::cout << " Done." << std::endl;

    std::cout << "Prepearing stream for writing...";
    
    std::ofstream os;
    seqan3::debug_stream_type<> deos{os};
    os.open(args.outTSV);
    deos << "Human Count" <<"\t"<< "Dog Count"<<"\t"<< "Path"<< std::endl;
    
    std::cout << " Done." << std::endl;


    std::cout << "Processing...";
    {   
        breadth_cur queue(human_idx.cursor());
        while(!queue.empty()){
           cursor_t cur(queue.next());
           // print_cur(cur,texts);
           auto label = cur.path_label(texts);
           {
            if(!has_node(chimp_idx,label) && !has_node(wolf_idx,label)){
                auto dog_count{has_node(dog_idx,label)};
                if(dog_count){
                    deos << cur.count()<< "\t" << dog_count <<"\t" << label<< std::endl;   
                }    
            }
            
           }
           
        }
    }        
    std::cout << " Done." << std::endl;
    std::cout << "Closing stream...";
    os.close();
    std::cout << " Done." << std::endl;

    std::cout << " All Complete. GZ GZ" << std::endl;
}