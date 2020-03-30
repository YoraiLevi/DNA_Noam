#include <string>
#include <fstream>
#include <iostream>
// #include <iostream>
#include <vector>
#include <algorithm>

#include <seqan3/argument_parser/all.hpp> 
#include <seqan3/search/fm_index/fm_index.hpp>
#include <common.cpp>
#define PROGRAM_NAME "Count"
#define VERSION "0.0.1"
#include <seqan3/core/detail/all.hpp>
struct cmd_arguments
{
    // std::vector<std::filesystem::path> file_paths{};
    // std::filesystem::path all{INDEX_FILE_NAME};
    std::vector<std::filesystem::path> file_paths{};


    std::filesystem::path outTSV{"count.tsv"};

};
void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
{
    parser.info.author = "Yorai Levi";
    parser.info.short_description = "Count the number of occurences of all subsrtings";
    parser.info.version = VERSION;
    // parser.add_option(args.file_paths, 'f', "file", "fasta input file to index");
    parser.add_option(args.file_paths, 'f', "file", "original fasta input files that were indexed");

    parser.add_option(args.outTSV, '\0', "out-file-name", "name of file to output, as for this version please include .tsv as a file extension");
    // parser.add_option(args.include, 'i', "include", "index meta to include");

    // parser.add_flag(args.verify, 'v', "verify", "verifies that saved index is accurate");
    // parser.add_flag(args.meta, 'm', "meta-only", "save meta data only");
    // parser.add_option(args.index_dir, 'd', "dir", "The output directory for storing the files.",
                    //   seqan3::option_spec::REQUIRED, seqan3::output_directory_validator{});
    //.fa, .fasta, .fna, .ffn, .ffa, .frn.fa, .fasta, .fna, .ffn, .ffa, .frn
    //parser.add_positional_option(args.file_paths,"fasta input file to index");//,seqan3::option_spec::REQUIRED,seqan3::input_file_validator{{"fa","fasta,fna"}});
}
int main(int argc, char **argv)
{
    cmd_arguments args = parse_args<cmd_arguments>(PROGRAM_NAME,argc, argv);
    // fm_index_t idx;
    // fm_index_t dog_idx,wolf_idx,chimp_idx,human_idx;
    // load_fm(idx);
    std::cout << "Loading...";
    std::vector<sequence_type> texts{};
    std::vector<id_type> ids{};
    load_records(args.file_paths,texts,ids,false);

    std::cout << "There are: "<< texts.size()<<"texts";
    std::cout << " Done." << std::endl;
    std::cout << "Creating Indexes...";  
    fm_index_t index(texts);
    std::vector<fm_index_t> v_idx;
    for(auto &text : texts){
        std::vector<sequence_type> t;
        t.push_back(text);
        v_idx.emplace_back(t);
    }
    std::cout << "There are: "<< v_idx.size()<<"indexes";
    std::cout << " Done." << std::endl;
    std::cout << "Prepearing stream for writing...";
    
    std::ofstream os;
    seqan3::debug_stream_type<> deos{os};
    os.open(args.outTSV);
    for(auto &id : ids){
        deos << id << "\t";
    }
    deos<<"Path"<< std::endl;
    
    std::cout << "Done." << std::endl;


    std::cout << "Processing...";
    {   
        breadth_cur queue(index.cursor());
        while(!queue.empty()){
           cursor_t cur(queue.next());
           // print_cur(cur,texts);
           auto label = cur.path_label(texts);
           {
            for(auto &idx : v_idx){
                deos << has_node(idx,label) << "\t";
            }
            deos <<label <<std::endl;
            // if(!has_node(chimp_idx,label) && !has_node(wolf_idx,label)){
            //     auto dog_count{has_node(dog_idx,label)};
            //     if(dog_count){
            //         deos << cur.count()<< "\t" << dog_count <<"\t" << label<< std::endl;   
            //     }    
            // }
            
           }
           
        }
    }        
    std::cout << "Done." << std::endl;
    std::cout << "Closing stream...";
    os.close();
    std::cout << "Done." << std::endl;

    std::cout << "All Complete. GZ GZ" << std::endl;
}