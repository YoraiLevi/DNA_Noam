#include <string>
#include <fstream>
#include <iostream>

#include <seqan3/argument_parser/all.hpp>    // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>      // for debug_stream
#include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
#include <seqan3/search/fm_index/fm_index.hpp>

#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
// #include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input

// #include "utils/records.h"

#define PROGRAM_NAME "FM-Indexer"
struct cmd_arguments
{
    std::vector<std::filesystem::path> file_paths{};
    std::filesystem::path index_dir{};
    bool verify{};
    bool meta{};
};
void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
{
    parser.info.author = "Yorai Levi";
    parser.info.short_description = "FM-INDEXES fasta formatted files";
    parser.info.version = "0.0.1";
    parser.add_option(args.file_paths, 'f', "file", "fasta input file to index");
    parser.add_flag(args.verify, 'v', "verify", "verifies that saved index is accurate");
    parser.add_flag(args.meta, 'm', "meta-only", "save meta data only");
    // parser.add_option(args.index_dir, 'd', "dir", "The output directory for storing the files.",
                    //   seqan3::option_spec::REQUIRED, seqan3::output_directory_validator{});
    //.fa, .fasta, .fna, .ffn, .ffa, .frn.fa, .fasta, .fna, .ffn, .ffa, .frn
    //parser.add_positional_option(args.file_paths,"fasta input file to index");//,seqan3::option_spec::REQUIRED,seqan3::input_file_validator{{"fa","fasta,fna"}});
}
cmd_arguments parse_args(int argc, char **argv)
{
    //throws -1 if failed
    seqan3::argument_parser myparser{PROGRAM_NAME, argc, argv}; // initialise myparser
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
auto load_records(const std::filesystem::path &file_path)
{
    throw "not implemented";
}
auto store_meta(const auto &ids)
{
    std::filesystem::path path = "index.meta.json";
    std::ofstream os{path};
    cereal::JSONOutputArchive archive{os};
    archive(CEREAL_NVP(ids));
    // archive(CEREAL_NVP(alphabet));
}
auto store_fm(const auto &index,const bool is_verify)
{
    
    std::filesystem::path path = "index.bin";
    {
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
    }
    if(is_verify)
    {
        seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> index2;
    {
        std::ifstream is{path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index2);
    }
 
    if (index == index2)
        std::cout << "The indices are identical!\n";
    else
        std::cout << "The indices differ!\n";
    }
    
}
auto load_fm()
{
    throw "not implemented";
}

int main(int argc, char **argv)
{
    cmd_arguments args = parse_args(argc, argv);
#ifndef NDEBUG
    seqan3::debug_stream << "Printing args" << std::endl;
    seqan3::debug_stream << args.file_paths << std::endl;
    seqan3::debug_stream << args.index_dir << std::endl;
#endif
    // using record_type = seqan3::sequence_file_input<>::record_type;
    using id_type = seqan3::sequence_file_input<>::id_type;
    using sequence_type = seqan3::sequence_file_input<>::sequence_type;
    // std::vector<record_type> records();
    std::vector<sequence_type> sequences{};
    std::vector<id_type> ids{};

    for (auto &file_path : args.file_paths)
    {
        seqan3::sequence_file_input fin{file_path};
        for (auto &[seq, id, qual] : fin)
        {
            #ifndef NDEBUG
            seqan3::debug_stream <<"seq found:"<< id << std::endl;
            #endif
            if(!args.meta){
            sequences.push_back(std::move(seq));    
            }
            ids.push_back(std::move(id));
        }
    }
#ifndef NDEBUG
    seqan3::debug_stream << ids << std::endl;
#endif
    store_meta(ids);
    if(!args.meta){
    seqan3::fm_index index(sequences);
    store_fm(index,args.verify);    
    }
}
