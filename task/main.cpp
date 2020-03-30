#include <string>
#include <vector>
#include <algorithm>

#include <seqan3/core/debug_stream.hpp>      // for debug_stream
#include <seqan3/argument_parser/all.hpp>    // includes all necessary headers
#include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>

//exclude //include //iterator //ostream save to file
void print_cur(const seqan3::fm_index_cursor<auto> &cur, const std::vector<auto> &genomes)
{
    seqan3::debug_stream
        /* << " last rank: " << cur.last_rank()
        << " count: " << cur.count()
        << " len: " << cur.query_length()
        << " locate: " << cur.locate()
        << " str: "*/
        << cur.path_label(genomes)
        << '\n';
}
void deep_cur(const seqan3::fm_index_cursor<auto> &cur, const std::vector<auto> &genomes)
{
    //the correct thing to do is probably writing an iterator and ++ing the iterator
    print_cur(cur, genomes);
    auto temp_cur(cur); //hopefully copy constructor
    if (!temp_cur.extend_right())
    { //if deepest then done.
        return;
    }
    do
    {
        deep_cur(temp_cur, genomes);
    } while (temp_cur.cycle_back());
}

// template <typename CursorType1, typename CursorType2, typename T1, typename T2>
// void double_deep_cur(const CursorType1 &cur_desired, const CursorType2 &cur_disliked, const T1 &genomes_desired, const T2 &genomes_disliked, std::vector<CursorType1> &out)
template <typename T1, typename T2,typename T3, typename T4>
void double_deep_cur(const seqan3::fm_index_cursor<T1> &cur_desired,const seqan3::fm_index_cursor<T2> &cur_disliked, const std::vector<T3> &genomes_desired,const std::vector<T4> &genomes_disliked, std::vector<seqan3::fm_index_cursor<T1>> &out)
{
    //hopefully this doesn't stackoverflow.
    //the correct thing to do is probably writing an iterator and ++ing the iterator
    // print_cur(cur, genomes);
    auto temp_cur_desired(cur_desired); //hopefully copy constructor
    auto temp_cur_disliked(cur_disliked);

    auto is_desired_extended = temp_cur_desired.extend_right();
    auto is_disliked_extended = temp_cur_disliked.extend_right();
    //desired 0 disliked 0
    //this means both branches died, who cares. stop
    //desired 0 disliked 1
    //this means desired branch died, who cares. stop.
    //=> desired is false, stop.
    // seqan3::debug_stream << "visiting:"<< is_desired_extended<< " "; print_cur(temp_cur_desired,genomes_desired);
    // seqan3::debug_stream << "visiting:"<< is_desired_extended<< " "; print_cur(temp_cur_disliked,genomes_disliked);
    if (!is_desired_extended)
    {
        seqan3::debug_stream << "death to nodes";
        return;
    }

    //desired 1 disliked 0
    //this means disliked branch died, and desired continues! good!!
    //desired 1 disliked 1
    //this means desired branch continues and disliked branch continues,
    //need to check if paths are the same.
    //different paths? good!!
    bool is_desired_cycle_back, is_disliked_cycle_back;
    do
    {
        if (!is_disliked_extended)
        {
            // seqan3::debug_stream << "accepted node1";
            print_cur(temp_cur_desired,genomes_desired);
            //out.push_back(temp_cur_desired);
            return;
        }
        else
        {
            auto desired_path = temp_cur_desired.path_label(genomes_desired);
            auto disliked_path = temp_cur_disliked.path_label(genomes_disliked);
            auto booooo = std::equal(desired_path.begin(),desired_path.end(),disliked_path.begin());
            if (!booooo)
            {
                // seqan3::debug_stream << "accepted node2";
                print_cur(temp_cur_desired,genomes_desired);
                //out.push_back(temp_cur_desired);
                return;
            }
        }
        //continue to add letter
        double_deep_cur(temp_cur_desired, temp_cur_disliked, genomes_desired, genomes_disliked,out);
        //We're changing letter!
        temp_cur_disliked.cycle_back();
        //desired 0 disliked 0
        //desired 0 disliked 1
        //stop.
        //desired 1 disliked 0
        //desired 1 disliked 1
        //continue to check nodes:
    } while (temp_cur_desired.cycle_back());
}
// template <typename CursorType1, typename CursorType2, typename T1, typename T2>
// std::vector<CursorType1> double_deep_cur(const CursorType1 &cur_desired, const CursorType2 &cur_disliked, const T1 &genomes_desired, const T2 &genomes_disliked)
template <typename T1, typename T2,typename T3, typename T4>
std::vector<seqan3::fm_index_cursor<T1>> double_deep_cur(const seqan3::fm_index_cursor<T1> &cur_desired,const seqan3::fm_index_cursor<T2> &cur_disliked, const std::vector<T3> &genomes_desired,const std::vector<T4> &genomes_disliked)
{
    std::vector<seqan3::fm_index_cursor<T1>> out(10);
    // out.reserve(10 ^ 6);
    // auto out = null;
    double_deep_cur(cur_desired, cur_disliked, genomes_desired, genomes_disliked,out);
    return out;
}
void print_suffixes(const auto &index, const std::vector<auto> &genomes)
{
    auto cur = index.cursor();
    deep_cur(cur, genomes);
}

struct cmd_arguments
{
    std::filesystem::path chimp_seq{};
    std::filesystem::path human_seq{};
    // std::filesystem::path dog_seq{};
    // std::filesystem::path wolf_seq{};
};
void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
{
    parser.info.author = "Yorai Levi";
    parser.info.short_description = "No fucking clue";
    parser.info.version = "0.0.1";

    parser.add_option(args.chimp_seq, '\0', "chimp", "dna sequence of chimp", seqan3::option_spec::REQUIRED, seqan3::input_file_validator{{"fa", "fna", "fasta"}});
    parser.add_option(args.human_seq, '\0', "human", "dna sequence of human", seqan3::option_spec::REQUIRED, seqan3::input_file_validator{{"fa", "fna", "fasta"}});
    // parser.add_option(args.dog_seq,'\0',"dog","dna sequence of household dog",seqan3::option_spec::REQUIRED,seqan3::input_file_validator{{"fa","fasta"}});
    // parser.add_option(args.wolf_seq,'\0',"wolf","dna sequence of wolf",seqan3::option_spec::REQUIRED,seqan3::input_file_validator{{"fa","fasta"}});
}
auto load_records(const std::filesystem::path &file_path)
{
    seqan3::sequence_file_input fin{file_path};
    using record_type = decltype(fin)::record_type;
    std::vector<record_type> records{};
    std::ranges::copy(fin, std::ranges::back_inserter(records));
    return records;
}
void program(const cmd_arguments &args)
{
    auto chimp_records = load_records(args.chimp_seq);
    auto human_records = load_records(args.human_seq);
    // seqan3::debug_stream << seqan3::get<seqan3::field::id>(chimp_records[0]) << '\n';
    auto chimp_seq = seqan3::get<seqan3::field::seq>(chimp_records[0]);
    auto human_seq = seqan3::get<seqan3::field::seq>(human_records[0]);
    seqan3::fm_index chimp_index{chimp_seq};
    seqan3::fm_index human_index{human_seq};

    auto chimp_cur = chimp_index.cursor();
    auto human_cur = human_index.cursor();
    // deep_cur(chimp_cur,chimp_seq);
    // deep_cur(human_cur,human_seq);
    
    auto out = double_deep_cur(human_cur, chimp_cur, human_seq, chimp_seq);
    // std::cout << out.size()<<std::endl;
    // for(auto cur : out){
    //     print_cur(cur,human_seq);
    // }
    // print_suffixes(human_cur,seqan3::get<seqan3::field::seq>(human_records[0]));
    
}
int main(int argc, char **argv)
{
    seqan3::argument_parser myparser{"ProgramName", argc, argv}; // initialise myparser
    cmd_arguments args{};

    initialise_argument_parser(myparser, args);

    try
    {
        myparser.parse(); // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const &ext) // catch user errors
    {
        seqan3::debug_stream << "[Error Message:] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
    program(args);
}

// int main()
// {
//     using seqan3::operator""_dna4;

//     std::vector<seqan3::dna4> genomes{"AATGAATGAACAGATAATAATAACAATGAAC"_dna4};
//     seqan3::fm_index index{genomes}; // build the index
//     print_suffixes(index, genomes);
// }
