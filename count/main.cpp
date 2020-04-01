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
        std::filesystem::path head_index_path{};
        std::vector<std::filesystem::path> file_paths{};
        std::vector<std::filesystem::path> index_paths{};


        //search configuration
        size_type min_depth{0};
        size_type max_depth{0};

        //output configuration
        std::filesystem::path outTSV{"count.tsv"};

    };
    void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
    {
        parser.info.author = "Yorai Levi";
        parser.info.short_description = "Count the number of occurences of all subsrtings";
        parser.info.version = VERSION;

        parser.add_option(args.file_paths, 'f', "file", "original fasta input files that were indexed(head)");
        parser.add_option(args.head_index_path, 'x', "head-index", "index of fasta input"); 

        parser.add_option(args.index_paths, 'i', "index", "index of fasta input to query");

        parser.add_option(args.outTSV, '\0', "out-file-name", "name of file to output, as for this version please include .tsv as a file extension");
        

        parser.add_option(args.max_depth, '\0', "max", "max depth to search");
        parser.add_option(args.min_depth, '\0', "min", "min depth to search");
    }
    int main(int argc, char **argv)
    {
        cmd_arguments args = parse_args<cmd_arguments>(PROGRAM_NAME,argc, argv);
        
        std::cout << "Loading..."<<std::endl;
        //loading head
        fm_index_t head_idx;
        std::vector<sequence_type> texts{};
        {
            std::vector<id_type> ids{};
            load_fm(head_idx,args.head_index_path);
            load_records(args.file_paths,texts,ids,false);
        }
        //loading indexes
        std::vector<fm_index_t> v_idx(args.index_paths.size());
        {
            std::vector<fm_index_t>::size_type i=0;
            for(auto &idx : v_idx)
            {
                load_fm(idx,args.index_paths[i++]);
            }
            //for(std::vector<fm_index_t>::size_type i=0;i<args.index_paths.size();i++)
            // {
            //     fm_index_t idx;
            //     load_fm(idx,args.index_paths[i]);
            //     v_idx.push_back(idx);
            // }
        }

        std::cout << "There are: "<< v_idx.size()<<" indexes"<<std::endl;
        std::cout << "Done." << std::endl;

        std::cout << "Prepearing stream for writing...";
        
        std::ofstream os;
        seqan3::debug_stream_type<> deos{os};
        os.open(args.outTSV);
        for(auto &id : args.index_paths)
        {
            deos << id << "\t";
        }
        deos<<"Path"<< std::endl;
        
        std::cout << "Done." << std::endl;


        std::cout << "Processing...";
        {   
            breadth_cur queue(head_idx.cursor());
            if(args.min_depth>0)
                queue.fast_forward(args.min_depth);

            fm_index_t::size_type total,count;
            // breadth_cur queue(human_idx.cursor());
            //fast forward to desired min depth;

            //scheduled printing:
            std::chrono::time_point<std::chrono::system_clock> past = std::chrono::system_clock::now();
            
            //main loop
            #pragma omp parallel
            while( (args.max_depth==0 || queue.depth()<args.max_depth) && !queue.empty())
            {
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
                    if(milliseconds.count()>30*1000)
                    {
                        past = std::chrono::system_clock::now();
                        auto givemetime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                        seqan3::debug_stream << ctime(&givemetime);
                        seqan3::debug_stream << "total: "<< total <<" count: "<<count<< " label: \n" << label;
                        seqan3::debug_stream <<std::endl;
                        count=0;
                    }
                    
                }
                    //main program
                {        
                    std::vector<size_type> counts;
                    for(auto &idx : v_idx)
                    {
                        counts.push_back(has_node(idx,label));
                    }
                    #pragma omp critical
                    {
                        for(auto &count : counts)
                        {
                        deos << count << "\t";        
                        }
                        deos << label<< std::endl;
                    }
                }
                   
            } //main loop end
        }
        std::cout << "Done." << std::endl;
        std::cout << "Closing stream...";
        os.close();
        std::cout << "Done." << std::endl;

        std::cout << "All Complete. GZ GZ" << std::endl;
    }