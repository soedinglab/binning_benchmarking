#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <dirent.h>
#include <unordered_set>
#include <unordered_map>


// Function to create a directory if it doesn't exist
void create_directory(const std::string& directory) {
    DIR* dir = opendir(directory.c_str());
    if (dir) {
        // Directory exists
        closedir(dir);
    } else {
        // Directory does not exist, create it
        if (mkdir(directory.c_str(), 0777) == -1) {
            std::cerr << "Error: Unable to create directory " << directory << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

// Function to load IDs from input file into a set
std::unordered_set<std::string> load_ids(const std::string& id_file) {
    std::unordered_set<std::string> ids;
    std::ifstream infile(id_file);
    std::string line;
    while (std::getline(infile, line)) {
        ids.insert(line);
    }
    return ids;
}

// Function to process the FASTQ file and split reads by IDs
void split_fastq_by_sampleids(
    const std::string& fastq_file,
    const std::unordered_set<std::string>& ids,
    const std::string& outdir, bool paired) {
    std::ifstream infile(fastq_file);
    std::string line;
    std::unordered_map<std::string, std::ofstream> outfiles_1;
    std::unordered_map<std::string, std::ofstream> outfiles_2;

    while (std::getline(infile, line)) {

        for (const auto& id: ids) {
            if (line.find(id) != std::string::npos && line[0] == '@') {
                bool is_read1 = line.find("/1") != std::string::npos;
                bool is_read2 = line.find("/2") != std::string::npos;
                if (paired && !(is_read1 || is_read2)) {
                    std::cout << "skipping unpaired reads\n";
                    break;
                }

                if (paired) {
                    if (is_read1) {
                        if (outfiles_1.find(id) == outfiles_1.end()) {
                            outfiles_1[id].open(outdir + "/" + id + "_1.fastq", std::ios_base::app);
                        }
                        outfiles_1[id] << line << '\n';
                    } else if (is_read2) {
                        if (outfiles_2.find(id) == outfiles_2.end()) {
                            outfiles_2[id].open(outdir + "/" + id + "_2.fastq", std::ios_base::app);
                        }
                        outfiles_2[id] << line << '\n';
                    }
                } else {
                    if (outfiles_1.find(id) == outfiles_1.end()) {
                        outfiles_1[id].open(outdir + "/" + id + ".fastq", std::ios_base::app);
                    }
                    outfiles_1[id] << line << '\n';
                }
                
                // Write the next 3 lines (sequence, plus sign, and quality score)
                for (int i = 0; i < 3; ++i) {
                    std::getline(infile, line);
                    if (paired) {
                        if (is_read1) {
                            outfiles_1[id] << line << '\n';
                        } else if (is_read2) {
                            outfiles_2[id] << line << '\n';
                        }
                    } else {
                        outfiles_1[id] << line << '\n';
                    }
                }
                break;
            }
        }
    }

    // Close all opened output files
    for (auto& pair : outfiles_1) {
        pair.second.close();
    }
    for (auto& pair : outfiles_2) {
        pair.second.close();
    }
}

int main(int argc, char* argv[]) {
    
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <id_file> <fastq_file> <outdir> [--paired]" << std::endl;
        return 1;
    }
    std::string id_file = argv[1];
    std::string fastq_file = argv[2];
    std::string outdir = argv[3];
    bool paired = false;
    if (argc > 4) {
        for (int i = 1; i < argc; ++i) {
            std::string flag = argv[i];
            if (flag == "--paired") {
                paired = true;
                ++i;
            }
        }
    }
    create_directory(outdir);

    // Load the IDs
    std::unordered_set<std::string> ids = load_ids(id_file);

    // Split the FASTQ file by IDs
    split_fastq_by_sampleids(fastq_file, ids, outdir, paired);

    std::cout << "FASTQ splitting completed." << std::endl;

    return 0;
}