#include <iostream>
#include <fstream>
#include <string>
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
void split_fastq_by_sampleids(const std::string& fastq_file, const std::unordered_set<std::string>& ids, const std::string& outdir) {
    std::ifstream infile(fastq_file);
    std::string line;
    std::unordered_map<std::string, std::ofstream> outfiles; // Map to store file streams by ID

    while (std::getline(infile, line)) {

        for (const auto& id: ids) {
            if (line.find(id) != std::string::npos && line[0] == '@') {
                // Open the output file for this ID if not already open
                if (outfiles.find(id) == outfiles.end()) {
                    outfiles[id].open(outdir + "/" + id + ".fastq", std::ios_base::app); // Open in append mode
                }

                // Write the header and the next 3 lines (sequence, plus sign, and quality score)
                outfiles[id] << line << '\n';
                for (int i = 0; i < 3; ++i) {
                    std::getline(infile, line);
                    outfiles[id] << line << '\n';
                }
                // move to next line
                break;
            }
        }
    }

    // Close all opened output files
    for (auto& pair : outfiles) {
        pair.second.close();
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <id_file> <fastq_file> <outdir>" << std::endl;
        return 1;
    }

    std::string id_file = argv[1];
    std::string fastq_file = argv[2];
    std::string outdir = argv[3];

    create_directory(outdir);

    // Load the IDs
    std::unordered_set<std::string> ids = load_ids(id_file);

    // Split the FASTQ file by IDs
    split_fastq_by_sampleids(fastq_file, ids, outdir);

    std::cout << "FASTQ splitting completed." << std::endl;

    return 0;
}