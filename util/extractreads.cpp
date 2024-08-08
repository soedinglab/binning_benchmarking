#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <dirent.h>

using IDsList = std::unordered_set<std::string>;
using ContigBinPairs = std::unordered_map<std::string, std::string>;
using ReadBinPairs = std::unordered_map<std::string, std::string>;

// get contig ids from a bin file
IDsList extractContigIDsFromFasta(const std::string& binfasta) {
    std::unordered_set<std::string> contig_ids;
    std::ifstream file(binfasta);
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] == '>') {
            // remove > character
            std::string contig_id = line.substr(1);
            contig_ids.insert(contig_id);
        }
    }
    return contig_ids;
}


// Function to list all files with a specific extension in a directory
std::vector<std::string> getfastalist(const std::string& dirPath, const std::string& extension) {
    std::vector<std::string> files;
    DIR *dir = opendir(dirPath.c_str());
    struct dirent *entry;

    if (dir == nullptr) {
        std::cerr << "Could not open directory: " << dirPath << "\n";
        return files;
    }

    while ((entry = readdir(dir)) != nullptr) {
        std::string filename = entry->d_name;
        if (filename.length() > extension.length() && 
            filename.substr(filename.length() - extension.length()) == extension) {
            files.push_back(dirPath + "/" + filename);
        }
    }

    closedir(dir);
    return files;
}

// collect contig ids and file name from all bin fasta files
ContigBinPairs processfastafiles(std::string inputdir, std::string& format){
    ContigBinPairs contigbinpairs;
    std::vector<std::string> fastafiles = getfastalist(inputdir, format);
    for (const auto& filepath : fastafiles) {
        IDsList contig_ids = extractContigIDsFromFasta(filepath);
        // Replace .fasta with .fastq
        std::string fastqfilename = filepath.substr(0, filepath.find_last_of('.')) + ".fastq";

        for (const auto& contig_id : contig_ids) {
            contigbinpairs[contig_id] = fastqfilename;
            // std::cout << contig_id << " " << fastqfilename << '\n';
        }
    }

    return contigbinpairs;

}
// Function to process map file that contains list of reads mapped to contigs
ReadBinPairs processMapFile(ContigBinPairs& contigbinpairs, const std::string& mapfile) {
    std::ifstream file(mapfile);
    std::string line;
    ReadBinPairs readbinpairs;
    while (std::getline(file, line)) {    
        size_t pos = line.find(' ');
        if (pos != std::string::npos) {
            std::string read_id = line.substr(0, pos);
            std::string contig_id = line.substr(pos + 1);
            pos = contig_id.find(' ');
            contig_id = contig_id.substr(0, pos);
            auto it = contigbinpairs.find(contig_id);
            if (it != contigbinpairs.end()) {
                readbinpairs[read_id] = it->second;
            }
        }
    }
    return readbinpairs;
}

// write fastq file only for mapped reads for each bin
void extractreads(const std::string& infastq, ReadBinPairs& readbinpairs) {
    std::ifstream infile(infastq);
    std::unordered_map<std::string, std::ofstream> outfiles;
    if (!infile.is_open()) {
        std::cerr << "Error opening infile\n";
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] != '@') continue;
        size_t endpos = line.find('/');
        // if condition single end read (eg., @S20R8323008)
        // else condition paired end reads (eg., @S20R8323008/1)
        std::string read_id = (endpos == std::string::npos) ? line.substr(1) : line.substr(1, endpos - 1);

        auto it = readbinpairs.find(read_id);
        if (it != readbinpairs.end()) {
            const std::string& fastqfilename = it->second;
            if (outfiles.find(fastqfilename) == outfiles.end()) {
                outfiles[fastqfilename].open(fastqfilename, std::ios::app);
                if (!outfiles[fastqfilename].is_open()) {
                    std::cerr << "Error opening outfile: " << fastqfilename << "\n";
                    continue;
                }
            }
            std::ofstream& outfile = outfiles[fastqfilename];
            outfile << line << '\n';  // Write the read ID line
            for (int i = 0; i < 3; ++i) {
                std::getline(infile, line);
                outfile << line << '\n';
            }
        } else {
            // Skip the next 3 lines (sequence, +, quality score), if read not matched
            for (int i = 0; i < 3; ++i) {
                std::getline(infile, line);
            }
        }
    }

    infile.close();
    for (auto& pair : outfiles) {
        pair.second.close();
    }
}

int main(int argc, char *argv[]) {

    if (argc == 1 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        std::cerr << "extractreads inputdir mapfile inputreads.fastq -f (optional, default=fasta)" << "\n";
        return 1;
    }
    if (argc > 6) {
        std::cerr << "you gave more input arguments!\n";
        return 1;
    }
    
    std::string format = "fasta";
    for (int i = 3; i < argc; ++i) {
        if (strcmp(argv[i], "-f") == 0) {
            if (i+1 < argc) {
                format = argv[i+1];
                ++i;
            }
        }
    }

    // fasta file for bins with full absolute path required
    std::string bindir = argv[1];
    // read and contig pair from all mapping files (using aligner2counts script)
    std::string mapfile = argv[2];
    // input fastq file to extract read (pooled reads, recommended)
    std::string infastq = argv[3];
    
    // Load contig IDs from all bin files
    ContigBinPairs contigbinpairs = processfastafiles(bindir, format);
    // Load read IDs from all each count file
    ReadBinPairs readbinpairs = processMapFile(contigbinpairs, mapfile);


    // Extract reads by IDs
    extractreads(infastq, readbinpairs);

    return 0;
}

