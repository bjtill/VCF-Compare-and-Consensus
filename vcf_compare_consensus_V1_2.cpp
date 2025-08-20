//BT June 13, 2025
//Create a consensus VCF from 2 or 3 input VCFs compressed or otherwise, and then creat a summar table that can be used to make a venn or other plot. 

//g++ -o intersect2 vcf_compare_consensus_V1_2.cpp -std=c++11
//# Basic usage (uncompressed)
//./vcf_compare file1.vcf file2.vcf

//# With compression
//./vcf_compare --compress file1.vcf.gz file2.vcf.gz file3.vcf.gz

//# Mixed input formats
//./vcf_compare --compress file1.vcf file2.vcf.gz file3.vcf

//Output Files:

//consensus_table.txt - Tab-delimited table perfect for downstream analysis, Venn diagrams, etc.
//consensus.vcf or consensus.vcf.gz - Consensus variants in VCF format

//Requirements:

//For compressed input: zcat (usually available with gzip)
//For compressed output: bgzip (part of samtools/htslib)


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

struct Variant {
    std::string chrom;
    int pos;
    std::string ref;
    std::string alt;
    std::string id;
    std::string qual;
    std::string filter;
    std::string info;
    std::string format;
    std::vector<std::string> samples;
    
    // Create a unique key for comparison
    std::string getKey() const {
        return chrom + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
    }
    
    // Equality operator for comparison
    bool operator==(const Variant& other) const {
        return getKey() == other.getKey();
    }
};

class VCFReader {
private:
    std::vector<std::string> header_lines;
    std::vector<std::string> sample_names;
    
    bool isGzipped(const std::string& filename) {
        return filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz";
    }
    
public:
    std::vector<Variant> readVCF(const std::string& filename) {
        std::vector<Variant> variants;
        std::string line;
        FILE* file_ptr = nullptr;
        std::ifstream regular_file;
        bool is_compressed = isGzipped(filename);
        
        if (is_compressed) {
            // Use popen to read compressed file
            std::string command = "zcat \"" + filename + "\"";
            file_ptr = popen(command.c_str(), "r");
            if (!file_ptr) {
                std::cerr << "Error: Cannot open compressed file " << filename << std::endl;
                std::cerr << "Make sure zcat is available in your PATH" << std::endl;
                return variants;
            }
        } else {
            // Regular file
            regular_file.open(filename);
            if (!regular_file.is_open()) {
                std::cerr << "Error: Cannot open file " << filename << std::endl;
                return variants;
            }
        }
        
        // Read lines
        char buffer[65536]; // Large buffer for VCF lines
        while (true) {
            bool line_read = false;
            if (is_compressed) {
                if (fgets(buffer, sizeof(buffer), file_ptr)) {
                    line = std::string(buffer);
                    // Remove trailing newline
                    if (!line.empty() && line.back() == '\n') {
                        line.pop_back();
                    }
                    line_read = true;
                }
            } else {
                if (std::getline(regular_file, line)) {
                    line_read = true;
                }
            }
            
            if (!line_read) break;
            if (line.empty()) continue;
            
            if (line[0] == '#') {
                header_lines.push_back(line);
                if (line.substr(0, 6) == "#CHROM") {
                    // Parse sample names from header
                    std::istringstream iss(line);
                    std::string token;
                    int col = 0;
                    while (iss >> token) {
                        if (col >= 9) { // Sample names start from column 10 (0-indexed 9)
                            sample_names.push_back(token);
                        }
                        col++;
                    }
                }
            } else {
                // Parse variant line
                Variant var = parseVariantLine(line);
                if (!var.chrom.empty()) {
                    variants.push_back(var);
                }
            }
        }
        
        // Close file
        if (is_compressed) {
            pclose(file_ptr);
        } else {
            regular_file.close();
        }
        
        std::cout << "Read " << variants.size() << " variants from " << filename << std::endl;
        return variants;
    }
    
    const std::vector<std::string>& getHeaderLines() const { return header_lines; }
    const std::vector<std::string>& getSampleNames() const { return sample_names; }
    
private:
    Variant parseVariantLine(const std::string& line) {
        Variant var;
        std::istringstream iss(line);
        std::string token;
        int col = 0;
        
        while (std::getline(iss, token, '\t')) {
            switch (col) {
                case 0: var.chrom = token; break;
                case 1: 
                    try {
                        var.pos = std::stoi(token);
                    } catch (const std::exception& e) {
                        std::cerr << "Warning: Invalid position " << token << std::endl;
                        var.pos = 0;
                    }
                    break;
                case 2: var.id = token; break;
                case 3: var.ref = token; break;
                case 4: var.alt = token; break;
                case 5: var.qual = token; break;
                case 6: var.filter = token; break;
                case 7: var.info = token; break;
                case 8: var.format = token; break;
                default:
                    if (col >= 9) {
                        var.samples.push_back(token);
                    }
                    break;
            }
            col++;
        }
        
        return var;
    }
};

class VCFComparator {
private:
    std::vector<std::vector<Variant> > vcf_data;
    std::vector<std::string> vcf_names;
    std::vector<std::string> header_lines;
    std::vector<std::string> sample_names;
    bool compress_output;
    
public:
    VCFComparator(bool compress = false) : compress_output(compress) {}
    
    void addVCF(const std::vector<Variant>& variants, const std::string& name,
                const std::vector<std::string>& headers, const std::vector<std::string>& samples) {
        vcf_data.push_back(variants);
        vcf_names.push_back(name);
        if (header_lines.empty()) {
            header_lines = headers;
            sample_names = samples;
        }
    }
    
    void compare() {
        if (vcf_data.size() < 2 || vcf_data.size() > 3) {
            std::cerr << "Error: Need 2 or 3 VCF files for comparison" << std::endl;
            return;
        }
        
        // Create sets for easier comparison
        std::vector<std::set<std::string> > variant_sets;
        std::vector<std::map<std::string, Variant> > variant_maps;
        
        for (size_t i = 0; i < vcf_data.size(); i++) {
            std::set<std::string> var_set;
            std::map<std::string, Variant> var_map;
            
            for (size_t j = 0; j < vcf_data[i].size(); j++) {
                const Variant& variant = vcf_data[i][j];
                std::string key = variant.getKey();
                var_set.insert(key);
                var_map[key] = variant;
            }
            
            variant_sets.push_back(var_set);
            variant_maps.push_back(var_map);
        }
        
        // Print summary table
        printSummaryTable(variant_sets);
        
        // Write summary table to file
        writeSummaryTable(variant_sets);
        
        // Generate consensus VCF
        generateConsensusVCF(variant_sets, variant_maps);
    }
    
private:
    void writeSummaryTable(const std::vector<std::set<std::string> >& variant_sets) {
        std::ofstream table_file("consensus_table.txt");
        if (!table_file.is_open()) {
            std::cerr << "Error: Cannot create consensus_table.txt" << std::endl;
            return;
        }
        
        // Write header
        table_file << "Comparison\tComparison_Files\tCount" << std::endl;
        
        // Individual VCF counts
        for (size_t i = 0; i < variant_sets.size(); i++) {
            table_file << (i + 1) << "\t" << vcf_names[i] << "\t" 
                      << variant_sets[i].size() << std::endl;
        }
        
        if (variant_sets.size() == 2) {
            // Two-way comparison
            std::set<std::string> intersection;
            std::set_intersection(variant_sets[0].begin(), variant_sets[0].end(),
                                variant_sets[1].begin(), variant_sets[1].end(),
                                std::inserter(intersection, intersection.begin()));
            
            table_file << "1,2\t" << vcf_names[0] << "," << vcf_names[1] << "\t"
                      << intersection.size() << std::endl;
                      
        } else if (variant_sets.size() == 3) {
            // Three-way comparison
            std::set<std::string> int_01, int_02, int_12, int_all;
            
            // Pairwise intersections
            std::set_intersection(variant_sets[0].begin(), variant_sets[0].end(),
                                variant_sets[1].begin(), variant_sets[1].end(),
                                std::inserter(int_01, int_01.begin()));
            
            std::set_intersection(variant_sets[0].begin(), variant_sets[0].end(),
                                variant_sets[2].begin(), variant_sets[2].end(),
                                std::inserter(int_02, int_02.begin()));
            
            std::set_intersection(variant_sets[1].begin(), variant_sets[1].end(),
                                variant_sets[2].begin(), variant_sets[2].end(),
                                std::inserter(int_12, int_12.begin()));
            
            // Three-way intersection
            std::set_intersection(int_01.begin(), int_01.end(),
                                variant_sets[2].begin(), variant_sets[2].end(),
                                std::inserter(int_all, int_all.begin()));
            
            table_file << "1,2\t" << vcf_names[0] << "," << vcf_names[1] << "\t"
                      << int_01.size() << std::endl;
            table_file << "1,3\t" << vcf_names[0] << "," << vcf_names[2] << "\t"
                      << int_02.size() << std::endl;
            table_file << "2,3\t" << vcf_names[1] << "," << vcf_names[2] << "\t"
                      << int_12.size() << std::endl;
            table_file << "1,2,3\t" << vcf_names[0] << "," << vcf_names[1] << "," << vcf_names[2] << "\t"
                      << int_all.size() << std::endl;
        }
        
        table_file.close();
        std::cout << "Summary table written to: consensus_table.txt" << std::endl;
    }

    void printSummaryTable(const std::vector<std::set<std::string> >& variant_sets) {
        std::cout << "\n=== VCF Comparison Summary ===" << std::endl;
        std::cout << std::setw(30) << "Comparison" << std::setw(15) << "Count" << std::endl;
        std::cout << std::string(45, '-') << std::endl;
        
        // Individual VCF counts
        for (size_t i = 0; i < variant_sets.size(); i++) {
            std::cout << std::setw(30) << ("Total in " + vcf_names[i]) 
                      << std::setw(15) << variant_sets[i].size() << std::endl;
        }
        
        if (variant_sets.size() == 2) {
            // Two-way comparison
            std::set<std::string> intersection;
            std::set_intersection(variant_sets[0].begin(), variant_sets[0].end(),
                                variant_sets[1].begin(), variant_sets[1].end(),
                                std::inserter(intersection, intersection.begin()));
            
            std::cout << std::setw(30) << ("Common in " + vcf_names[0] + " & " + vcf_names[1])
                      << std::setw(15) << intersection.size() << std::endl;
                      
        } else if (variant_sets.size() == 3) {
            // Three-way comparison
            std::set<std::string> int_01, int_02, int_12, int_all;
            
            // Pairwise intersections
            std::set_intersection(variant_sets[0].begin(), variant_sets[0].end(),
                                variant_sets[1].begin(), variant_sets[1].end(),
                                std::inserter(int_01, int_01.begin()));
            
            std::set_intersection(variant_sets[0].begin(), variant_sets[0].end(),
                                variant_sets[2].begin(), variant_sets[2].end(),
                                std::inserter(int_02, int_02.begin()));
            
            std::set_intersection(variant_sets[1].begin(), variant_sets[1].end(),
                                variant_sets[2].begin(), variant_sets[2].end(),
                                std::inserter(int_12, int_12.begin()));
            
            // Three-way intersection
            std::set_intersection(int_01.begin(), int_01.end(),
                                variant_sets[2].begin(), variant_sets[2].end(),
                                std::inserter(int_all, int_all.begin()));
            
            std::cout << std::setw(30) << ("Common in " + vcf_names[0] + " & " + vcf_names[1])
                      << std::setw(15) << int_01.size() << std::endl;
            std::cout << std::setw(30) << ("Common in " + vcf_names[0] + " & " + vcf_names[2])
                      << std::setw(15) << int_02.size() << std::endl;
            std::cout << std::setw(30) << ("Common in " + vcf_names[1] + " & " + vcf_names[2])
                      << std::setw(15) << int_12.size() << std::endl;
            std::cout << std::setw(30) << "Common in all three"
                      << std::setw(15) << int_all.size() << std::endl;
        }
        
        std::cout << std::string(45, '-') << std::endl;
    }
    
    void generateConsensusVCF(const std::vector<std::set<std::string> >& variant_sets,
                             const std::vector<std::map<std::string, Variant> >& variant_maps) {
        
        std::set<std::string> consensus_variants;
        
        if (variant_sets.size() == 2) {
            // Intersection of two sets
            std::set_intersection(variant_sets[0].begin(), variant_sets[0].end(),
                                variant_sets[1].begin(), variant_sets[1].end(),
                                std::inserter(consensus_variants, consensus_variants.begin()));
        } else if (variant_sets.size() == 3) {
            // Intersection of all three sets
            std::set<std::string> temp_intersection;
            std::set_intersection(variant_sets[0].begin(), variant_sets[0].end(),
                                variant_sets[1].begin(), variant_sets[1].end(),
                                std::inserter(temp_intersection, temp_intersection.begin()));
            
            std::set_intersection(temp_intersection.begin(), temp_intersection.end(),
                                variant_sets[2].begin(), variant_sets[2].end(),
                                std::inserter(consensus_variants, consensus_variants.begin()));
        }
        
        // Convert consensus variants to vector and sort
        std::vector<Variant> consensus_vars;
        for (std::set<std::string>::const_iterator it = consensus_variants.begin(); 
             it != consensus_variants.end(); ++it) {
            const std::string& key = *it;
            // Use variant from first VCF file (they should be identical)
            std::map<std::string, Variant>::const_iterator variant_it = variant_maps[0].find(key);
            if (variant_it != variant_maps[0].end()) {
                consensus_vars.push_back(variant_it->second);
            }
        }
        
        // Sort variants by chromosome and position
        std::sort(consensus_vars.begin(), consensus_vars.end(),
                 VariantComparator());
        
        // Write consensus VCF
        std::string output_filename = compress_output ? "consensus.vcf.gz" : "consensus.vcf";
        
        if (compress_output) {
            // Write to temporary file first, then compress
            std::string temp_filename = "consensus_temp.vcf";
            std::ofstream temp_file(temp_filename.c_str());
            
            if (!temp_file.is_open()) {
                std::cerr << "Error: Cannot create temporary file " << temp_filename << std::endl;
                return;
            }
            
            writeVCFContent(temp_file, consensus_vars);
            temp_file.close();
            
            // Compress using bgzip
            std::string compress_cmd = "bgzip -c \"" + temp_filename + "\" > \"" + output_filename + "\"";
            int result = system(compress_cmd.c_str());
            
            if (result != 0) {
                std::cerr << "Error: bgzip compression failed. Make sure bgzip is installed." << std::endl;
                std::cerr << "Writing uncompressed file instead." << std::endl;
                std::string mv_cmd = "mv \"" + temp_filename + "\" consensus.vcf";
                system(mv_cmd.c_str());
                output_filename = "consensus.vcf";
            } else {
                std::remove(temp_filename.c_str());
            }
        } else {
            std::ofstream outfile(output_filename.c_str());
            if (!outfile.is_open()) {
                std::cerr << "Error: Cannot create output file " << output_filename << std::endl;
                return;
            }
            writeVCFContent(outfile, consensus_vars);
            outfile.close();
        }
        
        std::cout << "\nConsensus VCF written to: " << output_filename << std::endl;
        std::cout << "Number of consensus variants: " << consensus_vars.size() << std::endl;
    }
    
    struct VariantComparator {
        bool operator()(const Variant& a, const Variant& b) const {
            if (a.chrom != b.chrom) return a.chrom < b.chrom;
            return a.pos < b.pos;
        }
    };
    
    void writeVCFContent(std::ofstream& outfile, const std::vector<Variant>& variants) {
        // Write header
        for (size_t i = 0; i < header_lines.size(); i++) {
            outfile << header_lines[i] << std::endl;
        }
        
        // Write variants
        for (size_t i = 0; i < variants.size(); i++) {
            const Variant& var = variants[i];
            outfile << var.chrom << "\t" << var.pos << "\t" << var.id << "\t"
                    << var.ref << "\t" << var.alt << "\t" << var.qual << "\t"
                    << var.filter << "\t" << var.info << "\t" << var.format;
            
            for (size_t j = 0; j < var.samples.size(); j++) {
                outfile << "\t" << var.samples[j];
            }
            outfile << std::endl;
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " [--compress] <vcf1> <vcf2> [vcf3]" << std::endl;
        std::cerr << "Compare 2 or 3 VCF files and generate consensus variants" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  --compress    Compress output VCF with bgzip" << std::endl;
        std::cerr << "  Input files can be .vcf or .vcf.gz (bgzip compressed)" << std::endl;
        return 1;
    }
    
    bool compress_output = false;
    int vcf_start_idx = 1;
    
    // Check for compress flag
    if (argc > 3 && std::string(argv[1]) == "--compress") {
        compress_output = true;
        vcf_start_idx = 2;
    }
    
    // Validate we have the right number of VCF files
    int num_vcfs = argc - vcf_start_idx;
    if (num_vcfs < 2 || num_vcfs > 3) {
        std::cerr << "Error: Need 2 or 3 VCF files for comparison" << std::endl;
        return 1;
    }
    
    VCFReader reader;
    VCFComparator comparator(compress_output);
    
    // Read VCF files
    for (int i = vcf_start_idx; i < argc; i++) {
        std::string filename = argv[i];
        std::cout << "Reading " << filename << "..." << std::endl;
        
        std::vector<Variant> variants = reader.readVCF(filename);
        if (variants.empty()) {
            std::cerr << "Error reading " << filename << std::endl;
            return 1;
        }
        
        comparator.addVCF(variants, filename, reader.getHeaderLines(), reader.getSampleNames());
    }
    
    // Perform comparison
    comparator.compare();
    
    return 0;
}
