/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file CLASSIFY.cpp
 * @author Francisco Javier Ortiz Molinero <javierom@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include <iostream>
#include "KmerCounter.h"

using namespace std;

/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(ostream& outputStream) {
    outputStream << "ERROR in CLASSIFY parameters" << endl;
    outputStream << "Run with the following parameters:" << endl;
    outputStream << "CLASSIFY [-k kValue] [-n nucleotidesSet] <file.dna> <profile1.prf> [<profile2.prf> <profile3.prf> ....]" << endl;
    outputStream << endl;
    outputStream << "Parameters:" << endl;
    outputStream << "-k kValue: number of nucleotides in a kmer (5 by default)" << endl;
    outputStream << "-n nucletiodesSet: set of possible nucleotides in a kmer (ACGT by default). "
                 << "It is used when learning a model for <file.dna>. " 
                 << "Note that the characters should be provided in uppercase" << endl;
    outputStream << "<profile1.prf> [<profile2.prf> <profile3.prf> ....] ....: "
                 << "names of the Profile models (at least one is mandatory)" << endl;
    outputStream << endl;
    outputStream << "This program obtains the identifier of the closest profile to the input DNA file" << endl;
    outputStream << endl;
}

bool min (double a, double b) {
    return a < b;
}
bool max (double a, double b) {
    return a > b;
}

int PosMinMax (double distances[], int num_args, bool(*MinMax)(double,double)) {
    double minmax = distances[0];
    int pos_minmax = 0;

    for (int i = 1; i < num_args; i++) {
        if ((*MinMax)(distances[i], minmax)) {
            minmax = distances[i];
            pos_minmax = i;
        }
    }
    
    return pos_minmax;
}

/**
 * This program prints the profile identifier of the closest profile model
 * for an input DNA file (<file.dna>) among the set of provided models:
 * <profile1.prf>, <profile2.prf>, ...
 * The program uses the KmerCounter class to obtain a Profile for the input
 * file <file.dna>. That Profile should be zipped, to eliminate kmers with
 * any missing nucleotide, and sorted in decreasing order of frequency of
 * kmers. After that, the program compares the learned Profile with the ones
 * provided by the arguments <profile1.prf> [<profile2.prf> <profile3.prf> ....]
 * It classifies the input DNA file with the identifier of the Profile with
 * a minor distance.
 * 
 * This program assumes that the profile files are already normalized and 
 * sorted by frequency. This is not checked in this program. Unexpected results
 * will be obtained if those conditions are not met.
 * 
 * Running sintax:
 * > CLASSIFY [-k kValue] [-n nucleotidesSet] <file.dna> <profile1.prf> [<profile2.prf> <profile3.prf> ....]
 * 
 * Running example: 
 * > CLASSIFY ../Genomes/human_chr6_s60000_l500000.dna ../Genomes/brewers_yeast_chrVII.s1_l500000.prf ../Genomes/chimpanzee_chr9_s1_l500000.prf ../Genomes/covidFullGenomeDNA.prf ../Genomes/drosophila_chr2L_s1_l500000.prf ../Genomes/ebolaFullGenomeDNA.prf ../Genomes/human_chr9_s10000_l500000.prf ../Genomes/monkeypoxFullGenomeDNA.prf ../Genomes/mouse_chr6_s3050050_l500000.prf ../Genomes/nematode_chrI_s1l500000.prf ../Genomes/rat_chr6_s1l500000.prf ../Genomes/zebrafish_chr6_s1l500000.prf
Distance to ../Genomes/brewers_yeast_chrVII.s1_l500000.prf (saccharomyces cerevisiae): 0.20294
Distance to ../Genomes/chimpanzee_chr9_s1_l500000.prf (pan troglodytes): 0.0643864
Distance to ../Genomes/covidFullGenomeDNA.prf (severe acute respiratory syndrome coronavirus 2): 0.194633
Distance to ../Genomes/drosophila_chr2L_s1_l500000.prf (drosophila melanogaster): 0.189238
Distance to ../Genomes/ebolaFullGenomeDNA.prf (ebolavirus zaire): 0.179686
Distance to ../Genomes/human_chr9_s10000_l500000.prf (homo sapiens): 0.0557804
Distance to ../Genomes/monkeypoxFullGenomeDNA.prf (monkey pox virus): 0.262987
Distance to ../Genomes/mouse_chr6_s3050050_l500000.prf (mus musculus): 0.088129
Distance to ../Genomes/nematode_chrI_s1l500000.prf (caenorhabditis elegans): 0.221075
Distance to ../Genomes/rat_chr6_s1l500000.prf (rattus norvegicus): 0.111126
Distance to ../Genomes/zebrafish_chr6_s1l500000.prf (danio rerio): 0.145231

Final decision: homo sapiens with a distance of 0.0557804
 * 
 * @param argc The number of command line parameters
 * @param argv The vector of command line parameters (cstrings)
 * @return 0 If there is no error; a value > 0 if error
 */

#include <regex>
#include <cmath>

const string VALID_NUCLEOTIDES_ADN = "ACGT";
const string COMPLEMENTARY_NUCLEOTIDES_ADN = "TGCA";
const string VALID_NUCLEOTIDES_ARN = "ACGU";
const string COMPLEMENTARY_NUCLEOTIDES_ARN = "UGCA";
const int DIM_VECTOR_KMER_FREQ = 2000;
#define ENDL "\n"

int main(int argc, char *argv[]) {
    // Process the main() arguments
    if (argc < 3) {
        showEnglishHelp(cerr);
        return 1;
    }
    
    int num_args, first_arg;
    int k = 5;
    string n = KmerCounter::DEFAULT_VALID_NUCLEOTIDES;
    
    bool sigo = true;
    int i = 1;
    while (sigo && i < argc-1) {
        if (string(argv[i]).at(0) == '-') {
            if (string(argv[i]) == "-k") {
                k = stoi(argv[i+1]);
                i += 2;
            }
            else if (string(argv[i]) == "-n") {
                n = argv[i+1];
                i += 2;
            }
            else {
                showEnglishHelp(cerr);
                return 1;
            }
        }
        else {
            sigo = false;
        }
    }
    
    num_args = argc - i - 1;
    first_arg = i + 1;
    
    if (string(argv[argc-1]).at(0)== '-' || string(argv[argc-2]).at(0) == '-') {
        showEnglishHelp(cerr);
        return 1;
    }
    
    // Calculate the kmer frecuencies of the input genome file using 
    //    a KmerCounter object
    KmerCounter kc(k,n);
    kc.calculateFrequencies(argv[first_arg-1]);
    
    // Obtain a Profile object for the input genome from the KmerCounter object
    Profile prf = kc.toProfile();
    
    // Zip the for the input genome Profile object
    prf.zip();
    
    // Sort the for the input genome Profile object
    prf.sort();
    
    // Use a loop to print the distance from the input genome to 
    //   each one of the provided profile models
    Profile* arrayProfiles;
    double* distances;
    arrayProfiles = new Profile[num_args];
    distances = new double[num_args];
    
    for (int k = 0; k < num_args; k++) {
        arrayProfiles[k].load(argv[k+ first_arg]);
    }
    
    for (int j = 0; j < num_args; j++) {
        distances[j] = prf.getDistance(arrayProfiles[j]);
        cout << "Distance to " << argv[j + first_arg] << " ("
             << arrayProfiles[j].getProfileId() << "): " 
             << distances[j] << endl;
    }
    
    // Print the identifier and distance to the closest profile
    int pos_min = PosMinMax(distances, num_args, min);    
    cout << "Final decision: "  << arrayProfiles[pos_min].getProfileId() 
         << "with a distance of " << distances[pos_min] << endl;;
    
    delete[] arrayProfiles;
    delete[] distances;
    
    return 0;
}

