/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerCounter.cpp
 * @author F. Javier Ortiz Molinero <javierom@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include <fstream>

#include "KmerCounter.h"

using namespace std;

/**
 * DEFAULT_VALID_NUCLEOTIDES is a c-string that contains the set of characters
 * that will be considered as valid nucleotides. 

 * The constructor of the class KmerCounter uses this c-string as a 
 * default parameter. It is possible to use a different c-string if that
 * constructor is used with a different c-string
 */
const char* const KmerCounter::DEFAULT_VALID_NUCLEOTIDES="ACGT";

KmerCounter::KmerCounter(int k, const std::string& validNucleotides) {
    _k = k;
    _validNucleotides = validNucleotides;
    _allNucleotides = Kmer::MISSING_NUCLEOTIDE + validNucleotides;
    _frequency = new int*[getNumRows()];
    for (int i = 0; i < getNumRows(); i++) {
        _frequency[i] = new int[getNumCols()];
    }
     
    initFrequencies();
}

KmerCounter::KmerCounter(const KmerCounter& orig) {
    _k = orig._k;
    _validNucleotides = orig._validNucleotides;
    _allNucleotides = orig._allNucleotides;
    _frequency = new int*[orig.getNumRows()];
    for (int i = 0; i < orig.getNumRows(); i++) {
        _frequency[i] = new int [orig.getNumCols()];
    }
    
    for (int i = 0; i < getNumRows(); i++) {
        for (int j = 0; j < getNumCols(); j++) {
            _frequency[i][j] = orig._frequency[i][j];
        }
    }
}

KmerCounter::~KmerCounter() {
    deallocate();
    _k = 0;
    _validNucleotides = DEFAULT_VALID_NUCLEOTIDES;
    _allNucleotides = Kmer::MISSING_NUCLEOTIDE + _validNucleotides;
}

int KmerCounter::getNumNucleotides() const {
    return _allNucleotides.length();
}

int KmerCounter::getK() const {
    return _k;
}

int KmerCounter::getNumKmers() const {
    return pow(_allNucleotides.length(),_k);
}

int KmerCounter::getNumberActiveKmers() const{
    int activeKmers = 0;
    
    for (int i = 0; i < getNumRows(); i++) {
        for (int j = 0; j < getNumCols(); j++) {
            if (_frequency[i][j] != 0)
                activeKmers++;
        }
    }
    
    return activeKmers;
}

std::string KmerCounter::toString() const{
    string outputString = _allNucleotides + " " + to_string(_k) + "\n";
    
    for(int row=0; row<this->getNumRows(); row++){
        for(int col=0; col<this->getNumCols(); col++){
            outputString += to_string((*this)(row,col)) + " ";
        }
        outputString += "\n";
    }
    
    return outputString;
}

void KmerCounter::increaseFrequency(const Kmer& kmer, int frequency) {
    for (int i = 0; i < int(kmer.toString().length()); i++) {
        if (!IsValidNucleotide(kmer.toString().at(i), _allNucleotides))
            throw std::invalid_argument(string("void "
                    "KmerCounter::increaseFrequency(const Kmer& kmer, int "
                    "frequency = 1): the given kmer contains invalid "
                    "nucleotides"));
    }
    
    int row, column;
    getRowColumn(kmer, row, column);
    _frequency[row][column] += frequency;
}

KmerCounter& KmerCounter::operator=(const KmerCounter& orig) {
    if (this != &orig) {
        deallocate();
        _k = orig._k;
        _validNucleotides = orig._validNucleotides;
        _allNucleotides = orig._allNucleotides;
        _frequency = new int*[orig.getNumRows()];
        for (int i = 0; i < orig.getNumRows(); i++)
            _frequency[i] = new int [orig.getNumCols()];
    
        for (int i = 0; i < getNumRows(); i++) {
            for (int j = 0; j < getNumCols(); j++) {
                _frequency[i][j] = orig._frequency[i][j];
            }
        }
    }
    
    return *this;
}

KmerCounter& KmerCounter::operator+=(const KmerCounter& kc) {
    if(_k != kc._k)
        throw std::invalid_argument(string("KmerCounter& "
                "KmerCounter::operator+=(const KmerCounter& kc): the number "
                "of nucleotides is different in each KmerCounter"));
    
    for (int i = 0; i < getNumRows(); i++) {
        for (int j = 0; j < getNumCols(); j++) {
            _frequency[i][j] += kc._frequency[i][j];
        }
    }
    
    return *this;
}

void KmerCounter::calculateFrequencies(const char* fileName) {
    ifstream input;
    input.open(fileName);
    
    if (input) {
        initFrequencies();
        string kmer;
        input >> kmer;
        int kmer_l = kmer.length();
        int nkmers;
        
        
        nkmers = kmer_l - _k + 1;
        for (int i = 0; i < nkmers; i++) {
            string subkmer = kmer.substr(i,_k);
            Kmer k(subkmer);
            k.normalize(_validNucleotides);
            increaseFrequency(k);
        }
        
        input.close();
    }
    else
        throw ios_base::failure(string("void "
                "KmerCounter::calculateFrequencies(const char* fileName): "
                "the given file cannot be opened\n"));
}

Profile KmerCounter::toProfile() const {
    Profile p;
    KmerFreq kf;
    for (int i = 0; i < getNumRows(); i++) {
        for (int j = 0; j < getNumCols(); j++) {
            if (_frequency[i][j] > 0) {
                Kmer k(getKmer(i,j));
                kf.setFrequency(_frequency[i][j]);
                kf.setKmer(k);
                p += kf;
            }
        }
    }
    
    return p;
}

int KmerCounter::getNumRows() const {
    return (pow(_allNucleotides.length(), (_k+1)/2));
}

int KmerCounter::getNumCols() const {
    return (pow(_allNucleotides.length(), (_k -((_k+1)/2))));
}

int KmerCounter::getIndex(const std::string& kmer) const{
    int index = 0;
    int base = 1;

    for (size_t i = 0; i < kmer.size(); i++) {
        size_t pos = _allNucleotides.find(kmer[kmer.size()-i-1]);
        if (pos == string::npos)
            return -1;
        index += pos * base;
        base *= _allNucleotides.size();
    }
    return index;
}

string KmerCounter::getInvertedIndex(int index, int nCharacters) const {
    string result(nCharacters, Kmer::MISSING_NUCLEOTIDE);

    for (int i = result.size(); i > 0; i--) {
        result[i - 1] = _allNucleotides[index % _allNucleotides.size()];
        index = index / _allNucleotides.size();
    }
    return result;
}

void KmerCounter::getRowColumn(const Kmer& kmer, int& row, int& column) const {
    string k = kmer.toString();
    string k1 = k.substr(0,(_k+1)/2);
    string k2 = k.substr((_k+1)/2,k.length());
    row = getIndex(k1);
    column = getIndex(k2);
}

Kmer KmerCounter::getKmer(int row, int column) const {
    if(row >= getNumRows() || column >= getNumCols() || row < 0 || column < 0)
        throw std::invalid_argument(string("Kmer KmerCounter::getKmer(int row, "
                "int column) const: the row or the column is not "
                "in the range"));
    
    string k1 = getInvertedIndex(row, (_k+1)/2);
    string k2 = getInvertedIndex(column, (_k+1)/2);
    string k = k1 + k2;
    Kmer kmer(k);
    return kmer;
}

void KmerCounter::initFrequencies() {
    for (int i = 0; i < getNumRows(); i++) {
        for (int j = 0; j < getNumCols(); j++) {
            _frequency[i][j] = 0;
        }
    }
}

const int& KmerCounter::operator()(int row, int column) const {
    return _frequency[row][column];
}

int& KmerCounter::operator()(int row, int column) {
    return _frequency[row][column];
}

void KmerCounter::deallocate() {
    for (int i = 0; i < getNumRows(); i++) {
        delete[] _frequency[i];
    }
    delete[] _frequency;
    _frequency = nullptr;
}