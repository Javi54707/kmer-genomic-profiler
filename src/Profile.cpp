/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "Profile.h"
#include <fstream>

using namespace std;

const string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";
const string Profile::MAGIC_STRING_B="MP-KMER-B-1.0";

Profile::Profile(): _profileId("unknown"), _size(0),
    _capacity(INITIAL_CAPACITY){
    _vectorKmerFreq = new KmerFreq[_capacity];
}

Profile::Profile(int size) {
    if (size < 0)
        throw out_of_range(string("Profile::Profile(int size): "
                "size must be at least 0"));
    
    _profileId = "unknown";
    _size = size;
    _capacity = size;
    _vectorKmerFreq = new KmerFreq[_size];
}

Profile::Profile(const Profile& orig) {
    _profileId = orig._profileId;
    _size = orig._size;
    _capacity = orig._capacity;
    _vectorKmerFreq = new KmerFreq[_capacity];
    
    for (int i = 0; i < _size; i++) {
        _vectorKmerFreq[i] = orig._vectorKmerFreq[i];
    }
}

Profile::~Profile() {
    deallocate();
    _size = 0;
    _capacity = INITIAL_CAPACITY;
}

Profile& Profile::operator=(const Profile& orig) {
    if (this != &orig) {
        deallocate();
        _profileId = orig._profileId;
        _size = orig._size;
        _capacity = orig._capacity;
        _vectorKmerFreq = new KmerFreq[_capacity];
        for (int i = 0; i < _size; i++) {
            _vectorKmerFreq[i] = orig._vectorKmerFreq[i];
        }
    }
    
    return *this;
}

const string& Profile::getProfileId() const {
    return _profileId;
}

void Profile::setProfileId(const std::string& id) {
    _profileId = id;
}

const KmerFreq& Profile::at(int index) const {
    if (index < 0 || index >= _size)
        throw out_of_range(string("const KmerFreq& Profile::at(int index) "
                "const: index must be between 0 and _size"));
    
    return _vectorKmerFreq[index];
}

KmerFreq& Profile::at(int index) {
    if (index < 0 || index >= _size)
        throw out_of_range(string("KmerFreq& Profile::at(int index): "
                "index must be between 0 and _size"));
    
    return _vectorKmerFreq[index];
}

int Profile::getSize() const {
    return _size;
}

int Profile::getCapacity() const {
    return _capacity;
}

double Profile::getDistance(const Profile& otherProfile) const {
    if (otherProfile._size <= 0 || _size <= 0)
        throw invalid_argument(string("double Profile::getDistance(const "
                "Profile& otherProfile) const: one of the Profiles (or both) "
                "is empty"));
    
    double distance = 0;
    
    for (int i = 0; i < _size; i++) {
        int rank_2 = otherProfile.findKmer(_vectorKmerFreq[i].getKmer());
        
        if (rank_2 == -1)
            rank_2 = otherProfile._size;
        
        distance += abs(i - rank_2);
    }
    
    distance /= _size*otherProfile._size;
    
    return distance;
}

int Profile::findKmer(const Kmer& kmer, int initialPos, int finalPos) const {
    int pos_kmer = -1;
    int i = initialPos;
    bool sigo = true;
    
    while (i <= finalPos && sigo) {
        if (kmer.toString() == _vectorKmerFreq[i].getKmer().toString()) {
            pos_kmer = i - initialPos;
            sigo = false;
        }
        
        i++;
    }
    
    return pos_kmer;
}

int Profile::findKmer(const Kmer& kmer) const {
    return findKmer(kmer, 0, _size - 1);
}

string Profile::toString() const {
    string profile = _profileId + '\n' + to_string(_size) + '\n';
    for (int i = 0; i < _size; i++) {
        profile += _vectorKmerFreq[i].toString() + '\n';
    }
    
    return profile;
}

void Profile::sort() {
    for (int i = 0; i < _size; i++) {
        int pos_max = i;
        for (int j = i + 1; j < _size; j++) {
            if (_vectorKmerFreq[j] > _vectorKmerFreq[pos_max]) {
                pos_max = j;
            }
        }
        
        KmerFreq temp = _vectorKmerFreq[i];
        _vectorKmerFreq[i] = _vectorKmerFreq[pos_max];
        _vectorKmerFreq[pos_max] = temp;
    }
}

void Profile::save(const char fileName[], char mode) const {
    ofstream output;
    output.open(fileName);
    
    if (output) {
        if (mode == 't') {
            output << MAGIC_STRING_T << endl;
            output << *this;
        }
        else if (mode == 'b') {
            output << MAGIC_STRING_B << endl;
            output << _profileId << endl;
            output << _size << endl;
            
            for (int i = 0; i < _size; i++) {
                _vectorKmerFreq[i].write(output);
                output << endl;
            }
        }
        else {
            output.close();
            throw std::invalid_argument(string("void Profile::save(const char* "
                    "fileName, char mode) const: mode is not valid"));
        }

        
        if (!output) {
            output.close();
            throw ios_base::failure(string("void Profile::save(const char*"
                    " fileName) const: an error ocurred while writing in "
                    "the file\n"));
        }
        
        output.close();
    }
    else {
        throw ios_base::failure(string("void Profile::save(const char*"
                " fileName[]) const: the given file cannot be opened\n"));
    }
}

void Profile::load(const char fileName[]) {
    ifstream input;
    input.open(fileName);
    
    if (!input) {
        throw ios_base::failure(string("void Profile::load(const char* "
                "fileName[]): the given file cannot be opened\n"));
    }
    
    deallocate();
    _size = 0;
    _capacity = INITIAL_CAPACITY;
    string magic_string;
    input >> magic_string;
    input.get();

    if (magic_string == MAGIC_STRING_T) {
        input >> *this;
    }
    else if (magic_string == MAGIC_STRING_B) {
        deallocate();
        string id;
        int size;
        getline(input,id);
        _profileId = id;
        input >> size;
        allocate(size);
        KmerFreq kf;

        for (int i = 0; i < size; i++) {
            kf.read(input);
            append(kf);
        }
    }
    else {
        input.close();
        throw invalid_argument(string("void Profile::load(const char* "
                "fileName): an invalid magic string has been found in "
                "the given file"));
    }

    input.close(); 
}

void Profile::append(const KmerFreq& kmerFreq) {  
    int pos_kmer = findKmer(kmerFreq.getKmer());
    if (pos_kmer != -1) {
        int freq = kmerFreq.getFrequency() + 
                   _vectorKmerFreq[pos_kmer].getFrequency();
        _vectorKmerFreq[pos_kmer].setFrequency(freq);
    }
    else {
        if (_size ==_capacity) {
            _capacity += BLOCK_SIZE;
            reallocate();
        }
        _vectorKmerFreq[_size] = kmerFreq;
        _size++;
    }
}

void Profile::normalize(const string& validNucleotides) {
    // Loop to traverse and normalize each one of the kmers in array
    for (int i = 0; i < _size; i++) {
        // Normalize kmer i
        Kmer aux = _vectorKmerFreq[i].getKmer();
        aux.normalize(validNucleotides);
        _vectorKmerFreq[i].setKmer(aux);
    }
    
    int i = 1;
    while (i < _size) {
        int pos = findKmer(_vectorKmerFreq[i].getKmer(), 0, i - 1);
        if (pos != -1) {
            int freq = _vectorKmerFreq[i].getFrequency() + 
                    _vectorKmerFreq[pos].getFrequency();
            _vectorKmerFreq[pos].setFrequency(freq);
            deletePos(i);
        }
        else
            i++;
    }
}

void Profile::deletePos(int pos) {
    if (pos >= _size || pos < 0) {
        throw out_of_range(
                string("void Profile::deletePos(int pos): pos is not on "
                "the range from 0 to _size-1 (both included)"));
    }
    
    for (int i = pos; i < _size - 1; i++) {
        _vectorKmerFreq[i] = _vectorKmerFreq[i + 1];
    }
    
    _size--;
}

void Profile::zip(const bool deleteMissing, int lowerBound) {
    int i = 0;
    while (i < _size) {
        if ((deleteMissing &&
          _vectorKmerFreq[i].getKmer().toString().find(Kmer::MISSING_NUCLEOTIDE)
            != std::string::npos) || 
            (_vectorKmerFreq[i].getFrequency() <= lowerBound)) {
            deletePos(i);
        }
        else {
            // No se elimina este elemento, avanzamos al siguiente.
            i++;
        }
    }
}

void Profile::join(const Profile& profile) {
    for (int i = 0; i < profile.getSize(); i++) {
        append(profile.at(i));
    }
}

const KmerFreq& Profile::operator[](int index) const {
    return _vectorKmerFreq[index];
}

KmerFreq& Profile::operator[](int index) {
    return _vectorKmerFreq[index];
}

Profile& Profile::operator+=(const KmerFreq& kmerFreq) {
    append(kmerFreq);
    return *this;
}

Profile& Profile::operator+=(const Profile& profile) {
    join(profile);
    return *this;
}

void Profile::allocate(int size) {
    _vectorKmerFreq = new KmerFreq[size];
}

void Profile::deallocate() {
    delete[] _vectorKmerFreq;
    _vectorKmerFreq = nullptr;
}

void Profile::reallocate() {
    KmerFreq* temp;
    temp = new KmerFreq[_capacity];
    copy(temp);
    deallocate();
    _vectorKmerFreq = temp;
}

void Profile::copy(KmerFreq copy[]) {
    for (int i = 0; i < _size; i++) {
        copy[i] = _vectorKmerFreq[i];
    }
}

std::ostream& operator<<(std::ostream& os, const Profile& profile) {
    os << profile.toString();
    return os;
}

std::istream& operator>>(std::istream& is, Profile& profile) {
    profile.~Profile();
    string id;
    int size;
    is >> id >> size;
    if (size < 0) {
        throw out_of_range(string("std::istream& operator>>(std::istream& is, "
                "Profile& profile): the size given must be positive"));
    }
    
    profile.setProfileId(id);
    
    for (int i = 0; i < size; i++) {
        KmerFreq kf;
        string kmer;
        int freq;
        is >> kmer >> freq;
        Kmer k(kmer);
        kf.setKmer(k);
        kf.setFrequency(freq);
        profile.append(kf);
    }
    
    return is;
}