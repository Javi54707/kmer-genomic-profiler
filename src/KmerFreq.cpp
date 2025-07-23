/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author F. Javier Ortiz Molinero <javierom@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "KmerFreq.h"

using namespace std;

#include <iostream>
#include <string>

KmerFreq::KmerFreq(): _kmer(Kmer()), _frequency(0) {}

const Kmer& KmerFreq::getKmer() const{
    return _kmer;
}

int KmerFreq::getFrequency() const {
    return _frequency;
}

void KmerFreq::setKmer(const Kmer& kmer) {
    _kmer = kmer;
}

void KmerFreq::setFrequency(const int frequency) {
    if (frequency < 0)
        throw std::out_of_range(
                std::string("setFrequency(const int& frequency): ") + 
                "frequency is negative");
    
    _frequency = frequency;
}

std::string KmerFreq::toString() const {
    return _kmer.toString() + " " + std::to_string(_frequency);
}

void KmerFreq::write(std::ostream& outputStream) const {
    outputStream << *this;
}

void KmerFreq::read(std::istream& inputStream) {
    _kmer.read(inputStream);
    inputStream >> _frequency;
}

std::ostream& operator<<(std::ostream& os, const KmerFreq& kmerFreq) {
    os << kmerFreq.toString();
    return os;
}

std::istream& operator>>(std::istream& is, KmerFreq& kmerFreq) {
    int f;
    Kmer k;
    is >> k >> f;
    kmerFreq.setKmer(k);
    kmerFreq.setFrequency(f);
    
    return is;
}

bool operator>(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    bool gt = false;
    if (kmerFreq1.getFrequency() > kmerFreq2.getFrequency())
        gt = true;
    else
        if (kmerFreq1.getFrequency() == kmerFreq2.getFrequency() && 
                kmerFreq1.getKmer().toString() < kmerFreq2.getKmer().toString())
            gt = true;
        
    return gt;
}

bool operator<(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    return !(kmerFreq1 >= kmerFreq2);
}

bool operator==(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    bool eq = false;
    if (kmerFreq1.getFrequency() == kmerFreq2.getFrequency() && 
            kmerFreq1.getKmer().toString() == kmerFreq2.getKmer().toString())
        eq = true;
    
    return eq;
}

bool operator!=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    return !(kmerFreq1 == kmerFreq2);
}

bool operator<=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    return !(kmerFreq1 > kmerFreq2);
}

bool operator>=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {
    bool get = false;
    if (kmerFreq1.getFrequency() > kmerFreq2.getFrequency())
        get = true;
    else
        if( kmerFreq1.getFrequency() == kmerFreq2.getFrequency() && 
               kmerFreq1.getKmer().toString() <= kmerFreq2.getKmer().toString())
            get = true;
        
    return get;
}