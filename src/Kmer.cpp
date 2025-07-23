/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file Kmer.cpp
 * @author Francisco Javier Ortiz Molinero <javierom@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "Kmer.h"

using namespace std;

#include <iostream>
#include <string>

Kmer::Kmer(int k) {
    if (k <= 0) {
        throw std::invalid_argument(
                std::string("Kmer(int k=1): ") + 
                "k must be a number greater than zero");
    }
    
    _text = std::string(k, MISSING_NUCLEOTIDE);
    
}

Kmer::Kmer(const std::string& text) {
    if (text.length() == 0) {
        throw std::invalid_argument(
                std::string("Kmer(const std::string& text): ") + 
                "text is an empty string");
    }
    else {
        _text = text;
    }
}

int Kmer::getK() const {
    return _text.length();
}

int Kmer::size() const {
    return _text.length();
}

std::string Kmer::toString() const {
    return _text;
}

const char& Kmer::at(int index) const {
    if (index < 0 || index >= getK()) {
        throw std::out_of_range(
                std::string("const char& Kmer::at(int index) const: ")
                + "invalid position" + std::to_string(index));
    }

    return _text.at(index);

}

char& Kmer::at(int index) {
    if (index < 0 || index >= getK()) {
        throw std::out_of_range(
                std::string("char& Kmer::at(int index): ") +
                "invalid position" + std::to_string(index));
    }
    
    return _text.at(index);
}

void Kmer::toLower() {
    for (int i = 0; i < getK(); i++) {
        at(i) = std::tolower(at(i));
    }
}

void Kmer::toUpper() {
    for (int i = 0; i < getK(); i++) {
        at(i) = std::toupper(at(i));
    }
}

void Kmer::normalize(const std::string& validNucleotides) {
    toUpper();
    
    for (int i = 0; i < getK(); i++) {
        if (!IsValidNucleotide(_text.at(i), validNucleotides))
            at(i) = MISSING_NUCLEOTIDE;
    }
}

Kmer Kmer::complementary(const std::string& nucleotides, 
        const std::string& complementaryNucleotides) const {
    if (nucleotides.length() != complementaryNucleotides.length()) {
        throw std::invalid_argument(
                std::string("Kmer Kmer::complementary(const std::string& "
                "nucleotides, const std::string& "
                "complementaryNucleotides) const: ")
                + "the sizes of nucleotides and complementaryNucelotides " +
                "are not the same");
    }
    
    Kmer CKmer(_text);
    
    for (int i = 0; i < (int)_text.length(); i++) {
        if (IsValidNucleotide(_text.at(i), nucleotides)) {
            CKmer.at(i) = 
                complementaryNucleotides.at(nucleotides.find(_text.at(i)));
        }
    }

    return CKmer;
    
}

void Kmer::write(std::ostream& outputStream) const {
    outputStream << *this;
}

void Kmer::read(std::istream& inputStream) {
    inputStream >> _text;
}

bool IsValidNucleotide(char nucleotide, const std::string& validNucleotides) {
    return (int(validNucleotides.find(nucleotide)) != -1);
}

void ToLower(Kmer& kmer) {
    for (int i = 0; i < kmer.getK(); i++) {
        kmer.at(i) = std::tolower(kmer.at(i));
    }
}

void ToUpper(Kmer& kmer) {
    for (int i = 0; i < kmer.getK(); i++) {
        kmer.at(i) = std::toupper(kmer.at(i));
    }
}

std::ostream& operator<<(std::ostream& os, const Kmer& kmer) {
    os << kmer.toString();
    return os;
}

std::istream& operator>>(std::istream& is, Kmer& kmer) {
    string text;
    is >> text;
    Kmer k(text);
    kmer = k;
    
    return is;
}
