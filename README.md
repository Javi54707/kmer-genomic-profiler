# Kmer Genomic Profiler

A C++ program for analyzing DNA sequences to generate k-mer frequency profiles for known species, and classify unknown samples by similarity. Developed as part of the “Metodología de la Programación” (MP) course at the University of Granada.

## 🧬 Features

- Reads genome files and extracts frequency counts of all k-length substrings (k-mers)
- Builds a profile for each known species
- Compares unknown genome sequences against profiles using distance metrics
- Designed to demonstrate modularity, abstraction and basic algorithmic efficiency

## 🛠️ Technologies

- C++
- Standard Library (file I/O, strings, maps)
- CLI-based input/output

## 📂 Structure

- `src/` – Source code (.cpp and .h files)
- `data/` – Input genome files for known species
- `output/` – Resulting classification or frequency tables

## ▶️ How to Compile
  ```bash
  g++ -std=c++11 -o kmer src/*.cpp
  ./kmer [input_file]
```

## 🎓 Academic Context
This project was developed during the 2023/24 academic year, developed as part of the “Metodología de la Programación” (MP) course of the Bachelor's Degree in Computer Engineering and Mathematics (UGR).

## 📜 License
Provided for educational purposes only.
Contact: javier.ortmol@gmail.com
