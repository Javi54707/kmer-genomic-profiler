touch tests//.timeout
CMD=" /home/usuario/Documentos/Copia/NetBeansProjects/Kmer5/./dist/LEARN/GNU-Linux/LEARN  -b -p 'homo sapiens' -o tests/output/human_chr9_s10000_l500000.prf ../Genomes/human_chr9_s10000_l500000.dna; dist/CLASSIFY/GNU-Linux/CLASSIFY ../Genomes/human_chr6_s60000_l500000.dna ../Genomes/brewers_yeast_chrVII.s1_l500000.prf ../Genomes/chimpanzee_chr9_s1_l500000.prf ../Genomes/covidFullGenomeDNA.prf ../Genomes/drosophila_chr2L_s1_l500000.prf ../Genomes/ebolaFullGenomeDNA.prf tests/output/human_chr9_s10000_l500000.prf ../Genomes/monkeypoxFullGenomeDNA.prf ../Genomes/mouse_chr6_s3050050_l500000.prf ../Genomes/nematode_chrI_s1l500000.prf ../Genomes/rat_chr6_s1l500000.prf ../Genomes/zebrafish_chr6_s1l500000.prf  1> tests//.out26 2>&1"
eval $CMD
rm tests//.timeout
