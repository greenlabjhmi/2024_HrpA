# 2024_HrpA

Analysis of ribosome profiling data from three strains derived from E. coli MG1655: wild-type, delta-hrpA, and the double knockout delta-hrpA delta-smrB

Developed by Allen Buskirk using code originally by Fuad Mohammad and Nick Guydosh

Dependencies:

python 2.7
cutadapt 2.1
bowtie 1.2.2
BCBio 
Biopython
Jupyter notebook


The raw sequencing data are available as FASTQ files and processed WIG files at the GEO.

Processing of the FASTQ files, mapping with bowtie, and storing ribosome density are all described in the iPython notebook. The GFF file for annotation of version 2 of the MG1655 genome is given as coli2.gff. All strains also express the SecM reporter (IRAGP) which sequence is given as fasta file. One notebook describes the processing of the reads mapping to the genome, and the other the reads mapping to the plasmid. 

