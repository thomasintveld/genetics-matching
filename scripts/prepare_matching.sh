#!/bin/bash
#
# Prepare BLAST matching
# Thomas in't Veld 2016-03-24
#
# Sets up workspace for matching process
# Downloads species files from NCBI (protein sequences, FASTA)
# Build blast database
#
#
#

mkdir -p ~/genetics-matching
cd ~/genetics-matching
mkdir -p ~/genetics-matching/blast-db
mkdir -p ~/genetics-matching/species-fasta

cd ~/genetics-matching/species-fasta

species=$1
#species='Monodelphis_domestica' # test

echo "starting for $species" `date`

rm -f $species.fa.gz
rm -rf $species.fa
rm -rf ~/genetics-matching/blast-db/$species"_ncbi"*

# grab species from NCBI

wget -q "ftp://ftp.ncbi.nih.gov/genomes/$species/protein/protein.fa.gz"
mv -f protein.fa.gz $species.fa.gz
gunzip $species.fa.gz

# create blast database from species fasta
makeblastdb -in $species.fa -dbtype 'prot' -out ~/genetics-matching/blast-db/$species"_ncbi"

echo "done for $species" `date`
