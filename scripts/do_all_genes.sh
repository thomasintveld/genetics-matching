#!/bin/bash

#
#
#  do_all_genes.sh
#  written by Thomas in't Veld, 2016-02-29
#  idea is to loop over all NCBI genes and calculate GC contents
#


# loop over gene list
while read species; do

  echo ""
  echo "STARTING "$species

  # make sure directory is clean
  rm -rf tmp
  rm -rf rna.gbk.gz

  # get genebank file from NCBI server, rename and unzip
  wget -q "ftp://ftp.ncbi.nih.gov/genomes/$species/RNA/rna.gbk.gz"
  # for RNA fasta:
  # wget -q "ftp://ftp.ncbi.nih.gov/genomes/$species/RNA/rna.fa.gz"
  # for Protein fasta:
  # wget -q "ftp://ftp.ncbi.nih.gov/genomes/$species/protein/protein.fa.gz"
  mv -f rna.gbk.gz $species.gbk.gz
  gunzip $species.gbk.gz

  # run the process_genbank script on top
  ./process_genbank.sh $species

  cp -f "$species"_gc.txt ~/Dropbox/Landscapes\ of\ CG\ content\ of\ vertebrate\ genomes\ 2016/species/lijst/
  rm -f $species.gbk

  echo "finished "$species
  echo ""
done < genelist.txt
