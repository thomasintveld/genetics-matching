#!/bin/bash

#
# Grab genbank protein file from
# NCBI website.
#
# Thomas in't Veld 2016-03-26
#
#
#

species=$1
LOC=~/genetics-matching/aggregate-protein/

cd $LOC

echo 'grabbing' $species

rm -rf $species.gbk.gz
rm -rf $species.gbk
rm -rf protein.gbk.gz

wget -q "ftp://ftp.ncbi.nih.gov/genomes/$species/protein/protein.gbk.gz"
mv -f protein.gbk.gz $species.gbk.gz
gunzip $species.gbk.gz

echo 'done'
