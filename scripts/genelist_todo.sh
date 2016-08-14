#!/bin/bash
#
# # do big blast matching for all species in this list
#
# for genome in `cat ./genes_todo.txt`
# do
#   echo "**************************************"
#   echo "Starting analysis for" $genome "at" `date`
#   Rscript do_blast_match.R $genome Homo_sapiens
#   echo "DONE reverse analysis for" $genome "at" `date`
#   echo "**************************************"
#   echo " "
# done

# do GC aggregation for all genes in list

for genome in `cat ./genes_todo.txt`
do
  echo "**************************************"
  echo "Starting analysis for" $genome "at" `date`
  ./get_genbank_files_rna.sh $genome
  cd ~/genetics-matching/aggregate-gc
  ../scripts/process_genbank_rna.sh $genome
  cd ~/genetics-matching/scripts
  echo "DONE reverse analysis for" $genome "at" `date`
  echo "**************************************"
  echo " "
done
