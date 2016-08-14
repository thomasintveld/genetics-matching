#!/bin/bash

#  process_genbank.sh
#
#
#  Created by Lieven Thorrez on 12/02/16.
#  Appended and extended by Thomas in't Veld on 2016-02-27
#  Appended and extended by Thomas in't Veld on 2016-05-16 to incorporate 'longest transcript length' aggregations
#   usage: ./process_genbank_rna [species]

# run in the aggregate-protein or aggregate-gc folder

# make sure to exit (and shout) on error
set -e
#
# echo -n "Which species (file has format [species].gbk) ? > "
# read species
# echo "Species is $species"
species=$1

# make a temp directory and cd in there
mkdir -p tmp
cd tmp

# assume infoseq is installed systemwide, otherwise you should add its directory to $PATH
infoseq ../"$species".gbk -nousa -nodatabase -noname -notype -noorganism -nodescription -outfile "$species"_infoseq.txt

# clean up infoseq, switch around the order of the columns to get accession || gc % || length
awk -F" " '{print $1}' "$species"_infoseq.txt > temp1_acc.txt
awk -F" " '{print $2}' "$species"_infoseq.txt > temp2_length.txt
awk -F" " '{print $3}' "$species"_infoseq.txt > temp3_gc.txt

# egrep "translation=\"([A-Z]+)" ../"$species".gbk > temp4_protein.txt

paste temp2_length.txt temp1_acc.txt temp3_gc.txt > "$species"_infoseq.txt

grep "VERSION" ../"$species".gbk > "$species"_accession.txt

python ../../scripts/gb2tab.py -f 'gene' ../"$species".gbk > "$species"_tabular.txt

cut -f4 "$species"_tabular.txt > "$species"_descr.txt

# grab the gene column from this description file
awk -F"/" '{print $2}' "$species"_descr.txt | awk -F"\"" '{print $2}' > "$species"_gene.txt.t

# make all gene symbols uppercase
tr a-z A-Z < "$species"_gene.txt.t >  "$species"_gene.txt

# add line with gene header to this file
{ echo "gene"; cat "$species"_gene.txt; } > "$species"_gene.txt.t
mv -f "$species"_gene.txt.t "$species"_gene.txt

# paste columns together in new file combining gc content with gene symbols and accession numbers
paste "$species"_gene.txt "$species"_infoseq.txt > "$species"_all.txt

# remove non-coding genes (XR and NR)
grep -v 'XR_' "$species"_all.txt | grep -v 'NR_' > "$species"_filtered.txt

# now first sort this filtered list (without sorting the header)
# sort either on transcript name (sort -bf -k1 -k2) or on transcript length (sort -bf -k1 -k4)
head -n1 "$species"_filtered.txt > "$species"_sorted.txt
tail -n+2 "$species"_filtered.txt | sort -bfdr -k1 -k2 >> "$species"_sorted.txt

# now per gene count occurence of transcript factors (after removing first line)
tail -n+2 "$species"_sorted.txt | cut -f1 | uniq -c | sed -e 's/^[ \t]*//' | cut -d' ' -f1 > "$species"_counts.txt
#       remove first line      | first col | count  | remove leading spaces from result | grab first column (the numbers) > output

# add header flag back again
{ echo "NB transcripts"; cat "$species"_counts.txt; } > "$species"_counts.txt.t
mv -f "$species"_counts.txt.t "$species"_counts.txt

# now sort based on gene symbol and accession number (case insensitive)
# first dump header so we don't sort that weirdly, then sort the rest of the file as we want it
# and only keep one entry per gene symbol (alphabetically first accession number)
head -n1 "$species"_sorted.txt > "$species"_almost.txt
tail -n+2 "$species"_sorted.txt | sort -fur -k1,1 >> "$species"_almost.txt

# paste together with counts and reorder back to normal
paste "$species"_almost.txt "$species"_counts.txt  | sort -fub -k1 > "$species"_final.txt



# move this final output to parent directory
mv "$species"_final.txt ../"$species"_gc.txt

cd ..
rm -rf tmp
