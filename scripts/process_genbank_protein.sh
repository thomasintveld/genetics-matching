#!/bin/bash

#  process_genbank.sh
#
#
#  Created by Lieven Thorrez on 12/02/16.
#  Appended and extended by Thomas in't Veld on 2016-02-27
#   usage: ./process_genbank [species]

# make sure to exit (and shout) on error
set -e
#
# echo -n "Which species (file has format [species].gbk) ? > "
# read species
# echo "Species is $species"
species=$1

# make a temp directory and cd in there
mkdir -p tmp
mkdir -p result
cd tmp


# assume infoseq is installed systemwide, otherwise you should add its directory to $PATH
infoseq ../genbank/"$species".gbk -nousa -noname -nodatabase -notype -noorganism -nogi -nodescription -delimiter '\t' -outfile "$species"_infoseq.txt
egrep -o "gene=\"([A-Za-z0-9\_\+]+)" ../genbank/"$species".gbk > "$species"_gene.txt
cat "$species"_gene.txt | sed -e $'s/gene\=\"//g' > temp_line.txt

(echo "gene"; cat temp_line.txt) > temp_all.txt

# replace spaces by tabs for better format

paste "$species"_infoseq.txt temp_all.txt > temp_merge.txt
#cat temp_merge.txt | sed -E $'s/ +/\t/g' > temp_bla.txt


cat temp_merge.txt | sed -E $'s/ +/\t/g' > ../result/"$species"_protein_info.txt


cd ..
rm -rf tmp
