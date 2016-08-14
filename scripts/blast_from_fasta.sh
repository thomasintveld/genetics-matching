#!/bin/bash
#
# Blast From Fasta
# Thomas in't Veld, 20160324
#
# Executes a blast matching (best bitscore wins) between a reference species (e.g. homo)
# and a subject species (e.g. gallus). Results are displayed as one line per homo protein.
# col1 : col2 is reference protein name and description.
# Best match in subject (according to bitscore) is added in cols 3 and 4 (name and desc)
# Identity scores (percentage and bitscore) are added in cols 5 and 6.
#
# Required:
# blastp (from blast+)
# run the script prepare_matching first for all species you want to analyse.
#

PREFIX='/Users/thomasintveld'


SPECIESFOLDER=$PREFIX'/genetics-matching/species-fasta/'
BLASTDBFOLDER=$PREFIX'/genetics-matching/blast-db/'

REFERENCE_SP=$1
SUBJECT_SP=$2

REFERENCE="$SPECIESFOLDER"$REFERENCE_SP".fa_seqs.csv"
SUBJECT="$BLASTDBFOLDER"$SUBJECT_SP"_ncbi"

# SPECIES='gallus_gallus'
# REFERENCE='~/blast-db/Homo_sapiens_ncbi.fa'

# EXECUTE FOR TEST RUN
# cat $REFERENCE | head -n 50 > $REFERENCE"_test"
# REFERENCE=$REFERENCE"_test"

# initialise outputfile
mkdir -p ~/genetics-matching/match-output
OUTPUTFILE="$PREFIX'/genetics-matching/match-output/r_"$REFERENCE_SP"__s_"$SUBJECT_SP"_out.txt"
echo "reference_name, reference_desc, subject_name, subject_desc, subject_identity, subject_bitscore" > $OUTPUTFILE

IFS=$'\n'
for j in `tail -n+2 $REFERENCE`
do

  name=`echo $j | awk -F"," '{print $1}'`
  name="${name%\"}"
  name="${name#\"}"


  desc=`echo $j | awk -F"," '{print $2}'`
  desc="${desc%\"}"
  desc="${desc#\"}"

  # remove comma's
  name=${name//,/}
  desc=${desc//,/}

  seq=`echo $j | egrep -o "[A-Z]{15,}"`

  topmatch=`echo $seq | blastp -db $SUBJECT -outfmt "10 std stitle" | head -n 1`
  # blast result
  # ('query_id', 'subject_id', 'identity_p', 'alignment_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score', 'desc'
  matchname=`echo $topmatch | awk -F"," '{print $2}'`
  matchident=`echo $topmatch | awk -F"," '{print $3}'`
  matchbitscore=`echo $topmatch | awk -F"," '{print $12}'`
  matchdesc=`echo $topmatch | awk -F"," '{print $13}'`

  echo "$name, $desc, $matchname, $matchdesc, $matchident, $matchbitscore" >> $OUTPUTFILE

  echo $name `date`

done

msg='DONE'
echo "================================"
echo $msg `date`
