#
# Blast Matching
# Thomas in't Veld 2016-03-24
#
# High performance script running over multiple cores matching all proteins from a reference genome
# (e.g. homo_sapiens) to a blast database of a subject genome (i.e. gallus_gallus).
# Result is a file with one row per protein for the reference genome and its matching protein in subject.
# Winner is decided with BLAST bitscore (with an EValue threshold configurable)
#
#


args<-commandArgs(TRUE)
library(seqinr)
library(snow)
source('genetics.helperfunctions.R')

#REFERENCE <- 'Homo_sapiens'
#SUBJECT <- 'Gallus_gallus'

REFERENCE <- args[1]
SUBJECT <- args[2]

TESTING <- FALSE

THRESHOLD <- 45
EVALUE <- 1E-5
EVALUETHRESHOLD <- 0.0001
BLASTPATH <- ''
BLASTDBPATH <- '~/genetics-matching/blast-db/'
SPECIESPATH <- '~/genetics-matching/species-fasta/'
OUTPUTPATH <- '~/genetics-matching/match-output/'


fastafile <- paste(SPECIESPATH, REFERENCE, '.fa', sep='')
logFile <- paste(OUTPUTPATH, 'r_', REFERENCE, '_s_', SUBJECT, '.txt', sep='')
outFile <- paste(OUTPUTPATH, 'r_', REFERENCE, '_s_', SUBJECT, '.RData', sep='')

referencePath <- paste(BLASTDBPATH, SUBJECT, '_ncbi', sep='')

speciesFasta <- loadFasta(fastafile)


cl <- makeCluster(8)

# initialise logfile
write('ref_name, ref_desc, match_name, match_desc, match_idnty, match_bitscore', file=logFile, append=FALSE)


getMatchedResult <- function(speciesFastaEntry){
  qname <-  attr(speciesFastaEntry[[1]], "name")
  qdesc <- attr(speciesFastaEntry[[1]], "Annot")
  qseq <- speciesFastaEntry[[1]][1]

  cat(paste('Blasting', qname, 'at', Sys.time()))

  blastresult <- blastThis(qseq)

  resulting <- c(name=qname, desc=gsub(',', '', qdesc), best_match=blastresult[[2]], best_match_desc=gsub(',', '', blastresult[[13]]), best_match_identity=blastresult[[3]], bitscore=blastresult[[12]])

  # write to disk
  write(paste(resulting, collapse=','), append=TRUE, file=logFile )

  return(resulting)

}

prodList <- speciesFasta

if(TESTING){
  prodList <- speciesFasta[sample(1:length(speciesFasta), 50)]
}

clusterExport(cl, c('prodList','getMatchedResult', 'blastThis', 'BLASTPATH', 'OUTPUTPATH', 'referencePath',  'BLASTDBPATH', 'EVALUE', 'deleteFile', 'log.INFO', 'logFile'))

resframe <- clusterApplyLB(cl, 1:length(prodList), function(i) return(getMatchedResult(prodList[i])))

save(resframe, file=outFile)



stopCluster(cl)
