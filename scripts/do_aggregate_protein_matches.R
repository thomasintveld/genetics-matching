args<-commandArgs(TRUE)
SPECIES <- args[1]

THRESHOLD <- 45
KEEPONLYDOUBLESIDEDMATCH <- TRUE

setwd('../aggregate-protein/')
#SPECIES <- 'Pseudopodoces_humilis'

LEFTFILE <- paste('../match-output/r_Homo_sapiens_s_', SPECIES, '.RData', sep='')
RIGHTFILE <- paste('../match-output/r_', SPECIES, '_s_Homo_sapiens.RData', sep='')


species_gc <- paste('../species-gc/', SPECIES, '_gc.txt', sep='')
homo_gc <- paste('../species-gc/', 'Homo_sapiens', '_gc.txt', sep='')

cat('Running for ')
cat(SPECIES)
cat(' ... Loading source files ...')

load(LEFTFILE) # to variable 'resframe' which is a list

# bind together into dataframe
resvec <- do.call('rbind', resframe)
left_df <- as.data.frame(resvec, stringsAsFactors = FALSE)


load(RIGHTFILE) # to variable 'resframe' which is a list

# bind together into dataframe
resvec <- do.call('rbind', resframe)
right_df <- as.data.frame(resvec, stringsAsFactors = FALSE)

cat(' loaded.\n')
cat('Calculating identities...')

left_df$best_match_identity <- as.numeric(left_df$best_match_identity)
right_df$best_match_identity <- as.numeric(right_df$best_match_identity)


# grab protein info for Homo and for Species
prot_homo <- read.table('result/Homo_sapiens_protein_info.txt', stringsAsFactors = FALSE, header=TRUE, sep="")
prot_spe <- read.table(paste('result/', SPECIES, '_protein_info.txt', sep=''), stringsAsFactors = FALSE, fill=TRUE, header=TRUE, sep='')

##### HELPER FUNCTIONS

testOverLap <- function(dfA, dfB){

  # reduce
  dfA <- unique(data.frame(gene.a=dfA$gene, match.a=dfA$matched_gene,stringsAsFactors=FALSE))
  dfB <- unique(data.frame(gene.b=dfB$gene, match.b=dfB$matched_gene,stringsAsFactors=FALSE))

  dfA_L <- merge(dfA, dfB, by.x='gene.a', by.y='match.b')
  dfB_R <- merge(dfB, dfA, by.x='gene.b', by.y='match.a')

  # count match between col 1 and col 3
  dfA_L['match'] <- (dfA_L$gene.a == dfA_L$gene.b)
  dfB_R['match'] <- (dfB_R$gene.b == dfB_R$gene.a)

  # compare
  table(dfA_L$match)
  table(dfB_R$match)

  # remove entries with FALSE backmatch
  dfA_L <- dfA_L[dfA_L$match,]

  return(dfA_L)

}


#convert to protein accession names
convertProteinNames <- function(xvec){

  matches <- regmatches(xvec, regexec("([A-Z]+)\\_([0-9]+)", xvec))
  matches_frame <- do.call(rbind, matches)
  matches_res <- matches_frame[,1]

  return(matches_res)
}

perGeneGetLongestTranscript <- function(dfx, gn){

  dfx <- dfx[dfx$gene==gn,]

  # throw away where identity < 45
  dfx_f <- dfx[dfx$best_match_identity>=THRESHOLD,]

  # check size, otherwise return 0
  if(nrow(dfx_f)<1){

    return(NULL)
  }

  maxL <- max(dfx_f$Length)

  dfx_f <- dfx_f[dfx_f$Length==maxL,]
  dfx_f <- dfx_f[1,]

  return(dfx_f)

}



#####

cat(' done.\n')
cat('Reduce protein level down to genes...')

left_df['Accession'] <- convertProteinNames(left_df$name)
right_df['Accession'] <- convertProteinNames(right_df$name)

left_df <- merge(left_df, prot_homo, by='Accession', all.x=TRUE)
right_df <- merge(right_df, prot_spe, by='Accession', all.x=TRUE)


left_df['match_accession'] <- convertProteinNames(left_df$best_match)
right_df['match_accession'] <- convertProteinNames(right_df$best_match)

cat('converting done.\n Getting longest transcripts...')
# reduce to longest transcripts
Lgenes <- unique(left_df$gene)
Ldf <- lapply(Lgenes, function(xn) return(perGeneGetLongestTranscript(left_df, xn)))
Ldf <- do.call('rbind', Ldf)

Rgenes <- unique(right_df$gene)
Rdf <- lapply(Rgenes, function(xn) return(perGeneGetLongestTranscript(right_df, xn)))
Rdf <- do.call('rbind', Rdf)

cat('done.\n')
cat('Doing aggregation steps now...')

Ldf['matched_gene'] <- sapply(Ldf$match_accession, function(mac) {
  return(prot_spe[prot_spe$Accession==mac,'gene'][[1]])
  })

Rdf['matched_gene'] <- sapply(Rdf$match_accession, function(mac) {
  return(prot_homo[prot_homo$Accession==mac,'gene'][[1]])
})

# test overlaps and remove matches with only one-sided match
if(KEEPONLYDOUBLESIDEDMATCH){
testFrame <- testOverLap(Ldf, Rdf)
Ldf <- Ldf[Ldf$gene %in% testFrame$gene.a,]
}

## GC step

sGC <- read.table(file=species_gc, sep='\t', header=TRUE, stringsAsFactors = FALSE)
hGC <- read.table(file=homo_gc, sep='\t', header=TRUE, stringsAsFactors = FALSE)


Ldf['homo_gc'] <- sapply(Ldf$gene, function(x){
  return(hGC[hGC$gene==x, 'X.GC'][1])
})

Ldf['matched_gc'] <- sapply(Ldf$matched_gene, function(x){
  return(sGC[sGC$gene==x, 'X.GC'][1])
})

Ldf['homo_gene_accession'] <- sapply(Ldf$gene, function(x){
  return(hGC[hGC$gene==x, 'Accession'][1])
})

Ldf['matched_gene_accession'] <- sapply(Ldf$matched_gene, function(x){
  return(sGC[sGC$gene==x, 'Accession'][1])
})

cat(' done. Saving.\n')



### FINAL FRAME

result <- data.frame(
  homo.gene=Ldf$gene,
  homo.accession=Ldf$homo_gene_accession,
  homo.gc=Ldf$homo_gc,
  matched.gene=Ldf$matched_gene,
  matched.accession=Ldf$matched_gene_accession,
  matched.gc=Ldf$matched_gc,
  matched.identity=Ldf$best_match_identity,
  stringsAsFactors = FALSE
)

result <- unique(result)


# put in context information
homofile <- '../aggregate-gc/Homo_sapiens_info.txt'
infoH <- read.csv(homofile, stringsAsFactors = FALSE, header=TRUE)

infoH <- unique(data.frame(homo.gene=infoH$Associated.Gene.Name, homo.chr=infoH$Chromosome.Name, homo.start=infoH$Gene.Start..bp., homo.end=infoH$Gene.End..bp., homo.gc=infoH$X..GC.content, stringsAsFactors = FALSE))

dfR <- merge(infoH, result, by.x='homo.gene', by.y='homo.gene', all.x=TRUE, all.y=FALSE)
# 
# result['homo.chromosome'] <- sapply(result$homo.gene, function(x) return( infoH[infoH$homo.gene==x,'homo.chr'][1] ))
# result['homo.start'] <- sapply(result$homo.gene, function(x) return( infoH[infoH$homo.gene==x,'homo.start'][1] ))
# result['homo.end'] <- sapply(result$homo.gene, function(x) return( infoH[infoH$homo.gene==x,'homo.end'][1] ))
# 
# 

dfR <- dfR[!is.na(dfR$homo.accession),]

dfR[dfR$homo.chr==1,'homo.chr'] <- '01'
dfR[dfR$homo.chr==2,'homo.chr'] <- '02'
dfR[dfR$homo.chr==3,'homo.chr'] <- '03'
dfR[dfR$homo.chr==4,'homo.chr'] <- '04'
dfR[dfR$homo.chr==5,'homo.chr'] <- '05'
dfR[dfR$homo.chr==6,'homo.chr'] <- '06'
dfR[dfR$homo.chr==7,'homo.chr'] <- '07'
dfR[dfR$homo.chr==8,'homo.chr'] <- '08'
dfR[dfR$homo.chr==9,'homo.chr'] <- '09'



dfR <- dfR[order(dfR$homo.chr, dfR$homo.start),]

# make sliding windows


dfR['homo.gc.w30'] <- c(rep(NA, 29), sapply(30:nrow(dfR), function(i) return( mean(dfR[i-29:i, 'homo.gc.y'] ) )  )  )
dfR['match.gc.w30'] <- c(rep(NA, 29), sapply(30:nrow(dfR), function(i) return( mean(dfR[i-29:i, 'matched.gc'] ) )  )  )

##
## TODO: buffer tussen chromosomen
## 

# bind to reference

ref_homo <- read.csv('../aggregate-gc/homo_reference.csv', stringsAsFactors = FALSE, header=TRUE)

## remove entries where no gene name available
# ref_homo <- ref_homo[!ref_homo$gene.symbol=='',]
ref_homo$homo_gc <- as.numeric(ref_homo$homo_gc)

species_ref <- data.frame(homo.gene=dfR$homo.gene, match.gene=dfR$matched.gene,match.gc=dfR$matched.gc, match.accession=dfR$matched.accession,match.identity=dfR$matched.identity , stringsAsFactors = FALSE)

ref_all <- merge(ref_homo, species_ref, by.x='gene.symbol', by.y='homo.gene',all.x=TRUE)

# order to homo order
ref_all <- ref_all[order(ref_all$order_homo),]

# sliding window
ref_all['homo.gc.w30'] <- c(rep(NA, 29), sapply(30:nrow(ref_all), function(i) return( mean(ref_all[i-29:i, 'homo_gc'] , na.rm = TRUE) )  )  )
ref_all['match.gc.w30'] <- c(rep(NA, 29), sapply(30:nrow(ref_all), function(i) return( mean(ref_all[i-29:i, 'match.gc'], na.rm=TRUE ) )  )  )

# make sure the buffers don't have a GC calculation
ref_all[ref_all$fransnummer=='buffer','homo.gc.w30'] <- NA
ref_all[ref_all$fransnummer=='buffer','match.gc.w30'] <- NA


library(ggplot2)
p <- ggplot(data=ref_all) + geom_line(aes(x=order_homo, y=homo.gc.w30)) + geom_line(aes(x=order_homo, y=match.gc.w30))
p

write.csv(ref_all, file=paste('../aggregate-protein/linked_files/Homo_sapiens_to_', SPECIES, '.csv', sep=''), row.names=FALSE)


cat('All done.')
