library(seqinr)
# makeblastdb -in Homo_sapiens_ncbi.fa -parse_seqids -dbtype prot
# http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/

blastFind <- function(ensid, fasta, referenceSpecies){

	code <- fasta[[ensid]][[1]]

	log.INFO(paste('Blasting', ensid, 'with ref', referenceSpecies, 'at', Sys.time()))

	blastResult <- blastThis(code, referenceSpecies)

	
		result <- interpretBlast(blastResult)
		
		result <- cbind('ENSID_REF'=ensid, result)

	if(WRITEINTERMEDIATETOFILE){

		if(result[1,2]=='MISSING'){
			writeOut(ensid, result[1,2], 'TRUE')
		}else{
			writeOut(ensid, result[1,2], 'FALSE')
		}

	}


	return(result)

}
	
blastThis <- function(seq){
  
	
	# blast against reference species
	# blastformula <- paste("blastp -query ",codeFile," -db ", 
	#	BLASTDBPATH, referenceSpecies, " -out ", outFile, " -outfmt 6 -evalue ", EVALUE, sep="")
	# blastformula <- paste(BLASTPATH, "blastp -query ",codeFile," -db ", 
	#	BLASTDBPATH, referenceSpecies, " -out ", outFile, ' -outfmt "6 std stitle" -evalue ', EVALUE, sep="")


	blastformula <- paste('echo "', seq,'" | blastp -db ', referencePath, ' -outfmt "10 std stitle" -evalue ', EVALUE , ' | head -n 1 ', sep='')
	
	
	
	res <- system(blastformula, intern=TRUE)
	# # Fields: query id, subject id, % identity, alignment length, 
	#   mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, desc
	
	# check
	if(length(res)==0 || is.na(res) || is.null(res)){
	  return(NULL)
	}
	
  res <- strsplit(res, ',')
  res <- res[[1]]
	  
	blastResult <- res
  
	names(blastResult) <- c('query_id', 'subject_id', 'identity_p', 'alignment_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score', 'desc')
	
	
  
	return(blastResult)


}

deleteFile <- function(file){

	#log.INFO(paste('Deleting file', file))
	try(system(paste("rm", file)), silent=TRUE)

}

interpretBlast <- function(blastResult){

	if(is.null(dim(blastResult))){
		return(MISSINGSTRING)
	}

	blastResult <- blastResult[(blastResult$V3 >= THRESHOLD) & (blastResult$V11 <= EVALUETHRESHOLD),]

	if(dim(blastResult)[[1]] >= 1){

		result <- data.frame('ENSID_FOUND'=blastResult$V2,
							'IDENTITY'=blastResult$V3,
							'EVALUE'=blastResult$V11)

		return(result)
	}else{
		return(MISSINGSTRING)
	}	

}

writeOut <- function(ensid, ensid_ref, result){

	write(paste(ensid, ensid_ref, result, sep=','), file=outputFile, append=TRUE)

}



log.INFO <- function(message){
	write(message, logFile, append=TRUE)

	print(message)

}

loadFasta <- function(file){

	df <- read.fasta(file, as.string=TRUE, forceDNAtolower=FALSE)

	return(df)

}