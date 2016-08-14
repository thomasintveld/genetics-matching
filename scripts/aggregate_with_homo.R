args<-commandArgs(TRUE)
SPECIES <- args[1]


# matched file location
matchfile <- paste('../aggregate-protein/linked_files/Homo_sapiens_to_', SPECIES, '.csv', sep='')
homofile <- '../aggregate-gc/Homo_sapiens_info.txt'


# load
df <- read.csv(matchfile, stringsAsFactors = FALSE, header=TRUE)
infoH <- read.csv(homofile, stringsAsFactors = FALSE, header=TRUE)

table(df$homo.gene %in% infoH$Associated.Gene.Name)


df['homo.chromosome'] <- sapply(df$homo.gene, function(x) return( infoH[infoH$Associated.Gene.Name==x,'Chromosome.Name'][1] ))
df['homo.start'] <- sapply(df$homo.gene, function(x) return( infoH[infoH$Associated.Gene.Name==x,'Gene.Start..bp.'][1] ))
df['homo.end'] <- sapply(df$homo.gene, function(x) return( infoH[infoH$Associated.Gene.Name==x,'Gene.End..bp.'][1] ))
