#!/SCAR_env/bin/Rscript

# Sorts a chromosome list with a single information entry (a two column file)
# Intended input is a tab-delimited list of chromosomes and their sizes, the path to which is entered at line 10
# Output is a sorted list of the input info

library('data.table') 
library('gtools')

chrs <- read.delim('/opt/conda/process/genome/yeast_chr_size.txt', header = FALSE, sep = '\t')
names(chrs) <- c('chrom', 'length')
sorted_chrs <- mixedsort(chrs[['chrom']], numeric.type = 'roman', roman.case = 'upper')
order_chrs <- mixedorder(chrs[['chrom']], numeric.type = 'roman', roman.case = 'upper')

sort_chr <- chrs[match(sorted_chrs, chrs$chrom),]

write.table(sort_chr, 'genome/sacCer_chr_sorted.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)