data_full = read.table("2012-07-12_SF_10-count_otu_table2_forR.txt",sep="\t",header=T,row.names=1,check.names=F)
data <- subset(data_full, select = -c(A23S2,A23S3,taxonomy))

correlationpvalue=c()
names_cor_pvalue=c()
matrix = t(data)
for (i in 1:2601){
cat(i,"\t")
  for (j in 1:2601){
    a<-cor.test(matrix[,i],matrix[,j],method="pearson")$p.value
    correlationpvalue<-c(correlationpvalue,a)
    b<-paste((paste("OTU-",colnames(matrix)[i], sep="")),(paste("OTU-",colnames(matrix)[j], sep="")),sep="  ")
    names_cor_pvalue<-c(names_cor_pvalue,b)
  }
}

    
  head(correlation)
head(data_full)
