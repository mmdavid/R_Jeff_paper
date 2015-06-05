data_full = read.table("2012-07-12_SF_10-count_otu_table2_forR.txt",sep="\t",header=T,row.names=1,check.names=F)
data <- subset(data_full, select = -c(A23S2,A23S3,taxonomy))

matrix = t(data)
correlation<-cor(matrix,method="pearson")
head(correlation)
head(data_full)
