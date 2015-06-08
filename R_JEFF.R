#Importing data
data_full = read.table("2012-07-12_SF_10-count_otu_table2_forR.txt",sep="\t",header=T,row.names=1,check.names=F)
data <- subset(data_full, select = -c(A23S2,A23S3,taxonomy))
matrix = t(data)

#filtering out the OTU not detected in at least 5 samples, quick and dirty
pb_5=c()
for (i in 1:(dim(matrix)[2])){
 a<-((length(subset(matrix[,i], matrix[,i] > 0)))>=5)
 pb_5<-c(pb_5,a)
 cat(i,"\t")
}

matrix_5samplesmin<-matrix[,pb_5]

#calculating the pvalue
correlationpvalue=c()
names_cor_pvalue=c()
for (i in 1:dim(matrix_5samplesmin)[2]){
cat(i,"\t")
  for (j in 1:dim(matrix_5samplesmin)[2]){
    a<-cor.test(matrix_5samplesmin[,i],matrix_5samplesmin[,j],method="pearson")$p.value
    correlationpvalue<-c(correlationpvalue,a)
    b<-paste((paste("OTU-",colnames(matrix_5samplesmin)[i], sep="")),(paste("OTU-",colnames(matrix_5samplesmin)[j], sep="")),sep="  ")
    names_cor_pvalue<-c(names_cor_pvalue,b)
  }
}

names(correlationpvalue)<-names_cor_pvalue

save(correlationpvalue, file="correlationpvalue.Rda")
save(names_cor_pvalue,file="names_cor_pvalue.Rda")

ok_correlationpvalue<-correlationpvalue
names(ok_correlationpvalue)<-names_cor_pvalue

#remove the one done twice and the one did with themselves
load("~/Documents/Jeff_paper/R_Jeff_paper/correlationpvalue.Rda")
load("~/Documents/Jeff_paper/R_Jeff_paper/names_cor_pvalue.Rda")

first_name<-lapply(names_cor_pvalue, function(x) strsplit(x,"  ")[[1]][1])
second_name<-lapply(names_cor_pvalue, function(x) strsplit(x,"  ")[[1]][2])

first_name_unlist<-unlist(first_name)
second_name_unlist<-unlist(second_name)
#ok the one with themselve are when first_name[i]==second_name[i]
#create a df


names_pb<-as.data.frame(first_name_unlist)
names_pb$second_name<-second_name_unlist

head(names_pb)
row.names(names_pb)<-names_cor_pvalue
name_pb_save<-names_pb


#with 2 millions values any loop is too slow, let's forget that: for i in 1:lenght() first_name[i]==second_name[i]
#I order them : good idea but complicated sortedpvalue <- with(names_pb,  names_pb[order(first_name_unlist, second_name) , ])
#how about I just rewrite the column inverting them !

head(names_pb)
names_pb$invert1<-names_pb$second_name
names_pb$invert2<-names_pb$first_name_unlist

#now merge with function


#thenI grep the legnth for each cycle: 
length(grep("OTU-0",sortedpvalue$first_name_unlist))
#1443  => ok so every 1443 line should be the same let's try
sortedpvalue[0, (1444 + 1), ]
sortedpvalue[(1444*2+ 2),]

newdata <- mtcars[order(mpg, cyl),]

#final correction
p.adjust(ok_correlationpvalue,method=c("fdr")
