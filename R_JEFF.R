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
save(matrix_5samplesmin, file="matrix_5samplesmin.Rda")
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

#now merge 
inverted_names<-paste (names_pb$invert1,names_pb$invert2, sep="  ")

names_pb$inverted_names<-inverted_names
save(names_pb, file="names_pb.Rda")

#now: remove duplicates with same column names
load("names_pb.Rda")

#I'm having problem with level and factors 
names_pb$first_name_unlist<-as.character(names_pb$first_name_unlist)
names_pb$invert2<-as.character(names_pb$invert2)

#ok I tried without looping but 2 millions = too big for a function. Remove the one done with themselves
samename<-names_pb$first_name_unlist == names_pb$second_name
names_pb_nothem<-names_pb[!samename,]

#other optnio but too slow
#themselves=c()
#for (i in 1:dim(names_pb)[1]){
#  a<-(names_pb$first_name_unlist[i] == names_pb$second_name[i])
#  themselves<-c(themselves,a)
#  cat(i,"\t")
#}
#save(themselves,file="themselves.Rda")
#names_pb_nothem<-names_pb[!themselves,]
save(names_pb_nothem,file="names_pb_nothem.Rda")

#next step: avoid stupid loop and just look at the one with name first < second
#let's try to add a column

names_pb_nothem$normalnames<-row.names(names_pb_nothem)

#other strategy: secdon number is alway bigger than the first one
head(names_pb_nothem)
first_name_num<-lapply(names_pb_nothem$first_name_unlist, function(x) strsplit(x,"-")[[1]][2])
second_name_num<-lapply(names_pb_nothem$second_name, function(x) strsplit(x,"-")[[1]][2])
first_name_num<-unlist(first_name_num)
second_name_num<-unlist(second_name_num)
names_pb_nothem$first_name_num<-first_name_num
names_pb_nothem$second_name_num<-second_name_num

oktoremove<-(first_name_num <  second_name_num)
okfinalJEffnames<-names_pb_nothem[oktoremove,]


save(okfinalJEffnames,file="okfinalJEffnames.Rda")
  
#final correction
Jeff_adjusted<-p.adjust(okfinalJEffnames$pvalue,method=c("fdr")
names(Jeff_adjusted)<-okfinalJEffnames$normalnames                       
save(Jeff_adjusted,file="Jeff_adjusted.Rda")            

passed_jeff<-(Jeff_adjusted<0.05)
passedok<-Jeff_adjusted[passed_jeff]

#let's compare with Jeff's file
used<-read.table("edges.txt")
used$V1<-as.character(used$V1)
used$V2<-as.character(used$V2)

names_used_new<-paste(used$V1,used$V2,sep="  ")
pb_check_order<-setdiff(names_used_new,names(passedok))
#ok some of them may have the same ont he other oder
jeff_network_firstpart<-sapply(pb_check_order, function(x) strsplit(x,"  ")[[1]][1])
jeff_network_secondpart<-sapply(pb_check_order, function(x) strsplit(x,"  ")[[1]][2])
jeff_names_inv<-paste(jeff_network_secondpart,jeff_network_firstpart,sep="  ")
pb_network<-setdiff(jeff_names_inv,names(passedok))

#NONE = everypone passed
#and with the old file 
used<-read.table("old_edge_list.txt")
used$V1<-as.character(used$V1)
used$V2<-as.character(used$V2)

names_used_new<-paste(used$V1,used$V2,sep="  ")
pb_check_order<-setdiff(names_used_new,names(passedok))
#ok some of them may have the same ont he other oder
jeff_network_firstpart<-sapply(pb_check_order, function(x) strsplit(x,"  ")[[1]][1])
jeff_network_secondpart<-sapply(pb_check_order, function(x) strsplit(x,"  ")[[1]][2])
jeff_names_inv<-paste(jeff_network_secondpart,jeff_network_firstpart,sep="  ")
pb_network<-setdiff(jeff_names_inv,names(passedok))

write.csv(pb_network,file="pb_network.oldfile.csv")
