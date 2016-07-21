# Identify MZs that have RNA-seq but need their MZs UK10K sequence data duplicated

T='skin'

# All EUROBATs tissue IDs
tissue.ids=read.table('/Users/Craig/Dropbox/skin.ids',head=T)

# IDs with RNA-seq AND UK10K data
tissue.uk10k=read.table('/Users/Craig/Dropbox/skin_UK10K.ids',head=T)

IDs=NULL

for(i in 1:nrow(tissue.uk10k)){
	
	if(tissue.uk10k$Zygosity[i] =="MZ"){

		index = which(tissue.uk10k$Study_Number[i] == tissue.ids$Study_Number)
		if(index > 0){

		fam = which(tissue.uk10k$Family[i] == tissue.ids$Family)

		if(length(fam) > 1){ 

		if(substring(tissue.uk10k$Study_Number[i],nchar(tissue.uk10k$Study_Number[i])) == 1){
		IDs=append(IDs,paste(substring(tissue.uk10k$Study_Number[i],1,nchar(tissue.uk10k$Study_Number[i])-1),"2",sep="")) } else{
		IDs=append(IDs,paste(substring(tissue.uk10k$Study_Number[i],1,nchar(tissue.uk10k$Study_Number[i])-1),"1",sep=""))
		}

		}


	}
}}

# list of MZs that have RNA-seq but need to be put in UK10K data
IDs

all_IDs=append(IDs,tissue.uk10k[,2])
all_IDs=all_IDs[mixedorder(all_IDs)]

file_out=paste("~/Dropbox/all_",T,".ids",sep="")
write.table(all_IDs,file_out,col.names=F,row.names=F,sep="\t",quote=F)
