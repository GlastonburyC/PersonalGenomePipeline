
args <- as.character(commandArgs(trailingOnly = TRUE))

SAMPLE_ID=as.character(args[1])
file1=args[2]
file2=args[3]

mat=read.table(file1,head=T)
pat=read.table(file2,head=T)

mat[,7]=mat[,7]+pat[,7]

colnames(mat)[7]=SAMPLE_ID

file_out=paste(SAMPLE_ID,".GeneCount.Final.txt",sep="")
write.table(mat,file_out,col.names=T,row.names=F,sep="\t",quote=F))
