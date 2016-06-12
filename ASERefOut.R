# Takes the maternal and paternal allele counts that are now relative to the universal reference alleles and calculated REF/ALT count

args <- as.character(commandArgs(trailingOnly = TRUE))

file1=args[1]
file2=args[2]
ase_mat = read.table(file1,head=T,stringsAsFactors=F)
ase_pat = read.table(file2,head=T,stringsAsFactors=F)

file_out= args[3]

merged_ase=data.frame()

for(i in 1:dim(ase_mat)[1]){

	index=which(ase_pat$variantID == ase_mat$variantID[i])
	if (length(index) < 1){

		merged_ase=rbind(merged_ase,ase_mat[i,])

	} else{

		c_refCount      = ase_mat$refCount[i]+ase_pat$refCount[index]
		c_altCount      = ase_mat$altCount[i]+ase_pat$altCount[index]
		c_totalCount    = ase_mat$totalCount[i]+ase_pat$totalCount[index]
		c_lowMAPQDepth  = ase_mat$lowMAPQDepth[i]+ase_pat$lowMAPQDepth[index]
		c_lowBaseQDepth = ase_mat$lowBaseQDepth[i]+ase_pat$lowBaseQDepth[index]
		c_rawDepth      = ase_mat$rawDepth[i]+ase_pat$rawDepth[index]
		c_otherBases    = ase_mat$otherBases[i]+ase_pat$otherBases[index]
		c_improperPairs = ase_mat$improperPairs[i]+ase_pat$improperPairs[index]
		mod=data.frame(contig=ase_mat$contig[i],position=ase_mat$position[i],variantID=ase_mat$variantID[i],
			       refAllele=ase_mat$refAllele[i],altAllele=ase_mat$altAllele[i],refCount=c_refCount,
			       altCount=c_altCount,totalCount=c_totalCount,lowMAPQDepth=c_lowMAPQDepth,
			       lowBaseQDepth=c_lowBaseQDepth,rawDepth=c_rawDepth,otherBases=c_otherBases,
			       improperPairs=c_improperPairs)
		merged_ase      = rbind(merged_ase,mod)
	}

}

for(i in 1:dim(ase_pat)[1]){

	index=which(ase_mat$variantID == ase_pat$variantID[i])
	if (length(index) < 1){

		merged_ase=rbind(merged_ase,ase_pat[i,])
	}
}

write.table(merged_ase,file_out,col.names=T,row.names=F,sep="\t",quote=F)
