# Converting whatshap's haplotype output file to a more concise form

require(stringr,quietly=T)

# Main function
PhaseEx <-function(chromosome){
	phase=c()
	block=c()
	hap0=c()
	hap1=c()
	vcf=hapdata[hapdata$V1 ==chromosome,]
	vcf2=vcf[str_detect(vcf$V9,"PS"),]
	for (i in c(1:num)){
	
		v10_vec=strsplit(as.character(vcf2$V10[i]),":")
	
		phase=append(phase,v10_vec[[1]][1])
		block=append(block,v10_vec[[1]][5])
		if (v10_vec[[1]][1]=="1|0"){
			hap0=append(hap0,as.character(vcf2$V5[i]))
			hap1=append(hap1,as.character(vcf2$V4[i]))
		}else if (v10_vec[[1]][1]=="0|1"){
			hap0=append(hap0,as.character(vcf2$V4[i]))
			hap1=append(hap1,as.character(vcf2$V5[i]))
		}
	data=cbind(vcf2[,1:5],phase)
	data=cbind(data,block)
	data=cbind(data,hap0)
	data=cbind(data,hap1)
	}
	return(data)
}


# Arguments
args<-commandArgs(TRUE)
if(length(args)!=2){
  stop("Usage: Rscript Phase_Extract.R
       <1. Input haplotype block file from whatshap. >
       <2. chromosome Parameter. For one chromosome only. >
       ",call=F)
}

hapfile<-args[1]
chromosome<-args[2]
phasefile=paste(chromosome,".phase.txt",sep="")
hapdata <- read.table(hapfile,header=F)

phasedata <- PhaseEx(chromosome)
# Save
write.table(phasedata,file=phasefile,quote=F,col.names=T,row.names=F)
