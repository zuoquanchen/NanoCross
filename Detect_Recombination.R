#Detection of recombinant molecules 
#under development by Chen Zuoquan  (China Jxlab)

require(Rsamtools,quietly=T)
require(foreach,quietly=T)
require(doParallel,quietly=T)

#Judging odd numbers
jishu <- function(x){-
  ifelse(x%%2 ==0,F,T)
}
#Confirmation of major haplotypes
mainHap <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#Generating haplotype block information
blockgap=function (vec){
	block_left=c()
	block_right=c()
	for (b in vec){
	dat_tmp=hap[hap$block==b,]
	tmpnum=nrow(dat_tmp)
	b_left=dat_tmp$V2[1]
	b_right=dat_tmp$V2[tmpnum]
	block_left=append(block_left,b_left)
	block_right=append(block_right,b_right)
	}
	blockdata=data.frame(block_left,block_right)
	return(blockdata)
}
#Main program
FindRecReads<- function(i,seq_num=0,m=0,gapin="false"){
	rec=c()
	pos_strati=as.numeric(pos[i])
	seq_loc= pos_strati
	seq_vec=strsplit(seq[i],"")
	cigar_vec=strsplit(cigar[i],split ="(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)
	#Iterate through CIGAR Get location information
	for (ci in cigar_vec[[1]]){
		m=m+1
		if (jishu(m)){
			num=as.numeric(ci)
		}else {
			if (ci=="M"){
				seq_loc=seq_loc+num
				seq_num=seq_num+num
			}else if(ci=="I"){
				seq_loc=seq_loc
				seq_num=seq_num+num
			}else if(ci=="D"){
				seq_loc=seq_loc+num
			}else if(ci=="S"){
			seq_num=seq_num+num
			}else if(ci=="H"){}
		}	
	}
	#Filtering reference span short 
	pos_end=seq_loc-1
	pos_vec=hapPos[which(hapPos >pos_strati & hapPos < pos_end)]
	blockdata3=blockdata[blockdata$block_left<=pos_end & blockdata$block_right >=pos_strati,]

	if (nrow(blockdata3)==1){
		block_vec=c(blockdata3$block_left[1]:blockdata3$block_right[1])
		pec_vec=intersect(pos_vec,block_vec)
		npos=length(pec_vec)
		hap_part=hap[hap$V2 %in% pec_vec,]
	}else {
		gapin="true"
		npos=0
	}
	vcf_have="true"
	if (npos < 13){
	vcf_have="false"
	}
    #phase molecules 
	if (vcf_have=="true" & gapin=="false"){
		seq_num=0;m=0;reads_geno=c()
		seq_loc=as.numeric(pos[i])
		for (loc in pec_vec){
			while (seq_loc < loc) {
				m=m+1
				if (jishu(m)){
					ci=cigar_vec[[1]][m]
					num=as.numeric(ci)
				}else {
					ci=cigar_vec[[1]][m]
					if (ci=="M"){
					seq_loc=seq_loc+num
					seq_num=seq_num+num
					}else if(ci=="I"){
					seq_num=seq_num+num
					}else if(ci=="D"){
					seq_loc=seq_loc+num
					}else if(ci=="S"){ 
					seq_num=seq_num+num 
					}else if(ci=="H"){}
				}
			}
			if(ci=="M"){
				n=seq_num-(seq_loc-loc-1)
				geno=seq_vec[[1]][n]
			}else if (ci=="D"){
				geno="N"
			}else{
				geno="N"
			}
			reads_geno=append(reads_geno,geno)
		}
		hap_A=as.character(hap_part$hap0)
		hap_B=as.character(hap_part$hap1)
		hap_vec=c()
		for (p in 1:npos){
			if (as.character(reads_geno[p]) ==as.character(hap_A[p])){
				hap_p="A"
				hap_vec=append(hap_vec,hap_p)
			}
			else if (as.character(reads_geno[p]) ==as.character(hap_B[p])){
				hap_p="B"
				hap_vec=append(hap_vec,hap_p)
			}
			else{
				hap_p="NA"
				hap_vec=append(hap_vec,hap_p)
			}
		}
		#print(hap_vec)
		hap_f=hap_vec[1:7]
		hap_f=hap_f[hap_f !="NA"]
		hap_f=append(hap_f,"Hap")
		hap_e=hap_vec[(npos-6):npos]
		hap_e=hap_e[hap_e !="NA"]
		hap_e=append(hap_e,"Hap")
		#Using the sliding window method to determine if reads are recombinant molecules
		if (mainHap(hap_f)!=mainHap(hap_e)){
			phase_vec=c()
			phase_pos=c()
			for (n in 4:(npos-3)){
			
			win_vec=hap_vec[(n-3):(n+3)]
			win_vec=win_vec[win_vec != "NA"]
			win_data=as.data.frame(table(win_vec))
			rownames(win_data)=win_data$win_vec
			if (nrow(win_data)==1){
				if (win_data$win_vec[1]=="A"){
					Anum=as.numeric(win_data$Freq[1])
					Bnum=0
				}else {
					Bnum=as.numeric(win_data$Freq[1])
					Anum=0
				}
			}else if (nrow(win_data)==2){
				Anum=as.numeric(win_data$Freq[1])
				Bnum=as.numeric(win_data$Freq[2])
			}else {
				Anum=0
				Bnum=0
			}
			
			if (length(win_vec)==0){
				phase="NA"
			}
			else{
				if ((Anum/length(win_vec))>=(5/6)){
					phase="A"
				}else if ((Bnum/length(win_vec))>=(5/6)){
					phase="B"
				}else {
					phase="NA"
				}
			}
			phase_pos=append(phase_pos,pos_vec[n])
			phase_vec=append(phase_vec,phase)
	        }
	        phase_data=data.frame(phase_vec,phase_pos)
			phase_data=phase_data[!phase_data$phase_vec %in% c("NA") ,]
			phase_num=nrow(phase_data)
			if (phase_num >=3 & mainHap(hap_f)==as.character(phase_data$phase_vec[1])){
				for(n in 2:phase_num){
					a=phase_data$phase_vec[n-1]
					b=phase_data$phase_vec[n]
					if (phase_data$phase_vec[1]!= phase_data$phase_vec[phase_num] ){
						if (a!=b){
							rec=c(i,qname[i],as.character(phase_data$phase_pos[n-1]),as.character(phase_data$phase_pos[n]))
							#print(rec)
						}
					} 
				}
			}
		}
	}
	return(rec)
}
args<-commandArgs(TRUE)
if(length(args)!=5){
  stop("Usage: Rscript Detect_Recombination.R
       <1. Input Bam files (Same as used for call variation files)  >
       <2. Input phase block file (output file of whatsphase) >
       <3. Number of threads >
	   <4. Chromsome (For one chromosome only) >
       <5. output filestem (output with columns Index, Reads, Pos_left, pos_right. For one chromosome only.) > ",call=F)   
}
bamFile <- args[1]
phaseFile <- args[2]
threads <- as.numeric(args[3])
chromosome <- args[4]
outFile <- args[5]

#Read data & Data preparation
data<-scanBam(bamFile)
hap=read.table(phaseFile,header=T)
hapPos=as.numeric(hap$V2)
# Do each chromosome
bam=data[[1]][data[[1]]$rname == chromosome]
cigar=as.character(bam$cigar)
seq=as.character(bam$seq)
qname=as.character(bam$qname)
pos=as.numeric(bam$pos)
blockdata=blockgap(unique(hap$block))
#run FindRecReads function parallel
cl <- makeCluster(threads) #not to overload your computer
registerDoParallel(cl)
recdata <- foreach(i=1:length(qname), .combine=rbind) %dopar% {
	recom=FindRecReads(i)
	recom
}
stopCluster(cl) #stop cluster
# Save
write.table(recdata,file=outFile,quote=F,col.names=F,row.names=F)
