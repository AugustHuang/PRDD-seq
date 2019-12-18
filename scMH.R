#!/bin/env Rscript
library(yyxMosaicHunter)

args<-commandArgs(TRUE)
verbose=TRUE
working_folder=args[1]

denovo_rate=10^-8
mosaic_rate=10^-2

MH=function(x) {
  ref=x[1]
  alt=x[2]
  math_mean=nchar(alt)/(nchar(ref)+nchar(alt))
  math_var=1.5/sqrt(nchar(ref)+nchar(alt))
  min_ci=math_mean-math_var
  max_ci=math_mean+math_var
  if (min_ci<0) {
    min_ci=0
  }
  if (max_ci>1) {
    max_ci=1
  }
  mh_result=yyx_wrapped_mosaic_hunter_for_one_site(ref,alt)
  cred_int=yyx_get_credible_interval(mh_result$likelihood_fun,c(min_ci,max_ci),0.95)
  all=c(mh_result$ref_het_alt_mosaic_posterior,cred_int$MLE,cred_int$CI)
  return (all)
}

sc_list=read.delim("1465_single-cells/1465_single-cells.list",header=T)[,1]
sc_dropout=read.delim("1465_single-cells/1465_single-cells.list",header=T)[,2]
sc_error=read.delim("1465_single-cells/1465_single-cells.list",header=T)[,3]

setwd(working_folder)

GMfile="1465_GM_BAQextraction.vcf"
GM=read.table(GMfile,header=T,fill=TRUE,sep='\t',quote='',comment.char='')
index=which((GM$ref_count+GM$alt_count)==0)
if (length(index)>0) {
	GM=GM[-index,]
}
BAQresult=t(apply(as.matrix(GM[,c(7,8)]),1,MH))
colnames(BAQresult)=c("Lref","Lhet","Lalt","Lmosaic","MLE","CIlow","CIhigh")
GM_MH=cbind(GM,BAQresult)
posterior=t(t(cbind(GM_MH$Lref,GM_MH$Lhet,GM_MH$Lalt,GM_MH$Lmosaic))*c(1-2*denovo_rate-denovo_rate^2-mosaic_rate,2*denovo_rate,denovo_rate^2,mosaic_rate))
GM_MH$Pref=posterior[,1]/rowSums(posterior)
GM_MH$Phet=posterior[,2]/rowSums(posterior)
GM_MH$Palt=posterior[,3]/rowSums(posterior)
GM_MH$Pmosaic=posterior[,4]/rowSums(posterior)
GM_MH$MAF=GM_MH$alt_count/(GM_MH$ref_count+GM_MH$alt_count)

template=cbind(GM_MH[,1:3],1:nrow(GM_MH))
colnames(template)[4]="index"

for (file in sc_list) {
	original=paste(file,"_BAQextraction.vcf",sep="")
	finalfile=paste(file,"_MHcall.txt",sep="")
	BAQ=read.table(original,header=T,fill=TRUE,sep='\t',quote='',comment.char='')
	index=which((BAQ$ref_count+BAQ$alt_count)==0)
	if (length(index)>0) {
		BAQ=BAQ[-index,]
	}
	BAQresult=t(apply(as.matrix(BAQ[,c(7,8)]),1,MH))
	colnames(BAQresult)=c("Lref","Lhet","Lalt","Lmosaic","MLE","CIlow","CIhigh")
	BAQfinal=cbind(BAQ,BAQresult)
	write.table(BAQfinal,file=finalfile,sep='\t',quote=F,col.names=T,row.names=F)
}

for (file in sc_list) {
	original=paste(file,"_MHcall.txt",sep="")
	finalfile=paste(file,"_MHcall_sort.txt",sep="")
	BAQ=read.table(original,header=T,fill=T,sep='\t',quote='',comment.char='')
	BAQtemp=merge(template,BAQ,all=T)
	BAQtemp=unique(BAQtemp)
	tempindex=which(is.na(BAQtemp$index))
	if (length(tempindex)) {
		BAQtemp=BAQtemp[-tempindex,]
	}
	BAQtemp=BAQtemp[order(BAQtemp$index),]
	write.table(BAQtemp,file=finalfile,sep='\t',quote=F,col.names=T,row.names=F)
}

posData=data.frame(GM_MH$chr,GM_MH$pos,GM_MH$ref,GM_MH$alt)
colnames(posData)=c("chr","pos","ref","alt")
bulkData=cbind(GM_MH$Pref,GM_MH$Phet,GM_MH$Palt,GM_MH$Pmosaic,GM_MH$MAF)
colnames(bulkData)=c("Bref","Bhet","Balt","Bmosaic","Bmaf")

for (file in sc_list) {
	original=paste(file,"_MHcall_sort.txt",sep="")
	after=paste(file,"_MHcall_PrepData.txt",sep="")
	SC=read.table(original,sep='\t',quote='',fill=T,comment.char='',header=T)
	SCData=cbind(SC$Lref,SC$Lhet,SC$Lalt,SC$Lmosaic)
	colnames(SCData)=c("Sref","Shet","Salt","Smosaic")
	allData=cbind(posData,bulkData,SCData)
	write.table(allData,file=after,sep='\t',quote=F,col.names=T,row.names=F)
}

##################################################################################################
# needed input:
# "chr" "pos" "ref" "alt" "Bref" "Bhet" "Balt"  "Bmosaic" "Bmaf" "Sref" "Shet" "Salt" "Smosaic"
##################################################################################################

PriorBulk = function (bulk) {
	pre_ref=as.numeric(bulk[1])
	pre_het=as.numeric(bulk[2])
	pre_alt=as.numeric(bulk[3])
	pre_mosaic=as.numeric(bulk[4])
	pre_sum=pre_ref+pre_het+pre_alt+pre_mosaic
	pre_ref=pre_ref/pre_sum
	pre_het=pre_het/pre_sum
	pre_alt=pre_alt/pre_sum
	pre_mosaic=pre_mosaic/pre_sum
	MAF=as.numeric(bulk[5])
	if (MAF<0.5 | MAF==0.5) {
		post_ref=pre_ref+(1-2*MAF)*pre_mosaic
		post_het=pre_het+2*MAF*pre_mosaic
		post_alt=pre_alt
	} else {
		post_ref=pre_ref
		post_het=pre_het+2*(1-MAF)*pre_mosaic
		post_alt=pre_alt+(2*MAF-1)*pre_mosaic
	}
	return (c(post_ref,post_het,post_alt))
}

ErrorDropoutAdjust=function(SC,error,dropout) {
	pre_ref=as.numeric(SC[1])
	pre_het=as.numeric(SC[2])+as.numeric(SC[4])
	pre_alt=as.numeric(SC[3])
	post_ref=(1-2*error)*pre_ref+dropout*pre_het
	post_het=2*error*(pre_ref+pre_alt)+(1-2*dropout)*pre_het
	post_alt=(1-2*error)*pre_alt+dropout*pre_het
	return (c(post_ref,post_het,post_alt))
}

Posterior=function(data,error,dropout) {
	bulk=data[5:9]
	SC=data[10:13]
	prior=PriorBulk(bulk)
	orilh=ErrorDropoutAdjust(SC,error,dropout)
	posterior=prior*orilh
	post_sum=posterior[1]+posterior[2]+posterior[3]
	posterior=posterior/post_sum
	names(posterior)=c("post_adj_ref","post_adj_het","post_adj_alt")
	return (posterior)
}

singlecell_post=function(file) {
	original=paste(file,"_MHcall_PrepData.txt",sep="")
	after=paste(file,"_MHcall_Posterior.txt",sep="")
	DATA=read.table(original,sep='\t',header=T)
	dropout=sc_dropout[which(sc_list==file)]
	error=sc_error[which(sc_list==file)]
	temp=t(apply(DATA,1,Posterior,error=error,dropout=dropout))
	temp=cbind(DATA,temp)
	write.table(temp,after,sep='\t',quote=F,col.names=T,row.names=F)
}

for (file in sc_list) {
  singlecell_post(file)
}

#####################

template=posData
for (file in sc_list) {
  original=paste(file,"_MHcall_Posterior.txt",sep="")
  DATA=read.table(original,sep='\t', header=T)
  template=cbind(template,DATA$post_adj_het)
}
colnames(template)[5:ncol(template)]=as.character(sc_list)
#hist(unlist(template[,5:ncol(template)]),breaks=100)
write.table(template,file="all.MHcall_Posterior_sum.txt",sep='\t',quote=F,col.names=T,row.names=F)

#####################

template=posData
x=template[,1:3]
y=template[,1:3]
for (file in sc_list) {
  original=paste(file,"_MHcall_sort.txt",sep="")
  SC=read.table(original,sep='\t',quote='',fill=T,comment.char='',header=T)
  temp=merge(SC,template[,1:3],all=T)
  temp=temp[order(temp$index),]
  x=cbind(x,temp$alt_count/(temp$ref_count+temp$alt_count))
  y=cbind(y,temp$alt_count+temp$ref_count)
}
x=cbind(template,x[,4:ncol(x)])
y=cbind(template,y[,4:ncol(y)])
colnames(x)[5:ncol(x)]=as.character(sc_list)
colnames(y)[5:ncol(y)]=as.character(sc_list)
#hist(unlist(x[,5:ncol(x)]),breaks=100)
#hist(unlist(y[,5:ncol(y)]),breaks=100)
write.table(x,file="all.MHcall_SCVAF_sum.txt",sep='\t',quote=F,col.names=T,row.names=F)
write.table(y,file="all.MHcall_depth_sum.txt",sep='\t',quote=F,col.names=T,row.names=F)

#########################################################

binary_0=function(cutoff,mindepth,minVAF,maxVAF,data,depth,vaf) {
	x=as.matrix(data[,5:ncol(data)])
	y=as.matrix(depth[,5:ncol(depth)])
	z=as.matrix(vaf[,5:ncol(vaf)])
	#x=x[,c(1,9,16,5,21,10,11,2,6,25,26,17,7,22,4,12,15,24,13,14,3,8,18,23)]
	#y=y[,c(1,9,16,5,21,10,11,2,6,25,26,17,7,22,4,12,15,24,13,14,3,8,18,23)]
	#z=z[,c(1,9,16,5,21,10,11,2,6,25,26,17,7,22,4,12,15,24,13,14,3,8,18,23)]
	x_cp=x
	x[which(x_cp<cutoff | x_cp==cutoff | y<mindepth | z<minVAF| z==minVAF | z>maxVAF | z==maxVAF)]=0
	x[which(x_cp>cutoff & y>(mindepth-1) & z>minVAF & z<maxVAF)]=1
	x[which(is.na(x_cp))]=0
	temp=apply(x,1,sum)
	x[which(is.na(x_cp))]=-1
	x=cbind(data[,1:4],x,temp)
	return (x)
}

binary=function(cutoff,mindepth,minVAF,maxVAF,data,depth,vaf) {
	x=as.matrix(data[,5:ncol(data)])
	y=as.matrix(depth[,5:ncol(depth)])
	z=as.matrix(vaf[,5:ncol(vaf)])
	#x=x[,c(1,9,16,5,19,21,10,11,2,6,25,26,17,7,22,4,12,15,24,13,20,14,3,8,18,23)]
	#y=y[,c(1,9,16,5,19,21,10,11,2,6,25,26,17,7,22,4,12,15,24,13,20,14,3,8,18,23)]
	#z=z[,c(1,9,16,5,19,21,10,11,2,6,25,26,17,7,22,4,12,15,24,13,20,14,3,8,18,23)]
	x_cp=x
	x[which(x_cp<cutoff | x_cp==cutoff | y<mindepth | z<minVAF| z==minVAF | z>maxVAF | z==maxVAF)]=0
	x[which(x_cp>cutoff & y>(mindepth-1) & z>minVAF & z<maxVAF)]=1
	x[which(is.na(x_cp))]=0
	temp=apply(x,1,sum)
	x[which(is.na(x_cp))]=-1
	x=cbind(data[,1:4],x,temp)
	return (x)
}

Phet_all=read.table("all.MHcall_Posterior_sum.txt",header=T,sep="\t",check.names=F)
DEP_all=read.table("all.MHcall_depth_sum.txt",header=T,sep="\t",check.names=F)
SCVAF_all=read.table("all.MHcall_SCVAF_sum.txt",header=T,sep="\t",check.names=F)

#cutoff=0.4
#mindepth=10
#minVAF=0.2
#maxVAF=1.01
#output_all=binary_0(cutoff,mindepth,minVAF,maxVAF,Phet_all,DEP_all,SCVAF_all)
#index_shared=which(output_all$temp>1)
#Phet_all=Phet_all[index_shared,]
#DEP_all=DEP_all[index_shared,]
#SCVAF_all=SCVAF_all[index_shared,]

cutoff=0.4
mindepth=1
minVAF=0.0
maxVAF=1.01
output_all=binary(cutoff,mindepth,minVAF,maxVAF,Phet_all,DEP_all,SCVAF_all)
write.table(output_all,file="all.MHcall_candidates.txt",sep='\t',quote=F,col.names=T,row.names=F)

type=matrix(NA,ncol=4,nrow=nrow(Phet_all))
for (i in 1:nrow(Phet_all)) {
	type1=length(which(as.numeric(Phet_all[i,5:(ncol(Phet_all))])<0.4))
	type3=length(which(as.numeric(Phet_all[i,5:(ncol(Phet_all))])>0.6))
	type2=length(which(as.numeric(Phet_all[i,5:(ncol(Phet_all))])<0.6 & as.numeric(Phet_all[i,5:(ncol(Phet_all))])>0.4))
	type[i,]=c(type1,type2,type3,ncol(Phet_all)-4-type1-type2-type3)
}
index=(which(type[,2]<5 & type[,4]<4 & type[,1]<(ncol(output_all)-5) & type[,1]>0.5*(ncol(output_all)-5)))
output_all_s1=output_all[index,]
SCVAF_all_s1=SCVAF_all[index,]
write.table(output_all_s1,file="all.MHcall_candidates_s1.txt",sep='\t',quote=F,col.names=T,row.names=F)

VAF_type=matrix(NA,ncol=4,nrow=nrow(SCVAF_all_s1))
for (i in 1:nrow(SCVAF_all_s1)) {
	type1=length(which(as.numeric(SCVAF_all_s1[i,5:ncol(SCVAF_all_s1)])<0.05))
	type2=length(which(as.numeric(SCVAF_all_s1[i,5:ncol(SCVAF_all_s1)])<0.3))-type1
	type4=length(which(is.na(SCVAF_all_s1[i,5:ncol(SCVAF_all_s1)])))
	type3=ncol(SCVAF_all_s1)-4-type1-type4-type2
	VAF_type[i,]=c(type1,type2,type3,type4)
}
fp=VAF_type[,2]/(VAF_type[,2]+VAF_type[,3])
VAF_type=cbind(VAF_type,fp)
#hist(fp,breaks=20)
#hist(unlist(SCVAF_all_s1[,5:20]),breaks=100)

index=which(fp<0.51 & VAF_type[,2]<(0.15*(ncol(SCVAF_all)-5)) & (VAF_type[,3]+VAF_type[,2])<(0.5*(ncol(SCVAF_all)-5)))
output_all_s2=output_all_s1[index,]
SCVAF_all_s2=SCVAF_all_s1[index,]
#hist(unlist(SCVAF_all_s2[,5:20]),breaks=100)
write.table(output_all_s2,file="all.MHcall_candidates_s2.txt",sep='\t',quote=F,col.names=T,row.names=F)
