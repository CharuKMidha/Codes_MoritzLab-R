#!/usr/bin/env Rscript

message()
message("MetaMSD (Meta Analysis for Mass Spectrometry Data)")
message()

#	Using BiocManager to install any missing required packages.

list.of.packages <- c("ROCR","stats","gridExtra","gtable","optparse","qvalue")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
	if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install(new.packages)
}


#	http://tuxette.nathalievilla.org/?p=1696
library("optparse")

option_list = list(
	make_option(c("-m","--metaanalysis"), type="character", default="Stouffer",
		help="Specify a meta-analysis test. Either Stouffer or Pearson. [ default = %default ]", metavar="character"),
	make_option(c("-c","--cutoff"), type="double", default=0.05,
		help="Specify a q-value cut off. [ default = %default ]", metavar="number"),
	make_option(c("-t","--top"), type="integer", default=15,
		help="Specify the number of proteins in Top-N differential protein list. [ default = %default ]", metavar="number"),
	make_option(c("-i", "--input"), type="character", default="input",
		help="Specify the input folder name. [ default = %default ]", metavar="character"),
	make_option(c("-o", "--output"), type="character", default="output",
		help="Specify the output folder name. [ default = %default ]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if ( ! opt$metaanalysis %in% c("Stouffer","Pearson") ){
	stop("MetaAnalysis must be either Stouffer or Pearson\n", call.=FALSE)
}
MetaAnalysis   = opt$metaanalysis
q.value.cutoff = opt$cutoff
Top.n          = opt$top
input          = opt$input
output         = opt$output
dir.create( output, showWarnings = FALSE, recursive = TRUE )

require(ROCR)
require(stats)
require(qvalue)
require(gridExtra)
require(gtable)


######################################################################


################
#Misc Functions#
################

left <- function(p, stat){left.p = ifelse(stat<0, p/2, 1-p/2); return(unlist(left.p))}

#############
#Fisher.test#
#############
na.perc <- function(x){mean(is.na(x))}
Fisher.test <- function(DEprotein, union=TRUE){
	if (union==FALSE){
		index=which(apply(DEprotein$p, 1, na.perc)>0)
		DEprotein$p[index,]=NA
	}
	return(apply(DEprotein$p, 1, fisher.test))
}

fisher.test <- function(p){
	if(length(which(p==0))>0){
		p[p==0]=1e-20
	}
	if (sum(!is.na(p))>1){
		if (sum(p, na.rm=TRUE)==0){
			fisher.p=1e-20
		}else{
			p=p[!is.na(p)]
			fisher.p=pchisq(-2*sum(log(p)), 2*length(p), lower.tail=FALSE)
		}
	}else{
		if (sum(!is.na(p))==0){fisher.p=NA}
		if (sum(!is.na(p))==1){fisher.p=p[!is.na(p)]}
	}	
	return(fisher.p)
}

##############
#Pearson.test#
##############
Pearson <- function(DEprotein, union=TRUE){
	if (union==FALSE){
		index=which(apply(DEprotein$p, 1, na.perc)>0)
		DEprotein$p[index,]=NA
	}
	p = array(c(NA), dim=nrow(DEprotein$p))
	side = array(c(NA), dim=nrow(p))
	
	for (i in 1:nrow(DEprotein$p)){
		tmp=pearson(DEprotein$p[i,], DEprotein$stat[i,])
		p[i]=tmp$p
		side[i]=tmp$side
	}
	return(list(p=p, side=side))
}

pearson <- function(p.list, stat.list){
	left.p = left(p.list, stat.list)
	right.p = 1-left.p
	fisher.left=fisher.test(left.p)	
	fisher.right=fisher.test(right.p)
	pearson = 2 * min(fisher.left, fisher.right)
	pearson = min(c(1, pearson))
	pearson.side = ifelse(fisher.left < fisher.right, -1,1)
	return(list(p=pearson, side=pearson.side))
}

##########
#Stouffer#
##########
Stouffer <- function(DEprotein, union=TRUE, wt=NULL){
	if (union==FALSE){
		index=which(apply(DEprotein$p, 1, na.perc)>0)
		DEprotein$p[index,]=NA
	}
	p = array(c(NA), dim=nrow(DEprotein$p))
	side = array(c(NA), dim=nrow(p))
	if (length(wt)==0){wt=rep(1,ncol(DEprotein$p))}	
	for (i in 1:nrow(DEprotein$p)){
		tmp=stouffer.one.side(DEprotein$p[i,], DEprotein$stat[i,], wt)
		p[i]=tmp$p
		side[i]=tmp$side
	}
	return(list(p=p, side=side))
}

stouffer.one.side <- function(p.list, stat.list, wt){
	left.p = left(p.list, stat.list)
	right.p = 1-left.p
	stouffer.left=stouffer.left.side(left.p, wt)	
	stouffer.right=stouffer.right.side(right.p, wt)
	stouffer = 2 * min(stouffer.left, stouffer.right)
	stouffer = min(c(1, stouffer))
	stouffer.side = ifelse(stouffer.left < stouffer.right, -1,1)
	return(list(p=stouffer, side=stouffer.side))
}

stouffer.left.side <- function(x, wt){
	if (length(!is.na(x))>0){
		wt=wt[!is.na(x)]; x=x[!is.na(x)]
		z=qnorm(x, lower.tail=FALSE)
		#com.z=sum(z)/sqrt(length(z))
		com.z=sum(wt*z)/sqrt(sum(wt^2))
		p=pnorm(com.z, lower.tail=FALSE)
	}else{
		p=NA
	}
	return(p)
}

stouffer.right.side <- function(x, wt){
	if (length(!is.na(x))>0){
		wt=wt[!is.na(x)]; x=x[!is.na(x)]
		z=qnorm(x, lower.tail=TRUE)
		com.z=sum(wt*z)/sqrt(sum(wt^2))
		p=pnorm(com.z, lower.tail=TRUE)
	}else{
		p=NA
	}
	return(p)
}

######################################
#no. of detected proteins at q-values#
######################################
q.no <- function(z){
	a <- c(0,sort(unique(z), decreasing=FALSE))
	b=NULL
	for (i in 1:length(a)){
		b[i]=sum(z <= a[i], na.rm=TRUE)
	}
	return(list(x=c(0,a), y=c(0,b)))
}


#############################################
#measure performance of integrative approach#
#############################################
IDR.fun <- function(q.int, q.ind, cut){
	num1 <- which(q.int<cut)
	num2 <- union(which(q.ind[,1]>cut), which(is.na(q.ind[,1])))
	for (i in 2:ncol(q.ind)){
		tmp = union(which(q.ind[,i]>cut), which(is.na(q.ind[,i])))
		num2 <- intersect(num2, tmp)
	}
	num=length(intersect(num1, num2))
	den=length(num1)
	return(num/den)
}

IRR.fun <- function(q.int, q.ind, cut){
	num1 <- union(which(q.int > cut), which(is.na(q.int)))
	num2 <- which(q.ind[,1]<cut)
	for (i in 2:ncol(q.ind)){
		tmp = which(q.ind[,2]<cut)
		num2 <- union(num2, tmp)
	}
	num=length(intersect(num1, num2))
	den=length(num2)
	return(num/den)
}

######################################################################

input_files = list.files(input, pattern=".txt$",full.names=TRUE)
if( length(input_files) == 0 )
	stop("No files found in input dir\n", call.=FALSE)

message()
message(paste0("Found ",length(input_files)," in the ",input," directory. Processing..."))
message()

i = 1
for( f in input_files ) {
	message(f)

	this_data = read.csv(file=f, header=TRUE, sep="")
	#	Set each data set to have unique names for Sign and Pvalue by adding numeric suffix (Sign.1, Pvalue.1)
	names(this_data) = c("Protein",paste("Sign",i,sep="."),paste("Pvalue",i,sep="."))

	if( exists( "userdata" ) ){
		message("Merging data frame ...")
		userdata = merge( userdata, this_data, by="Protein", all=TRUE )
	} else {
		message("Creating initial data frame ...")
		userdata = this_data
	}
	rm(this_data)

	i = i + 1
}
rm(i)

message()
message("Merged data.frame contains the following columns.")
message(paste(colnames(userdata), collapse=", "))
message()

N.dataset = length( input_files )


#	userdata should have columns ...
#	Protein, Sign.1, Pvalue.1, Sign.2, Pvalue.2, ...
#	Use grepl to select those matching.
DEprotein=list(Protein=userdata$Protein,
	p=cbind(userdata[names(userdata)[grepl("Pvalue", names(userdata))]]),
	stat=cbind(userdata[names(userdata)[grepl("Sign", names(userdata))]]))
str(DEprotein)
#calculate q-value for each dataset
ind.q = NULL
for (i in 1:N.dataset){
	ind.q[[i]]=qvalue(DEprotein$p[,i],pi0 = 1)$qvalues
}

#Meta Analysis
if (MetaAnalysis=="Stouffer"){
	result=Stouffer(DEprotein)
} else if (MetaAnalysis=="Pearson"){
	result=Pearson(DEprotein)
} else {
	#if ((MetaAnalysis!="Stouffer") & (MetaAnalysis!="Pearson")){
	#	should never have gotten here
	cat("Choose between Stouffer and Pearson.")
	result=NULL
}

if (length(result)>0){
	q.meta=qvalue(result$p)$qvalues
	meta.result = data.frame(
		Protein     = DEprotein$Protein )

#	separated simply for output comparison
#	could put these other 3 columns here, but would change the output order
#		Protein     = DEprotein$Protein,
#		Sign.meta   = result$side,
#		pvalue.meta = result$p,
#		qvalue.meta = qvalue(result$p)$qvalues)

	#	these 3 columns are from each dataset
	#	still not sure if/when I need to add the [,]
	for (i in 1:N.dataset){
		meta.result[paste0("Sign.set",i)]   = ifelse( userdata[paste("Sign",i,sep=".")] > 0, 1, -1)
		meta.result[paste0("pvalue.set",i)] = userdata[paste("Pvalue",i,sep=".")]
		meta.result[paste0("qvalue.set",i)] = ind.q[[i]]
	}

	meta.result$Sign.meta   = result$side
	meta.result$pvalue.meta = result$p
	meta.result$qvalue.meta = qvalue(result$p)$qvalues

}

######################################################################

# 1. Output meta analysis results
write.table(meta.result, file=paste(output,"MetaAnalysisResult.txt", sep="/"),
	col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

######################################################################

# 2. a qplot (q-value threshold vs. # of detected proteins)
meta.no <- q.no(meta.result$qvalue.meta)

set.no = c()	#	list()
for (i in 1:N.dataset){
	set.no[[i]] <- q.no(meta.result[paste0("qvalue.set",i)][,])
}

sums = c()
sums[1] = sum(meta.result$qvalue.meta<q.value.cutoff, na.rm=TRUE)
for (i in 1:N.dataset){
	sums[length(sums)+1] = sum( meta.result[paste0("qvalue.set",i)] < q.value.cutoff, na.rm=TRUE )
}
max.value = max( sums )
	
max.value = round(max.value+5)
pdf(paste(output,"qplot.pdf", sep="/"))
par(mar=c(5,5,3,3))
plot(meta.no$x, meta.no$y, type="l", lwd=2,
	xlim=c(0,q.value.cutoff), ylim=c(0, max.value),
	xlab="q-value threshold", ylab="# of differential proteins",
	cex.lab=1.5, cex.axis=1.5, lty=1, col=1)

for (i in 1:N.dataset){
	points(set.no[[i]]$x, set.no[[i]]$y, type="l", lwd=2, lty=1, col=i+1)
}

legend_labels = c("Meta Analysis")
for (i in 1:N.dataset){
	legend_labels[i+1] = paste0("Single Analysis (Set ",i,")")
}
legend("topleft", legend_labels, lty=c(1),col=c(1,2:(N.dataset+1)), lwd=c(2), cex=1.3)

#	This prints "null device\n 1". Seems to work when removed.
#dev.off()

######################################################################

# 3. Top differential proteins detected by Meta Analysis
meta.result.top = meta.result[order(meta.result$qvalue.meta),][1:Top.n,]
#sorted by q-value and take the Top.n protein with the smallest q-values

Sign = ifelse(meta.result.top$Sign.meta==-1, "-", "+")

top.result = data.frame(  Rank = c(1:Top.n),
	Protein=meta.result.top$Protein,
	Sign = ifelse(meta.result.top$Sign.meta==-1, "-", "+"),
	Pvalue = meta.result.top$pvalue.meta,
	Qvalue = meta.result.top$qvalue.meta)
rownames(top.result)=NULL
pdf(paste(output,"TopDifferentialProteins.pdf",sep="/"))
tt3 <- ttheme_minimal(
	core=list(bg_params = list(fill = rep(blues9[2:3],Top.n)[1:Top.n], col=NA),
	fg_params=list(fontface=1)),
	colhead=list(fg_params=list(col="navyblue", fontface=1L)))
grid.arrange(
	tableGrob(top.result, theme=tt3, rows=NULL), nrow=1)

write.table(top.result, file=paste(output,"TopDifferentialProteins.txt", sep="/"),
	col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

#	This prints "null device\n 1". Seems to work when removed.
#dev.off()

######################################################################

# 4. Summarize Meta Analysis Results
proteins.datasets = c()
for (i in 1:N.dataset){
	proteins.datasets[[i]] = which(meta.result[paste0("qvalue.set",i)][,] < q.value.cutoff )
}

for (i in 1:N.dataset){
	if( exists("proteins.intersect") ){
		proteins.intersect = intersect(proteins.intersect, proteins.datasets[[i]])
		proteins.union     = union(proteins.union, proteins.datasets[[i]])
	} else {
		proteins.intersect = proteins.datasets[[i]]
		proteins.union     = proteins.datasets[[i]]
	}
}

meta.result.qvalue.sets = cbind( meta.result[names(meta.result)[grepl("qvalue.set", names(meta.result))]])

IDR = IDR.fun(meta.result$qvalue.meta,
	meta.result.qvalue.sets,
	q.value.cutoff)

IRR = IRR.fun(meta.result$qvalue.meta,
	meta.result.qvalue.sets,
	q.value.cutoff)

summary = rbind( sum(meta.result$qvalue.meta<q.value.cutoff, na.rm=TRUE))
for (i in 1:N.dataset){
	summary = rbind( summary, length(proteins.datasets[[i]]) )
}
summary = rbind( summary, length(proteins.intersect), length(proteins.union))

diagnosis = rbind(paste(round(IDR*100,2),"%"),
	paste(round(IRR*100,2), "%"))

rownames_summary = c("Meta Analysis")
for (i in 1:N.dataset){
	rownames_summary = append(rownames_summary,c(paste0("Single Analysis ",i)))
}
rownames_summary = append(rownames_summary,
	c("Intersection among Single Analyses","Union among Single Analyses"))
rownames(summary) = rownames_summary

colnames(summary)=c("# of detected proteins")

rownames(diagnosis)=c("Integration-driven Discovery Rate (IDR)",
	"Integration-driven Revision Rate (IRR)")

colnames(diagnosis)=c("Meta-Analysis Evaluation")

pdf(paste(output,"Summary.pdf",sep="/"), height=4, width=6)

tt1 <- ttheme_minimal(
	core=list(bg_params = list(fill = rep(blues9[2:3],nrow(summary))[1:nrow(summary)], col=NA),
	fg_params=list(fontface=1)),
	colhead=list(fg_params=list(col="navyblue", fontface=1L)))

tt2 <- ttheme_minimal(
	core=list(bg_params = list(fill = rep(blues9[2:3],nrow(diagnosis))[1:nrow(diagnosis)], col=NA),
	fg_params=list(fontface=1)),
	colhead=list(fg_params=list(col="navyblue", fontface=1L)))

grid.arrange(
	tableGrob(summary, theme=tt1),
	tableGrob(diagnosis, theme=tt2))

write.table(summary, file=paste(output,"Summary.txt", sep="/"),
	col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
write.table(diagnosis, file=paste(output,"Diagnosis.txt", sep="/"),
	col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")


#	This prints "null device\n 1". Seems to work when removed.
#dev.off()