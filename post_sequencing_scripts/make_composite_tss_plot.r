## Read gene/TSS profile files

blocks = seq(from=-3000, to=2975,by=25)

f = read.table("Eugene_KR_H3K9me2_CD4_2x.profile.6000.240.txt", header=TRUE,sep="\t",quote="")
b = read.table("Eugene_WT_H3K9me2_CD4_2x.profile.6000.240.txt", header=TRUE,sep="\t",quote="")

make_composite_plot(f,b,c("KR", "WT"),"Eugene_H3K9me2_CD4_2x",blocks)

f = read.table("Eugene_KR_H3K4me3_CD4.profile.6000.240.txt", header=TRUE,sep="\t",quote="")
b = read.table("Eugene_WT_H3K4me3_CD4.profile.6000.240.txt", header=TRUE,sep="\t",quote="")

make_composite_plot(f,b,c("KR", "WT"),"Eugene_H3K4me2_CD4_2x",blocks)

base = "/media/bigdisk/sequencing/profiles/Eugene"

f = read.table(paste(sep="/", base,"Eugene_KR_H3K9me2_CD4_3x_no_chrX.profile.6000.240.txt"), header=TRUE,sep="\t",quote="")
b = read.table(paste(sep="/", base,"Eugene_WT_H3K9me2_CD4_3x_no_chrX.profile.6000.240.txt"), header=TRUE,sep="\t",quote="")

make_composite_plot(f,b,c("KR", "WT"),"Eugene_H3K9me2_CD4_3x_no_chrX",blocks)

f = read.table(paste(sep="/", base,"Eugene_KR_H3K4me3_CD4_no_chrX.profile.6000.240.txt"), header=TRUE,sep="\t",quote="")
b = read.table(paste(sep="/", base,"Eugene_WT_H3K4me3_CD4_no_chrX.profile.6000.240.txt"), header=TRUE,sep="\t",quote="")

make_composite_plot(f,b,c("KR", "WT"),"Eugene_H3K4me2_CD4_3x_no_chrX",blocks)


make_composite_plot = function(f,b,legend,output_basename,blocks) {
	non_numeric_cols = c("symbol","chr","start","end","strand")
	f = f[-match(non_numeric_cols,names(f))]
	b = b[-match(non_numeric_cols,names(b))]

	f = as.vector(unlist(mean(f)))
	b = as.vector(unlist(mean(b)))

	filename = paste(output_basename,"tss.pdf", sep=".")
	pdf(filename)
	ymax = max(f,b)
	plot(x=blocks,y=f, col="red", pch=20, cex=0.5, xlab="Distance from TSS", ylab="Mean(RPKM)", ylim=c(0,ymax))
	points(x=blocks,y=b, pch=20, cex=0.5)

	legend_v = legend
	legend("topright", legend_v, pch=20, cex=0.5, col=c("red","black"))
	dev.off()
}