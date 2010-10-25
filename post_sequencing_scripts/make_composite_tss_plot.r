f = as.list(read.table("scores_f.txt"))$V1
b = as.list(read.table("scores_b.txt"))$V1

x = as.list(seq(-1000,1000))
filename = paste(tail(strsplit(getwd(),"/")[[1]],1),"tss.pdf", sep=" ")
pdf(filename)
ymax = max(f,b)
plot(x=x,y=f, col="red", pch=20, cex=0.5, xlab="Distance from TSS", ylab="Cumulative Score", ylim=c(0,ymax))
points(x=x,y=b, pch=20, cex=0.5)

legend_v = c("ChIP", "Control")
legend("topright", legend_v, pch=20, cex=0.5, col=c("red","black"))
