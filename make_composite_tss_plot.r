
f = as.list(read.table("scores_f.txt"))$V1
b = as.list(read.table("scores_b.txt"))$V1

x = as.list(seq(-1000,1000))
pdf("tss.pdf")
plot(x=x,y=f, col="red")
points(x=x,y=b)

legend_v = c("ChIP", "Control")
legend("topright", legend_v, pch=1, col=c("red","black"))
