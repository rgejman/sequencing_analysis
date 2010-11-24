num_blocks = 80
filename = "Eugene_WT_H3K9me2_CD4_all_profile.txt"
t = read.table(filename, header=TRUE,sep="\t",quote="")
col_names = unlist(lapply(0:(num_blocks-1), function(x) paste("b",x,sep="")))
data = t[col_names]
row.names(data) = t$symbol


# This plot helps us find the optimal # of clusters
wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:50) wss[i] <- sum(kmeans(data, centers=i)$withinss)
pdf("wss.pdf")
plot(1:50, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()

clusters = c(15) #16,6,10,20,

for(cluster in clusters) {
	fit <- kmeans(data, cluster, iter.max=100,nstart=50)
	# get cluster means 
	#aggregate(data,by=list(fit$cluster),FUN=mean)
	# append cluster assignment
	d <- data.frame(data, t$symbol)
	d <- data.frame(d, fit$cluster)

	d = d[order(d$fit.cluster),]
	
	color_limits = c(0,10)

	pdf(paste("heatmap-",cluster,".pdf",sep=""))
	heatmap(data.matrix(d[col_names]), labRow=NA, Rowv=NA, Colv=NA, col=rev(gray(1:256/256)),scale="none",zlim=color_limits)
	dev.off()
}

