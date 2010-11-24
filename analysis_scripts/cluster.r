num_blocks = 20
filename = "file"
t = read.table(filename, header=TRUE,sep="\t",quote="")
col_names = unlist(lapply(0:(num_blocks-1), function(x) paste("b",x,sep="")))
data = t[col_names]

# This plot helps us find the optimal # of clusters
wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:50) wss[i] <- sum(kmeans(data, centers=i)$withinss)
pdf("wss.pdf")
plot(1:50, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()

clusters = 15

fit <- kmeans(data, clusters, iter.max=100,nstart=50)
# get cluster means 
#aggregate(data,by=list(fit$cluster),FUN=mean)
# append cluster assignment
data <- data.frame(data, t$symbol)
data <- data.frame(data, fit$cluster)

data = data[order(data$fit.cluster),]
#row.names(data) = data$symbol

pdf(paste("heatmap-",clusters,".pdf",sep=""))
heatmap(data.matrix(data[col_names]),labRow=NA,Rowv=NA, Colv=NA)
dev.off()