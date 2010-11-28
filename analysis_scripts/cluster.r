num_blocks = 80
base_col_names = unlist(lapply(0:(num_blocks-1), function(x) paste("b",x,sep="")))
non_numeric_col_names = c("symbol","chr","start","end","strand")

f1 = "Eugene_WT_H3K4me3_CD4.profile.6000.80.txt"
f2 = "Eugene_KR_H3K4me3_CD4.profile.6000.80.txt"

t1 = read.table(f1, header=TRUE,sep="\t",quote="")
t2 = read.table(f2, header=TRUE,sep="\t",quote="")

data1 = t1[-match(non_numeric_col_names,names(t1))]
data2 = t2[-match(non_numeric_col_names,names(t2))]

row.names(data1) = t1$symbol
row.names(data2) = t2$symbol

data = merge(data1,data2, all.x=FALSE,all.y=TRUE,sort=FALSE, by="row.names")
# Remove the "row.names" column which was added automatically on merge
data = data[,-(1:1)]

# We will now log-transform the data.
# To do this we should set all RPKM values <1 to some constant (I chose 1 because log(1) = 0)
less_than_one = (data < 1)
data[less_than_one] = 1
data = log(data)


# This plot helps us find the optimal # of clusters
wss <- (nrow(data)-1)*sum(apply(data,2,var))
for (i in 2:25) {
	wss[i] <- sum(kmeans(data, iter.max=50,nstart=5, centers=i)$withinss)
	print(paste("Cluster for WSS", i,"done",sep=" "))
}
pdf("wss.pdf")
plot(1:25, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()

clusters = c(10,15) #16,6,10,20,

for(cluster in clusters) {
	
	######## DETERMINE # OF BREAKS ########
	
	# K9 columns have mean=~2 and sd=~2 (without log)
	# K4 columns have mean=~2 and sd that ranges from 7-25
	
	# FOR K9 Columns:
	
	max_mean_plus_sd = max(mean(data)+sd(data))
	#breaks = c(seq(from=0,to=(max(data)-1),by=(max(data)-1)/255),10000)
	breaks = c(seq(from=0,to=max_mean_plus_sd,by=0.5),10000)
	colors = rev(gray(0:(length(breaks)-2)/(length(breaks)-2)))
	
	# For K4 Columns:
	#breaks = c(seq(from=0,to=25,by=1),10000)
	#colors = rev(gray(0:(length(breaks)-2)/(length(breaks)-2)))
	
	######## CLUSTERING ########
	
	fit <- kmeans(data, cluster, iter.max=50,nstart=5)
	# get cluster means 
	#aggregate(data,by=list(fit$cluster),FUN=mean)
	# append cluster assignment
	d <- data.frame(data, t1$symbol)
	names(d)[names(d)=="t1.symbol"] = "symbol"
	d <- data.frame(d, fit$cluster)

	d = d[order(d$fit.cluster),]
	
	png(paste("heatmap-",cluster,".png",sep=""), res=300, height=44,width=34,units="in")
	#pdf(paste("heatmap-",cluster,".pdf",sep=""))
	
	par(mfrow=c(1,2))
		
	data_plot = t(d[,colnames(data1)]) # transpose the array
	image(1:nrow(data_plot),1:ncol(data_plot),data.matrix(data_plot),breaks=breaks,col=colors,axes=FALSE,xlab="WT H4K9me2",ylab="")
	data_plot = t(d[,colnames(data2)])
	image(1:nrow(data_plot),1:ncol(data_plot),data.matrix(data_plot),breaks=breaks,col=colors,axes=FALSE,xlab="KR H4K9me2",ylab="")
	dev.off()
	
	#

	#png(paste("heatmap-",cluster,".png",sep=""), res=300, height=44,width=34,units="in")
	#par(mfrow=c(1,2))
	#heatmap(data.matrix(d[colnames(data1),]),labRow=NA,Rowv=NA,Colv=NA,col=colors,scale="none",zlim=color_limits)
	#heatmap(data.matrix(d[colnames(data2),]),labRow=NA,Rowv=NA,Colv=NA,col=colors,scale="none",zlim=color_limits)
	#dev.off()
	
	## Output the symbols in each group.
	write.table(d,file=paste("cluster-",cluster,".txt",sep=""),sep="\t",quote=FALSE)
	
	print(paste("Finished cluster", cluster,sep=" "))
}