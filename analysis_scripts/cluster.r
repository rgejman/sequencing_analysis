max_mean_plus_sd = function(d,stdevs=2) {
	return(max(mean(d, trim=0.05)+(stdevs * sd(d))));
}

remove_less_one_and_log = function(d) {
	# Changes values <1 to 1 and log-transforms the data
	# To do this we should set all RPKM values <1 to some constant (I chose 1 because log(1) = 0)
	less_than_one = (d < 1)
	d[less_than_one] = 1
	d = log(d)
	return(d)
}

sizes_to_breaks = function(sizes) {
	breaks = c(0)
	last_index = 0
	for(size in sizes) {
		last_index = size+last_index
		breaks = c(breaks,last_index)
	}
	return(breaks)
}

make_wss_plot = function(data,name,n=25) {
	# This plot helps us find the optimal # of clusters
	wss <- (nrow(data)-1)*sum(apply(data,2,var))
	for (i in 2:n) {
		wss[i] <- sum(kmeans(data, iter.max=50,nstart=5, centers=i)$withinss)
		print(paste("Cluster for WSS", i,"done",sep=" "))
	}
	pdf(paste("wss",".",name,".pdf",sep=""))
	plot(1:n, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
	dev.off()
}

#data = data.frame of all values
#cols_per_heatmap = data.frame like this:
# WTK9	KRK9	WTK4	KRK4
# b0	b0		b0		b0
# b1	b1		b1		b1
# etc	etc		etc		etc
# data_to_add is useful if you want to cluster on one set of data, but put another cluster aside it.
make_heatmaps = function(data, breaks, symbol, cols_per_heatmap, clusters, postfix, cols_to_cluster=c()) {
	if(length(cols_to_cluster) == 0) {
		cols_to_cluster = colnames(data)
	}
	for(cluster in clusters) {
		######## CLUSTERING ########
		fit <- kmeans(subset(data, select=cols_to_cluster), cluster, iter.max=50,nstart=10)

		# get cluster means 
		# append cluster assignment
		d <- data.frame(data, symbol)
		d <- data.frame(d, fit$cluster)
		

		d = d[order(d$fit.cluster),]

		######## MAKE FIGURES ########

		png(paste("heatmap-",cluster,".",postfix,".png",sep=""), res=300, height=44,width=34,units="in")
		#pdf(paste("heatmap-",cluster,".pdf",sep=""))
		
		nheatmaps = ncol(cols_per_heatmap)

		par(mfrow=c(1,nheatmaps+1))
		
		for(n in 1:nheatmaps) {
			brk = unlist(breaks[n])
			colors = rev(gray(0:(length(brk)-2)/(length(brk)-2)))
			data_plot = t(d[,as.vector(cols_per_heatmap[,n])]) # transpose the array
			nr = nrow(data_plot)
			nc = ncol(data_plot)
			image(1:nr,1:nc,data.matrix(data_plot),breaks=brk,col=colors,axes=FALSE,xlab=colnames(cols_per_heatmap)[n],ylab="")
		}
		## Now we add the cluster labels
				
		label_breaks = c(0,seq(from=1,to=length(fit$size),by=1))
		colors = gray(1:(length(label_breaks)-1)/(length(label_breaks)))
		data_plot = t(data.matrix(d$fit.cluster))
		image(1:nrow(data_plot),1:ncol(data_plot),data_plot,breaks=label_breaks,col=colors,axes=FALSE,ylab="",xlab="Clusters")
		dev.off()

		######## OUTPUT SYMBOLS ########

		## Output the symbols in each group.
		write.table(d,file=paste("cluster-",cluster,".",postfix,".txt",sep=""),sep="\t",quote=FALSE)

		print(paste("Finished cluster", cluster,sep=" "))
	}
}

num_blocks = 80
#base_col_names = unlist(lapply(0:(num_blocks-1), function(x) paste("b",x,sep="")))
non_numeric_col_names = c("symbol","chr","start","end","strand")

f1 = "Eugene_WT_H3K4me3_CD4.profile.6000.80.txt"
f2 = "Eugene_KR_H3K4me3_CD4.profile.6000.80.txt"
f3 = "Eugene_WT_H3K9me2_CD4_all.profile.6000.80.txt"
f4 = "Eugene_KR_H3K9me2_CD4_all.profile.6000.80.txt"

t1 = read.table(f1, header=TRUE,sep="\t",quote="")
t2 = read.table(f2, header=TRUE,sep="\t",quote="")
t3 = read.table(f3, header=TRUE,sep="\t",quote="")
t4 = read.table(f4, header=TRUE,sep="\t",quote="")


WTK4 = t1[-match(non_numeric_col_names,names(t1))]
KRK4 = t2[-match(non_numeric_col_names,names(t2))]
WTK9 = t3[-match(non_numeric_col_names,names(t3))]
KRK9 = t4[-match(non_numeric_col_names,names(t4))]


row.names(WTK4) = t1$symbol
row.names(KRK4) = t2$symbol
row.names(WTK9) = t3$symbol
row.names(KRK9) = t4$symbol

WTK4_KRK4 			= merge(WTK4,KRK4, all.x=FALSE,all.y=TRUE,sort=FALSE, by="row.names")
WTK4_KRK4 			= WTK4_KRK4[,-(1:1)] # Remove the "row.names" column which was added automatically on merge

WTK9_KRK9 			= merge(WTK9,KRK9, all.x=FALSE,all.y=TRUE,sort=FALSE, by="row.names")
WTK9_KRK9 			= WTK9_KRK9[,-(1:1)] # Remove the "row.names" column which was added automatically on merge

# Again, we have to assume that the rows are ordered by symbol name and that the symbol names are the same, etc.
WTK4_KRK4_WTK9_KRK9 = merge(WTK4_KRK4,WTK9_KRK9, all.x=FALSE,all.y=TRUE,sort=FALSE, by="row.names")
WTK4_KRK4_WTK9_KRK9 = WTK4_KRK4_WTK9_KRK9[,-(1:1)] # Remove the "row.names" column which was added automatically on merge


### DATA NORMALIZATION ###
## OPTIONS:
## 1. Change values <1 to 1 and log-transform the data. This smoothes everything and is nice to look for big trends.
## 2. Don't normalize. The values get capped automatically in the heatmap making stage.
# Changes values <1 to 1 and log-transform the data
WTK4_KRK4 			= remove_less_one_and_log(WTK4_KRK4)
WTK9_KRK9 			= remove_less_one_and_log(WTK9_KRK9)
WTK4_KRK4_WTK9_KRK9 = remove_less_one_and_log(WTK4_KRK4_WTK9_KRK9)

clusters = c(8,10,15) #16,6,10,20,

make_wss_plot(WTK4_KRK4, "WTK4_KRK4")
make_wss_plot(WTK9_KRK9, "WTK9_KRK9")
make_wss_plot(WTK4_KRK4_WTK9_KRK9, "WTK4_KRK4_WTK9_KRK9")

bk = c(seq(from=0,to=max_mean_plus_sd(WTK4_KRK4),by=0.5),10000)
breaks = list(bk,bk)
cols_per_heatmap = data.frame(WTK4=colnames(WTK4), KRK4=colnames(KRK4))
make_heatmaps(WTK4_KRK4, breaks, t1$symbol, cols_per_heatmap, clusters, "WTK4_KRK4")

bk = c(seq(from=0,to=max_mean_plus_sd(WTK9_KRK9),by=0.5),10000)
breaks = list(bk,bk)
cols_per_heatmap = data.frame(WTK9=colnames(WTK9), KRK9=colnames(KRK9))
make_heatmaps(WTK9_KRK9, breaks, t1$symbol, cols_per_heatmap, clusters, "WTK9_KRK9")

k4bk = c(seq(from=0,to=max_mean_plus_sd(WTK4_KRK4),by=0.5),10000)
k9bk = c(seq(from=0,to=max_mean_plus_sd(WTK9_KRK9),by=0.5),10000)
breaks = list(k4bk,k4bk,k9bk,k9bk)
cols_per_heatmap = data.frame(WTK4=colnames(WTK4), KRK4=colnames(KRK4), WTK9=colnames(WTK9), KRK9=colnames(KRK9))
make_heatmaps(WTK4_KRK4_WTK9_KRK9, breaks, t1$symbol, cols_per_heatmap, clusters, "WTK4_KRK4_WTK9_KRK9")

## Cluster only on K9
k4bk = c(seq(from=0,to=max_mean_plus_sd(WTK4_KRK4),by=0.5),10000)
k9bk = c(seq(from=0,to=max_mean_plus_sd(WTK9_KRK9),by=0.5),10000)
breaks = list(k4bk,k4bk,k9bk,k9bk)
cols_per_heatmap = data.frame(WTK4=colnames(WTK4), KRK4=colnames(KRK4), WTK9=colnames(WTK9), KRK9=colnames(KRK9))
cols_to_cluster = c(colnames(WTK9), colnames(KRK9))
make_heatmaps(WTK4_KRK4_WTK9_KRK9, breaks, t1$symbol, cols_per_heatmap, clusters, "WTK4_KRK4_WTK9_KRK9_clustered_on_K9", cols_to_cluster)
