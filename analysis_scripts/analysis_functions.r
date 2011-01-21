normalize_data = function(data, inf_to=1,na_to=1,lt=0.1, gt=NA) {
	data[data == Inf] = inf_to
	data[data == -Inf] = inf_to
	data[is.na(data)] = na_to
	if(is.numeric(lt)) {
		data[data < lt] = lt
	}
	if(is.numeric(gt)) {
		data[data > gt] = gt
	}
	return(data);
}

remove_rows_with_inf_or_na = function(data) {
	
	return(data[apply(data,1,none_are_na),])
	return(data[apply(data,1,none_are_inf),])
	return(data[apply(data,1,none_are_nan),])
}

none_are_na = function(l) {
	return(!any(is.na(l)))
}

none_are_inf = function(l) {
	return(!any(is.inf(l)))
}

none_are_nan = function(l) {
	return(!any(is.nan(l)))
}

remove_rows_with_mean_lt = function(data, lt) {
	return(data[apply(data,1,mean) > lt,])
}

make_colors = function(data,by=0.5,max=10000) {
	colors = c(seq(from=0,to=max_mean_plus_sd(data),by=by),max)
	return(colors)
}

read_table_remove_cols_set_row_names = function(file,cols,sep="	",quote="",header=TRUE) {
	t 				= read.table(file,header=header,sep=sep,quote=quote)
	row.names(t) 	= t$symbol
	t 				= t[-match(cols,names(t))]
	return(t)
}

merge_all = function(data_frames, all.x=FALSE,all.y=TRUE,sort=FALSE,by="row.names"){
	data = data.frame()
	last_names = c()
	for(item in data_frames) {
		if(length(last_names) == 0) {
			last_names = row.names(item)
		}
		if(last_names != row.names(item)) {
			print("Error! Row names do not all match in merge_all()")
		}
	}
	for(item in data_frames) {
		data = merge(data,item, all.x=all.x,all.y=all.y,sort=sort, by=by)
		data = data[,-(1:1)] # Remove the "row.names" column which was added automatically on merge
		row.names(data) = row.names(item) # this will set row.names multiple times, but that's OK
	}
	return(data)
}

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
make_heatmaps = function(data, breaks, cols_per_heatmap, clusters, postfix, color_function="gray", cols_to_cluster=c()) {
	symbol = row.names(data)
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
			if(color_function == "gray") {
				colors = rev(gray(0:(length(brk)-2)/(length(brk)-2)))
			}
			else if(color_function == "redgreen") {
				colors = redgreen((length(brk)-2)/(length(brk)-2)+1)
			}
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