library("marray")
library("geneplotter")
even = function(x) {
	if(x %% 2) return(F) else return(T)
}

odd = function(x) {
	if(x %% 2) return(T) else return(F)
}

smoothscatter_cor = function(x,y,mark,xlab,ylab,x_low_lim,x_high_lim,y_low_lim=NA,y_high_lim=NA,line=T,bw=NA) {
	if(is.na(y_low_lim) | is.na(y_high_lim)) {
		y_low_lim = x_low_lim
		y_high_lim = x_high_lim
	}
	main = paste(sep=" ", mark, xlab, "v", ylab)
	png(paste(sep=".",main,"png"), width=1000, height=1000)
	c = cor(x, y, method = "pearson")
	#cr =colorRampPalette(c("white","#FFFFCC","#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#800026"))
	smoothScatter(x, y, bandwidth=bw, main=main,  xlab=xlab, ylab=ylab, ylim=c(y_low_lim,y_high_lim), xlim=c(x_low_lim,x_high_lim)) #colramp=cr
	if(line) {
		abline(0,1, col="black") # regression line (y~x) 
	}
	mtext(paste(sep="", "r=",c,",  r^2=",c*c),side=3)
	dev.off();
}

scatter_cor = function(x,y,mark,xlab,ylab,x_low_lim,x_high_lim,y_low_lim=NA,y_high_lim=NA,line=T) {
	if(is.na(y_low_lim) | is.na(y_high_lim)) {
		y_low_lim = x_low_lim
		y_high_lim = x_high_lim
	}
	main = paste(sep=" ", mark, xlab, "v", ylab)
	png(paste(sep=".",main,"png"), width=1000, height=1000)
	c = cor(x, y, method = "pearson")
	plot(x, y, main=main,  xlab=xlab, ylab=ylab, ylim=c(y_low_lim,y_high_lim), xlim=c(x_low_lim,x_high_lim))
	if(line) {
		abline(0,1, col="red") # regression line (y~x) 
	}
	mtext(paste(sep="", "r=",c,",  r^2=",c*c),side=3)
	dev.off();
}

scatter_abline_cor = function(x,y,mark,xlab,ylab,x_low_lim,x_high_lim,y_low_lim=NA,y_high_lim=NA,line=T) {
	if(is.na(y_low_lim) | is.na(y_high_lim)) {
		y_low_lim = x_low_lim
		y_high_lim = x_high_lim
	}
	main = paste(sep=" ", mark, xlab, "v", ylab)
	png(paste(sep=".",main,"png"), width=1000, height=1000)
	c = cor(x, y, method = "pearson")
	plot(x, y, main=main,  xlab=xlab, ylab=ylab, ylim=c(y_low_lim,y_high_lim), xlim=c(x_low_lim,x_high_lim))
	if(line) {
		abline(lm(y~x), col="red") # regression line (y~x) 
	}
	mtext(paste(sep="", "r=",c,",  r^2=",c*c),side=3)
	dev.off();
}

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

lappend <- function(lst, obj) {
    lst[[length(lst)+1]] <- obj
    return(lst)
}

remove_rows_with_inf_or_na = function(data) {
	data = data[!apply(data,1,function(y){any(is.na(y))}),]
	data = data[!apply(data,1,function(y){any(!is.finite(y))}),]
	data = data[!apply(data,1,function(y){any(is.nan(y))}),]
	return(data)
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


# Bidirectional colors only apply to fold change values. e.g. vector(a) / vector(b)


make_bidirectional_colors = function(data,ncolors=100) {
	range = color_range(data,ncolors=ncolors)
	return(range)
}

read_table_remove_cols = function(file,cols,sep="	",quote="",header=TRUE) {
	t 				= read.table(file,header=header,sep=sep,quote=quote)
	t 				= t[-match(cols,names(t))]
	return(t)
}

read_table_remove_cols_set_row_names = function(file,cols,sep="	",quote="",header=TRUE) {
	t 				= read.table(file,header=header,sep=sep,quote=quote)
	row.names(t) 	= t$symbol
	t 				= t[-match(cols,names(t))]
	return(t)
}

merge_all = function(data_frames, all.x=F,all.y=F,sort=F,by="row.names"){
	data = merge(data_frames[[1]], data_frames[[2]], all.x=all.x,all.y=all.y,sort=sort, by=by)
	row.names(data) = data[,1]
	data = data[,-(1:1)]
	data_frames = data_frames[-(1:2)]
	for(item in data_frames) {
		data = merge(data,item, all.x=all.x,all.y=all.y,sort=sort, by=by)
		row.names(data) = data[,1]
		data = data[,-(1:1)] # Remove the "row.names" column which was added automatically on merge
	}
	return(data)
}

#Use "color_range" instead of max_mean_plus_sd or min_mean_minus_sd
color_range = function(d,stdevs=3,na.rm=T,ncolors=100) {
	d[d == -Inf] 	= NA
	d[d == Inf] 	= NA
	m = mean(unlist(d),na.rm=na.rm)
	sd = sd(unlist(d),na.rm=na.rm) * stdevs
	max = m+sd
	min = m-sd
	if(max > max(unlist(d))) {
		max = max(unlist(d))
	}
	if(min < min(unlist(d))) {
		min = min(unlist(d))
	}
	print(paste("Color Range Min: ",min, " max: ", max," mean: ", m, " sd: ", sd))
	by = (max-min) / ncolors
	breaks = c(-1000000,seq(from=min,to=max,by=by),1000000)
	return(breaks)
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
	na_to = breaks[[1]][floor(length(breaks[[1]])/2)]
	min_inf_to = breaks[[1]][1]
	max_inf_to = breaks[[1]][length(breaks[[1]])]
	print(paste(sep=" ", "Setting NA/NaN to", na_to))
	print(paste(sep=" ", "Setting -Inf to", min_inf_to))
	print(paste(sep=" ", "Setting Inf to", max_inf_to))
	data[is.na(data)] = na_to #changes NaN also.
	data[data == -Inf] 	= min_inf_to
	data[data == Inf]	= max_inf_to
	
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
				if(even(length(brk))) {
					step = brk[length(brk)] - brk[length(brk)-1]
					brk = c(brk,brk[length(brk)]+step)
				}
				colors = maPalette(low="green", high="red",mid="white",k=(length(brk)-1))
			}
			else if(color_function == "redblue") {
				if(even(length(brk))) {
					step = brk[length(brk)] - brk[length(brk)-1]
					brk = c(brk,brk[length(brk)]+step)
				}
				colors = maPalette(low="blue", high="red",mid="white",k=(length(brk)-1))
			}
			else if(color_function == "whitered") {
				if(even(length(brk))) {
					step = brk[length(brk)] - brk[length(brk)-1]
					brk = c(brk,brk[length(brk)]+step)
				}
				colors = maPalette(low="white", high="red",mid="pink",k=(length(brk)-1))
			}
			else if(color_function == "bluewhite") {
				if(even(length(brk))) {
					step = brk[length(brk)] - brk[length(brk)-1]
					brk = c(brk,brk[length(brk)]+step)
				}
				colors = maPalette(low="blue", high="white",k=(length(brk)-1))
			}
			else if(color_function == "whiteyellowred") {
				if(even(length(brk))) {
					step = brk[length(brk)] - brk[length(brk)-1]
					brk = c(brk,brk[length(brk)]+step)
				}
				colors = maPalette(low="white", mid="yellow", high="red",k=(length(brk)-1))
			}
			data_plot = t(d[,as.vector(cols_per_heatmap[,n])]) # transpose the array
			nr = nrow(data_plot)
			nc = ncol(data_plot)
			image(1:nr,1:nc,data.matrix(data_plot),breaks=brk,col=colors,axes=F,xlab=colnames(cols_per_heatmap)[n],ylab="")
		}
		## Now we add the cluster labels
				
		label_breaks = c(0,seq(from=1,to=length(fit$size),by=1))
		colors = gray(1:(length(label_breaks)-1)/(length(label_breaks)))
		data_plot = t(data.matrix(d$fit.cluster))
		image(1:nrow(data_plot),1:ncol(data_plot),data_plot,breaks=label_breaks,col=colors,axes=F,ylab="",xlab="Clusters")
		dev.off()

		######## OUTPUT SYMBOLS ########

		## Output the symbols in each group.
		write.table(d,file=paste("cluster-",cluster,".",postfix,".txt",sep=""),sep="\t",quote=F)

		print(paste("Finished cluster", cluster,sep=" "))
	}
}

make_hclust_heatmap = function(data, main, colors, breaks) {
	data		= data.matrix(data)
	distance 	= dist(data)
	cluster		= hclust(distance, method="ward")
	dendrogram	= as.dendrogram(cluster)
	Rowv <- rowMeans(data, na.rm = T)
	dendrogram <- reorder(dendrogram, Rowv)
	
	## Produce the heatmap from the calculated dendrogram.
	## Don't allow it to re-order rows because we have already re-ordered them above.
	
	reorderfun = function(d,w) { d }
	png(paste(main,".png",sep=""), res=150, height=22,width=17,units="in")
	
	heatmap(data,col=colors,breaks=breaks,scale="none",Colv=NA,Rowv=dendrogram,labRow=NA, reorderfun=reorderfun)

	dev.off()
	
	
	## Re-order the original data using the computed dendrogram
	rowInd = rev(order.dendrogram(dendrogram))
	di = dim(data)
	nc = di[2L]
	nr = di[1L]
	colInd = 1L:nc
	data_ordered <- data[rowInd, colInd]
	write.table(data_ordered, paste(main,".txt",sep=""),quote=F, sep="\t",row.names=T, col.names=T)
}


heatmap.2.custom = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, hclust_method="ward",
	cols_to_cluster=colnames(x),
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
		if(length(cols_to_cluster) != length(colnames(x))) {
			hcr <- hclustfun(distfun(x[,cols_to_cluster]), method=hclust_method)
		}
		else {
			hcr <- hclustfun(distfun(x), method=hclust_method)
		}
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
		if(length(cols_to_cluster) != length(colnames(x))) {
			hcr <- hclustfun(distfun(x[,cols_to_cluster]), method=hclust_method)
		}
		else {
			hcr <- hclustfun(distfun(x), method=hclust_method)
		}
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)), method=hclust_method)
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)), method=hclust_method)
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x))
            tmpbreaks[length(tmpbreaks)] <- max(abs(x))
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}