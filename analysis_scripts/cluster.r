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
#WTK4_KRK4 			= remove_less_one_and_log(WTK4_KRK4)
#WTK9_KRK9 			= remove_less_one_and_log(WTK9_KRK9)
#WTK4_KRK4_WTK9_KRK9 = remove_less_one_and_log(WTK4_KRK4_WTK9_KRK9)

clusters = c(8,10,15) #16,6,10,20,

#make_wss_plot(WTK4_KRK4, "WTK4_KRK4")
#make_wss_plot(WTK9_KRK9, "WTK9_KRK9")
#make_wss_plot(WTK4_KRK4_WTK9_KRK9, "WTK4_KRK4_WTK9_KRK9")

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
