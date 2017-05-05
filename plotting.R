require(jsonlite)
require(ggplot2)
require(reshape2)
require(ggdendro)
require(cowplot)

boolean_heatmap = function(pathIn,pathOut){
	Data = do.call(rbind,fromJSON(pathIn))
	pData = apply(Data,2,as.logical)
	Dendro = hclust(dist(t(pData)),'average')
	pDendro = ggdendrogram(Dendro) + theme(axis.text.y=element_blank())
	pData = pData[,Dendro$order]
	pData = melt(pData); colnames(pData) = c('Insertion','Sample','Exists')
	pHeatmap = ggplot(pData,aes(x=Sample,y=Insertion,fill=Exists))+
		geom_tile()+scale_fill_brewer(palette="RdYlBu")+
		theme(
			axis.ticks.x=element_blank(),
			axis.text.x=element_blank(),
			axis.line=element_blank()
		)+guides(fill=F)+
		labs(x=NULL,y=NULL)
	Plot = plot_grid(pDendro,pHeatmap,ncol=1,nrow=2,align='v',rel_heights=c(1,2.5))
	ggsave(Plot,filename=pathOut,device='pdf')
}

if(!interactive()){
	Args = as.list(commandArgs(trailingOnly=T))
	do.call(Args[[1]],Args[2:length(Args)])
}
