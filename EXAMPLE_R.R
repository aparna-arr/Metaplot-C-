library(ggplot2)
library(reshape2)
pdf(file="metaplot_outfile.pdf", width=12, height=8)
plot<-read.table("metaplot_outfile.txt", header=T)
plot.melt<-melt(plot[,c('bp', 'cleangenes')], id.vars=1)
ggplot(plot.melt, aes(x=bp, y=value, colour=variable, group=variable)) +
geom_line()+
geom_smooth() +
#geom_point()+
#theme_bw() +
ggtitle("Metaplot") 
#theme(panel.grid.minor=element_blank()) +
#scale_colour_brewer(palette="Set1", name="Bed")
