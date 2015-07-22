
#######################################################
# plot relavant expression matrix
plotEset<-function(eset, df){
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  if(is.data.frame(df)){
    gene_ids <- df$Gene
  }else{
    gene_ids <- df
  }
  print(gene_ids)
  if(length(gene_ids)>10) stop("too many genes")
  mt_df <- as.data.frame(exprs(eset)[gene_ids,])
  mt_df$Gene <- rownames(mt_df)
  # convert mt_df to long table
  mt_df <- melt(mt_df, id.var="Gene", variable.name="ID", value.name="count")
  mt_df$logFC <- log(mt_df$count + 0.5, 2)
  # convert dose/time variable to factor
  mt_df <- merge(mt_df, pData(eset), by="ID") %>% mutate(., dose = factor(dose, levels=c(1, 3, 10)), time=factor(time, levels=sort(unique(time))))
  
  ## plot control over time
  p1 <- ggplot(subset(mt_df, drug=="D"), aes(time, logFC, color=Gene))
  p1 <- p1 + geom_jitter(size=3, position = position_jitter(width = .1), alpha=0.7, aes(type=Gene))
  p1 <- p1 + stat_summary(fun.y=mean, geom="line", aes(group=Gene))
  p1 <- p1 + scale_color_manual(values=colorRampPalette(brewer.pal(10, "Paired"))(length(unique(mt_df$Gene))), guide=guide_legend(keywidth=0.5, keyheight=0.5))
  p1 <- p1 + ggtitle("Control over time")
  #p1
  ## plot drug over time, at intermediate concentration
  p2 <- ggplot(subset(mt_df, drug!="D" & dose==3), aes(time, logFC, color=Gene))
  p2 <- p2 + facet_grid(.~ drug)
  p2 <- p2 + geom_jitter(size=3, position = position_jitter(width = .1), alpha=0.7, aes(type=Gene))
  p2 <- p2 + stat_summary(fun.y=mean, geom="line", aes(group=Gene))
  p2 <- p2 + scale_color_manual(values=colorRampPalette(brewer.pal(10, "Paired"))(length(unique(mt_df$Gene))), guide=guide_legend(keywidth=0.5, keyheight=0.5))
  p2 <- p2 + ggtitle("Drug median dose over time")
  #p2
  ## plot drug concentration at 
  p3 <- ggplot(subset(mt_df, drug!="D" & time==24), aes(dose, logFC, color=Gene))
  p3 <- p3 + facet_grid(.~ drug)
  p3 <- p3 + geom_jitter(size=3, position = position_jitter(width = .1), alpha=0.7, aes(type=Gene))
  p3 <- p3 + stat_summary(fun.y=mean, geom="line", aes(group=Gene))
  p3 <- p3 + scale_color_manual(values=colorRampPalette(brewer.pal(10, "Paired"))(length(unique(mt_df$Gene))), guide=guide_legend(keywidth=0.5, keyheight=0.5))
  p3 <- p3 + ggtitle("Drug dose series at 1day")
  #p3
  deco <- function(p, yrange){
    p <- p + ylab("Expression level") 
    p <- p + scale_y_continuous(limits=yrange)
    p <- p + theme(plot.margin = unit(c(0.2,0.1,0.2,0.1), "cm"))
    p <- p + theme( panel.background = element_rect(fill='white',colour='black'))
    p <- p + theme( panel.grid.major.y = element_line(linetype=3, color="darkgray"),  panel.grid.major.x = element_line(linetype=1, color="lightgray"))
    p <- p + theme( legend.text = element_text(size = rel(.7)), axis.text = element_text(size=rel(.7)), axis.title= element_text(size=rel(.8)))
    return(p)
  }
  yrange <- range(mt_df$logFC)
  p1 <- deco(p1, yrange)
  p2 <- deco(p2, yrange)
  p3 <- deco(p3, yrange)
  p1 <- p1 + theme(legend.position="none")
  grid.newpage() # Open a new page on grid device
  pushViewport(viewport(layout = grid.layout(2, 5)))
  print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2:5)) 
  print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 2:5))
}

# 
# #######################################################
# # plot relavant expression matrix
# plotEset<-function(eset, df){
#   library(reshape2)
#   library(ggplot2)
#   library(RColorBrewer)
#   library(gridExtra)
#   if(is.data.frame(df)){
#     gene_ids <- df$Gene
#   }else{
#     gene_ids <- df
#   }
#   print(gene_ids)
#   mt_df <- as.data.frame(exprs(eset)[gene_ids,])
#   mt_df$Gene <- rownames(mt_df)
#   # convert mt_df to long table
#   mt_df <- melt(mt_df, id.var="Gene", variable.name="ID", value.name="logFC")
#   # convert dose/time variable to factor
#   mt_df <- merge(mt_df, pData(eset), by="ID") %>% mutate(., dose = factor(dose, levels=c(1, 3, 10)), time=factor(time, levels=sort(unique(time))))
#   
#   ## plot control over time
#   p1 <- ggplot(subset(mt_df, drug=="D"), aes(time, logFC, color=Gene))
#   p1 <- p1 + geom_jitter(size=3, position = position_jitter(width = .1), alpha=0.7, aes(type=Gene))
#   p1 <- p1 + stat_summary(fun.y=mean, geom="line", aes(group=Gene))
#   p1 <- p1 + scale_color_manual(values=colorRampPalette(brewer.pal(10, "Paired"))(length(unique(mt_df$Gene))), guide=guide_legend(keywidth=0.5, keyheight=0.5))
#   p1 <- p1 + ggtitle("Control over time")
#   #p1
#   ## plot drug over time, at intermediate concentration
#   p2 <- ggplot(subset(mt_df, drug!="D" & dose==3), aes(time, logFC, color=Gene))
#   p2 <- p2 + geom_jitter(size=3, position = position_jitter(width = .1), alpha=0.7, aes(type=Gene))
#   p2 <- p2 + stat_summary(fun.y=mean, geom="line", aes(group=Gene))
#   p2 <- p2 + scale_color_manual(values=colorRampPalette(brewer.pal(10, "Paired"))(length(unique(mt_df$Gene))), guide=guide_legend(keywidth=0.5, keyheight=0.5))
#   p2 <- p2 + ggtitle("Drug median dose over time")
#   #p2
#   ## plot drug concentration at 
#   p3 <- ggplot(subset(mt_df, drug!="D" & time==24), aes(dose, logFC, color=Gene))
#   p3 <- p3 + geom_jitter(size=3, position = position_jitter(width = .1), alpha=0.7, aes(type=Gene))
#   p3 <- p3 + stat_summary(fun.y=mean, geom="line", aes(group=Gene))
#   p3 <- p3 + scale_color_manual(values=colorRampPalette(brewer.pal(10, "Paired"))(length(unique(mt_df$Gene))), guide=guide_legend(keywidth=0.5, keyheight=0.5))
#   p3 <- p3 + ggtitle("Drug dose series at 1day")
#   #p3
#   deco <- function(p, yrange){
#     p <- p + ylab("Expression level") 
#     p <- p + scale_y_continuous(limits=yrange)
#     p <- p + theme(plot.margin = unit(c(0.2,0.1,0.2,0.1), "cm"))
#     p <- p + theme( panel.background = element_rect(fill='white',colour='black'))
#     p <- p + theme( panel.grid.major.y = element_line(linetype=3, color="darkgray"),  panel.grid.major.x = element_line(linetype=1, color="lightgray"))
#     p <- p + theme( legend.text = element_text(size = rel(.7)), axis.text = element_text(size=rel(.7)), axis.title= element_text(size=rel(.8)))
#     return(p)
#   }
#   yrange <- range(mt_df$logFC)
#   p1 <- deco(p1, yrange)
#   p2 <- deco(p2, yrange)
#   p3 <- deco(p3, yrange)
#   #pl <- list(p1, p2, p3)
#   grid.arrange(p1, p2, p3, ncol=3)
# }