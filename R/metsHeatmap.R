

metsHeatmap <- function(dat,compid,covariates=NULL,biochem=NULL,method="pearson",
                        hclust=T,title,subtitle=NULL,plot_colors=NULL,
                        xaxis_size=7, yaxis_size=7,
                        title_size=12, subtitle_size=8,plot_theme=NULL){

     require(ggplot2, quietly=T)


     # Set the color scheme
     if (is.null(plot_colors)) colors <- c("blue","white","red")
     if (!is.null(plot_colors)) colors <- plot_colors



     # Set the theme - users can provide their own
     if (is.null(plot_theme)){
     heatTheme <- theme(
        axis.text.x.bottom=element_text(size=xaxis_size,hjust=.5),
        axis.text.y.left=element_text(size=yaxis_size),
        legend.position="right",
        legend.direction = "vertical",
        #legend.key.height=unit(2,"cm"),
        legend.title=element_blank(),
        plot.margin=unit(c(-1,0.5,-1,0.5),"cm"),
        plot.title=element_text(hjust=0,face="bold",size=title_size,vjust=0),
        plot.caption = element_text(hjust = 0,size=subtitle_size))
     } else {
          heatTheme <- plot_theme
     }

     # Capitalize the correlation method (i.e., "pearson -> "Pearson)
     capMethod <- paste0(toupper(substr(method,1,1)),
                         tolower(substr(method,2,nchar(method))))

     # Create a subtitle describing the correlations
     # I.e. "Simple pairwise Pearson correlation coefficients"
     # I.e. "Partial Pearson correlation coefficients adjusted for Age_int and Race

     # Create a covariate vector
     if (!is.null(covariates)){
          l <- length(covariates)
          if (l==1) covar_string<- covariates
          if (l > 1){
               last <- covariates[l]
               rest <- covariates[covariates != last]
               covar_string <- paste0(paste(rest, collapse=", "),
                                      ", and ", last)
          }
     }
     # Construct the subtitle
     # Users can supply their own, or I'll make one for them
     if (is.null(subtitle)){
          if (is.null(covariates)) {
               sub <- paste0("Simple pairwise ",capMethod," correlation coefficients")
          } else {
               sub <- paste0("Partial ",capMethod, " correlation coefficients adjusted for\n",
                             covar_string)
          }
     } else sub <- subtitle

     # Create the plots

     # Step 1 - generate the correlation matrix
     suppressMessages(
     mat <- dplyr::select(corMatrix(dat,compid,compid,covariates,method),-COMP_ID)
     )
     row.names(mat) <- names(mat)


     # Step 2 - generate the heatmap
     g <- ggcorrplot::ggcorrplot(mat,hc.order=hclust,
                                 outline.color="white",title=title,
                                 legend.title="Correlations",
                                 colors = colors)

     # Step 3 - add some titles and stuff
     g <- g + labs(title=title,caption=sub)


     # By default, the x- and y-axis labels will just be the COMP_IDs
     # If the user supplies a biochem option, BIOCHEMICAL names will be added
     suppressWarnings(
     if (!is.null(biochem)){
     v <- biochem[biochem$COMP_ID %in% g$data$Var2,c("COMP_ID","BIOCHEMICAL")] # subset to just the COMP_IDs I need
     order <- data.frame(COMP_ID=unique(g$data$Var2),order=1:nrow(v),stringsAsFactors=F)
     v <- left_join(v,order,"COMP_ID")
     v <- v[order(v$order),]
     v$Label <- paste0(v$order,". ",v$BIOCHEMICAL)
     g$data$Var2 <- factor(g$data$Var2,v$COMP_ID,v$Label)
     g$data$Var1 <- factor(g$data$Var1,v$COMP_ID,v$order)
     }
     )
}

