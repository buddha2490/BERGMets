metsBoxplot <- function(dat, classvar, metabolite, ylog=F, boxcolors=NULL,
                        ref_line=NULL, title=NULL, subtitle=NULL,xlab=NULL,
                        ylab=NULL, xaxis_size =12,yaxis_size =12,titlesize=14,
                        plot_theme=NULL){


     # DEFINE DATA
     boxdf <- dat
     boxdf$CAT <- droplevels(boxdf[[classvar]])
     boxdf$METABOLITE <- boxdf[[metabolite]]

     # Labels
     if (!is.null(ylab)){
          ylabel <- ylab
     } else { ylabel <- metabolite}
     if (ylog==T) ylabel <- paste0("log(",ylabel,")")
     if (!is.null(xlab)) {
          xlabel <- xlab
     } else { xlabel <- classvar }

     # Theme
     if (is.null(boxcolors)) mycolors <- "grey"
     if (!is.null(boxcolors)) mycolors <- boxcolors
     if (!is.null(plot_theme)) {
          boxtheme <- theme
     } else {
          boxtheme <- theme(axis.text=element_text(size=10),
               axis.title.x=element_text(size=xaxis_size, face="bold"),
               axis.title.y=element_text(size=yaxis_size, face="bold"),
               axis.text.x=element_text(size=xaxis_size-1, face="plain"),
               axis.text.y=element_text(size=yaxis_size-1, face="plain"),
               title = element_text(size=titlesize, face="bold"),
               plot.subtitle = element_text(size=titlesize-2, face="plain"))
     }





     # Draw the initial plot
     g <-  ggplot(boxdf, aes(x=CAT, y=METABOLITE, fill=CAT)) +
          geom_boxplot(fill=mycolors, outlier.shape = NA) +
          geom_jitter(width = 0.3, show.legend=F, size=0.5)


     # Add the labels
     if (!is.null(subtitle)){
          g <- g + ggtitle(title, subtitle=subtitle) + xlab(xlabel) + ylab(ylabel)
     } else {
          g <- g + ggtitle(title) + xlab(xlabel) + ylab(ylabel)
     }

     # Reference line
     if (!is.null(ref_line)) {
          g <- g + geom_hline(yintercept=ref_line, linetype="dashed")
     }

     # Plot Y on a log scale?
     if (ylog==T) {g <- g + scale_y_log10()}

     # print the figure
     g + boxtheme



}




