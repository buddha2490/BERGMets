metsROC <- function(dat,
                    outcome,
                    compid,
                    covariates,
                    normalize=T,
                    xlabel = "1 - Specificity",
                    ylabel = "Sensitivity",
                    title = NULL,
                    colors = NULL,
                    title_size = 10,
                    subtitle_size = 10,
                    xaxis_size = 9,
                    yaxis_size = 9,
                    plot_theme = NULL) {
     require(ggplot2)


# data checks

# outcome
if (anyNA(dat[[outcome]])) {
     n <- nrow(dat[!is.na(dat[[outcome]]),])
     message(paste0("Outcome includes missing values - only ",n," observations will be used"))
}
if (sum(unique(dat[[outcome]]),na.rm=T) != 1){
     stop("Outcome variable must be coded as [0,1] for a logistic regression")
}
if (!is.numeric(dat[[outcome]])){
     stop("Outcome variable must be coded as a numeric variable")
}

# subset data based on outcome
dat <- dat[!is.na(dat[[outcome]]),]

# normalize mets
if (normalize==T){
   message("Metabolite values will be log transformed, mean centered, and normalized prior to modeling")
     dat <- normalizeMets(dat,compid)
}

# covariates
if (anyNA(dat[,covariates])) {
     warning("You have missing values in your covariates.  These observations will be dropped from the models")
}



     # Build the theme
     if (is.null(plot_theme)) {
          ROCTheme <- theme(
               axis.text.x = element_text(size = xaxis_size, hjust = .5),
               axis.text.y = element_text(size = yaxis_size),
               axis.title.x = element_text(size = xaxis_size + 1),
               axis.title.y = element_text(size = yaxis_size + 1),
               legend.position = "right",
               legend.direction = "vertical",
               legend.text = element_text(size = subtitle_size),
               legend.title = element_text(face = "bold", size = subtitle_size +
                                                1),
               plot.title = element_text(
                    hjust = 0,
                    face = "bold",
                    size = title_size,
                    vjust = 0
               ),
               plot.caption = element_text(hjust = 0, size = subtitle_size),
               plot.subtitle = element_text(size = title_size)
          )
     } else {
          ROCTheme = plot_theme
     }

     # Set all the text

     # Title
     if (is.null(title)) {
          if (length(compid) == 1) {
               ROCtitle <- paste0("ROC Curve: Classifying ",
                                  outcome,
                                  " with ",
                                  compid)
          } else {
               ROCtitle <- paste0("ROC Curves: Classifying ",
                                  outcome,
                                  " with selected metabolites")
          }
     } else {
          ROCtitle <- title
     }



     # Subtitle
     if (is.null(covariates)) {
          ROCsubtitle <- NULL
     } else {
          if (length(covariates) == 1)
               covar.string <- covariates
          if (length(covariates) > 1) {
               last <- covariates[length(covariates)]
               rest <- covariates[covariates != last]
               covar.string <- paste0(paste(rest, collapse = ", "),
                                      " and ", last)
          }
          ROCsubtitle <-
               paste0("Logistic regression models adjusted for ",
                      covar.string)
     }




     # Step 1 - fit the models for each compid
     models <- lapply(compid, function(x) {
          y <- formula(paste0(outcome, "~", paste0(c(
               x, covariates
          ), collapse = "+")))
          fit <- glm(y, data = dat, family = "binomial")
          fitstats <- data.frame(
               fit$model,
               values = fit$fitted.values,
               stringsAsFactors = F
          )
          fitstats$Model <- x
          fitstats <- fitstats[, c(outcome, "values", "Model")]
          names(fitstats) <- c("OUTCOME", "values", "Model")
          return(fitstats)
     }) %>% do.call("rbind", .)
     models$Model <- factor(models$Model)



     # Step 2 - fit the plot
     g <- ggplot(models, aes(d = OUTCOME, m = values, color = Model)) +
          plotROC::geom_roc(n.cuts = 0, size = 1.00) +
          geom_abline(
               slope = 1,
               intercept = 0,
               linetype = "dashed",
               size = 1.00
          )

     auc <- calc_auc(g)
     auc$group <- factor(auc$group, auc$group, compid)
     auc$label <- paste0(auc$group, " (",
                         format(round(auc$AUC, 3), nsmall = 3),
                         ")")

     if (is.null(colors)) {
          g <- g + scale_color_hue(name = "Metabolite (AUC)",
                                   labels = auc$label)
     } else {
          g <- g + scale_color_manual(values = colors,
                                      name = "Metabolite (AUC)",
                                      labels = auc$label)
     }
     g <- g + labs(title = ROCtitle,
                   subtitle = ROCsubtitle) +
          xlab(xlabel) + ylab(ylabel) +
          ROCTheme

     if (length(compid) == 1) {
          AUC <- format(round(auc$AUC, 3), nsmall = 3)
          g <- g +
               ggtitle(paste0(ROCtitle, "\nAUC = ", AUC)) +
               theme_update(
                    axis.text.x = element_text(size = xaxis_size, hjust = .5),
                    axis.text.y = element_text(size = yaxis_size),
                    axis.title.x = element_text(size = xaxis_size + 1),
                    axis.title.y = element_text(size = yaxis_size + 1),
                    plot.title = element_text(face = "bold", size = title_size),
                    legend.position = "none",
                    plot.subtitle = element_text(size = subtitle_size),
               )
     }
     return(g)
}
