library(plotrix)
library(RColorBrewer)

CI_plotter <- function(dd, ci, col, x_offset) {
    ## Purpose: Helper function for gradient_CI_plot; adds points with
    ##    error bars to an existing gradient model component plot.
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   dd (data frame): NOAA sites in rows, model component
    ##      drawdowns (pptv) in columns
    ##   ci: (data frame): NOAA sites in rows, model component drawdowns
    ##      standard errors (pptv) in columns
    ##   col: color sequence for plotting; one color per model component
    ##   x_offset (float): horizontal offset for the points/error bars.
    ## ----------------------------------------------------------------------
    ## Author: Timothy W. Hilton, Date: 16 Feb 2016, 17:27

    n_sites <- length(dd)
    plotCI(1:n_sites + x_offset, dd, ci, col=col, add=TRUE)
}

gradient_CI_plot <- function(dd, ci)
{
  ## Purpose: produce a scatter plot with error bars of STEM model
  ##    component drawdowns at NOAA sites.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##    dd (data frame): NOAA sites in rows, model component
  ##       drawdowns (pptv) in columns
  ##   ci: (data frame): NOAA sites in rows, model component drawdowns
  ##      standard errors (pptv) in columns
  ## ----------------------------------------------------------------------
  ## Author: Timothy W. Hilton, Date: 16 Feb 2016, 17:27

    n_sites <- ncol(dd)
    n_models <- nrow(dd)
    pal <- brewer.pal(n_models, 'Dark2')

    plotCI(1:n_models, dd[[1]], ci[[1]], xaxt='n',
           main="drawdown with error bars",
           ylab='Drawdown',  xlab='site',  col=pal[[1]])
    ## points(x, y, col=pal[[1]], type='l')
    axis(1, at=1:n_models, labels=site_names)

    x_offset <- seq(from=0.0, by=0.1, length.out=n_models)
    mapply(CI_plotter, dd, ci, pal[2:n_models], x_offset)

    legend(x='right', legend=models, pch=marker_sequence, col=pal)
}

models <- c('CASA-GFED3', 'Can-IBIS', 'SiB mech')
site_names <- c('NHA', 'CMA', 'SCA', 'WBI', 'THD')
marker_sequence <- c('o', 'x', 5, '+', '*')
marker_sequence <- c(0, 1, 3, 4, 5)

dd <- as.data.frame(matrix(runif(length(models) * length(site_names)),
                           ncol=length(models)),
                    row.names=models)
ci <- as.data.frame(matrix(runif(length(models) * length(site_names)),
                           ncol=length(models)),
                    row.names=models)
names(dd) <- models
names(ci) <- models
gradient_CI_plot(dd, ci)
