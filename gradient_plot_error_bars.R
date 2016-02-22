library(plotrix)
library(RColorBrewer)
library(tidyr)  # could also use reshape2

row_normalizer <- function(dd, norm_site='NHA') {
    norm_data <- dd[norm_site, ]
    for (i in seq(1, nrow(dd))) {
        dd[i, ] <- dd[i, ] / norm_data
    }
    return(dd)
}

ci_normalizer <- function(dd, ci, norm_site='NHA') {
    dd_hi <- dd + ci
    dd_lo <- dd - ci

    ddnorm <- row_normalizer(dd, norm_site)
    ddnorm_hi <- row_normalizer(dd_hi, norm_site)
    ddnorm_lo <- row_normalizer(dd_lo, norm_site)

    ci_hi <- ddnorm_hi - ddnorm
    ci_lo <- ddnorm - ddnorm_lo

    return(list(dd=ddnorm, ci_hi=ci_hi, ci_lo=ci_lo))
}


dummy_data <- function() {
    ## return a dummy set of random "observations" and "confidence
    ## intervals" with the same row and column labels that the real
    ## data will have.  Useful for testing the plotting code and
    ## whether plotrix is able to produce the plot I want.
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
    return(list(dd=dd, ci=ci))
}

##' .. Helper function for gradient_CI_plot
##'
##' adds points with error bars to an existing gradient model
##' component plot.
##' @title
##' @param dd (data frame): NOAA sites in rows, model component
##' drawdowns (pptv) in columns
##' @param ci (data frame): NOAA sites in rows, model component
##' drawdowns standard errors (pptv) in columns
##' @param col color sequence for plotting; one color per model
##' component
##' @param x_offset (float): horizontal offset for the points/error
##' bars.
##' @return
##' @author Timothy W. Hilton
##' @export
CI_plotter <- function(dd, ci_hi, ci_lo, col, x_offset, marker) {
    n_sites <- length(dd)
    plotCI(1:n_sites + x_offset, dd, uiw=ci_hi, liw=ci_lo,
           col=col, pch=marker, add=TRUE)
}

##' .. produce a scatter plot with error bars of STEM model component
##' drawdowns at NOAA sites.
##'
##' .. content for \details{} ..
##' @title
##' @param df (data frame): NOAA sites in rows, model component
##' drawdowns (pptv) in columns
##' @param dd_col (string): name of the column in df containing
##' drawdown values
##' @param se_col (string): name of the column in df containing
##' standard errors
##' @return
##' @author Timothy W. Hilton
##' @export
gradient_CI_plot <- function(df, dd_col='dd', se_col='dd_se_neff', norm=FALSE) {
    n_sites <- nlevels(df[['site_code']])
    n_models <- nlevels(df[['model']])
    site_names <- levels(df[['site_code']])
    models <- levels(df[['model']])
    pal <- brewer.pal(n_models, 'Dark2')
    marker_sequence <- c(0, 1, 3, 4, 5)

    ## dfw_dd: "data frame, wide; drawdown"
    dfw_dd <- spread(df[, c('model', 'site_code', 'dd')],
                     model, dd)
    dfw_dd <- dfw_dd[, !names(dfw_dd) %in% "site_code"]
    row.names(dfw_dd) <- site_names
    dfw_se <- spread(df[, c('model', 'site_code', 'dd_se_neff')],
                     model, dd_se_neff)
    dfw_se <- dfw_se[, !names(dfw_se) %in% "site_code"]
    row.names(dfw_se) <- site_names
    ## dfw_dd: "data frame, wide; confidence interval"
    dfw_ci <- dfw_se * 1.96 ## 95% confidence interval

    if (norm) {
        dd_list <- ci_normalizer(dfw_dd, dfw_ci, 'NHA')
        t_str <- "normalized drawdown with 95% confidence intervals"
    } else {
        dd_list <- list(dd=dfw_dd, ci_hi=dfw_ci, ci_lo=dfw_ci)
        t_str <- "drawdown with 95% confidence intervals"
    }

    plotCI(1:n_sites,
           dd_list[['dd']][[1]],
           uiw=dfw_ci[['ci_hi']][[1]],
           liw=dfw_ci[['ci_lo']][[1]],
           xaxt='n',
           main=t_str,
           ylab='Drawdown  (pptv)',
           xlab='site',
           col=pal[[1]],
           ylim=range(cbind(dd_list[['dd']] + dd_list[['ci_hi']],
                            dd_list[['dd']] - dd_list[['ci_lo']])),
           xlim=c(1, n_sites + 1),
           pch=marker_sequence[[1]])
    ## points(x, y, col=pal[[1]], type='l')
    axis(1, at=1:n_sites, labels=site_names)

    x_offset <- seq(from=0.0, by=0.05, length.out=n_models)
    mapply(CI_plotter, dd_list[['dd']], dd_list[['ci_hi']], dd_list[['ci_lo']],
           pal, x_offset, marker_sequence[1:n_models])

    legend(x='right', legend=models, pch=marker_sequence, col=pal)
    return(list(dd=dfw_dd, ci=dfw_ci))
}

df <- read.csv('./model_components_18Feb.csv')
fig3a_data <- df[df[['site_code']] %in% c('NHA', 'CMA', 'SCA') &
                     df[['model']] %in% c('casa_gfed_161', 'SiB_calc',
                                          'SiB_mech', 'canibis_C4pctLRU'),
                 c('model', 'site_code', 'dd', 'dd_se_neff')]
fig3a_data <- droplevels(fig3a_data)
pdf('ECoast_gradient_with_std_error_norm.pdf')
data <- gradient_CI_plot(fig3a_data, norm=TRUE)
dev.off()
