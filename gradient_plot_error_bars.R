library(plotrix)
library(RColorBrewer)
library(tidyr)  # could also use reshape2
library(boot)

##' offsets are calculated from an arbitrary center
##'
##' @title calculate horizontal offsets for N points with specified
##' interval
##' @param npoints (int): number of points in the sequence
##' @param gapwidth (float): space between consecutive points
##' @return numeric array containing npoints horizontal offset values
##' @author Timothy W. Hilton
##' @export
calculate_hoffset <- function(npoints, gapwidth) {
    width_total <- gapwidth * (npoints - 1)
    offsets <- seq(from=-(width_total / 2.0),
                   to=(width_total / 2.0),
                   length.out=npoints)
    return(offsets)
}
##' Normalize all rows in a data frame to the row with label
##' norm_site.  Helper function for ci_normalizer -- should not be
##' called by user.
##'
##'
##' @title normalize data frame rows (helper function)
##' @param dd (data frame): data frame containing drawdown
##' observations, labeled by site code (rows) and model (columns)
##' @param norm_site (string): site code (also the row label) of the
##' site to normalize against.
##' @return dd, normalized to site_code
##' @author Timothy W. Hilton
##' @export
row_normalizer <- function(dd, norm_site='NHA') {
    norm_data <- dd[norm_site, ]
    for (i in seq(1, nrow(dd))) {
        dd[i, ] <- dd[i, ] / norm_data
    }
    return(dd)
}
##' Normalization is calculated relative to the specified site.  The
##' normalized confidence intervals are calculated by finding the
##' large span within dd and ci, and nomalizing against the specified
##' site.
##'
##' @title Normalize a data frame containing confidence intervals.
##' @param dd (data frame): data frame containing drawdown
##' observations, labeled by site code (rows) and model (columns)
##' @param ci (data frame): data frame containing drawdown confidence
##' intervals, labeled by site code (rows) and model (columns)
##' @param norm_site (string): site code (also the row label) of the
##' site to normalize against.
##' @return list labeled "dd", "ci_hi", and "ci_lo".  Each element
##' contains a data frame labeled by site code (rows) and model
##' (columns) containing normalized drawdowns, upper confidence
##' intervals, and lower confidence intervals, respectively.
##' @author Timothy W. Hilton
##' @export
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

##' return a dummy set of random "observations" and "confidence
##' intervals" with the same row and column labels that the real data
##' will have.  Useful for testing the plotting code and whether
##' plotrix is able to produce the plot I want.
##'
##' @title generate a set of random "observations" and "confidence
##' intervals"
##' @return list with labels "dd", "ci".  Each element contains a data
##' frame labeled by site code (rows) and model (columns) containing
##' normalized drawdowns and confidence intervals, respectively.
##' @author Timothy W. Hilton
##' @export
dummy_data <- function() {
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
gradient_CI_plot <- function(df, dd_col='dd', se_col='dd_se_neff', norm=FALSE,
                             site_names=list()) {
    n_sites <- nlevels(df[['site_code']])
    n_models <- nlevels(df[['model']])
    models <- levels(df[['model']])
    pal <- brewer.pal(n_models, 'Dark2')
    marker_sequence <- c(0, 1, 3, 4, 5)

    if (length(site_names) == 0) {
        site_names <- levels(df[['site_code']])
    }

    ## dfw_dd: "data frame, wide; drawdown"
    dfw_dd <- spread(df[, c('model', 'site_code', 'dd')],
                     model, dd)
    row.names(dfw_dd) <- dfw_dd[['site_code']]
    dfw_dd <- dfw_dd[, !names(dfw_dd) %in% "site_code"]
    ## make sure sites are in the requested order
    dfw_dd <- dfw_dd[site_names, ]

    dfw_se <- spread(df[, c('model', 'site_code', 'dd_se_neff')],
                     model, dd_se_neff)
    row.names(dfw_se) <- dfw_se[['site_code']]
    dfw_se <- dfw_se[, !names(dfw_se) %in% "site_code"]
    ## make sure sites are in the requested order
    dfw_se <- dfw_se[site_names, ]
    ## dfw_dd: "data frame, wide; confidence interval"
    dfw_ci <- dfw_se * 1.96 ## 95% confidence interval

    if (norm) {
        dd_list <- ci_normalizer(dfw_dd, dfw_ci, 'NHA')
        t_str <- "normalized drawdown with 95% confidence intervals"
    } else {
        dd_list <- list(dd=dfw_dd, ci_hi=dfw_ci, ci_lo=dfw_ci)
        t_str <- "drawdown with 95% confidence intervals"
    }

    x_offset <- calculate_hoffset(n_models, 0.01)
    plotCI(1:n_sites + x_offset[[1]],
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
           xlim=c(1 + min(x_offset), n_sites + max(x_offset)),
           pch=marker_sequence[[1]])
    ## points(x, y, col=pal[[1]], type='l')
    axis(1, at=1:n_sites, labels=site_names)

    mapply(CI_plotter, dd_list[['dd']], dd_list[['ci_hi']], dd_list[['ci_lo']],
           pal, x_offset, marker_sequence[1:n_models])

    ## legend(x='right', legend=models, pch=marker_sequence, col=pal)
    return(list(dd=dfw_dd, ci=dfw_ci))
}

myboot <- function(x) {
    boot(x[['dd']],
         statistic=function(data, ind) return(mean(data[ind])),
         R=5000)
}

df <- read.csv('./model_components_25Feb.csv')
df[['clim_bounds']] <- 'CONST'
df[['clim_bounds']][grepl('climatological', df[['model']])] <- 'CLIM'
df[['clim_bounds']] <- as.factor(df[['clim_bounds']])

dfl <- split(df, f=df[, c('site_code', 'clim_bounds')], drop=TRUE)
## dfl <- lapply(dfl, function(x) x[sort(x[['model']]), ])
## bar <- cbind(dfl[['SCA.CLIM']][, c('model', 'dd')],
##              dfl[['SCA.CONST']][c('model', 'dd')])

boot_results <- lapply(dfl[c('NHA.CONST', 'NHA.CLIM',
                             'CMA.CONST', 'CMA.CLIM',
                             'SCA.CONST', 'SCA.CLIM')],
                       FUN=myboot)
boot_ci_results <- lapply(boot_results, boot.ci,
                          type=c("norm","basic", "perc", "bca"))
dfboot <- data.frame(row.names=names(boot_results),
                     dd=rep(NA, length(boot_results)),
                     ci_lo=rep(NA, length(boot_results)),
                     ci_hi=rep(NA, length(boot_results)))
for (this_set in names(boot_results)) {
    dfboot[[this_set, 'dd']] <- boot_results[[this_set]][['t0']]
    dfboot[[this_set, 'ci_lo']] <- boot_ci_results[[this_set]][['basic']][[4]]
    dfboot[[this_set, 'ci_hi']] <- boot_ci_results[[this_set]][['basic']][[5]]
}


fig3a_data <- df[df[['site_code']] %in% c('NHA', 'CMA', 'SCA'),
                 c('model', 'site_code', 'dd', 'dd_se_neff')]
fig3a_data <- droplevels(fig3a_data)
data <- gradient_CI_plot(fig3a_data, norm=FALSE, site_names=c('NHA', 'CMA', 'SCA'))
