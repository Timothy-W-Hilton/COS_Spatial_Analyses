library(plotrix)
library(RColorBrewer)
library(tidyr)  # could also use reshape2
library(boot)

merge_obs <- function(dfboot) {
    obs <- read.csv('./ocs_dd_renamed.csv')
    obs <- obs[, c('sample_site_code', 'NOAA.obs')]
    names(obs) <- c('site', 'obs')
    result <- merge(dfboot, obs)
    return(result)
}


##' human-readable names for COS Fplant models
##'
##' The data frame column names are more machine-oriented: no spaces, caps, etc.  These are nicer-looking strings for e.g. plot labels.
##' @title
##' @return list of strings. The data frame column labels are the list
##' names and the human-readble strings are the list elements.
##' @author Timothy W. Hilton
##' @export
human_readable_model_names <- function() {
    return(list(canibis_161='Can-IBIS, LRU=1.61',
                canibis_C4pctLRU= 'Can-IBIS, LRU=C3/C4',
                casa_gfed_161='CASA=GFED3, LRU=1.61',
                casa_gfed_C4pctLRU='CASA=GFED3, LRU=C3/C4',
                SiB_calc='SiB, LRU=1.61',
                SiB_mech='SiB, mechanistic'))
}

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
    num_idx = unlist(lapply(dd, is.numeric))
    norm_data <- dd[norm_site, num_idx]
    for (i in seq(1, nrow(dd))) {
        dd[i, num_idx] <- dd[i, num_idx] / norm_data
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

##' .. produce a scatter plot with error bars of STEM model component
##' drawdowns at NOAA sites.
##'
##' .. content for \details{} ..
##' @title
##' @param df (data frame): NOAA sites in rows, model component
##' drawdowns (pptv) in columns
##' @param dd_col (string): name of the column in df containing
##' drawdown values
##' @param ci_hi_col (string): name of the column in df containing
##' upper confidence interval widths
##' @param ci_lo_col (string): name of the column in df containing
##' lower confidence interval widths
##' @param t_str (string): Main title string for the plot
##' @param site_names (vector of strings): Three letter site codes;
##' the gradient will be plotted in this order from left to right.
##' @param norm_site (string): row label of a site to normalize the
##' data against.  If unspecified (default), no normalization is
##' performed.
##' @return
##' @author Timothy W. Hilton
##' @export
gradient_CI_plot <- function(df,
                             dd_col='dd',
                             ci_hi_col='ci_hi',
                             ci_lo_col='ci_lo',
                             t_str='gradient plot',
                             site_names=list(),
                             norm_site='') {

    n_sites <- length(site_names)
    models <- rev(sort(unique(df[['Fplant']])))
    n_models <- length(models)
    pal <- brewer.pal(n_models, 'Paired')
    marker_sequence <- seq(0, n_models - 1)
    x_offset <- calculate_hoffset(n_models, 0.075)
    ylab_str <- 'Drawdown  (pptv)'

    if (nchar(norm_site) > 0) {
        ylab_str <- paste("Drawdown normalized to", norm_site)
        df_norm <- by(df, df[['Fplant']], function(x) {
            orig_row_names <- row.names(x)
            row.names(x) <- x[['site']]
            x <- row_normalizer(x, 'NHA')
            row.names(x) <- orig_row_names
            return(x)})
        df_norm <- do.call(rbind, df_norm)
        df <- df_norm
    }
    ylim <- range(df[df[['site']] %in% site_names, c('ci_lo', 'ci_hi', 'obs')])

    idx = (df[['site']] %in% site_names) & (df[['Fplant']]==models[[1]])
    this_df <- df[idx, ]
    row.names(this_df) <- this_df[['site']]
    this_df <- this_df[site_names, ]

    with(this_df,
         plotCI(1:n_sites + x_offset[[1]],
                dd, uiw=(ci_hi - dd), liw=(dd - ci_lo),
                xaxt='n',
                main=t_str,
                ylab=ylab_str,
                xlab=NA,
                col=pal[[1]],
                ylim=ylim,
                xlim=c(1 + min(x_offset), n_sites * 1.6),
                pch=marker_sequence[[1]]))
    axis(1, at=1:n_sites, labels=site_names)
    points(1:n_sites, this_df[['obs']], pch=8)

    for (i in 2:n_models) {
        cat(paste('plotting models[', models[[i]], ']\n'))
        idx = (df[['site']] %in% site_names) & (df[['Fplant']]==models[[i]])
        this_df <- df[idx, ]
        row.names(this_df) <- this_df[['site']]
        this_df <- this_df[site_names, ]
        with(this_df,
             plotCI(x=1:n_sites + x_offset[[i]],
                    y=dd, uiw=ci_hi - dd, liw=dd - ci_lo,
                    add=TRUE,
                    col=pal[[i]],
                    pch=marker_sequence[[i]]))
    }
    mod_strs <- unlist(human_readable_model_names()[models])
    mod_strs <- c('observed', mod_strs)
    marker_sequence <- c(8, marker_sequence)
    pal <- c('#000000', pal)
    legend(x='right', legend=mod_strs, pch=marker_sequence, col=pal, cex=0.7)
}

myboot <- function(x) {
    boot(x[['dd']],
         statistic=function(data, ind) return(mean(data[ind])),
         R=5000)
}

df <- read.csv('./model_components_25Feb.csv')
df[['Fbounds']] <- 'CONST'
df[['Fbounds']][grepl('climatological', df[['model']])] <- 'CLIM'
components <- strsplit(x=as.character(df[['model']]), split='-')
df[['Fplant']] <- unlist(lapply(components, function(x) x[[1]]))
df[['Fsoil']] <- unlist(lapply(components, function(x) x[[2]]))
df[['Fanthro']] <- unlist(lapply(components, function(x) x[[3]]))

dfl <- split(df[, c('site_code', 'Fplant', 'Fsoil', 'Fanthro',
                    'Fbounds', 'dd')],
             f=df[, c('site_code', 'Fplant')], drop=TRUE)

boot_results <- lapply(dfl,
                       FUN=myboot)
boot_ci_results <- lapply(boot_results, boot.ci,
                          type=c("norm","basic", "perc", "bca"))
dfboot <- data.frame(row.names=names(boot_results),
                     dd=rep(NA, length(boot_results)),
                     ci_lo=rep(NA, length(boot_results)),
                     ci_hi=rep(NA, length(boot_results)),
                     Fplant=rep(NA, length(boot_results)),
                     site=rep(NA, length(boot_results)))
for (this_set in names(boot_results)) {
    dfboot[[this_set, 'dd']] <- boot_results[[this_set]][['t0']]
    dfboot[[this_set, 'ci_lo']] <- boot_ci_results[[this_set]][['basic']][[4]]
    dfboot[[this_set, 'ci_hi']] <- boot_ci_results[[this_set]][['basic']][[5]]
    dfboot[[this_set, 'site']] <- unlist(strsplit(this_set, '\\.'))[[1]]
    dfboot[[this_set, 'Fplant']] <- unlist(strsplit(this_set, '\\.'))[[2]]
}

dfboot_orig <- dfboot
dfboot <- merge_obs(dfboot)

if (TRUE) {

    norm_site <- ''
    plot_width <- 6 # inches
    plot_height <- 10 # inches
    if (nchar(norm_site) == 0){
        pdf(file='gradients_bootstrapCIs.pdf',
            width=plot_width, height=plot_height)
        cat('writing gradients_bootstrapCIs.pdf\n')
    } else {
        pdf(file='gradients_bootstrapCIs_norm.pdf',
            width=plot_width, height=plot_height)
        cat('writing gradients_bootstrapCIs_norm.pdf\n')
    }
    oldpar<-panes(matrix(1:3,nrow=3,byrow=TRUE))
    # par(mar=c(bottom, left, top, right)’)
    par(mar=c(2,5,1.6,1))
    gradient_CI_plot(dfboot, t_str='East Coast',
                     site_names=c('NHA', 'CMA', 'SCA'),
                     norm_site=norm_site)

    gradient_CI_plot(dfboot, t_str='mid-continent West to East (Dry to Wet)',
                     site_names=c('CAR', 'WBI', 'AAO', 'HIL', 'CMA'),
                     norm_site=norm_site)

    gradient_CI_plot(dfboot, t_str='mid-continent North to South',
                     site_names=c('ETL', 'DND', 'LEF', 'WBI', 'BNE', 'SGP', 'TGC'),
                     norm_site=norm_site)
    dev.off()
    par(oldpar)

    ## write.csv(df[, c("site_code", "longitude",  "latitude",
    ##                  "dd", "Fbounds", "Fplant", "Fsoil", "Fanthro")],
    ##           file='model_components_26Feb.csv')
}
