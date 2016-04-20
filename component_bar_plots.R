library(lattice)
library(latticeExtra)
library(tidyr)
library(plotrix)  ## for std.error

parse_components <- function() {
    df <- read.csv('ocs_dd_renamed.csv')

    # remove boundaries sites
    df <- df[!(df[['sample_site_code']] %in%
             c('THD', 'EST', 'HIP', 'ACG', 'TGC')), ]

    df <- df[, c("sample_site_code",
                 "Kettle.Fsoil", "Hybrid.Fsoil",
                 "Anthropogenic..Kettle", "Anthropogenic..Zumkehr",
                 "CASA.GFED3..LRU.1.61", "CASA.GFED3..LRU.C3.C4",
                 "SiB..mechanistic.canopy", "SiB..prescribed.canopy",
                 "Can.IBIS..LRU.1.61", "Can.IBIS..LRU.C3.C4",
                 "climatological.boundaries")]

    ## use prettier name strings
    new_names <- c("site code",
                   "fSoil Kettle et al. (2002)", "fSoil Whelan et al. (2016)",
                   "fAnthro Kettle et al. (2002)", "fAnthro Campbell et al. (2015)",
                   "CASA-GFED3, LRU 1.61", "CASA-GFED3, LRU C3/C4",
                   "SiB, mechanistic", "SiB, LRU 1.61",
                   "Can-IBIS, LRU 1.61", "Can-IBIS, LRU C3/C4",
                   "clim bounds")
    names(df) <- new_names
    ## reorder columns to group Anthro, Soil, bounds, plant components together
    df <- df[, c("site code", "CASA-GFED3, LRU 1.61",
                 "CASA-GFED3, LRU C3/C4", "SiB, mechanistic", "SiB, LRU 1.61",
                 "Can-IBIS, LRU 1.61", "Can-IBIS, LRU C3/C4",
                 "fSoil Kettle et al. (2002)", "fSoil Whelan et al. (2016)",
                 "fAnthro Kettle et al. (2002)", "fAnthro Campbell et al. (2015)",
                 "clim bounds")]
    df[['const bounds']] <- 0

    return(df)
}

reformat_components <- function(df) {
    names0 <- names(df)
    dft <- gather_(df,
                   key_col='component',
                   value_col='drawdown',
                   gather_cols=names0[2:length(names0)],
                   factor_key=TRUE)
    dft[['type']] <- 'COS plant flux'
    dft[['type']][grepl('[sS]oil', dft[['component']])] <- 'COS soil flux'
    dft[['type']][grepl('[aA]nthro', dft[['component']])] <- 'COS anthropogenic flux'
    dft[['type']][grepl('[bB]ounds', dft[['component']])] <- 'boundary conditions'
    dft[['type']] <- as.factor(dft[['type']])

    dft[['component']] <- as.character(dft[['component']])
    dft[['component']] <- gsub('fSoil ', '', dft[['component']])
    dft[['component']] <- gsub('fAnthro ', '', dft[['component']])
    dft[['component']] <- gsub('clim bounds', 'climatological', dft[['component']])
    dft[['component']] <- gsub('const bounds', 'constant', dft[['component']])
    dft[['component']] <- as.factor(dft[['component']])

    dft <- dft[, c("site code", "component", "type", "drawdown")]
    return(dft)
}

class.subset <- function(df, class.name) {
    dfout <- df[df[['type']] == class.name, ]
    dfout <- droplevels(dfout)
    contains_kettle <- grepl('[Kk]ettle', levels(dfout[['component']]))
    if (any(contains_kettle)) {
        idx <- which(contains_kettle)
        if (idx == 1) {
            levels(dfout[['component']]) <- rev(levels(dfout[['component']]))
        }
    }
    return(dfout)
}

df0 <- parse_components()
df <- reformat_components(df0)

dd_mean <- aggregate(x=df['drawdown'],
                     by=df[, c('component', 'type')],
                     FUN=mean)
dd_se <- aggregate(x=df['drawdown'],
                   by=df[, c('component', 'type')],
                   FUN=std.error)
dd <- merge(dd_mean, dd_se, by=c('component', 'type'),
            suffixes=c('.mean', '.se'))
dd[['CI95']] <- dd[['drawdown.se']] * 1.96

## adapted from
## http://stackoverflow.com/questions/19975390/add-error-bars-seperately-in-lattice-line-plot
panel.errorbars <- function(x, y, lx, ux, ...) {
                        panel.xyplot(x, y, col='black', pch='x', cex=1.5, ...)
                        panel.arrows(x0 = lx, x1 = ux,
                                     y0 = y, y1 = y,
                                     col = "black",
                                     ends='both',
                                     angle=90,
                                     length=0.05,
                                     ...)
                    }
full_xlim <- c(with(dd, min(drawdown.mean - CI95)),
               with(dd, max(drawdown.mean + CI95)))

this_subset <- class.subset(dd, 'COS plant flux')
plt_plant <- xyplot(component ~ drawdown.mean | type,
                    data=this_subset,
                    ylab='',
                    scales = list(y=list(rot=0)),
                    xlim=c(-10, 60),
                    lx = this_subset$drawdown.mean - this_subset$CI95,
                    ux = this_subset$drawdown.mean + this_subset$CI95,
                    panel=panel.errorbars)
this_subset <- class.subset(dd, 'COS soil flux')
plt_soil <- xyplot(component ~ drawdown.mean | type,
                   data=this_subset,
                   ylab='',
                   scales = list(y=list(rot=0)),
                   xlim=c(-10, 60),
                   lx = this_subset$drawdown.mean - this_subset$CI95,
                   ux = this_subset$drawdown.mean + this_subset$CI95,
                   panel=panel.errorbars)
this_subset <- class.subset(dd, 'COS anthropogenic flux')
plt_anthro <- xyplot(component ~ drawdown.mean | type,
                     data=this_subset,
                     scales = list(y=list(rot=0)),
                     ylab='',
                     xlim=c(-10, 60),
                     lx = this_subset$drawdown.mean - this_subset$CI95,
                     ux = this_subset$drawdown.mean + this_subset$CI95,
                     panel=panel.errorbars)
this_subset <- class.subset(dd, 'boundary conditions')
plt_bounds <- xyplot(component ~ drawdown.mean | type,
                     data=this_subset,
                     xlab='COS drawdown (pptv)',
                     ylab='',
                     scales = list(y=list(rot=0)),
                     xlim=c(-10, 60),
                     lx = this_subset$drawdown.mean - this_subset$CI95,
                     ux = this_subset$drawdown.mean + this_subset$CI95,
                     panel=panel.errorbars)

pls <- c(plt_bounds, plt_anthro, plt_soil, plt_plant, layout=c(1, 4))

pdf('/tmp/model_components.pdf')
print(pls)
dev.off()
