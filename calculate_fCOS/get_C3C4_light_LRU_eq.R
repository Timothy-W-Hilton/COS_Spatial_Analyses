library(lattice)
library(RColorBrewer)

##parse and pre-process data
fpath <- file.path(Sys.getenv('HOME'), 'work', 'Data',
                   'Stimler_COS_exchange_data.csv')
all_data <- read.csv(fpath)
light_data <- all_data[-1, c('PAR_umol.m.2s.1', 'lru', 'plant')]
names(light_data) <- c('PAR', 'LRU', 'species')
light_experiments <- grep('-light', as.character(light_data[['species']]))
light_data <- light_data[light_experiments, ]
light_data[['species']] <- as.factor(sub('-light', '',
                                         as.character(light_data[['species']])))
isC4 <- light_data[['species']] %in% c('corn', 'sorghum', 'amaranthus')
light_data[['C3C4']][isC4] <- 'C4'
light_data[['C3C4']][!isC4] <- 'C3'
light_data[['C3C4']] <- as.factor(light_data[['C3C4']])


pal <- brewer.pal(length(levels(light_data[['species']])), 'Dark2')
mytheme <- standard.theme("pdf")
mytheme[['superpose.symbol']][['pch']] <- c(15,16,17,3,4)
mytheme[['superpose.symbol']][['col']] <- pal

xy <- xyplot(LRU~PAR|C3C4, data=light_data, groups=species,
             auto.key=TRUE,
             col.line=pal,
             par.settings=mytheme,
             panel.groups=function(x, y, ...){
                 panel.xyplot(x, y,  ...)
                 panel.loess(x, y, ...)
             },
             panel=function(x,y,...) {
                 ## custom panel function to add an overall loess line
                 panel.superpose(x,y,...)
                 panel.loess(x, y, lty='dashed', col='black', lwd=2)
             }
         )
pdf(file=file.path(Sys.getenv('HOME'), 'work', 'Plots', 'SpatialAnalysisPaper',
        'stimler_PAR_LRU.pdf'))
print(xy)
dev.off()

C3_lo <- loess(LRU~PAR, data=light_data[light_data[['C3C4']] == 'C3', ])
C4_lo <- loess(LRU~PAR, data=light_data[light_data[['C3C4']] == 'C4', ])

C3_pred <- predict(C3_lo, newdata=seq(0, max(light_data[['PAR']], na.rm=TRUE)))
C4_pred <- predict(C4_lo, newdata=seq(0, max(light_data[['PAR']], na.rm=TRUE)))


