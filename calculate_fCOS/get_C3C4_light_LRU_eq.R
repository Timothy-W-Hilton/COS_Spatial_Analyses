library(lattice)
library(RColorBrewer)

##parse and pre-process data
fpath <- file.path(Sys.getenv('HOME'), 'work', 'Data',
                   'Stimler_COS_exchange_data.csv')
all_data <- read.csv(fpath)
all_data <- all_data[-1, c('PAR_umol.m.2s.1', 'lru', 'plant')]
names(all_data) <- c('PAR', 'LRU', 'species')
light_experiments <- grep('-light', as.character(all_data[['species']]))
light_data <- all_data[light_experiments, ]

light_data[['species']] <- as.factor(sub('-light', '',
                                         as.character(light_data[['species']])))
isC4 <- light_data[['species']] %in% c('corn', 'sorghum', 'amaranthus')
light_data[['C3C4']][isC4] <- 'C4'
light_data[['C3C4']][!isC4] <- 'C3'
light_data[['C3C4']] <- as.factor(light_data[['C3C4']])

mytheme <- standard.theme("pdf")
mytheme[['superpose.symbol']][['pch']] <- c(15,16,17,3,4)
mytheme[['superpose.symbol']][['col']] <-
    brewer.pal(length(levels(light_data[['species']])), 'Dark2')
xyplot(LRU~PAR|C3C4, data=light_data, groups=species,
       auto.key = TRUE,
       par.settings = mytheme)
