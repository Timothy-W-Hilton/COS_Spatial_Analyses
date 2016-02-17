library(plotrix)
library(RColorBrewer)



models <- c('CASA-GFED3', 'Can-IBIS', 'SiB mech')
site_names <- c('NHA', 'CMA', 'SCA', 'WBI', 'THD')
marker_sequence <- c('o', 'x', 5, '+', '*')
marker_sequence <- c(0, 1, 3, 4, 5)

n_sites <- length(site_names)
x <- 1:n_sites
y <- runif(n_sites)
err <- runif(n_sites)

pal <- brewer.pal(n_sites, 'Dark2')

plotCI(x, y, err, xaxt='n',
       main="Basic plotCI",  ylab='Drawdown',  xlab='site',  col=pal[[1]])
## points(x, y, col=pal[[1]], type='l')
axis(1, at=1:5, labels=site_names)

x <- 1:n_sites - (0.1)
y < -runif(n_sites)
err <- runif(n_sites)
plotCI(x, y, err, col=pal[[2]], add=TRUE, pch='x')
## points(x, y, col=pal[[2]], type='l')

legend(x='right', legend=models, pch=marker_sequence, col=pal)
