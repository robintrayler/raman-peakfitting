## load required packages------------------------------------------------------
library(tidyverse)
## 128 - one peak
## 206 two peaks
## 464 two peaks
## subsetting OK 

## read in data----------------------------------------------------------------
data <- read_csv('~/Desktop/spectra.csv', col_names = T)
spectra <- colnames(data)[-1]
pdf(width = 10, height = 10)
## loop through the data-------------------------------------------------------
for(i in 1:length(spectra)){
  ## make plots of backgrounds-------------------------------------------------
  plot(x = data[['wavenumber']], 
       y = data[[spectra[i]]], 
       type = 'l',
       lwd = 3,
       col = 'grey',
       xlab = 'wavenumber',
       ylab = '',
       ylim = c(0, max(data[spectra[i]])),
       main = spectra[i],
       yaxt = 'n')
  abline(v = c(128, 206, 465), 
         col = rgb(0,0,0,0.25), lwd = 10)
  ## smooth the raw data with a spline-----------------------------------------
  data[[spectra[i]]] <- smooth.spline(x = data[['wavenumber']], 
                                      y = data[[spectra[i]]], 
                                      spar = 0.2)$y
  lines(x = data[['wavenumber']],
        y = data[[spectra[i]]], 
        col = 'tomato',
        lwd = 2,
        lty = 3)
  
  ## select backgrounds--------------------------------------------------------
  bgs <- data %>% filter((wavenumber > 100 & wavenumber < 115) | 
                           (wavenumber > 145 & wavenumber < 147) |
                           # (wavenumber > 180 & wavenumber < 185) |
                           (wavenumber > 260 & wavenumber < 300) |
                           (wavenumber > 425 & wavenumber < 445) |
                           # (wavenumber > 485 & wavenumber < 487) |
                           (wavenumber > 520 & wavenumber < 535) |
                           (wavenumber > 575 & wavenumber < 600))
  
  points(bgs[['wavenumber']], bgs[[spectra[i]]], pch = 19, cex = 0.5)
  
  ## fit a spline to backgrounds-----------------------------------------------
  smoothby = 0.8
  bg.spline <- smooth.spline(x = bgs[['wavenumber']], 
                             y = bgs[[spectra[i]]], 
                             spar = smoothby)
  lines(bg.spline, lwd = 2, lty = 2, col = 'green')
  
  ## remove the backgrounds----------------------------------------------------
  f <- approxfun(x = bg.spline$x, y = bg.spline$y) # interpolation function
  corrected <- data[[spectra[i]]] - f(data[['wavenumber']]) # subtract the spline from the data
  corrected[corrected < 0] = 0 # get rid of negitive values 
  lines(data$wavenumber, corrected, lwd = 2, col = 'steelblue') # add a line to the plot
  data[[spectra[i]]] <- corrected
  
  ## add a legend--------------------------------------------------------------
  legend('topleft',
         legend = c('raw', 
                    'smoothed', 
                    'baseline points', 
                    'baseline spline', 
                    'baseline subtracted', 
                    'peak of interest'),
         lwd = c(2,2,NA,2,2, 10),
         lty = c(1,3,NA,2,1, 1),
         pch = c(NA,NA,19,NA,NA, NA),
         col = c('grey', 
                 'tomato', 
                 'black', 
                 'black', 
                 'steelblue', 
                 rgb(0,0,0,.25)),
         cex = .8)
}
dev.off()

## subset 128 peak-------------------------------------------------------------
P128 <- data %>% filter(wavenumber > 106 & wavenumber < 150)
plot(x = P128[['wavenumber']], 
     y = P128[['X117']], 
     type = 'o',
     pch = 19, 
     col = 'grey')
source('~/Dropbox/Documents/Rstudio/Peak.Fit.R')
fit <- Peak.Fit(D = data.frame(x = P128[['wavenumber']], y = P128[['X117']]),
                centers = 128,
                widths = 5,
                scale = 200, 
                type = 'G')
lines(x = P128[['wavenumber']], 
      y = dnorm(x = P128[['wavenumber']], fit$centers, fit$widths) * fit$scales)

## MCMC peakfitting------------------------------------------------------------
npeaks = 1
itterations = 100000

## function for summing gaussians----------------------------------------------
gaussSum <- function(x, center, width, scale){
  temp <- matrix(ncol = length(center), nrow = length(x))
  for(i in 1:length(center)){
    temp[, i] <- dnorm(x, mean = center, sd = width) * scale
  }
  return(apply(X = temp, MARGIN = 1, FUN = sum))
}

## preallocate some storage----------------------------------------------------
centerStore <- matrix(ncol = npeaks, nrow = itterations)
widthStore <- matrix(ncol = npeaks, nrow = itterations)
scaleStore <- matrix(ncol = npeaks, nrow = itterations)
sdStore <- vector(length = itterations)
## make some initial guesses---------------------------------------------------
y <- P128[['X118']]
centerCurrent <- 129.5
widthCurrent <- 3
scaleCurrent <- max(y)*widthCurrent*2
sdCurrent <- rexp(1, rate = 0.1)
## sd prior--------------------------------------------------------------------
# sdx <- seq(0, 20, by = 0.01)
# p <- dexp(sdx, rate = 0.1)
# plot(sdx, p)
## propose new values----------------------------------------------------------

for(i in 1:itterations){
  centerProposed <- rnorm(1, centerCurrent, 2)
  widthProposed <- abs(rnorm(1, widthCurrent, 0.1))
  scaleProposed <- abs(rnorm(1, scaleCurrent, 5))
  # scaleProposed <- 160
  # sdProposed <- rexp(1, rate = 0.05)
  fitCurrent <- gaussSum(x = P128[['wavenumber']], 
                         center = centerCurrent, 
                         width = widthCurrent, 
                         scale = scaleCurrent)
  
  fitProposed <- gaussSum(x = P128[['wavenumber']], 
                          center = centerProposed, 
                          width = widthProposed, 
                          scale = scaleProposed)
  
  probProposed <- sum(log(dnorm(y, mean = fitProposed, sd = 2)))
  probProposed <- ifelse(is.infinite(probProposed), -100000, probProposed)
  
  probCurrent <- sum(log(dnorm(y, mean = fitCurrent, sd = 2)))
  probCurrent <- ifelse(is.infinite(probCurrent), -100000, probCurrent)
  
  if(exp(probProposed - probCurrent) > runif(1)){
    centerCurrent <- centerProposed
    widthCurrent <- widthProposed
    scaleCurrent <- scaleProposed
    # sdCurrent <- sdProposed 
  }
  centerStore[i, ] <- centerCurrent
  widthStore[i, ] <- widthCurrent
  scaleStore[i, ] <- scaleCurrent
  # sdStore[i] <- sdCurrent
  
}
layout(1)
plot(x = P128[['wavenumber']], 
     y = y, 
     type = 'n',
     pch = 19, 
     col = 'grey')
for(i in seq(2000, itterations, 200)){
  lines(P128[['wavenumber']], gaussSum(P128[['wavenumber']], 
                                       center = centerStore[i],
                                       width = widthStore[i],
                                       scale = scaleStore[i]),
        col = rgb(0,0,1,0.05))
}

lines(x = P128[['wavenumber']], 
    y = y, 
    type = 'o',
    pch = 19, 
    col = 'black')

# layout(c(1,2,3))
hist(centerStore[2000:itterations],breaks = 50)
hist(widthStore[2000:itterations], breaks = 50)
hist(scaleStore[2000:itterations], breaks = 50)

