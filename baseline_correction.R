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
       main = spectra[i])
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
                           # (wavenumber > 145 & wavenumber < 180) |
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