# Plot graph for stairwayplot2
rm(list=ls())
library(pracma)

# Function for minor ticks on x-axis
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
  
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  
  labels <- sapply(major.ticks,function(i)
    as.expression(bquote(10^ .(i)))
  )
  axis(ax,at=major.ticks,labels=labels,...)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  
  
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}

# Read in output file from stairwayplot2
full <- read.table("albo.final.summary", header=TRUE)
full$year <- log10(full$year)
full$Ne_median <- log10(full$Ne_median)
full$Ne_2.5. <- log10(full$Ne_2.5.) #numbers refer to confidence intervals
full$Ne_97.5. <- log10(full$Ne_97.5.)
full$Ne_12.5. <- log10(full$Ne_12.5.)
full$Ne_87.5. <- log10(full$Ne_87.5.)

pdf("albo_stairway_plot.pdf", width=8, height=6)
plot(u$year, u$Ne_median, xaxt="n", yaxt="n", xlim = c(0, 6), ylim = c(2.5, 8), col="white")
minor.ticks.axis(1,9,mn=0,mx=6) #x axis, lower limit, upper limit
minor.ticks.axis(2,9,mn=0,mx=8) #y axis, lower limit, upper limit

lines(full$year, full$Ne_median, pch=16, lwd=3)
lines(full$year, full$Ne_2.5., lty=3)
lines(full$year, full$Ne_97.5., lty=3)
lines(full$year, full$Ne_12.5., lty=5)
lines(full$year, full$Ne_87.5., lty=5)

dev.off()