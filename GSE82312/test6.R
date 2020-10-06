par(mfrow=c(1, 2), mar=rep(0.3, 4))
TernaryPlot(alab="Redder \u2192", blab="\u2190 Greener", clab="Bluer \u2192",
            lab.col=c('red', 'darkgreen', 'blue'),
            point='right', lab.cex=0.8, grid.minor.lines=0,
            grid.lty='solid', col=rgb(0.9, 0.9, 0.9), grid.col='white', 
            axis.col=rgb(0.6, 0.6, 0.6), ticks.col=rgb(0.6, 0.6, 0.6),
            padding=0.08)
data_points <- list(
  R = c(255, 0, 0), 
  O = c(240, 180, 52),
  Y = c(210, 222, 102),
  G = c(111, 222, 16),
  B = c(25, 160, 243),
  I = c(92, 12, 243),
  V = c(225, 24, 208)
)
AddToTernary(points, data_points, pch=21, cex=2.8, 
             bg=vapply(data_points, 
                       function (x) rgb(x[1], x[2], x[3], 128, maxColorValue=255),
                       character(1))
)

TernaryPlot(atip = "100% local", btip = "100% cis", ctip = "100% trans",
            lab.cex = 0.8,
            grid.lines = 0)
AddToTernary(points, input[, 7:9])
