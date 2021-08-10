labels <- c("very abrupt", "abrupt", "clear", "gradual", "diffuse")
codes <- c("V", "A", "C", "G", "D")
offsets = c(0.5, 2, 5, 15, 30) / 3
names(offsets) <- codes

gifski::save_gif(sapply(c(seq_along(offsets),
                          rev(seq_along(offsets))),
                        function(i, ylim = c(0, 1), xlim = c(-20, 20)) {
                          x <- rnorm(1e4, mean = 0, sd = offsets[i])
                          
                          dens <- density(x)
                          
                          plot(dens, xaxt = 'n', bty = 'n', 
                               ylim = ylim, xlim = xlim, 
                               cex.axis = 2, cex.main = 2,
                               main =  paste0("Half Boundary Thickness, cm\n",
                                              labels[i]))
                          
                          axis(1, pos = 0, cex.axis=1.5)
                          
                          polygon(c(-3*offsets[i],
                                    dens$x[abs(dens$x) < 3*offsets[i]], 
                                    3*offsets[i]), 
                                  c(0, 
                                    dens$y[abs(dens$x) < 3*offsets[i]], 
                                    0), 
                                  col = "forestgreen", lty = 2, lwd = 2)
                          
                          abline(v = 0, lty = 3, lwd = 2, ylim = c(0, 1))
                        }), "boundaries.gif")
