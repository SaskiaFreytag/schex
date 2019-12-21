transform_radial <- function(data, xoff = 0, yoff = 0) {
  phi <- (data$x * 60 - 30)*(pi/180)
  Y <- (data$y + yoff) * cos(phi) - xoff * sin(60*pi/360)
  X <- (data$y + yoff) * sin(phi) + 0.5 + xoff * cos(60*pi/360)
  tibble(x = X, y = Y)
}

colourfan_grob <- function(nrow, ncol, nmunch = 10) {
  # the trick is that we first make square polygons and then transform coordinates
  dx <- 1 / ncol
  dy <- 1 / nrow
  
  # grid of base points
  x <- rep((0:(ncol-1))/ncol, nrow)
  y <- rep(((nrow-1):0)/nrow, each = ncol)
  
  # turn into polygon boundaries
  x <- unlist(lapply(x, function(x) c(x+dx*(0:nmunch)/nmunch, x+dx*(nmunch:0)/nmunch)))
  y <- unlist(lapply(y, function(y) c(rep(y, nmunch + 1), rep(y+dy, nmunch + 1))))
  
  # now transform coordinates and make polygon
  data <- transform_radial(tibble(x, y))
}


grid.newpage()
pushViewport(viewport(width=0.6, height=0.6,
     xscale=c(0, 1), yscale=c(0, 1)))
dat <- colourfan_grob(nrow=4, ncol=8)
id <- rep(1:(4*8), each = 22)
colours <- c("#421964", "#3d4080", "#326389",
             "#30818c", "#399e8b", "#57ba7e",
             "#92d267", "#dde15c", "#6b588e",
             "#6b588e", "#638c9f", "#638c9f",
             "#72b99b", "#72b99b", "#c4db7d",
             "#c4db7d", "#8f95b1", "#8f95b1",
             "#8f95b1", "#8f95b1", "#acd3a8",
             "#acd3a8", "#acd3a8", "#acd3a8",
             "#b7cac9", "#b7cac9", "#b7cac9",
             "#b7cac9", "#b7cac9", "#b7cac9",
             "#b7cac9", "#b7cac9")
grid.polygon(dat$x, dat$y, id, gp = gpar(fill = colours, 
        col = colours, lwd = 1, lty = 1))

label.x.pos <- dat[dat$x==min(dat$x),]
label.y.pos <- dat[dat$y==min(dat$y),]

label.x.pos <- transform_radial(tibble(x = seq(0,1,1/4), y = 1), yoff = 0.04)
label.y.pos <- transform_radial(tibble(x = 0, y = seq(0,1,1/4)), xoff = -0.07,
                                yoff=-0.05)

grid.text(c("0", "1", "2", "3", "4"), label.x.pos$x, label.x.pos$y)
grid.text(c("0%", "25%", "50%", "75%", "100%"), label.y.pos$x, label.y.pos$y)

vp<-viewport(x=0.2, y=0.5,
             width=0.2, height=0.2)
pushViewport(vp)
grid.text("uncertainty", 0.1, 0.1, rot=300)
upViewport()

vp<-viewport(x=0.5, y=1.2,
             width=0.2, height=0.2)
pushViewport(vp)
grid.text("expression", 0.5, 0.2)


dev.off()
