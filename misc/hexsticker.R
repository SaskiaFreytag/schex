
library(fungible)
library(hexSticker)
library(hexbin)
library(ggplot2)
library(RColorBrewer)
library(gameofthrones)


da <- monte(seed = 153, nvar = 2, nclus = 3, clus.size = c(500, 500, 500),
      eta2 = c(0.919, 0.812, 0.761), cor.list = NULL,
      random.cor = FALSE, skew.list = NULL, kurt.list = NULL,
      secor = NULL, compactness = NULL, sortMeans = TRUE)

set.seed(51)
dat <- data.frame(
  x = da$data[,2],
  y = da$data[,3],
  col = as.factor(sample(1:5, 1500, replace=TRUE)))

ggplot(dat, aes(x=x, y=y, color=col)) + geom_point() + 
  theme_void()

hrdat <- hexbin(dat$x, dat$y, xbins=50, IDs=TRUE)

cID <- hrdat@cID

hrdat1 <- data.frame(
  x = hcell2xy(hrdat)$x,
  y = hcell2xy(hrdat)$y,
  col = factor(tapply(dat$col, cID, function(z) names(sort(table(z), decreasing=T))[1]))
)


cols <- brewer.pal(5, "Set2")
gg <- ggplot(hrdat1, aes(x=x, y=y, fill=col)) + geom_hex() + 
  theme_void() + guides(fill=FALSE) + scale_fill_got(discrete = TRUE, option  ="Daenerys")

font_add_google('Stylish', 'sans-serif')

sticker(gg, package="schex", p_size=8, s_x=1, s_y=.75, s_width=1.3, s_height=1,
        filename="~/Desktop/schex_hex.png", h_fill="#f1f1f1", h_color="#d2d2d2",
        p_color="#d2d2d2", p_family = "sans-serif")
