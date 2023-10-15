q()

load("vcfppR.RData")
library(hexSticker)
library(ggplot2)

library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Gochi Hand", "gochi")
## Automatically use showtext to render text for future devices
showtext_auto()

names(o)

d <- data.frame(x = unlist(o), 
                counts = rep(names(o),times = sapply(o,length)))

p <- ggplot(d, aes(x = counts, y = x)) + xlab("") + ylab("") +
  geom_boxplot()
p <- p + theme_classic()
p <- p + theme_transparent()

sticker(p, package="vcfppR", 
        p_size=20, s_x=0.9, s_y=.7, s_width=1.3, s_height=1.1,
        filename="vcfppR.png", h_fill = "red", h_color = "orange")

dev.off()
