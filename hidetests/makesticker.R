## 25 Mar 2025
## Make package sticker
library(hexSticker)
library(magick)

sticker("/Users/ssu/Documents/Main/GitHub/ratesci/man/figures/MNplot.png",
        package="ratesci",
        p_size=40,
        p_x = 0.9,
        p_y = 1.5,
        p_color = "#FFF8B4",
        s_x=1.05,
        s_y=0.95,
        s_width=1,
        filename="man/figures/imgfile.png",
        white_around_sticker = TRUE,
        h_color = "#F6AA47",
        dpi = 600)

#  fuzz <- 50
fuzz <- 20
p <- image_read("man/figures/imgfile.png")
pp <- p %>%
  image_fill(color = "transparent", refcolor = "white", fuzz = fuzz, point = "+1+1") %>%
  image_fill(color = "transparent", refcolor = "white", fuzz = fuzz, point = paste0("+", image_info(p)$width-1, "+1")) %>%
  image_fill(color = "transparent", refcolor = "white", fuzz = fuzz, point = paste0("+1", "+", image_info(p)$height-1)) %>%
  image_fill(color = "transparent", refcolor = "white", fuzz = fuzz, point = paste0("+", image_info(p)$width-1, "+", image_info(p)$height-1))
image_write(image = pp, path = "man/figures/hex_transp.png")

use_logo("man/figures/hex_transp.png")

