# Common setup of packages and fonts
# Install libraries
update.packages("ggplot2")
install.packages("patchwork")
install.packages("ggh4x")
install.packages("writexl")
install.packages("extrafont")

library(extrafont)

# Manual step: put Arial TTF files in ~/fonts
font_import(paths = '~/fonts')
