library(ggsci)
library(scales)

test_colour1 <- ggsci::pal_uchicago(alpha = 0.3)(9)
test_colour2 <- ggsci::pal_uchicago(alpha = 1)(9)
library(scales)

show_col(colours = c(test_colour1), borders = NA)
show_col(colours = c(test_colour2), borders = NA)

test_colour1
test_colour2



test_colour3 <- ggsci::pal_futurama(alpha = 1)(12)
show_col(colours = c(test_colour3), borders = NA)
test_colour3


test_colour4 <- ggsci::pal_futurama(alpha = 0.7)(12)
show_col(colours = c(test_colour4), borders = NA)
test_colour4


