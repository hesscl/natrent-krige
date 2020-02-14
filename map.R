#### Plot the Kriging Grid ----------------------------------------------------

library(ggplot2)
library(patchwork)

theme_krige <- function(...){
  theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(.8, "inch"),
          legend.key.height = unit(.1, "inch"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

#plot predicted values and standard errors
g1 <- ggplot() +
  geom_hex(data = GridBAUsPred, 
           aes(x = lng, y = lat, fill = predicted),
           stat = "identity", color = NA) +
  scale_fill_distiller(palette = "Spectral",
                       name = paste0("Predicted\n", bed_size, "B Rent"),
                       labels = scales::dollar) +
  geom_density_2d(data = as_tibble(query@coords), 
                  aes(x = lng, y = lat),
                  colour = "grey80", alpha = .5) + 
  theme_krige() +
  xlab("") + ylab("")

#Similar to above but with s.e.
g2 <- ggplot() + 
  geom_hex(data= GridBAUsPred,
           aes(x = lng, y = lat, fill = sqrt(variance)),
           stat = "identity") +
  scale_fill_distiller(palette = "BrBG",
                       name = "SE",
                       labels = scales::dollar,
                       guide = guide_legend(title = "SE")) +
  theme_krige() +
  theme(legend.key.size = unit(.2, "inch")) +
  xlab("") + ylab("")

if(!dir.exists(file.path(paste0("./output/maps/", metro_nickname)))){
  dir.create(file.path(paste0("./output/maps/", metro_nickname)))
}

g1 + g2 +
  ggsave(file.path(paste0("./output/maps/", metro_nickname, "/kriged", bed_size, "B.pdf")),
         width = 12, height = 8, dpi = 300)