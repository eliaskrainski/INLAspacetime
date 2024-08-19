
library(sf)
library(INLAspacetime)
library(ggplot2)

gg0 <- ggplot() + theme_minimal()

w <- worldMap()

E <- Earth_poly(100)

g <- world_grid(300)


o <- st_difference(E, st_union(w))
    
gE <- world_grid(300, domain = w)
gO <- world_grid(300, domain = o)

ggplot() + theme_minimal() +
    geom_sf(data = w, fill = 'gray')

ggplot() + theme_minimal() +
    geom_sf(data = E, fill = rgb(0.5,0.9,1,0.5)) +
    geom_sf(data = w, fill = 'gray') 

ggplot() + theme_minimal() +
    geom_sf(data = o, fill = rgb(.3,.7,1,.5)) + 
    geom_sf(data = w, fill = rgb(1, .5, .3, 0.5))

ggplot() + theme_minimal() +
    geom_sf(data = gE, color = rgb(1.0, 0.7, 0.3, 0.5)) +
    geom_sf(data = gO, color = rgb(0.3, 0.7, 1.0, 0.5))
