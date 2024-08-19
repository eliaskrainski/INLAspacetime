library(INLA)
library(inlabru)
library(INLAspacetime)
library(sf)
library(ggplot2)
library(bench)

s <- 30
pl <- cbind(c(0,1,1,0,0), c(0,0,1,1,0)) * s
pol <- st_sfc(st_polygon(list(pl)))

b1 <- cbind(c(1/7, 1/7, 2/7, 5/7, 6/7, 6/7, 5/7, 2/7, 1/7),
            c(3/7, 2/7, 1/7, 1/7, 2/7, 3/7, 2/7, 2/7, 3/7)) * s
b2 <- cbind(c(3/5, 4/5, 4/5, 3/5, 3/5),
            c(3/5, 3/5, 4/5, 4/5, 3/5)) * s
b3 <- cbind(c(1/5, 2/5, 2/5, 1/5, 1/5),
            c(3/5, 3/5, 4/5, 4/5, 3/5)) * s
bar <- st_sfc(st_multipolygon(list(
    st_polygon(list(b1)),
    st_polygon(list(b2)),
    st_polygon(list(b3)))))

ggplot() + theme_minimal() +
    geom_sf(data = pol) +
    geom_sf(data = bar, fill = "blue")

for(r in sqrt(c(30,20,15,10, 7, 5, 3:1))) {
    mesh <- inla.mesh.2d(
        loc.domain = pl,
        offset = c(1, 5),
        max.edge = c(r/3, r),
        cutoff = r/6)
    cat("n(mesh) =", mesh$n, "\n")
    print(mark(
        inla.mesh.fem(mesh),
        mesh2fem.barrier(mesh, unlist(fm_contains(bar, mesh))),
        check = FALSE
    ))
}

