library(sf)
library(fmesher)
library(INLAspacetime)
library(ggplot2)
library(bench)

s <- 100
pl <- cbind(c(0,1,1,0,0), c(0,0,1,1,0)) * s
pol <- st_sfc(st_polygon(list(pl)))

dataf <- data.frame(
    rr=c(78, 53, 36, 25.5, 18, 12.5, 8.8, 6.15)/10,
    n = NA, th=NA, tm=NA, itsec1=NA, itsec2=NA)
##dataf <- dataf[1:min(7,nrow(dataf)),]

e <- c("c0", "c1", "g1", "g2", "va", "ta")
for(k in 1:nrow(dataf)) {
    r <- dataf$rr[k]
    t0 <- Sys.time()
    hl <- fm_hexagon_lattice(
        st_buffer(pol, r), r/2)
    t1 <- Sys.time()
    mesh <- fm_mesh_2d(
        loc = hl, 
        offset = sqrt(s),
        max.edge = r*2, 
        cutoff = r/4)
    t2 <- Sys.time()
    cat("r =", r, "n. nodes =", mesh$n, "\n")
    a <- fm_fem(mesh)[e]
    b <- mesh2fem(mesh)[e]
    for(i in e[1:4])
        cat(i, all.equal(a[[i]]@x, b[[i]]@x), "")
    for(i in e[5:6])
        cat(i, all.equal(a[[i]], b[[i]]), "")
    cat("\n")
    bmk <- mark(
        fm_fem(mesh)[e],
        mesh2fem(mesh)[e], check = FALSE
    )
    print(bmk)
    dataf$n[k] <- mesh$n
    dataf$th[k] <- difftime(t1, t0, units = 'secs')
    dataf$tm[k] <- difftime(t2, t1, units = 'secs')
    dataf$itsec1[k] <- bmk$'itr/sec'[1]
    dataf$itsec2[k] <- bmk$'itr/sec'[2]
}

par(mfrow = c(1, 2), mar = c(4,4,0.5,0.5), mgp = c(3,2,0), bty = "n")
plot(dataf$n, dataf$itsec1, pch = 19, log = 'xy',
     ylim = range(dataf$itsec1, dataf$itsec2),
     xlab = "n", ylab = 'It. per sec.')
points(dataf$n, dataf$itsec2, pch = 8, col = 2)
legend("topright", c("fm", "R"), pch = c(19,8), col = 1:2, bty = 'n')
plot(data.frame(n=dataf$n, efr=dataf$itsec1/dataf$itsec2),
     pch = 19, bty = 'n', log = 'x',
     ylab = 'relative efficience fm/R')

if(FALSE) {
    plot(mesh)
    plot(pol, add = TRUE, lwd = 2)
}


## FEM matrices with barrier domain(s)

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

for(k in 1:nrow(dataf)) {
    r <- dataf$rr[k]
    hl <- fm_hexagon_lattice(
        st_buffer(pol, r), r/2)
    mesh <- fm_mesh_2d(
        loc = hl, 
        offset = sqrt(s),
        max.edge = r*2, 
        cutoff = r/4)
    cat("r =", r, "n. nodes =", mesh$n, "\n")
    itr <- unlist(fm_contains(bar, mesh))
    t0 <- Sys.time()
    a <- fm_fem(mesh)
    t1 <- Sys.time()
    print(t1-t0)
    b <- mesh2fem.barrier(mesh, barrier.triangles = itr)
    print(Sys.time()-t1)
    cat('c1 :', all.equal(a$c1@x, b$I@x),
        'g1 :', all.equal(a$g1@x, INLAtools::Sparse(Reduce("+", b$D))@x),
        'g2 :', all.equal(a$c0@x, INLAtools::Sparse(Reduce("+", b$C))@x), "\n")
}

