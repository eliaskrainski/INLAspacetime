# Illustrative code for Finite Element matrices of a mesh in 2d domain.

Illustrative code for Finite Element matrices of a mesh in 2d domain.

Illustrative code for Finite Element matrices when some triangles are in
a barrier domain.

## Usage

``` r
mesh2fem(mesh, order = 2, barrier.triangles = NULL)

mesh2fem.barrier(mesh, barrier.triangles = NULL)
```

## Arguments

- mesh:

  a 2d mesh object.

- order:

  the desired order.

- barrier.triangles:

  integer index to specify the triangles in the barrier domain

## Value

a list object containing the FE matrices.

a list object containing the FE matrices for the barrier problem.
