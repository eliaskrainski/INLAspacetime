# The SPDE model with transparent barriers

## The transparent barrier model

The transparent barrier model consider a domain $\Omega$ which is
partitioned into $k$ sub-domains, $\Omega_{d}$ for
$d \in \{ 1,\ldots,k\}$, where $\cup_{d = 1}^{k}\Omega_{d} = \Omega$.
The SPDE model assumes a common marginal variance and a particular to
each $\Omega_{d}$, $r_{d}$.

From Bakka et al. (2019), the precision matrix is
$$\mathbf{Q} = \frac{1}{\sigma^{2}}\mathbf{R}_{r}{\widetilde{\mathbf{C}}}^{- 1}\mathbf{R}_{r}{\mspace{6mu}\text{for}\mspace{6mu}}\mathbf{R}_{r} = \mathbf{C} + \frac{1}{8}\sum\limits_{d = 1}^{k}r_{d}^{2}\mathbf{G}_{d},\;\;\;{\widetilde{\mathbf{C}}}_{r} = \frac{\pi}{2}\sum\limits_{d = 1}^{k}r_{d}^{2}{\widetilde{\mathbf{C}}}_{d}$$
where $\sigma^{2}$ is the marginal variance. The Finite Element Method -
FEM matrices: $\mathbf{C}$, defined as
$$\mathbf{C}_{i,j} = \langle\psi_{i},\psi_{j}\rangle = \int_{\Omega}\psi_{i}(\mathbf{s})\psi_{j}(\mathbf{s})\partial\mathbf{s},$$
computed over the whole domain, while $\mathbf{G}_{d}$ and
${\widetilde{\mathbf{C}}}_{d}$ are defined as a pair of matrices for
each subdomain
$$\left( \mathbf{G}_{d} \right)_{i,j} = \langle 1_{\Omega_{d}}\nabla\psi_{i},\nabla\psi_{j}\rangle = \int_{\Omega_{d}}\nabla\psi_{i}(\mathbf{s})\nabla\psi_{j}(\mathbf{s})\partial\mathbf{s}\;{\mspace{6mu}\text{and}\mspace{6mu}}\;\left( {\widetilde{\mathbf{C}}}_{d} \right)_{i,i} = \langle 1_{\Omega_{d}}\psi_{i},1\rangle = \int_{\Omega_{d}}\psi_{i}(\mathbf{s})\partial\mathbf{s}.$$

In the case when $r = r_{1} = r_{2} = \ldots = r_{k}$ we have
$\mathbf{R}_{r} = \mathbf{C} + \frac{r^{2}}{8}\mathbf{G}$ and
${\widetilde{\mathbf{C}}}_{r} = \frac{\pi r^{2}}{2}\widetilde{\mathbf{C}}$
giving
$$\mathbf{Q} = \frac{2}{\pi\sigma^{2}}\left( \frac{1}{r^{2}}\mathbf{C}{\widetilde{\mathbf{C}}}^{- 1}\mathbf{C} + \frac{1}{8}\mathbf{C}{\widetilde{\mathbf{C}}}^{- 1}\mathbf{G} + \frac{1}{8}\mathbf{G}{\widetilde{\mathbf{C}}}^{- 1}\mathbf{C} + \frac{r^{2}}{64}\mathbf{G}{\widetilde{\mathbf{C}}}^{- 1}\mathbf{G} \right)$$
which coincides with the stationary case in Lindgren and Rue (2015),
when using $\widetilde{\mathbf{C}}$ in place of $\mathbf{C}$.

## Implementation

We assume $r_{d} = p_{d}r$ for known $p_{1},\ldots,p_{k}$ constants.
This gives
$${\widetilde{\mathbf{C}}}_{r} = \frac{\pi r^{2}}{2}\sum\limits_{d = 1}^{k}p_{d}^{2}{\widetilde{\mathbf{C}}}_{d} = \frac{\pi r^{2}}{2}{\widetilde{\mathbf{C}}}_{p_{1},\ldots,p_{k}}{\mspace{6mu}\text{and}\mspace{6mu}}\frac{1}{8}\sum\limits_{d = 1}^{k}r_{d}^{2}\mathbf{G}_{d} = \frac{r^{2}}{8}\sum\limits_{d = 1}^{k}p_{d}^{2}{\widetilde{\mathbf{G}}}_{d} = \frac{r^{2}}{8}{\widetilde{\mathbf{G}}}_{p_{1},\ldots,p_{k}}$$
where ${\widetilde{\mathbf{C}}}_{p_{1},\ldots,p_{k}}$ and
${\widetilde{\mathbf{G}}}_{p_{1},\ldots,p_{k}}$ are pre-computed
matrices.

## References

Bakka, H., J. Vanhatalo, J. Illian, D. Simpson, and H. Rue. 2019.
“Non-Stationary Gaussian Models with Physical Barriers.” *Spatial
Statistics* 29 (March): 268–88.
https://doi.org/<https://doi.org/10.1016/j.spasta.2019.01.002>.

Lindgren, Finn, and Havard Rue. 2015. “Bayesian Spatial Modelling with
R-INLA.” *Journal of Statistical Software* 63 (19): 1–25.
