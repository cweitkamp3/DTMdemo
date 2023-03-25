# DTMdemo

This packages implements tools for the calculation and analysis of Distance-to-Measure (DTM)-densities proposed in [Weitkamp et al. (2022)](https://arxiv.org/abs/2205.07689).
It is the basis for all simulation presented in said paper.

## Installation

In order to install this package, run:

```R
library(devtools)
install_github("cweitkamp3/DTMdemo")
```
## Usage
In the following, we illuatrate the usage of the main functions. First of all the note that the (r'th power of the) emprical DTM-function (see [Chazal et al. (2016)](https://www.jmlr.org/papers/v18/15-484.html) for a formal definition)
for all points in a sample can be calculated like this:
```R
square = matrix(0,n,2)
square[,1] = runif(n,0,1)
square[,2] = runif(n,0,1)


DTM.function(square, m=1)
```
The estimation and illustration of the DTM-densities of a list of samples can be done as presented below.
```R
disc = matrix(0,n,2)
radius = 0.5
r = radius*sqrt(runif(n,0,1))
theta = runif(n,0,1)*2*pi
disc[,1]= r*cos(theta)
disc[,2]= r*sin(theta)

sample = list(square, disc)

den.list=DTM.density(sample, m=1)
plot(den.list[[1]], ylim =c(0,5), main = "")
lines(den.list[[2]], col = "red")
```
Clearly, it is important to be able to compare different DTM-densities. As proposed in [Weitkamp et al. (2022)](https://arxiv.org/abs/2205.07689) this can be achieved
with the L1-distance. An approximated L1-distance matrix between all elements of a list of DTM-densities can be calculated like this:
```R
dist.mat = L1.density.dist(den.list, 0, 1)
```
## References

F. Chazal, B. Fasy, F. Lecci, B. Michel, A. Rinaldo and L. Wasserman. "Robust topological inference: Distance to a measure and kernel distance". The Journal of Machine Learning Research, 18(1), 5845-5884, 2016.

C. Weitkamp, K. Proksch, T. Staudt, B. Lelandais and C. Zimmer. "From small scales to large scales: Distance-to-Measure density based geometric analysis of complex data". Preprint arXiv:2205.07689, 2022.
