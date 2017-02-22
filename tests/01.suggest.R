## tests for auto suggest
library("shrt")
library("MultiPattern")

num.samples = 10

snames = p0("S", 1:num.samples)
aa = mtrx(rnorm(8*num.samples), col.names=letters[1:8], row.names=snames)
bb = mtrx(round(runif(8*num.samples)), col.names=LETTERS[1:8], row.names=snames)
ab = cbind(aa, bb)

my1 = MPnew(snames, data=list(all=ab))
MPeasyConfig(my1, data="all", type="subspace1")
MPeasyConfig(my1, data="all", type="subspaceR")
print(my1)


my2 = MPnew(snames, data=list(all=ab))
MPsuggestConfig(my2, data="all", verbose=T)
print(my2)


