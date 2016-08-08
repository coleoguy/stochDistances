source("stochDistances.R")
## do 100 itterates un correlated and correlated


## First we simulate some data
library(geiger)
library(phytools)
library(viridis)
tree <- sim.bdtree(b=0.1, d=0, stop="taxa", n=100, seed=1)
par <- rbind(c(-.025, .025), c(.075, -.075))
par2 <- rbind(c(-.025, .025), c(.075, -.075))
data1 <- sim.char(tree, par, model = "discrete")
data2 <- sim.char(tree, par2, model = "discrete")
map1 <- make.simmap(tree, x=data1[,,1],nsims=1, model="ARD")
map2 <- make.simmap(tree, x=data2[,,1],nsims=1, model="ARD")
cols<-setNames(viridis(2),1:2)
par(mfcol=c(1,2))
plotSimmap(map1, colors=cols, lwd=3)
plotSimmap(map2, colors=cols, lwd=3)

# and now we check this basic uncorrelated data
nocor.samerate <- testDistances(tree, 
                     trait1=data1[,,1], 
                     trait2=data2[,,1], 
                     n=100, model="ARD")





# so now we want to construct 
# 100 no correlation equal rates
par <- rbind(c(-.05, .05), c(.05, -.05))
par2 <- rbind(c(-.05, .05), c(.05, -.05))
data1 <- data2 <- list()
for(i in 1:100){
  tree <- sim.bdtree(b=0.1, d=0, stop="taxa", n=100)
  data1[[i]] <- sim.char(tree, par, model = "discrete")[,,1]
  data2[[i]] <- sim.char(tree, par2, model = "discrete")[,,1]
}

# 100 no correlation assymetric rates data sets: both traits bias 1 to 2
par <- rbind(c(-.025, .025), c(.075, -.075))
par2 <- rbind(c(-.025, .025), c(.075, -.075))
data1 <- data2 <- list()
for(i in 1:100){
  tree <- sim.bdtree(b=0.1, d=0, stop="taxa", n=100)
  data1[[i]] <- sim.char(tree, par, model = "discrete")[,,1]
  data2[[i]] <- sim.char(tree, par2, model = "discrete")[,,1]
}

## ANALYZE

# 100 classic correlation 1 to 2 high in state 2 but normal in other
a<-.025
b<-.050
par <- rbind(c(-2*a, a, a, 0), c(a, -2*a, 0, a), 
             c(a, 0, -(a+b), b), c(0, a, a, -2*a))
data1 <- data2 <- list()
for(i in 1:100){
  tree <- sim.bdtree(b=0.1, d=0, stop="taxa", n=100)
  x <- sim.char(tree, par, model = "discrete")[,,1]
  y <- z <- x
  y[y <= 2] <- 1
  y[y > 2] <- 2
  z[z == 1 | z == 3] <- 1
  z[z == 2 | z == 4] <- 2
  data1[[i]] <- y
  data2[[i]] <- z
  }

# 100 classic correlation 1 to 2 low in state 2 but normal in other
a<-.05
b<-.025
par <- rbind(c(-2*a, a, a, 0), c(a, -2*a, 0, a), 
             c(a, 0, -(a+b), b), c(0, a, a, -2*a))
data1 <- data2 <- list()
for(i in 1:100){
  tree <- sim.bdtree(b=0.1, d=0, stop="taxa", n=100)
  x <- sim.char(tree, par, model = "discrete")[,,1]
  y <- z <- x
  y[y <= 2] <- 1
  y[y > 2] <- 2
  z[z == 1 | z == 3] <- 1
  z[z == 2 | z == 4] <- 2
  data1[[i]] <- y
  data2[[i]] <- z
}
