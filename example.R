
# Example Code

## Packages and functions
library(cshapes)
library(sp)
source("R/functions.R")

## Cshapes
sud <- cspacefillrsud <- cshapes::cshp()
sud <- sud[sud$gwcode %in% 625:626,] 
sud <- sf::as_Spatial(sud)

plot(sud)
## Split 
sud.split <- get_poly_const_parts(poly = sud, 
                                  size.cutoff = .001, 
                                  base.shape = NULL)

## Make voronois
sud.cells <- sample_vorcells(spdf = sud.split, 
                             res = .01, 
                             size = 10000, 
                             size.km2 = T, seed = 1, 
                             iter.max = 100, sample.type = "nonaligned",
                             ncore = 1)
plot(sud.cells)
plot(sud.split, border = "red", add = T)