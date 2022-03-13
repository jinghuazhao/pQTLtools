#!/usr/bin/bash

library("GeneNet")

# A random network with 20 nodes and 10 percent (=19) edges
true.pcor <- ggm.simulate.pcor(20, 0.1)

# convert to edge list
test.results <- ggm.list.edges(true.pcor)[1:19,]

# Rgraphviz
nlab <- LETTERS[1:20]
gr <- network.make.graph( test.results, nlab)
gr
num.nodes(gr)
edge.info(gr)
gr2 <- network.make.graph( test.results, nlab, drop.singles=TRUE)
gr2
num.nodes(gr2)
edge.info(gr2)

# plot network
library("Rgraphviz")
plot(gr, "fdp")
plot(gr2, "fdp")
