rm(list=ls())

library(dagitty)
library(ggdag)
theme_set(theme_dag())


exposure <- "x"
outcome <- "y"
adjustmentSet <- c("c1","c2")
instrument <- NA

formula <-
  "x -> y
  x <- c1 ->y
  x <- c2 ->y
  x <- kuerp1 ->y"

dag.empirical <- dagitty(paste0("dag {",formula,"}"))

exposures(dag.empirical) <- exposure
outcomes(dag.empirical) <- outcome
latents(dag.empirical) <- 
plot(graphLayout(dag.empirical))
node.names.orig <- names(dag.empirical)
adjustmentSets.orig <- adjustmentSets(dag.empirical,"x","y")

potential.kuerps <- node.names.orig[!node.names.orig %in% exposure & !node.names.orig %in% outcome & !node.names.orig %in% adjustmentSet ]
kuerps <- potential.kuerps %in% adjustmentSets.orig[1]


# Assume that everything that isn't
if (is.na(adjustmentSet)){
  
}

setequal(adjustmentSet,adjustmentSets(dag.empirical,exposure,outcome)[1])






edges(dag.empirical)

dag.new <- dagitty(paste0(substr(dag.empirical[1],1,nchar(dag.empirical[1])-3),"\n","x<-u->y","}\n"))
