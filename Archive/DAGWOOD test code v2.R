# Current implementation assumes that only the minimal adjustment set is in the root DAG.
# Also assumes that there does not exist any nodes called "UUER"

# Setup
{
  rm(list=ls())
  
  library(dagitty)
  library(ggdag)
  theme_set(theme_dag())
}

# Default input parameters
{
  formula.DAG <- NA
  exposure <- NA
  outcome <- NA
  adjustmentSet <- NA
  instrument <- NA
  KUERs <- NA
}

# Input parameters
{
  # Input style has the user enter the root DAG separately from the KUERs
  formula.DAG <-
    "x -> y
    x <- c1 ->y
    x <- c2 ->y
    x -> m
    m -> y"
  formula.KUERs <- "x <- k1 -> y"
  exposure <- "x"
  outcome <- "y"
}

# Function for determining if branch candidate passes DAGWOOD rules
test.dag.branch.candidate <- function(dag.branch.candidate,restriction.type=NA,changes.made=NA) {
  # Rule 1: test if different than the root DAG
    rule.1 <- ifelse(dag.branch.candidate != dag.root,
                     "Passed","Failed, not different than root DAG")
  # Rule 2: test if forms a valid DAG
    adjustment.set.branch.candidate <- dagitty::adjustmentSets(dag.branch.candidate,effect="direct")
    rule.2 <- ifelse(!length(adjustment.set.branch.candidate)==0,
    "Passed","Failed, not valid identifiable causal DAG")
  # If both of these are passed, continue on, otherwise set rules 3a and 3b to NA and stop
    if (rule.1=="Passed" && rule.2=="Passed"){
      # Rule 3a: Test if requires a change in the adjustment set
        rule.3a <- ifelse(!setequal(unlist(adjustment.set.branch.candidate),unlist(adjustment.set.root)),
               "Passed","Failed, no required change in adjustment set")
      # Rule 3b: Test if changes the number of front door paths
        rule.3b <- ifelse(n.paths.frontdoor.root != length(dagitty::paths(dag.branch.candidate,directed=TRUE)),
                          "Passed","Failed, no change in number of frontdoor paths")
      # Determine verdict
        verdict <- ifelse(rule.3a=="Passed" | rule.3b=="Passed","Passed","Failed, not a valid branch DAG")
    } else {
      verdict <- "Failed, not a valid branch DAG"
      rule.3a <- NA
      rule.3b <- NA
    }
  output <- data.frame(verdict,rule.1,rule.2,rule.3a,rule.3b,restriction.type,changes.made,dag.branch.candidate[1])
}

# Main DAGWOOD function
{
  # Establish the root DAG and all relavant parameters of that root DAG
  {
    # Properties of the root DAG
      dag.root <- dagitty(paste0("dag {",formula.DAG,"}"))
      dagitty::exposures(dag.root) <- exposure
      dagitty::outcomes(dag.root) <- outcome
      adjustment.set.root <- dagitty::adjustmentSets(dag.root,effect="direct")
      nodes.root <- names(dag.root)
      n.paths.frontdoor.root <- length(dagitty::paths(dag.root,directed=TRUE))
      edges.root <- dagitty::edges(dag.root)
    # Properties of the KUERs (if any), set aside for later (not part of the root DAG)
      if (!is.na(formula.KUERs)){
        dag.root.KUERs <- dagitty(paste0("dag {",paste0(formula.DAG,"\n",formula.KUERs,"}")))
        nodes.root.KUERs <- names(dag.root.KUERs)
        nodes.KUERs <- nodes.root.KUERs[!nodes.root.KUERs %in% nodes.root]
      }
  }
  
  # Find every combinatorial of nodes
    combs <- combn(nodes.root,2)
    n.combs <- ncol(combs)
  
  # Function for testing if valid DAGWOOD object, taking a branch DAG (root DAG in memory)
    test.single.uuer.comb <- function(i){
      nodes.UUER.temp <- combs[,i]
      # Three possibilities: ->, <-, or <-UUER->
        # First, try -> and test
          addition <- paste0(nodes.UUER.temp[1],"->",nodes.UUER.temp[2])
          dag.branch.candidate <- dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
          dagitty::exposures(dag.branch.candidate) <- exposure
          dagitty::outcomes(dag.branch.candidate) <- outcome
          forward <- test.dag.branch.candidate(dag.branch.candidate = dag.branch.candidate,restriction.type = "ER",paste0("Added ",addition))
         # First, try <- and test
          addition <- paste0(nodes.UUER.temp[1],"<-",nodes.UUER.temp[2])
          dag.branch.candidate <- dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
          dagitty::exposures(dag.branch.candidate) <- exposure
          dagitty::outcomes(dag.branch.candidate) <- outcome
          backward <- test.dag.branch.candidate(dag.branch.candidate = dag.branch.candidate,restriction.type = "ER",paste0("Added ",addition))
        # Lastly, try <-UUER-> and test
          addition <- paste0(nodes.UUER.temp[1],"<-UUER->",nodes.UUER.temp[2])
          dag.branch.candidate <- dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
          dagitty::exposures(dag.branch.candidate) <- exposure
          dagitty::outcomes(dag.branch.candidate) <- outcome
          bidirectional <- test.dag.branch.candidate(dag.branch.candidate = dag.branch.candidate,restriction.type = "ER",paste0("Added ",addition))
        # Finally, merge into one data frame
          output <- rbind(forward,backward,bidirectional)
        return(output)
    }
    
    UUERs <- do.call("rbind",lapply(1:n.combs,function(x) test.single.uuer.comb(x)))
}
