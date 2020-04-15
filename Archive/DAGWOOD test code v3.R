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
  # formula.DAG <-
  #   "x -> y
  #   x <- c1 ->y
  #   x <- c2 ->y
  #   x -> m
  #   m -> y"
  formula.KUERs <- "x <- k1 -> y"
  exposure <- "x"
  outcome <- "y"
    formula.DAG <-
    "x -> y"
}

# Function for determining if branch candidate passes DAGWOOD rules
test.DAG.branch.candidate <- function(DAG.branch.candidate,restriction.type=NA,changes.made=NA) {
  # Rule 1: test if different than the root DAG
    rule.1 <- ifelse(DAG.branch.candidate != DAG.root,
                     "Passed","Failed, not different than root DAG")
  # Rule 2: test if forms a valid DAG (note: requires also testing flipped exposure/outcome)
    adjustment.set.branch.candidate <- dagitty::adjustmentSets(DAG.branch.candidate,effect="direct")
    rule.2 <- ifelse(!length(adjustment.set.branch.candidate)==0,
    "Passed","Failed, not valid identifiable causal DAG")
  # If both of these are passed, continue on, otherwise set rules 3a and 3b to NA and stop
    if (rule.1=="Passed" && rule.2=="Passed"){
      # Rule 3a: Test if requires a change in the adjustment set
        rule.3a <- ifelse(!setequal(unlist(adjustment.set.branch.candidate),unlist(adjustment.set.root)),
               "Passed","Failed, no required change in adjustment set")
      # Rule 3b: Test if changes the number of front door paths
        rule.3b <- ifelse(n.paths.frontdoor.root != length(dagitty::paths(DAG.branch.candidate,directed=TRUE)),
                          "Passed","Failed, no change in number of frontdoor paths")
      # Determine verdict
        verdict <- ifelse(rule.3a=="Passed" | rule.3b=="Passed","Passed","Failed, not a valid branch DAG")
    } else {
      verdict <- "Failed, not a valid branch DAG"
      rule.3a <- NA
      rule.3b <- NA
    }
    
  DAG.branch.candidate <- DAG.branch.candidate[1]
  output <- data.frame(verdict,rule.1,rule.2,rule.3a,rule.3b,restriction.type,changes.made,DAG.branch.candidate,
                       stringsAsFactors = FALSE)
}

# Main DAGWOOD function
{
  # Establish the root DAG and all relavant parameters of that root DAG
  # Keep the root DAG in the main functional environment
  {
    # Properties of the root DAG
      DAG.root <- dagitty(paste0("dag {",formula.DAG,"}"))
      dagitty::exposures(DAG.root) <- exposure
      dagitty::outcomes(DAG.root) <- outcome
      adjustment.set.root <- dagitty::adjustmentSets(DAG.root,effect="direct")
      nodes.root <- names(DAG.root)
      n.paths.frontdoor.root <- length(dagitty::paths(DAG.root,directed=TRUE))
      edges.root <- dagitty::edges(DAG.root)
    # Properties of the KUERs (if any), set aside for later (not part of the root DAG)
      if (!is.na(formula.KUERs)){
        DAG.root.KUERs <- dagitty(paste0("dag {",paste0(formula.DAG,"\n",formula.KUERs,"}")))
        nodes.root.KUERs <- names(DAG.root.KUERs)
        nodes.KUERs <- nodes.root.KUERs[!nodes.root.KUERs %in% nodes.root]
      }
  }
  
  # Identify exclusion restrictions
  {
    # Find every combinatorial of nodes
      combs <- combn(nodes.root,2)
      n.combs <- ncol(combs)
    
    # Function for testing if valid DAGWOOD object, taking a branch DAG (root DAG in memory)
      test.single.uuer.comb <- function(i){
        nodes.UUER.temp <- combs[,i]
        # Three possibilities: ->, <-, or <-UUER->
          # First, try -> and test
            addition <- paste0(nodes.UUER.temp[1],"->",nodes.UUER.temp[2])
            DAG.branch.candidate <- dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            forward <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,restriction.type = "ER",paste0("Added ",addition))
           # First, try <- and test
            addition <- paste0(nodes.UUER.temp[1],"<-",nodes.UUER.temp[2])
            DAG.branch.candidate <- dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            backward <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,restriction.type = "ER",paste0("Added ",addition))
          # Lastly, try <-UUER-> and test
            addition <- paste0(nodes.UUER.temp[1],"<-UUER->",nodes.UUER.temp[2])
            DAG.branch.candidate <- dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            bidirectional <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,restriction.type = "ER",paste0("Added ",addition))
          # Finally, merge into one data frame
            output <- rbind(forward,backward,bidirectional)
          return(output)
      }
      
      UUERs <- do.call("rbind",lapply(1:n.combs,function(x) test.single.uuer.comb(x)))
      UUERs <- UUERs[UUERs$verdict=="Passed",]
  }
  
  # Identify misdirection restrictions
  {
    # Start with master edges list, adding column for whether or not edge has been flipped
      edges.root.tracking <- edges.root[c("v","e","w")]
      edges.root.tracking <- edges.root[c("v","e","w")]
      edges.root.tracking$v <- as.character(edges.root.tracking$v)
      edges.root.tracking$e <- as.character(edges.root.tracking$e)
      edges.root.tracking$w <- as.character(edges.root.tracking$w)
      edges.root.tracking$e.orig <- edges.root.tracking$e
    # Generate columns to keep track of what has been flipped and what is available as a candidate to flip
      edges.root.tracking$flipped <- 0
      edges.root.tracking$flip.candidate <- 0

    # Function start
      MR.search.recursive <- function(edges.candidate.tracking,i){
        # TEMP FOR DEBUGGING
          edges.candidate.tracking <- edges.root.tracking
          i <- 1
        
        # First, reset the flip candidates to 0
          edges.candidate.tracking$flip.candidate <- 0
        # Put the edge down, flip it....
          if(edges.candidate.tracking$e[i] == "->"){
            edges.candidate.tracking$e[i] <- "<-"
            node.target.new <- edges.candidate.tracking$v[i]
          } else {
            edges.candidate.tracking$e[i] <- "->"
            node.target.new <- edges.candidate.tracking$w[i]
          }
          edges.candidate.tracking$flipped[i] <- 1
        # Convert to dagitty form
          DAG.branch.candidate <- dagitty::dagitty(paste0("dag {",paste(paste0(edges.candidate.tracking$v,edges.candidate.tracking$e,edges.candidate.tracking$w),collapse="\n"),"}"))
          dagitty::exposures(DAG.branch.candidate) <- exposure
          dagitty::outcomes(DAG.branch.candidate) <- outcome
        # Collect all changes for reporting
          recording.temp <- edges.candidate.tracking[edges.candidate.tracking$flipped==1,]
          changes.made <- paste(paste0(edges.candidate.tracking[edges.candidate.tracking$flipped==1,]$v,edges.candidate.tracking[edges.candidate.tracking$flipped==1,]$e,edges.candidate.tracking[edges.candidate.tracking$flipped==1,]$w),collapse=", ")
        # Test to see if this change results in a valid DAGWOOD branch DAG
          test.candidate <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,restriction.type = "MR",paste0("Flipped ",changes.made))
          if (test.candidate$verdict[1]== "Passed"){
            # If an MR is found, SUCCESS! Return it, game over
            return(test.candidate)
          } else {
            # Find all of edges connected to the new target edge which are unflipped
              edges.candidate.tracking$flip.candidate <- ifelse((edges.candidate.tracking$w==node.target.new | edges.candidate.tracking$v==node.target.new) & 
                                                                  edges.candidate.tracking$flipped==0,
                                                                1,0)
            # If nothing is available, return an empty set
              if(sum(edges.candidate.tracking$flip.candidate)==0){
                return(test.candidate[0,])
              } else {
                # Find the row indexes for the new candidates
                new.candidates <- as.numeric(rownames(edges.candidate.tracking[edges.candidate.tracking$flip.candidate==1,]))
                do.call("rbind",lapply(1:length(new.candidates),function(x) MR.search.recursive(edges.candidate.tracking,new.candidates[x])))
              }
          }
      }
    # Run for each edge
      MRs <- unique(do.call("rbind",lapply(1:nrow(edges.root.tracking),function(x) MR.search.recursive(edges.root.tracking,x))))
      
  }
  
  # Combine into one big set of branch DAGs
  DAGs.branch <- rbind(UUERs,MRs)
}
