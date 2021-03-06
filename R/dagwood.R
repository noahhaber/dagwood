#' @title DAGs with omitted objects displayed
#'
#' @description This function takes a root DAG and generates a DAGWOOD object from it. Details are available
#' from a preprint, here: https://arxiv.org/abs/2004.04251
#'
#' Current limitations to this package:
#' Instrumental variables support is experimental (conditional IVs not yet supported)
#' Does not currently support non-minimal adjustment sets from root DAGs
#' Does not currently match KUER branch DAG matching with the UUER branch DAGs
#'
#' Future features: improved graphical outputs, additional branch DAG types, possibly a GUI
#' Improve DAG import features and data checking (e.g. make sure IV is valid on entry)
#'
#' @param formula.DAG A string formula describing the DAG, in the style of DAGitty
#' @param exposure A character string identifying the exposure of interest (must match the formula above)
#' @param outcome A character string identifying the outcome of interest (must match the formula above)
#' @param formula.KUERs Not yet implemented
#' @param instrument The character string identifying which node is the instrumental variable of interest
#' @keywords DAG, causal inference
#' @export
#' @import dagitty
#' @importFrom utils combn
#' @examples
#' formula.DAG <-"Chocolate -> Alzheimers
#' Chocolate <- Education -> Alzheimers
#' Chocolate -> CV
#' CV -> Alzheimers"
#' exposure <- "Chocolate"
#' outcome <- "Alzheimers"
#' branch.DAGs <- dagwood(formula.DAG,exposure,outcome)$DAGs.branch

# Main DAGWOOD function
dagwood <- function(formula.DAG,exposure=NA,outcome=NA,formula.KUERs=NA,instrument=NA){
  # Clean up formula/dagitty objects
  {
    # Check if the formula provided is a dagitty object, and fill in appropriately
      if(dagitty::is.dagitty(formula.DAG)){
        DAG.root <- formula.DAG
        if (length(dagitty::exposures(formula.DAG))==0){
          dagitty::exposures(DAG.root) <- exposure
        } else {
          exposure <- dagitty::exposures(DAG.root)
        }
        if (length(dagitty::outcomes(formula.DAG))==0){
          dagitty::outcomes(DAG.root) <- outcomes
        } else {
          outcome <- dagitty::outcomes(DAG.root)
        }
        formula.DAG <- paste0(edges(DAG.root)$v,edges(DAG.root)$e,edges(DAG.root)$w,collapse=" \n ")
      } else {
        DAG.root <- dagitty(paste0("dag {",formula.DAG,"}"))
        dagitty::exposures(DAG.root) <- exposure
        dagitty::outcomes(DAG.root) <- outcome
      }
  }
  # Properties of the root DAG
    nodes.root <- names(DAG.root)
    edges.root <- dagitty::edges(DAG.root)
  # Properties of the KUERs (if any), set aside for later (not part of the root DAG)
    if (!is.na(formula.KUERs)){
      DAG.root.KUERs <- dagitty(paste0("dag {",paste0(formula.DAG,"\n",formula.KUERs,"}")))
      nodes.root.KUERs <- names(DAG.root.KUERs)
      nodes.KUERs <- nodes.root.KUERs[!nodes.root.KUERs %in% nodes.root]
    }

  # Identify exclusion restrictions
  {
    # Find every combinatorial of nodes
      combs <- utils::combn(nodes.root,2)
      n.combs <- ncol(combs)

    # Function for testing if valid DAGWOOD object, taking a branch DAG (root DAG in memory)
      test.single.uuer.comb <- function(i){
        nodes.UUER.temp <- combs[,i]
        # Three possibilities: ->, <-, or <-UUER->
          # First, try -> and test
            addition <- paste0(nodes.UUER.temp[1],"->",nodes.UUER.temp[2])
            DAG.branch.candidate <- dagitty::dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            forward <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root,instrument=instrument,restriction.type = "ER",changes.made=paste0("Added ",addition))
            forward <- forward[forward$verdict=="Passed",]
          # First, try <- and test
            addition <- paste0(nodes.UUER.temp[1],"<-",nodes.UUER.temp[2])
            DAG.branch.candidate <- dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            backward <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root,instrument=instrument,restriction.type = "ER",changes.made=paste0("Added ",addition))
            backward <- backward[backward$verdict=="Passed",]
          # Lastly, try <-UUER-> and test
            addition <- paste0(nodes.UUER.temp[1],"<-UUER->",nodes.UUER.temp[2])
            DAG.branch.candidate <- dagitty(paste0("dag {",formula.DAG,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            bidirectional <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root,instrument=instrument,restriction.type = "ER",changes.made=paste0("Added ",addition))
            bidirectional <- bidirectional[bidirectional$verdict=="Passed",]
          # Finally, merge into one data frame
            output <- rbind(forward,backward,bidirectional)
          return(output)
      }

      UUERs <- do.call("rbind",lapply(1:n.combs,function(x) test.single.uuer.comb(x)))
      #UUERs <- UUERs[UUERs$verdict=="Passed",]
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
      edges.root.tracking$distance <- 0

    # Recursive function for searching for MRs
      MR.search.recursive.outer <- function(edges.candidate.tracking,index.edge,depth=1,search.direction="downstream"){
        # First, attempt to find matches at current depth
          tests.at.level <- MR.search.recursive.at.depth(edges.candidate.tracking=edges.candidate.tracking,i=index.edge,depth.target=depth,search.direction=search.direction)
        # If it passes (or returns a blank), return the output
          if (any(tests.at.level$verdict=="Passed")|nrow(tests.at.level)==0){
            return(tests.at.level)
          } else {
            # If not, iterate and run the next level
            return(MR.search.recursive.outer(edges.candidate.tracking,index.edge,depth=depth+1))
          }
      }
    # Recursive function for searching for MRs, searching until it hits a target depth
      MR.search.recursive.at.depth <- function(edges.candidate.tracking,i,depth.target,search.direction="downstream"){
        # First, reset the flip candidates to 0
          edges.candidate.tracking$flip.candidate <- 0
        # Flip an edge
          if(edges.candidate.tracking$e[i] == "->"){
            edges.candidate.tracking$e[i] <- "<-"
            node.target.new <- edges.candidate.tracking$v[i]
            node.trailing.new <- edges.candidate.tracking$w[i]
          } else {
            edges.candidate.tracking$e[i] <- "->"
            node.target.new <- edges.candidate.tracking$w[i]
            node.trailing.new <- edges.candidate.tracking$v[i]
          }
          edges.candidate.tracking$flipped[i] <- 1
        # Convert to dagitty form
          DAG.branch.candidate <- dagitty::dagitty(paste0("dag {",paste(paste0(edges.candidate.tracking$v,edges.candidate.tracking$e,edges.candidate.tracking$w),collapse="\n"),"}"))
          dagitty::exposures(DAG.branch.candidate) <- exposure
          dagitty::outcomes(DAG.branch.candidate) <- outcome
        # Collect all changes for reporting
          recording.temp <- edges.candidate.tracking[edges.candidate.tracking$flipped==1,]
          changes.made <- paste(paste0(recording.temp$v,recording.temp$e,recording.temp$w),collapse=", ")
        # Test to see if this change results in a valid DAGWOOD branch DAG
        # Only test at target depth
          if (sum(edges.candidate.tracking$flipped)==depth.target){
            test.candidate <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root,instrument=instrument,restriction.type = "MR",changes.made=paste0("Flipped ",changes.made))
            return(test.candidate)
          } else {
            # If not yet at target depth, iterate to find more flip candidates
              if (search.direction == "downstream"){
                # Find all of edges connected to the new target edge which are unflipped
                  edges.candidate.tracking$flip.candidate <- ifelse((edges.candidate.tracking$w==node.target.new | edges.candidate.tracking$v==node.target.new) &
                                                                    edges.candidate.tracking$flipped==0,
                                                                  1,0)
              } else if (search.direction == "bidirectional") {
                # Find all of edges connected to either side of the newly flipped edge
                  edges.candidate.tracking$flip.candidate <- ifelse((edges.candidate.tracking$w==node.target.new | edges.candidate.tracking$v==node.target.new | edges.candidate.tracking$w==node.trailing.new | edges.candidate.tracking$v==node.trailing.new) &
                                                                      edges.candidate.tracking$flipped==0,
                                                                    1,0)
              }
            # If nothing is available, return an empty set
              if(sum(edges.candidate.tracking$flip.candidate)==0){
                test.candidate <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root)
                test.candidate$flips <- sum(edges.candidate.tracking$flipped)
                return(test.candidate[0,])
              } else {
            # Otherwise, iterate recursively over all possible pathways
                # Find the row indexes for the new candidates
                  new.candidates <- as.numeric(rownames(edges.candidate.tracking[edges.candidate.tracking$flip.candidate==1,]))
                # Run the function for all of the new candidates
                  output <- do.call("rbind",lapply(1:length(new.candidates),function(x) MR.search.recursive.at.depth(edges.candidate.tracking,new.candidates[x],depth.target=depth.target,search.direction=search.direction)))
                  return(output)
              }
          }
      }
    # Run for each edge
      MRs <- unique(do.call("rbind",lapply(1:nrow(edges.root.tracking),function(x) MR.search.recursive.outer(edges.root.tracking,x))))
  }

  # Combine into one big set of branch DAGs
    DAGs.tested <- rbind(UUERs,MRs)
    DAGs.branch <- DAGs.tested[DAGs.tested$verdict=="Passed",]

  # Export
    return(list("DAGs.branch" = DAGs.branch,"DAGs.tested"=DAGs.tested))
}

# Function for determining if branch candidate passes DAGWOOD rules
test.DAG.branch.candidate <- function(DAG.branch.candidate,DAG.root,instrument=NA,restriction.type=NA,changes.made=NA) {
  # Properties of the root DAG
    adjustment.set.root <- dagitty::adjustmentSets(DAG.root,effect="direct")
    adjustment.set.branch.candidate <- dagitty::adjustmentSets(DAG.branch.candidate,effect="direct")
    n.paths.frontdoor.root <- length(dagitty::paths(DAG.root,directed=TRUE)$paths)
    nodes.root <- names(DAG.root)
    edges.root <- dagitty::edges(DAG.root)

  # Create a reversed version for later use
    DAG.branch.candidate.reversed <- DAG.branch.candidate
    exposure <- exposures(DAG.root)
    outcome <- outcomes(DAG.root)
    dagitty::exposures(DAG.branch.candidate.reversed) <- outcome
    dagitty::outcomes(DAG.branch.candidate.reversed) <- exposure
    adjustment.set.branch.candidate.reverse <- dagitty::adjustmentSets(DAG.branch.candidate.reversed,effect="direct")

  # Rule 1: test if different than the root DAG
    rule.1 <- ifelse(DAG.branch.candidate != DAG.root,
                     "Passed","Failed, not different than root DAG")
  # Rule 2: Test if forms a valid, identifiable DAG.
    # Two parts: First check and make sure dagitty can identify an adjustment set (may be empty)
    # (note: requires also testing flipped exposure/outcome)
      rule.2 <- ifelse(!length(adjustment.set.branch.candidate)==0,
      "Passed","Failed, not valid identifiable causal DAG")
      if (rule.2!="Passed"){
        rule.2 <- ifelse(!length(adjustment.set.branch.candidate.reverse)==0,
          "Passed","Failed, not valid identifiable causal DAG")
      } else {}
    # Second: Check for bi-directional edges without nodes or non-directional edges
      if (rule.2=="Passed"){
        # Part 1: check the edges function for edge types
          edges.DAG.branch.candidate <- edges(DAG.branch.candidate)
          rule.2 <- ifelse(!any(edges.DAG.branch.candidate$e=="<->") | !any(edges.DAG.branch.candidate$e=="--"),
              "Passed","Failed, not valid identifiable causal DAG")
          if (rule.2=="Passed"){
            # Swap edges, and make sure each doesn't already exist
              edges.DAG.branch.candidate.swapped <- edges.DAG.branch.candidate
              edges.DAG.branch.candidate.swapped$w <- edges.DAG.branch.candidate$v
              edges.DAG.branch.candidate.swapped$v <- edges.DAG.branch.candidate$w
              rule.2 <- ifelse(nrow(unique(rbind(edges.DAG.branch.candidate,edges.DAG.branch.candidate.swapped)))==nrow(edges.DAG.branch.candidate)*2,
                     "Passed","Failed, not valid identifiable causal DAG")
          } else {}
      } else {}
  # If both of these are passed, continue on, otherwise set rules 3a and 3b to NA and stop
    if (rule.1=="Passed" && rule.2=="Passed"){
      # Rule 3a: Test if requires a change in the adjustment set
        if (is.na(instrument)){
          rule.3a <- ifelse(!setequal(unlist(adjustment.set.branch.candidate),unlist(adjustment.set.root)),
               "Passed","Failed, no required change in adjustment set")
          rule.3a <- ifelse(length(adjustment.set.branch.candidate)==0,
                 "Not determinable, no valid adjustment set",rule.3a)
        } else {
          IV.set.branch.candidate <- dagitty::instrumentalVariables(DAG.branch.candidate)
          ivs.branch.list <- apply(matrix(IV.set.branch.candidate, byrow = TRUE),1,function(x) paste0(unlist(x,recursive=FALSE),collapse="|"))
          rule.3a <- ifelse(!instrument %in% ivs.branch.list ,
                     "Passed","Failed, no required change in instrument adjustment set")
        }
      # Rule 3b: Test if changes the number of front door paths
        rule.3b <- ifelse(n.paths.frontdoor.root != length(dagitty::paths(DAG.branch.candidate,directed=TRUE)$paths),
                          "Passed","Failed, no change in number of frontdoor paths")
        # If this fails, flip exposure/outcome. If this passes but previous failed, this test passes
        if (rule.3b!="Passed"){
          rule.3b <- ifelse(length(adjustment.set.branch.candidate.reverse)!=0,
               "Passed",rule.3b)
        }
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
  return(output)
}
