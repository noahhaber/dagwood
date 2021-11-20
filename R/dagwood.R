#' @title DAGs with Omitted Objects Displayed (DAGWOOD)
#'
#' @description DAGs with omitted objects displayed (DAGWOOD) are a framework to help reveal key hidden
#' assumptions in a causal DAG. This package provides an implementation of the DAGWOOD algorithm.
#' Details on how DAGWOOD works, can be used, and should be interpreted are available
#' from our preprint, here: https://arxiv.org/abs/2004.04251
#' This package is implemented based on the DAGITTY package by Johannes Textor.
#' More information is available here: http://dagitty.net/learn/index.html
#' 
#' DAGWOODS take a root DAG and generates DAGWOOD branch DAGs from it.
#' At present, there are two types of branch DAGs:
#' 1) Exclusion branch DAGs represent nodes which, if exist, violate key causal model assumptions.
#' The most common example are confounders. There are two types of exclusion branch DAGs at present:
#' known exclusion branch DAGs (KEBDs) representing nodes that an analyst knows exist, but did not
#' include them in the analytical model (for example, due to lack of data), and unknown exclusion branch
#' DAGs (UEBDs) representing unknown or otherwise undeclared nodes which violate key model assumptions.
#' 2) Misdirection branch DAGs represent DAGs in which a minimal number of DAG arrows might be flipped
#' in order to violate key model assumptions. Common examples might be the case when something that is
#' adjusted-for as a confounder is a collider, or if missing time nodes exist such that there is "reverse"
#' causality between the exposure and outcome.
#' 
#' DAGWOOD objects currently output all branch DAGs as $DAGs.branch from the DAGWOOD object, which includes
#' details of why the branch DAG was tested, the type of branch DAG, what was changed from the original root
#' DAG, and the branch DAG itself as a dagitty object. The branch DAGs can be viewed, manipulated, etc as
#' any other dagitty-based DAG.
#' 
#' Full details of all DAGs tested as potential branch DAGs can be found in $DAGs.tested.
#' 
#' Current limitations to this package:
#' Instrumental variables support is experimental (conditional IVs not yet supported)
#' Does not currently support non-minimal adjustment sets from root DAGs
#' Does not currently match KEBD branch DAG matching with the UEBD branch DAGs
#' Does not check for errors for the fixed arrows (i.e. assumes that they match the main formula correctly)
#'
#' Future features: improved graphical outputs, additional branch DAG types, possibly a GUI
#' Improve DAG import features and data checking (e.g. make sure IV is valid on entry)
#' 
#' Abbreviations used:
#' BD: branch DAG
#' DAG: directed acyclic graph
#' DAGWOOD: directed acyclic graph with omitted objects displayed
#' EBD: exclusion branch DAG
#' IV: instrumental variable
#' KEBD: known exclusion branch DAG
#' MBD: misdirection branch DAG
#' UEBD: unknown exclusion branch DAG
#'
#' @param DAG.root A string formula describing the DAG, in the format style of DAGitty
#' @param exposure A character string identifying the exposure of interest (must match the formula above)
#' @param outcome A character string identifying the outcome of interest (must match the formula above)
#' @param KEBDs (not yet implemented)
#' @param instrument The character string identifying which node is the instrumental variable of interest
#' @param fixed.arrows (experimental) These arrows are prevented from flipping direction in the misdirection branch DAG algorithm. Nodes and arrows should be entered in the form of DAGitty.
#' @keywords DAG causal inference DAGWOOD
#' @export
#' @import dagitty
#' @importFrom utils combn
#' @examples
#' 
#' # Generate a DAGWOOD from an example root DAG:
#' DAG.root <-"Chocolate -> Alzheimers
#' Chocolate <- Education -> Alzheimers
#' Chocolate -> CV
#' CV -> Alzheimers"
#' 
#' # Identify the exposure and outcome of interest
#' exposure <- "Chocolate"
#' outcome <- "Alzheimers"
#' 
#' # Generate branch DAGs
#' branch.DAGs <- dagwood(DAG.root,exposure,outcome)$DAGs.branch
#' 
#' # Display the first branch DAG in the list
#' library(ggdag)
#' ggdag(branch.DAGs$DAG.branch.candidate[1]) + theme_dag()

# Main DAGWOOD function
dagwood <- function(DAG.root,exposure=NA,outcome=NA,KEBDs=NA,instrument=NA,fixed.arrows=NA){
  # Clean up formula/dagitty objects
  {
    # Check if the formula provided is a dagitty object, and fill in appropriately
      if(dagitty::is.dagitty(DAG.root)){
        DAG.root <- DAG.root
        if (length(dagitty::exposures(DAG.root))==0){
          dagitty::exposures(DAG.root) <- exposure
        } else {
          exposure <- dagitty::exposures(DAG.root)
        }
        if (length(dagitty::outcomes(DAG.root))==0){
          dagitty::outcomes(DAG.root) <- outcomes
        } else {
          outcome <- dagitty::outcomes(DAG.root)
        }
        DAG.root <- paste0(edges(DAG.root)$v,edges(DAG.root)$e,edges(DAG.root)$w,collapse=" \n ")
      } else {
        DAG.root <- dagitty(paste0("dag {",DAG.root,"}"))
        dagitty::exposures(DAG.root) <- exposure
        dagitty::outcomes(DAG.root) <- outcome
      }
  }
  # Properties of the root DAG
    nodes.root <- names(DAG.root)
    edges.root <- dagitty::edges(DAG.root)
  # Properties of the KEBDs (if any), set aside for later (not part of the root DAG)
    if (!is.na(KEBDs)){
      DAG.root.KEBDs <- dagitty(paste0("dag {",paste0(DAG.root,"\n",KEBDs,"}")))
      nodes.root.KEBDs <- names(DAG.root.KEBDs)
      nodes.KEBDs <- nodes.root.KEBDs[!nodes.root.KEBDs %in% nodes.root]
    }

  # Identify exclusion branch DAGs
  {
    # Find every combinatorial of nodes
      combs <- utils::combn(nodes.root,2)
      n.combs <- ncol(combs)

    # Function for testing if valid DAGWOOD object, taking a branch DAG (root DAG in memory)
      test.single.UEBD.comb <- function(i){
        nodes.UEBD.temp <- combs[,i]
        # Three possibilities: ->, <-, or <-UEBD->
          # First, try -> and test
            addition <- paste0(nodes.UEBD.temp[1],"->",nodes.UEBD.temp[2])
            DAG.branch.candidate <- dagitty::dagitty(paste0("dag {",DAG.root,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            forward <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root,instrument=instrument,BD.type = "ER",changes.made=paste0("Added ",addition))
            forward <- forward[forward$verdict=="Passed",]
          # First, try <- and test
            addition <- paste0(nodes.UEBD.temp[1],"<-",nodes.UEBD.temp[2])
            DAG.branch.candidate <- dagitty(paste0("dag {",DAG.root,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            backward <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root,instrument=instrument,BD.type = "ER",changes.made=paste0("Added ",addition))
            backward <- backward[backward$verdict=="Passed",]
          # Lastly, try <-UEBD-> and test
            addition <- paste0(nodes.UEBD.temp[1],"<-UEBD->",nodes.UEBD.temp[2])
            DAG.branch.candidate <- dagitty(paste0("dag {",DAG.root,"\n",addition,"}"))
            dagitty::exposures(DAG.branch.candidate) <- exposure
            dagitty::outcomes(DAG.branch.candidate) <- outcome
            bidirectional <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root,instrument=instrument,BD.type = "ER",changes.made=paste0("Added ",addition))
            bidirectional <- bidirectional[bidirectional$verdict=="Passed",]
          # Finally, merge into one data frame
            output <- rbind(forward,backward,bidirectional)
          return(output)
      }

      UEBDs <- do.call("rbind",lapply(1:n.combs,function(x) test.single.UEBD.comb(x)))
      #UEBDs <- UEBDs[UEBDs$verdict=="Passed",]
  }

  # Identify misdirection branch DAGs
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
      
    # Add information in for any fixed arrows
      if (!is.na(fixed.arrows)){
        edges.root.fixed <- edges(dagitty(dagitty(paste0("dag {",fixed.arrows,"}"))))[c("v","w","e")]
        edges.root.fixed$fixed.arrow <- 1
        edges.root.tracking <- merge(edges.root.tracking,edges.root.fixed,by=c("v","e","w"),all.x=TRUE,all.y=FALSE)
        edges.root.tracking$fixed.arrow <- ifelse(is.na(edges.root.tracking$fixed.arrow),0,edges.root.tracking$fixed.arrow)
      } else {
        edges.root.tracking$fixed.arrow <- 0
      }

    # Recursive function for searching for MBDs
      MBD.search.recursive.outer <- function(edges.candidate.tracking,index.edge,depth=1,search.direction="downstream"){
        # First, attempt to find matches at current depth
          tests.at.level <- MBD.search.recursive.at.depth(edges.candidate.tracking=edges.candidate.tracking,i=index.edge,depth.target=depth,search.direction=search.direction)
        # If it passes (or returns a blank), return the output
          if (any(tests.at.level$verdict=="Passed")|nrow(tests.at.level)==0){
            return(tests.at.level)
          } else {
            # If not, iterate and run the next level
            return(MBD.search.recursive.outer(edges.candidate.tracking,index.edge,depth=depth+1))
          }
      }
    # Recursive function for searching for MBDs, searching until it hits a target depth
      MBD.search.recursive.at.depth <- function(edges.candidate.tracking,i,depth.target,search.direction="downstream"){
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
            test.candidate <- test.DAG.branch.candidate(DAG.branch.candidate = DAG.branch.candidate,DAG.root=DAG.root,instrument=instrument,BD.type = "MBD",changes.made=paste0("Flipped ",changes.made))
            return(test.candidate)
          } else {
            # If not yet at target depth, iterate to find more flip candidates
              if (search.direction == "downstream"){
                # Find all of edges connected to the new target edge which are unflipped and not fixed
                  edges.candidate.tracking$flip.candidate <- ifelse((edges.candidate.tracking$w==node.target.new | edges.candidate.tracking$v==node.target.new) &
                                                                    edges.candidate.tracking$flipped==0 & edges.candidate.tracking$fixed.arrow == 0,
                                                                  1,0)
              } else if (search.direction == "bidirectional") {
                # Find all of edges connected to either side of the newly flipped edge which are unflipped and not fixed
                  edges.candidate.tracking$flip.candidate <- ifelse((edges.candidate.tracking$w==node.target.new | edges.candidate.tracking$v==node.target.new | edges.candidate.tracking$w==node.trailing.new | edges.candidate.tracking$v==node.trailing.new) &
                                                                      edges.candidate.tracking$flipped==0 & edges.candidate.tracking$fixed.arrow == 0,
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
                  output <- do.call("rbind",lapply(1:length(new.candidates),function(x) MBD.search.recursive.at.depth(edges.candidate.tracking,new.candidates[x],depth.target=depth.target,search.direction=search.direction)))
                  return(output)
              }
          }
      }
    # Run for each edge that isn't fixed
      unfixed.edges <- c(1:nrow(edges.root.tracking))[edges.root.tracking$fixed.arrow==0]
      MBDs <- unique(do.call("rbind",lapply(unfixed.edges,function(x) MBD.search.recursive.outer(edges.root.tracking,x))))
  }

  # Combine into one big set of branch DAGs
    DAGs.tested <- rbind(UEBDs,MBDs)
    DAGs.branch <- DAGs.tested[DAGs.tested$verdict=="Passed",]

  # Export
    return(list("DAGs.branch" = DAGs.branch,"DAGs.tested"=DAGs.tested))
}

# Function for determining if branch candidate passes DAGWOOD rules
test.DAG.branch.candidate <- function(DAG.branch.candidate,DAG.root,instrument=NA,BD.type=NA,changes.made=NA) {
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
  output <- data.frame(verdict,rule.1,rule.2,rule.3a,rule.3b,BD.type,changes.made,DAG.branch.candidate,
                       stringsAsFactors = FALSE)
  return(output)
}
