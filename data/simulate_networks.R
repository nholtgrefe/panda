library(SiPhyNetwork)
library(ape)
set.seed(42)

################### NEWICK CREATION ########################################
## Converts evonet object to an eNewick string with branch lengths,
## internal node names (of the form N[xx]), and each hybrid of the form N[xx]#H[yy]
############################################################################

# Adds internal labels to newick string, starting the labeling at 'label_counter'
add_internal_labels <- function(newick, label_counter) {
  locs <- gregexpr("\\)", newick, perl = TRUE)[[1]]
  
  if (locs[1] != -1) {
    parts <- regmatches(newick, list(locs))[[1]]
    replacements <- character(length(parts))
    
    for (i in seq_along(parts)) {
      replacements[i] <- paste0(")N", label_counter)
      label_counter <- label_counter + 1
    }
    
    regmatches(newick, list(locs)) <- list(replacements)
  }
  
  list(newick = newick, label_counter = label_counter)
}

# Adds internal labels to the hybrids of the newick string, starting the labeling at 'label_counter'
add_hybrid_labels <- function(newick, label_counter, hybrid_labels) {
  locs <- gregexpr("([\\w\\.]+)?(#H\\d+)", newick, perl = TRUE)[[1]]
  if (locs[1] != -1) {
    matches <- regmatches(newick, list(locs))[[1]]
    
    replacements <- sapply(matches, function(m) {
      label_part  <- sub("([\\w\\.]+)?(#H\\d+)", "\\1", m, perl = TRUE)
      hybrid_part <- sub("([\\w\\.]+)?(#H\\d+)", "\\2", m, perl = TRUE)
      
      if (hybrid_part %in% names(hybrid_labels)) {
        lbl <- hybrid_labels[[hybrid_part]]
      } else {
        lbl <- paste0("N", label_counter)
        hybrid_labels[[hybrid_part]] <<- lbl
        label_counter <<- label_counter + 1
      }
      
      paste0(lbl, hybrid_part)
    }, USE.NAMES = FALSE)
    
    regmatches(newick, list(locs)) <- list(replacements)
  }
  list(newick = newick, label_counter = label_counter, hybrid_labels = hybrid_labels)
}

# Main function that returns the eNewick string of the network
evonet_to_newick <- function(net) {
  
  if (getNetworkLevel(net) == 0) {
    newick <- write.tree(net, file = "", append = FALSE)
  } else {
    newick <- write.evonet(net, file = "", append = FALSE)
  }
  
  label_counter <- 1
  hybrid_labels <- list()
  
  # Step 1
  res1 <- add_internal_labels(newick, label_counter)
  newick <- res1$newick
  label_counter <- res1$label_counter
  
  # Step 2
  res2 <- add_hybrid_labels(newick, label_counter, hybrid_labels)
  newick <- res2$newick
  label_counter <- res2$label_counter
  hybrid_labels <- res2$hybrid_labels
  
  return(newick)
}

############################################################################
############################################################################




################### SIMULATOR PER LEVEL #####################################
## simulate 'num_per_level' networks with 'n_tips' taxa for each 'target_level'.
# Also returns the used nu, mu and lambda value for each network
############################################################################

sample_networks_by_level <- function(target_levels = 0:10,
                                     num_per_level = 100,
                                     n_tips = 100,
                                     nu_range = c(0.0, 0.2),
                                     lambda = 1.0,
                                     mu = 0.2,
                                     hyb_props = c(0.5, 0.25, 0.25),
                                     inheritance_fxn = make.beta.draw(10, 10),
                                     max_attempts = 1e6) {
  
  # Convert levels to character for named list
  target_levels <- sort(unique(target_levels))
  level_names <- as.character(target_levels)
  
  # Initialize storage: list for each target level
  networks_by_level <- setNames(vector("list", length(target_levels)), level_names)
  for (lvl in level_names) networks_by_level[[lvl]] <- list()
  
  # Counters
  counts <- setNames(rep(0, length(target_levels)), level_names)
  
  attempt <- 1
  
  while (any(counts < num_per_level) && attempt <= max_attempts) {
    # Random hybridization rate
    nu <- runif(1, nu_range[1], nu_range[2])
    
    # Simulate one network
    net <- sim.bdh.taxa.ssa(
      n = n_tips,
      numbsim = 1,
      lambda = lambda,
      mu = mu,
      nu = nu,
      hybprops = hyb_props,
      hyb.inher.fxn = inheritance_fxn,
      complete=FALSE #Removes all extinct taxa s.t. net has exactly n tips
    )[[1]]
    
    # Only valid networks
    if (!inherits(net, "numeric")) {
      # --- check for parallel edges using regex ---
      if (getNetworkLevel(net) == 0) {
        newick <- write.tree(net, file = "", append = FALSE)
      } else {
        newick <- write.evonet(net, file = "", append = FALSE)
      }
      pattern <- "#H([0-9]+):[0-9.]+,#H\\1:[0-9.]+"
      if (grepl(pattern, newick)) {
        cat("Parallel edges!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        next
      }
      # -------------------------------
      
      level <- getNetworkLevel(net)
      level_chr <- as.character(level)
      
      # Check if this level is in target_levels and needs more networks
      if (level_chr %in% level_names && counts[level_chr] < num_per_level) {
        # Store as a list with network + parameters
        networks_by_level[[level_chr]][[length(networks_by_level[[level_chr]]) + 1]] <- list(
          network = net,
          nu = nu,
          mu = mu,
          lambda = lambda
        )
        
        counts[level_chr] <- counts[level_chr] + 1
        cat("Level", level, "count:", counts[level_chr], " (nu =", round(nu,5), ")\n")
      }
    }
    
    attempt <- attempt + 1
    if (attempt %% 1000 == 0) cat("Attempt:", attempt, "\n")
  }
  
  if (any(counts < num_per_level)) {
    warning("Max attempts reached before filling all target levels.")
  }
  
  return(networks_by_level)
}


############################################################################
############################################################################





################### MAIN SIMULATOR #########################################
# Main function to perform a simulation that samples networks.
# Takes as input
# -number_of_taxa_range: a list of n_taxa to simulate
# -target_level_range: a list of levels to simulate
# -set_size: the number of networks to simulate for a given level and n_taxa
# -nu_range_fun: a function that returns the min and max nu-value for a number of taxa n
# Returns a string with all simulated networks in newick format and a header
############################################################################

perform_simulation <- function(number_of_taxa_range=c(25, 50, 100, 200),
                               target_level_range=0:10,
                               set_size = 100,
                               nu_range_fun = function(n) c(0.0, 2.0/n),
                               lambda = 1.0,
                               mu = 0.2,
                               hyb_props = c(0.5, 0.25, 0.25),
                               inheritance_fxn = make.beta.draw(10, 10),
                               max_attempts = 1e6) {
  
  # Initialize counter and storage
  net_id_counter <- 1
  all_lines <- character()
  all_lines <- c(all_lines, "id ntips level nu mu lambda newick")
  
  # Loop over each number of taxa
  for (n in number_of_taxa_range) {
    
    cat("Simulating networks with", n, "tips...\n")
    
    # Determine nu range for this n
    nu_range <- nu_range_fun(n)
    
    # Sample networks
    res <- sample_networks_by_level(n_tips = n,
                                    num_per_level = set_size,
                                    target_levels = target_level_range,
                                    nu_range = nu_range,
                                    lambda = lambda,
                                    mu = mu,
                                    hyb_props = hyb_props,
                                    inheritance_fxn = inheritance_fxn,
                                    max_attempts = max_attempts)
    
    # Loop over levels and networks
    for (lvl in names(res)) {
      nets_list <- res[[lvl]]
      
      for (net_entry in nets_list) {
        net <- net_entry$network
        nu <- net_entry$nu
        lambda_val <- net_entry$lambda
        mu_val <- net_entry$mu
        ntips <- ape::Ntip(net)

        # Network ID XXXXX
        net_id <- sprintf("%05d", net_id_counter)
        
        # Convert to Newick string
        #newick <- write.net(net, file = "", append = FALSE)
        newick <- evonet_to_newick(net)

        # Build line
        line <- paste(net_id, ntips, lvl, round(nu,5), round(mu_val,5),
                      round(lambda_val,5), newick)

                all_lines <- c(all_lines, line)
        
        net_id_counter <- net_id_counter + 1
      }
    }
  }
  
  cat("Simulation complete.\n")
  return(all_lines)
}


############################################################################
############################################################################





################### MAIN ##################################################
# Main function that calls other functions to do simulation and save output
############################################################################

main <- function() {
  
  # Output file
  folder <- "/home/nholtgreve/Documents/Projects/Phylogenetics/Panda MAPPD/SiPhySimulations/"
  output_file <- file.path(folder, "simulated_networks3.txt")
  
  # Remove output file if it exists
  #if (file.exists(output_file)) return()
  
  # String containing networks and information (including header)
  res <- perform_simulation(
            number_of_taxa_range=c(20, 50, 100, 200),
            target_level_range=0:15,
            set_size = 100,
            nu_range_fun = function(n) c(0.0, 2.0/n)
          )

  # Write all lines to file at once
  writeLines(res, con = output_file)
}

############################################################################
############################################################################


