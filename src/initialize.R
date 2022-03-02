## ******************************************************************
## ******************** load source functions ***********************
## ******************************************************************
library('RColorBrewer') ## for generating colour palettes for plotting
library('parallel') ## for using parallel processing
## library('cubature') ## for 2D integration

## ******************************************************************
## ******************** source R files ******************************
## ******************************************************************
source('src/base_prms.R') ## creates the parameter list
source('src/meta_functions.R') ## run different types of simulations
source('src/plotting.R') ## plotting code
source('src/population.R') ## all the core biological functions
source('src/variation_functions.R') ## generate resource landscapes
source('src/landscape.R') ## advanced landscape generation functions

## ******************************************************************
## ******************** initialize the population *******************
## ******************************************************************

initialize_pop <- function(prms, dispersal_evolution) {
  ## assign individuals random initial locations
  x <- runif(prms$N)*prms$s_dim
  y <- runif(prms$N)*prms$s_dim

  ## assign individuals uniform initial dispersal probability
  p1 <- rep.int(prms$init_p1, prms$N)
  
  ## if(dispersal_evolution) p2 <- runif(prms$N)
  ## else  	p2 <- 1
  p2 <- 1
  
  ## MC: new diapause-specific traits
  diap_prob <- rep.int(prms$init_diap, prms$N) #runif(prms$N) #probability of entering diapause in a time-step
  diap_length <- 1 #time-steps of diapause
  ## MC: To add later
  ## sense_tau <- #weight sensitivity to tau_i, or other local conditions
  
  current_diap <- 0 #how many time-steps you currently are from being active

  fec_cost <- matrix(0,prms$N)
  colnames(fec_cost) <- c('fec_cost')

  bins <- prms$o_space[bins(x, y, prms)]

  cbind(x=x, y=y, p1=p1, p2=p2, diap_prob=diap_prob, diap_length=diap_length, current_diap=current_diap, fec_cost, bin=bins)
}
