setwd('~/Documents/Github/Code-for-dormancy-dispersal-coevolution')
rm(list=ls())
source('src/initialize.R')
library(abind)
## ******************************************************************
ALT <- FALSE

#create parameters; other parameters to be modified can be explored in src/base_prms.R
prms <- base_prms(land_type="noisy", #landscape type: "striped" or "noisy"
					acl=0.001, #acl, only for noisy landscapes
                   f_max=75, 
                   plot=FALSE)

#generate landscape functions in prms
prms<-update_land_prms(prms)

#make a population
pop <- initialize_pop(prms)

#run a single iteration
out <- run_sim(prms,
                pop,
                dispersal_evolution=FALSE)

#Run multiple replicates for a given set of parameters (will save to /model_output)
out <- run_multiple(reps=2,
                    land_type="noisy",
                    stripes=1,
                    acl=0.01,
                    sigma_s=0.05,
                    gens=1e2,
                    init_p1=0.05,
                    init_diap=0,
                    mut_rate = 0.005,
                    mortality=0.05,
                    dispersal_evolution=FALSE,
                    diapause_evolution=TRUE,
                    diap_length_evolution=FALSE,
                    ncores=1,
                    print_every=10,
                    num_saves=10)
