## function to check if the population is fixed at a trait
fixation <- function(pop, trait)
  sum(abs(diff(pop[,trait])))==0

## iterate a single generation
next_gen <- function(prms,
                     pop,
                     dispersal_evolution,
                     diapause_evolution,
                     diap_length_evolution) {

  pop <- death(prms, pop)
  ## MC: remove diapausers
  active_pos <- pop[,'current_diap']==0
  active_pop <- pop[active_pos,]
  inactive_pop <- pop[!active_pos,]
  
  if(ALT) {
    distances <- dist_mat_alt(prms, active_pop)
  } else {
    distances <- dist_mat(prms, active_pop)
  }

  active_pop <- birth(prms, active_pop, distances)
  if(nrow(active_pop)==1) return(pop)

  if(dispersal_evolution){
    active_pop <- mutation(prms,active_pop,'p1')
    ## active_pop <- mutation(prms,active_pop,'p2')
  }
  if(diap_length_evolution){
    active_pop <- mutation(prms, active_pop, 'diap_length')
  }
  if(diapause_evolution){
    active_pop <- mutation(prms, active_pop, 'diap_prob')
  }
  
  ## MC: add diapause step
  pop <- rbind(active_pop, inactive_pop)
  pop <- diapause(prms, pop)
  active_pos <- pop[,'current_diap']==0
  active_pop <- pop[active_pos,]
  inactive_pop <- pop[!active_pos,]

  active_pop <- move(prms, active_pop)
  
  pop <- rbind(active_pop, inactive_pop)

  return(pop)
}

## run the model
## PARAMETERS:
##    prms - parameter list
##    plot - boolean value governing plotting
##    brk_if_fix - boolean value: stop sim if a type fixes
## OUTPUT:
##   - Saves parameters used for the simulation and population statistics
##    for each generation.
##   - Returns a list with the saved values

run_sim <- function(prms,
                    plot,
                    brk_if_fix=TRUE,
                    dispersal_evolution=FALSE,
                    diapause_evolution=TRUE,
                    diap_length_evolution=FALSE) {
  
  ## create save structure and specify generations to save
  saved_pops <- vector(mode='list', prms$num_saves)
  save_gens <- floor(seq(from=1,to=prms$num_gens,length=prms$num_saves))
  save_counter <- 1

  columnnames <- list(gen=save_gens,
                      char=c('pop_size',
                             'mean_trait_p1',
                             'mean_trait_p2',
                             'mean_diap_prob',
                             'mean_diap_length',
                             'mean_fec')) ##add diapause and environmental sensitivity reporting here
  pop_stats <- matrix(NA,
                      nrow=prms$num_saves,
                      ncol=length(columnnames[['char']]),
                      dimnames=columnnames)

  pop <- initialize_pop(prms, dispersal_evolution=dispersal_evolution)


  
  for(i in 1:prms$num_gens) {
    ## plot every 'print_every' generations
    if(!is.na(prms$print_every)) {
      if(i%%prms$print_every==0) {
        cat(sprintf('Generation %d, ', i))
        cat(sprintf('Pop Size %d, ', nrow(pop)))
        cat(sprintf('mean(p1)=%2.3f, ', mean(pop[,'p1'])))
    ##    cat(sprintf('mean(p2)=%2.3f\n', mean(pop[,'p2'])))
        cat(sprintf('mean(diap)=%2.3f\n', mean(pop[,'diap_prob'])))
        if(prms$plot) plot_hm(prms, i, pop)
      }
    }
    ## next generation
    pop <- next_gen(prms,
                    pop,
                    dispersal_evolution=dispersal_evolution,
                    diapause_evolution=diapause_evolution,
                    diap_length_evolution=diap_length_evolution)
    
    if(length(pop)==1) {
      cat(pop, '\n')
      break
    }

    ## save population
    if(i %in% save_gens) {
      saved_pops[[save_counter]] <- pop

      pop_size <- nrow(pop)
      mean_trait_p1 <- mean(pop[,'p1'])
      mean_trait_p2 <- mean(pop[,'p2'])
      mean_diap_prob <- mean(pop[,'diap_prob'])
      mean_diap_length <- mean(pop[,'diap_length'])
      mean_fec <- mean(pop[,'fec_cost']) 
      pop_char <- cbind(pop_size,
                        mean_trait_p1,
                        mean_trait_p2,
                        mean_diap_prob,
                        mean_diap_length,
                        mean_fec)
      pop_stats[save_counter,] <- pop_char
      
      save_counter <- save_counter+1
    }

    ## check for fixation
    if(brk_if_fix)
      if(fixation(pop, 1)) {
        save(saved_pops, file=file.path(prms$file_path, 'saved_pops.RData'))
        break
      }
  }
  ## save(saved_pops, pop_stats, prms,
  ##      file=file.path(prms$file_path, 'out.RData', fsep=''))

  out <- list(prms=prms,
              saved_pops=saved_pops,
              pop_stats=pop_stats)
  return(out)
}

## run multiple replicate simulations
run_multiple <- function(reps,
                         land_type,
                         stripes=1,
                         acl,
                         gens,
                         dispersal_evolution=TRUE,
                         diapause_evolution=TRUE,
                         diap_length_evolution=FALSE,
                         ncores=1, 
                         init_p1=0.05,
                         init_diap=0,
                         f_max=75, 
						 c=1/15,
						 sigma_s=0.05,
                         mortality=0,
                         mut_rate=0.01,
                         print_every=50,
                         num_saves=10){
  
  run.rep <- function(i) {
    if(land_type=="uniform" | land_type=="one_gauss"){
      if(land_type=="uniform") {
	prms.rep <- base_prms(land_type=land_type,
                              het=0)
      } else {
        prms.rep <- base_prms(land_type=land_type,
                              het=2)
      }      
    }
     else if(land_type=="noisy"){	
       prms.rep <- base_prms(land_type=land_type,
       							acl=acl,
       							het=2)

       mat <- cn_2D(acl=acl, n=128, amp=1)$cn
       prms.rep$land_mat <- abs(min(mat)) + mat
       prms.rep$land_max <- max(prms.rep$land_mat)
       
     }
    else if(land_type=="striped"){
      prms.rep <- base_prms(land_type=land_type,
                            stripes=stripes,
                            het=2)
    }
    else print("not a valid landscape type")
    ## adjust prms as needed
    
    prms.rep$num_gens <- gens
    prms.rep$f_max <- f_max
    prms.rep$c <- c
    prms.rep$sigma_s <- sigma_s
    prms.rep$mortality <- mortality
    prms.rep$mut_rate <- mut_rate
    prms.rep$init_p1 <- init_p1
    prms.rep$init_diap <- init_diap
    prms.rep$print_every <- print_every
    prms.rep$num_saves   <- num_saves
    
    prms.rep <- prms_update(prms.rep)
    prms.rep <- update_land_prms(prms.rep)
    ##plot_hm(prms.rep, "did I make the right landscape?")
    #if(!dispersal_evolution) {
     # prms$sigma_m_max<-0.05			
    #}

    rep_out.list <- run_sim(prms=prms.rep,
            dispersal_evolution=dispersal_evolution,
            diapause_evolution=diapause_evolution,
            diap_length_evolution=diap_length_evolution)
            
    ##########save it here
	rep_number <- 1 + length(list.files(path=paste("./model_output/", dir_name, "/", sep="")))

  	if(land_type=="striped"){
  		file_name <-
    		paste("model_output/", dir_name, sprintf("/%d_diapause_ev_%s_disp_ev_%s_diap_length_ev_%s_stripes_%d_init_disp%.2f_init_diap%.2f_%dx%d.Rdata",
            rep_number, diap_ev, disp_ev, diap_length_ev, stripes, init_p1, init_diap, reps, gens), sep="")
	}
	else   	file_name <-
    	paste("model_output/", dir_name, sprintf("/%d_diapause_ev_%s_disp_ev_%s_diap_length_ev_%s_%s_acl_%.3f_init_disp%.2f_init_diap%.2f_%dx%d.Rdata",
            rep_number, diap_ev, disp_ev, diap_length_ev, land_type, acl,init_p1, init_diap, reps, gens), sep="")
	save(rep_out.list, file=file_name)
  }

##define variables for naming convention
	disp_ev <- ifelse(dispersal_evolution, 'T', 'F')
  	diap_ev <- ifelse(diapause_evolution, 'T', 'F')
  	diap_length_ev <- ifelse(diap_length_evolution, 'T', 'F')
  	
##make directory to save reps
	if(land_type=="striped"){
		dir_name <-sprintf("diapause_ev_%s_disp_ev_%s_diap_length_ev_%s_stripes_%d_init_disp%.2f_init_diap%.2f_x%d", diap_ev, disp_ev, diap_length_ev, stripes, init_p1, init_diap, gens)
		}
	else if(land_type=="noisy"){
				dir_name <-sprintf("diapause_ev_%s_disp_ev_%s_diap_length_ev_%s_noisy_acl_%.3f_init_disp%.2f_init_diap%.2f_x%d", diap_ev, disp_ev, diap_length_ev, acl, init_p1, init_diap, gens)
	}
	
	if(dir_name %in% list.files("./model_output")==FALSE){
		dir.create(paste("./model_output/", dir_name, sep=""))
	}
	
  if(ncores==1)
    out.list <- lapply(1:reps, run.rep)
  if(ncores>1)
    out.list <- mclapply(1:reps,
                         run.rep,
                         mc.cores=ncores,
                         mc.preschedule=FALSE)
}
