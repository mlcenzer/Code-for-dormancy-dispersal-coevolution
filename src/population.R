## ******************************************************************
## ******************** population functions ************************
## ******************************************************************

## vector-based wrap-around distance function
dist_wa <- function(prms, x, ii, jj){
  f <- function(i, j) {
    dx <- abs(x[i,1] - x[j,1])
    dy <- abs(x[i,2] - x[j,2])
    gtr_x <- dx>prms$s_dim/2
    gtr_y <- dy>prms$s_dim/2
    dx[gtr_x] <- (prms$s_dim-dx)[gtr_x]
    dy[gtr_y] <- (prms$s_dim-dy)[gtr_y]
    out <- sqrt(dx^2 + dy^2)
  }
  cbind(i=ii, j=jj, distance=f(ii,jj))
}


## Compute relevent pairwise distances
##
## Version for use with binning optimization
dist_mat <- function(prms, pop) {
  ## individuals compared to individuals
  dists <- data.frame(i=integer(), j=integer(), distance=double())

  ## patch_k_dists: returns distances relevent to patch k
  patch_k_dists <- function(k) {
    locals <- bin_lookup_list[[k]]
    if(is.null(locals)) return(NULL)
    comp_inds <- unlist(bin_lookup_list[prms$adjacencies[[k]]])
    if (is.null(comp_inds)) comp_inds <- integer(0)


    n_within  <- length(locals)
    n_between <- length(comp_inds)

    ii_within <- locals[rep.int(1:n_within, 1:n_within)]
    jj_within <- locals[sequence(1:n_within)]

    ii_between <- c(rep.int(locals, n_between))
    jj_between <- c(rep.int(comp_inds, rep.int(n_within, n_between)))

    ii <- c(ii_within, ii_between)
    jj <- c(jj_within, jj_between)

    return(dist_wa(prms, pop[, c('x','y')], ii, jj))
  }

  ## construct list to quickly identify which individuals are in a
  ## particular bin
  bin_lookup_list <- vector('list', length=prms$o_dim^2)
  vals <- split(1:nrow(pop),pop[,'bin'])
  bin_lookup_list[as.numeric(names(vals))] <- vals

  dists_by_k <- lapply(1:(prms$o_dim^2), patch_k_dists)
  dists <- do.call(rbind, dists_by_k)
  keep <- which(dists[,'distance']<prms$crit_distance)
  dists[keep,]
}
##
## Original version
dist_mat_alt <- function(prms, pop) {
  dists <- dist_wa_alt(prms, pop[,c('x','y')])
  keep <- which(dists[,'distance']<prms$crit_distance)
  dists[keep,]
}

## compute tau_i (Equation 3)
compute_tau_i <- function(prms, pop, distances) {

##  if(prms$land_type=="uniform" | prms$land_type=="one_gauss") {
##    kxy.non.normalized <- function(x, y)
##      one_gauss(x, y, sigma=0.1, het=prms$het, s_dim=1)
## }
##  else if(prms$land_type=="noisy"){
##    kxy.non.normalized <- function(x, y) noisy(x, y, mat=prms$land_mat, res=128)
##  }
##  else if(prms$land_type=="striped"){
##    kxy.non.normalized <- function(x, y) {
##      striped(x, y, stripes=prms$stripes, s_dim=prms$s_dim)
##    }
##  }
  ## normalize the landscape
##  kxy <- normalize_f(kxy.non.normalized,
##                     prms$s_dim, prms$integrate_to)

  ## changed from dnorm b/c that gave -ve values

  denom <- 2 * pi * prms$sigma_s^2
  num <-  distances[,'distance'] / prms$sigma_s
  comp_coeffs <- exp(-0.5 * num^2 ) / denom

  ## comp_sums denotes the summations in Equation 3
  comp_sums <- sapply(split(c(comp_coeffs, comp_coeffs),
                            c(distances[,'i'], distances[,'j'])),
                      sum) - 1/denom

  ## the last subtracted term above is to account for double
  ## counting individual's effects on themselves.

  rho_i <- prms$kxy(pop[,'x'], pop[,'y']) / comp_sums ## Equation 3
  tau_i <- rho_i / (prms$c + rho_i) ## Equation 4 (Holling type II)
  tau_i
}

## birth
birth <- function(prms, pop, distances) {
  tau_i <- compute_tau_i(prms, pop, distances)

  f_i <- prms$f_max * tau_i # (Equation 7, without gamma)

  ## number of offspring per individual
  num <- rpois(nrow(pop), f_i)
  ## cat("Num offspring: ", sum(num), "\n")

  if(sum(num)<20) return("extinction")

  parent <- rep(1:length(num), num)
  pop[parent,]
}

## mutation
mutation <- function(prms, pop, trait_col){

  if(trait_col == "diap_length"){
    ## who mutates * sampling mutation amount; possible mutations are -3,-2,-1,1,2,3 generations. Probabilities decay in an approximately gaussian way. Because of limitations on the minimum value, not all mutants will show a phenotypic change.
    mut_array <- rbinom(nrow(pop), 1, prob=prms$mut_rate)*sample(c(-3,-2,-1,1,2,3), nrow(pop), prob=c(0.047,0.27,0.68, 0.68, 0.27,0.047), replace=TRUE)
    pop[,trait_col] <- pop[,trait_col] + mut_array
    pop[,trait_col] <- pmax(1, pop[,trait_col])
  }
	
  if(trait_col == 'p1' |
     trait_col == 'p2' |
     trait_col == 'diap_prob') {

    ## draw and add mutational effects
    mut_array <- rnorm(nrow(pop), mean=0, sd=prms$mut_rate)
    pop[,trait_col] <- pop[,trait_col] + mut_array

    if(trait_col == 'p2' |
       trait_col == 'diap_prob') {
      pop[,trait_col] <- pmin(1, pop[,trait_col])
    }

  }

  ## all traits need to be positive:
  pop[,trait_col] <- pmax(0, pop[,trait_col])    	  
  pop
}

## movement (standard wrap-around boundaries)
move <- function(prms, pop) {
  d <- rnorm(nrow(pop), mean=0, sd=pop[,'p1'])
  move <- matrix(rbinom(n=nrow(pop), size=1, prob=pop[,'p2']))
  moving_pop <- pop[which(move==1),]
  stay_pop <- pop[which(move==0),]
  d <- d[which(move==1)]

  a <- runif(nrow(moving_pop))*2*pi

  moving_pop[,'x'] <- (moving_pop[,'x'] + cos(a) * d)%%prms$s_dim
  moving_pop[,'y'] <- (moving_pop[,'y'] + sin(a) * d)%%prms$s_dim
  pop <- rbind(moving_pop, stay_pop)
  pop[,'bin'] <- prms$o_space[bins(pop[,'x'], pop[,'y'], prms)]

  pop
}

diapause <- function(prms, pop) {

  ## MC: first, split the pop by who is currently in diapause to avoid
  ## immediate re-entry of emerging diapausers into diapause (unless
  ## we decide we want that option)
  old_diap_vec<-pop[,'current_diap']>0
  old_diapausers <- pop[old_diap_vec,, drop=FALSE]
  active_diap_vec<-pop[,'current_diap']==0
  active <- pop[active_diap_vec,, drop=FALSE]

  ## MC: shift these up one time-step
  old_diapausers[, 'current_diap'] <- old_diapausers[,'current_diap']-1

  ## MC: identify new diapausers
  diapause <- rbinom(n=nrow(active), size=1, prob=active[,'diap_prob'])
  non_diapausers <- active[diapause==0,, drop=FALSE]
  new_diapausers <- active[diapause==1,, drop=FALSE]

  ## MC: how long will new diapausers diapause for? As a first pass,
  ## keep it strictly their diap_length with no variation.
  diapause_lengths <- new_diapausers[,'diap_length']

  ## MC: given that they diapause, must be at least one timestep; not
  ## necessary now but could be later if randomization is implemented.
  diapause_lengths[diapause_lengths<1]<-1

  ## MC: update current_diaps for new_diapausers
  new_diapausers[,'current_diap']<-diapause_lengths

  pop <- rbind(old_diapausers, non_diapausers, new_diapausers)

  pop
}

death <- function(prms, pop){
  if(prms$mortality==0) return(pop)
  else {
    deaths <- rbinom(n=nrow(pop), size=1, prob=(1-prms$mortality))
    survivors <- deaths==1
    return(pop[survivors,,drop=FALSE])
  }	
}
