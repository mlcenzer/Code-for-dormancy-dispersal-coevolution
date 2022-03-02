## ******************************************************************
## ******************** basic parameter set *************************
## ******************************************************************
base_prms <- function(land_type='striped',
                      stripes=2,
                      acl=0.1,
                      c=1/15,
                      ext_width=1.0,
                      file_path='model_output',
                      f_max=200,
                      het=2,
                      land_mat=matrix(0, 1, 1),
                      mut_rate=0.01,
                      N=50,
                      num_ext=5,
                      num_gens=1e4,
                      num_traits=4, ##MC: Update for diapause and environmental sensitivity
                      print_every=1e1,
                      num_saves=10,
                      s_dim=1,
                      sig=0.8,
                      sigma_s=0.05,
                      init_p1=0.05,
                      init_diap=0,
                      mortality=0,
                      plot=FALSE) {
  inputs <- as.list(environment())
  prms <- prms_update(inputs)
  prms[order(names(prms))] ## order for consistency
}

## DESCRIPTION OF PARAMETERS
##
## num_traits: number of traits
##
## sigma_s: spatial radius of competition
## crit_distance: critical interaction distance
## N: number of offspring
##
## c: Strength of competition for resources
## f_max: maximum fecundity
##
## *** environmental parameters ***
## s_dim: model space dimensions
## kxy: function specifying variation in carrying capacity
## het: heterogeneity of simulation environment
##
## *** simulation parameters ***
## num_gens: number of generations
## print_every: prints pop stats every print_every generations
## file_path: name of folder to save output to
## mut_rate: sd of gaussian curve from which mutations are sampled
##
## *** dependent parameters ***
## integrate_to: total amount of resources
## o_dim: overlay dimension
## o_space: overlay space - matrix-lookup for bin numbers
## adjacencies: matrix of adjacencies for each cell of overlay space
## crit_distance: critical interaction distance
## adjacencies: list of four adjacent bins for each bin in the o-space
## ******************************************************************

## generate dependent parameters
prms_update <- function(prms.in) {
  ## if(prms.in$x_len != prms.in$y_len)
  ##   cat('x_len != y_len! Bins won\'t work!\n')

  prms.in$integrate_to <- prms.in$s_dim^2

  ## individual with coordinates (x,y) is at grid
  ## cell floor((x,y)/crit_dist) in o-space (overlay space)
  prms.in$o_dim <- floor(prms.in$s_dim/(3*prms.in$sigma_s))
  prms.in$o_space <- matrix(1:(prms.in$o_dim^2), nrow=prms.in$o_dim)

  ## critical interaction distance
  prms.in$crit_distance <- prms.in$s_dim/prms.in$o_dim

  ## adjacencies returns a list of four adjacent bins for each bin in
  ## the o-space.
  adjacencies <- function(o_space=prms.in$o_space) {
    total <- prms.in$o_dim^2
    mod <- function(x, y)
      ifelse(x%%y==0, y, x%%y)
    get_below <- function(ii)
      floor((ii-1)/prms.in$o_dim)*prms.in$o_dim + ii%%(prms.in$o_dim) + 1
    get_beside1 <- function(ii)
      mod((ii+prms.in$o_dim)%%(total), total)
    get_beside2 <- function(ii)
     mod((ii-prms.in$o_dim)%%(total), total)

    get_adjacent <- function(x)
      c(below=get_below(x),
        beside=get_beside1(x),
        diag1=get_below(get_beside1(x)),
        diag2=get_below(get_beside2(x)))

    all_pairs <- lapply(o_space, get_adjacent)

    ## all pairs of patches to be "compared"
    ## pairs <- cbind(i=rep(o_space, each=3), j=unlist(all_pairs))
    return(all_pairs)
  }

  ## adjacency list
  prms.in$adjacencies <- adjacencies()
  
  return(prms.in)
}

update_land_prms<-function(prms.in){
	  ##trying to add back landscape normalization to update_prms
  
  if(prms.in$land_type=="uniform" | prms.in$land_type=="one_gauss") {
    kxy.non.normalized <- function(x, y)
      one_gauss(x, y, sigma=0.1, het=prms.in$het, s_dim=1)
  }
  else if(prms.in$land_type=="noisy"){
  	   mat <- cn_2D(acl=prms.in$acl, n=128, amp=1)$cn
       prms.in$land_mat <- abs(min(mat)) + mat
       prms.in$land_max <- max(prms.in$land_mat)
    kxy.non.normalized <- function(x, y) noisy(x, y, mat=prms.in$land_mat, res=128)
  }
  else if(prms.in$land_type=="striped"){
    kxy.non.normalized <- function(x, y) {
      striped(x, y, stripes=prms.in$stripes, s_dim=prms.in$s_dim)
    }
  }
  ## normalize the landscape
  prms.in$kxy <- normalize_f(kxy.non.normalized,
                     prms.in$s_dim, prms.in$integrate_to)

	return(prms.in)
}

## construct bin lookup
bins <- function(x, y, prms) {
  cbind(ceiling(x/prms$crit_distance), ceiling(y/prms$crit_distance))
}
