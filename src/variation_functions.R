## ******************************************************************
## ******************** normalize the carrying capacity *************
## ******************************************************************

## normalize carrying capacity using functional numeric integration
normalize_f <- function(kxy, s_dim, integrate_to=1) {
  ## create normalized carrying capacity
  integrand <- function(x) kxy(x[1], x[2])
  kxy_integral <- adaptIntegrate(integrand, lowerLimit=c(0,0), upperLimit=c(s_dim, s_dim))$integral
  function(x, y) (kxy(x, y)/kxy_integral) * integrate_to
}

## A single gaussian peak in the middle of the landscape
## changed from a copy of gaussian() to implement a more intuitive
## version of changing het
one_gauss <- function(x, y, sigma, het, s_dim) {
  if(het==0) return(rep(1,length(x)))

  base <- (1-2*(1+het)*exp(-1/(32*sigma^2)))/het
  center <- s_dim/2
  G <- base + exp((-1)*((x-center)^2+(y-center)^2)/(2*sigma^2))*het
  G
}

## acl : linear indication of autocorrelation
##      -> acl == 0.142 is max autocorrelation
## res : resolution of matrix
## amp : proportional to difference between high and low peaks
noisy <- function(x, y, mat, res) {
  ## land_mat <- cn_2D(acl=acl, n=res, amp=amp)$cn
  row <- ceiling(x*res)%%(res+1)
  col <- ceiling(y*res)%%(res+1)
  row[which(row==0)] <- res
  col[which(col==0)] <- res
  #browser()
  noisy_fn <- mat[cbind(row, col)]
  noisy_fn
}

## stripes : number of equal width stripes
striped <- function(x, y, stripes=1, s_dim){
  if(stripes==0) return(1)
  bounds <- s_dim/(stripes*2) #0--> <bounds = uninhabitable, bounds--> <2*bounds = habitable, etc
                                        
  k<-1+floor(x/bounds)%%2
  k
}


