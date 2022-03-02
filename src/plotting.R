## ************* Plot sigmas coloured in time ***************
plot_trait <- function(prms, pop, trait) {
  ## range of trait values observed
  plot_pop <- function(trait_vals, s, gen) {
    cols_dis <- colorRampPalette(brewer.pal(12, "Spectral"))
    plot(pop[,'x']/prms$s_dim, pop[,'y']/prms$s_dim,
         xlim=c(0,1), ylim=c(0,1),
         xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', pch=20,
         col=cols_dis(length(pop[,'x'])), cex = 0.75)
    axis(1, at=c(0:3)/3, yaxt = 'n', xaxt = 'n', labels=FALSE)
    axis(2, at=c(0:2)/2, yaxt = 'n', xaxt = 'n', labels=FALSE)
    title(main=s)
  }
  par(mar=c(.5, 1, 1, .5))
  plot_pop(pop[,"p1"], nrow(pop), i)
}
## *********************************************************

## heatmap landscape plot
plot_hm <- function(prms, i, pop){
	##generate prms$kxy; repeated code from compute_tau_i, but mysteriously blows up if it's in the actual prms list.
	
  if(prms$land_type=="uniform" | prms$land_type=="one_gauss") {
    kxy.non.normalized <- function(x, y)
      one_gauss(x, y, sigma=0.1, het=prms$het, s_dim=1)
  }
  else if(prms$land_type=="noisy"){
    kxy.non.normalized <- function(x, y) noisy(x, y, mat=prms$land_mat, res=128)
  }
  else if(prms$land_type=="striped"){
    kxy.non.normalized <- function(x, y) {
      striped(x, y, stripes=prms$stripes, s_dim=prms$s_dim)
    }
  }
  ## normalize the landscape
  kxy <- normalize_f(kxy.non.normalized,
                     prms$s_dim, prms$integrate_to)
	
  cols_dis <- colorRampPalette(c("white", "black"))

  x <- y <- seq(0, prms$s_dim, length=257)
  landscape_mat <- outer(x, y, kxy)
  colnames(landscape_mat) <- rownames(landscape_mat) <- x

  image(landscape_mat, col=cols_dis(12), xaxt='n', yaxt='n')

	if(!is.na(pop)){
  		cols_dis <- colorRampPalette(c("blue", "red"))
		points(pop[,'x'], pop[,'y'], pch=20,
    	    col=cols_dis(length(pop[,'x'])), cex=0.5)
    	  }
  title(main=i)
}

## points and lines in a specified colour
plot_points <- function(xx, yy, colour){
  lines(xx, yy, col=colour, lwd=2)
  points(xx, yy, pch=16, cex=0.3)
}

## NOT FUNCTIONAL
plot_val_comb <- function(directory, colour){
  files <- list.files(directory)
  files <- files[c(TRUE, FALSE)]

  for (i in 1:length(files)) {
    curr_file <- files[i]
    file <- paste(directory, curr_file, sep='/')
    load(file, verbose=TRUE)
    ## plot_points(pop_stats[,'mean_trait_p1'], colour=colour)
  }
}

##
plot_boxes <- function(directory, rows, num_gens, colour){
  files <- list.files(directory)
  files <- files[c(TRUE,FALSE)] ## only list pop_stat files

  get_vals <- function(ii) {
    curr_file <- files[ii]
    file <- paste(directory, curr_file, sep='/')
    load(file, verbose=TRUE)
    pop_stats[,'mean_trait_p1'][rows]
  }

  vals <- sapply(1:length(files), get_vals)
  means <- rowMeans(vals)
  stde <- apply(vals, 1, sd)
  make_errorbars(rows, means, stde, num_gens, col=colour)
}

## x - values run
## y - means
## stde - sd of means
make_errorbars <- function(x, y, stde, gens, axes=c(FALSE,FALSE), col, ...) {
  plot_points(x, y, col)
  if(axes[1])
    axis(1, at=x, labels=x, las=0)
  if(axes[2])
    axis(2, at=(0:10/10), labels=(0:10/10), las=2)

  xx <- c(x[1], x, rev(x[-1]))
  yy <- c((y-stde)[1], y+stde, rev((y-stde)[-1]))
  polygon(xx, yy, density=40, col=col)
  ## arrows(xx, yy+stde, xx, yy-stde, length=0.05, angle=90, code=3)
}

plot_reps <- function(directory, print_every, num_gens) {
  files <- list.files(directory)

  ## if run was completed, don't iterate over the val_list files
  if (file.exists("~/val_list.RData")) files <- files[-length(files)]

  num_files <- length(files)
  rows <- which((1:num_gens)%%print_every==0)

  cols_dis <- colorRampPalette(brewer.pal(8, "Dark2"))
  colour <- cols_dis(num_files)

  layout(matrix(c(rep(1,20), 2, rep(3,3), 2, rep(3,3)), nrow=4, ncol=7))
  plot(NA, pch=16, xlab='Generation', main=directory,
       ylab='Mean sigma_m', ylim=c(0, 2), xlim=c(1, num_gens), las=1)

  for (i in 0:(num_files-1))
    plot_boxes(paste(directory, i, sep='/'), rows, num_gens, colour[i+1])

  ## REMEMBER: change substr length when directory name changes
  if(grepl('het', directory)) {
    val <- as.numeric(substr(directory, start=28, stop=100000L))
    fn <- function(x,y) one_gauss(x, y, sigma=0.08, het=val, s_dim=1)
  } else if (grepl('sig', directory)) {
    val <- as.numeric(substr(directory, start=28, stop=100000L))
    fn <- function(x,y) one_gauss(x, y, sigma=val, het=2, s_dim=1)
  }

  oldpar <- par(mai=c(0.01,0.01,0.01,0.01))
  plot(0,0, type="n", ann=FALSE, axes=FALSE)
  legend(-0.7, 0, legend=(0:(num_files-1)), col=colour, cex=2, lwd=2, lty=1)
  par(oldpar)

  oldpar <- par(mai=c(0.01,0.3,0.2,0.2))
  plot_landscape_3D(fn)
  par(oldpar)

  ## legend("topright", inset=0.05, legend=(0:(num_files-1)),
  ##       col = colour, cex = 0.8, lwd = 1, lty = 1)

  ## reset layout
  layout(matrix(1, nrow=1, ncol=1))
}
## *********************************************************

plot_landscape_3D <- function(landscape, nlines=100) {
  ## landscape <- normalize_f(landscape, 1, 1)
  x <- y <- seq(0, 1, length=nlines)

  z <- outer(x, y, landscape)
  persp(x, y, z, zlim=c(0,10), theta=30, phi=5, expand=1,
        col="lightblue", ltheta=120, shade=0.75, border=NA,
        ticktype="detailed", nticks=2, xlab="X", ylab ="Y",
        zlab="Z")

  title('Landscape', line=-4, cex.main=2)
}
