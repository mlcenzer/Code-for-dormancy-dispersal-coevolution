####### All of these require that you have stored data summaries that are filed in a particular format (as they are in the included ./model_output folder). The folders must be labelled "noisy" and "striped" and the seven folders in each must maintain their names/order to properly generate these figures. However, a savvy tinkerer could make these work for other organizations by modifying the code here and in the plot_functions_local.R script.

setwd("~/Documents/GitHub/Code-for-dormancy-dispersal-coevolution")
source('src/initialize.R')
source('MC plotting/plot_functions_local.R')
library(viridis)
library(stringr)
library(autoimage)
library(gridBase)
library(grid)


####Figure 1
write_diap_summary_figure_four_panel <- function(write_to="./figures", load_from="./model_output", legend_text=c("0.01","0.05","0.50"), stripes=c(1:16), outlier=NA, zero=NA, width=F){

	if(!is.na(outlier)) {stripes_all <- c(stripes, outlier)}
  	if(!is.na(zero)) {stripes_all <- c(zero, stripes_all)}
  	
  	diap_list <- make_summary_multi(variable="diap", stripes=stripes_all, load_from=load_from)

	disp_list <- make_summary_multi(variable="disp", stripes=stripes_all, load_from=load_from)
	
	correlation_list <- make_summary_corr(land_var=stripes_all, load_from="./model_output/striped/disp_ev_T_0.05_init_diap_0.00_100k", acl=F, stripes=T)
	
  
  #write to pdf
  setwd(write_to)
  pdf(file="coevolution_4_panel_width.pdf", width=5.5, height=15)

 	split.screen(rbind(c(0.16, 0.44, 0.87, 1), c(0.44, 0.72, 0.87, 1), c(0.72, 1, 0.87, 1), c(0, 1, 0.57, 0.87), c(0, 1, 0.21, 0.57), c(0, 1, 0, 0.21)))
	
  if(width==F) screen(1)
  else if(width==T) screen(3)
 	par(mai=c(.2, .05, .3, .1), ps=18) #(bottom,left,top, right)
  	prms_1 <- base_prms(land_type = "striped", stripes=1, het=2)
	plot_hm(prms_1, i="", pop=NA)
	if(width==F) mtext("A. Landscapes", adj=0.01)

  screen(2)
	par(mai=c(.2, .05, .3, .1), ps=18) #(bottom,left,top, right)
  	prms_2 <- base_prms(land_type = "striped", stripes=8, het=2)   
  	plot_hm(prms_2, i="", pop=NA)

  if(width==F) screen(3)
  else if(width==T) screen(1)
	par(mai=c(.2, .05, .3, .1), ps=18) #(bottom,left,top, right)
  	prms_3 <- base_prms(land_type = "striped", stripes=16, het=2)
  	plot_hm(prms_3, i="", pop=NA) 
	if(width==T) mtext("A. Landscapes", adj=0.01)
	    	
 screen(5)
   par(mai=c(.8, .9, .55, .1), ps=18) #(bottom,left,top, right)
  plot_multi_stripes(summary=diap_list[[4]], ylim=c(0,0.55), ylab="Dormancy probability", header="C. Coevolution - dormancy", legend=F, outlier=outlier, zero=zero, width=width)
  lines(x=c(1,1), y=c(-0.5,2), lty=3)
  lines(x=c(9,9), y=c(-0.5,2), lty=3)
  lines(x=c(16,16), y=c(-0.5,2), lty=3)
  mtext("\\\\", side=1, line=-0.3, at=c(0))
  if(!is.na(zero)) lines(x=c(17.5,17.5), y=c(-0.5, 2), lty=1)
  if(width==F) axis(side=1, at=19, labels="32")
  else if(width==T) axis(side=1, at=c(-2, 19), labels= c("0.016", "H"))

  screen(4)
  par(mai=c(.3, .9, .2, .1), ps=18) #(bottom,left,top, right)    
   plot_multi_stripes(summary=disp_list[[4]], ylim=c(0,1.7), ylab="Dispersal distance", header="B. Coevolution - dispersal", legend=F, outlier=outlier, zero=zero, width=width)
  lines(x=c(1,1), y=c(-0.5,2), lty=3)
  lines(x=c(9,9), y=c(-0.5,2), lty=3)
  lines(x=c(16,16), y=c(-0.5,2), lty=3)
  mtext("\\\\", side=1, line=-0.3, at=c(0))
  if(!is.na(zero)) lines(x=c(17.5,17.5), y=c(-0.5, 2), lty=1)
  if(width==F) axis(side=1, at=19, labels="32")
  else if(width==T) axis(side=1, at=c(-2, 19), labels= c("0.016", "H"))

  
  screen(6)
  	par(mai=c(.8, .9, .2, .1), ps=18) #(bottom,left,top, right)
	color = c(rgb(0.5,0.5,0.5, alpha=0.4), viridis((length(stripes_all)-2), alpha=0.4), rgb(1,0,0,alpha=0.4))
  
  	plot(NA, xlim=c(0.05,2.45), yaxt='n', ylim=c(0,0.6),
       xlab="Dispersal distance",
       ylab="Dormancy probability", las=1)
       mtext("D. Coevolution trade-off", adj=0.01)
       axis(side=2, at=seq(0:0.6, by=0.2), las=1)
  	for(var in stripes_all) {
    	name <- as.character(var)
    	norm_disp<- correlation_list[name][[1]][,'dispersal']                               	
    	points(correlation_list[name][[1]][,'diapause']~norm_disp, col=color[which(stripes_all==var)], pch=19)
  	}
  
  close.screen(all.screens=TRUE)
  dev.off()

  }

write_diap_summary_figure_four_panel(stripes=1:16, outlier=32, zero=0, load_from="./model_output/striped", width=T)






####Figure 2
write_kin_comp_summary_figure <- function(write_to="./figures", load_from="./model_output/striped", stripes, outlier=F, zero=F){
  
  folders_list=list.files(path=load_from)
  tau_change_list=list()
  list_num <- 1
  for(folder in folders_list){
    all_files<-list.files(path=paste(load_from, "/", folder, sep=""))	
    tau_change_summary<-list()
    file_num <- 1
    for(file in all_files){
      load(paste(load_from, "/", folder, "/", file, sep=""))
#zero case has pop_size in the list
	  if(file_num==1)	tau_change_summary[[file_num]] <- summary_list[[4]]
      else tau_change_summary[[file_num]] <- summary_list[[3]]
      names(tau_change_summary)[[file_num]]<-strsplit(strsplit(file, "stripes_")[[1]][2], "_")[[1]][1]
      file_num <- file_num + 1
    }
    tau_change_list[[list_num]] <- tau_change_summary
    names(tau_change_list)[[list_num]] <- folder
    list_num <- list_num + 1
  }

  setwd(write_to)
  pdf(file="kin_comp_summary.pdf", width=12, height=10.5)
  split.screen(rbind(c(0,0.53,0.53,1), c(0.53,1,0.53,1), c(0,0.53,0,0.53), c(0.53,1,0,0.53)))
  
  par(mai=c(0.5, 1.1, 0.4, 0), ps=18)
  screen(1)
  plot_kin_comp_tau(summary_list=tau_change_list[[1]], stripes_all=stripes, disp_ev=F, header=expression(paste("A. Dispersal fixed (", sigma["m"], "=0.01)")), ylim=c(0,0.009), outlier=outlier, zero=zero, ylab=T)
  
  screen(2)
  par(mai=c(0.5, 0.4, 0.4, 0.1))
  plot_kin_comp_tau(summary_list=tau_change_list[[2]], stripes_all=stripes, disp_ev=F, header=expression(paste("B. Dispersal fixed (", sigma["m"], "=0.05)")), ylim=c(0,0.009), outlier=outlier, zero=zero, yaxt="n", ylab=F)
	axis(side=2, labels=F)
	  add.gray.scale <- function() {
    	vp <- baseViewports()
    	pushViewport(vp$inner,vp$figure,vp$plot)
    	## adjust x, y, and height to move location around, etc
    	pushViewport(viewport(x=0.91, y=0.85, width=0.04, height=.4,
        	                  just=c("left","top")))
    	op <- par(plt=gridPLT(),new=T)
    	on.exit(par(op))
    	cols <- c(viridis(length(stripes)))
   	 	if(!is.na(outlier)) cols <- c(viridis(length(stripes)-1), rgb(1,0,0))
   	 	if(!is.na(zero)) cols <- c(rgb(0.5,0.5,0.5), viridis(length(stripes)-2), rgb(1,0,0))
    	image(matrix(1:17, ncol=17),
        	  col=cols, xaxt="n", yaxt="n")
    	popViewport(4)
  	}    
  add.gray.scale()
  ## can then add text manually
  text(0.016, x=-1, y=1, pos=2, cex=0.7)
  text(0.5, x=-1, y=0.05, pos=2, cex=0.7)
  text(x=0, y=1.22, expression(frac(1,"2n"["p"])), cex=0.8)
  
  screen(3)
  par(mai=c(0.9, 1.1, 0.4, 0))
  plot_kin_comp_tau(summary_list=tau_change_list[[3]], stripes_all=stripes, disp_ev=F, header=expression(paste("C. Dispersal fixed (", sigma["m"], "=0.50)")), ylim=c(0,0.009), outlier=outlier, zero=zero, ylab=T)

  screen(4)
  par(mai=c(0.9, 0.4, 0.4, 0.1))
  plot_kin_comp_tau(summary_list=tau_change_list[[4]], stripes_all=stripes, disp_ev=T, header="D. Coevolution", ylim=c(0,0.009), outlier=outlier, zero=zero, yaxt="n", ylab=F)
  	axis(side=2, labels=F)
  close.screen(all.screens=TRUE)
  dev.off()
}

write_kin_comp_summary_figure(stripes=c(0:16, 32), outlier=32, zero=0)





####Figure 3
write_fixed_disp_figure <- function(write_to="./figures", load_from="./model_output/striped", legend_text_1, legend_text_2, stripes, outlier, zero, width){
	
	if(!is.na(outlier)) {stripes_all <- c(stripes, outlier)}
	if(!is.na(zero)) {stripes_all <- c(zero, stripes_all)}
  	
  	diap_list <- make_summary_multi(variable="diap", stripes=stripes_all, load_from=load_from)

	disp_list <- make_summary_multi(variable="disp", stripes=stripes_all, load_from=load_from)#browser()

	#write to pdf
  	setwd(write_to)
	#browser()
	pdf(file="fixed_stripes_two_panel_width.pdf", width=5.5, height=9.5)
	split.screen(rbind(c(0,1,0.53,1), c(0,1,0,0.53)))
	   screen(1)
	  par(mai=c(.3, .9, .5, .1), ps=18) #(bottom,left,top, right)
   			plot_multi_stripes(summary=disp_list[[5]], summaries_list=disp_list[6:7], ylim=c(0,1.8), ylab="Dispersal distance", header="", legend=T, legend_text=legend_text_2, outlier=outlier, zero=zero, legend_title=expression(paste(italic("q")["f"])), width=width)
  			 mtext("A. Dispersal with fixed dormancy", adj=0.01)	
  		lines(x=c(1,1), y=c(-0.5,2), lty=3)
  		lines(x=c(9,9), y=c(-0.5,2), lty=3)
  		lines(x=c(16,16), y=c(-0.5,2), lty=3)
 		if(!is.na(zero)) lines(x=c(17.5,17.5), y=c(-0.5, 2), lty=1)
 		if(width==F) {axis(side=1, at=19, labels="32")
  			mtext("\\\\", side=1, line=-0.3, at=17) 
  		}			
		else if(width==T){			
			 axis(side=1, at=c(-2, 19), labels= c("0.016", "H"))
 			 mtext("\\\\", side=1, line=-0.3, at=c(0))
 		}
  	
  		screen(2)
	   	par(mai=c(.8, .9, .5, .1), ps=18) #(bottom,left,top, right)
  		plot_multi_stripes(summary=diap_list[[1]], summaries_list=diap_list[2:3], ylim=c(0,0.55), ylab="Dormancy probability", header="", legend=T, legend_text=legend_text_1, outlier=outlier, zero=zero, legend_title = expression(paste(sigma["m"])), width=width)
  			mtext("B. Dormancy with fixed dispersal", adj=0.01)
  		lines(x=c(1,1), y=c(-0.5,2), lty=3)
  		lines(x=c(9,9), y=c(-0.5,2), lty=3)
  		lines(x=c(16,16), y=c(-0.5,2), lty=3)
 		if(!is.na(zero)) lines(x=c(17.5,17.5), y=c(-0.5, 2), lty=1)
 		if(width==F) {axis(side=1, at=19, labels="32")
  			mtext("\\\\", side=1, line=-0.3, at=17) 
  		}			
		else if(width==T){			
			 axis(side=1, at=c(-2, 19), labels= c("0.016", "H"))
 			 mtext("\\\\", side=1, line=-0.3, at=c(0))
 		}
  
	close.screen(all.screens = TRUE)
	dev.off()
}
write_fixed_disp_figure(stripes=1:16, outlier=32, zero=0, legend_text_1 = c("0.01", "0.05", "0.5"), legend_text_2=c("0.1","0.3", "0.5"), width=T)







####Figure 4
write_diap_summary_figure_four_panel_noisy <- function(write_to="./figures", load_from="./model_output/noisy", acl=c(0.001, 0.01, 0.1), land_type="noisy"){
  	#browser()
  diap_list <- make_summary_multi(variable="diap", acl=acl, load_from=load_from, stripes=NA)

 disp_list <- make_summary_multi(variable="disp", acl=acl, load_from=load_from, stripes=NA)
 
   correlation_list <- make_summary_corr(land_var=acl, load_from="./model_output/noisy/acl_runs", acl=T, stripes=F)
     acl <- sort(log10(acl), decreasing=T)  
      #browser()
  #write to pdf
  setwd(write_to)
  pdf(file="coevolution_noisy_4_panel.pdf", width=5.5, height=15)
 	split.screen(rbind(c(0.16, 0.44, 0.87, 1), c(0.44, 0.72, 0.87, 1), c(0.72, 1, 0.87, 1), c(0, 1, 0.57, 0.87), c(0, 1, 0.21, 0.57), c(0, 1, 0, 0.21)))
   
   screen(1)
	par(mai=c(.2, .05, .3, .1), ps=18) #(bottom,left,top, right)
  	prms_3 <- base_prms(land_type = land_type, acl=0.001)
  	if(land_type=="noisy"){
  		#browser()
  	   mat <- cn_2D(acl=prms_3$acl, n=128, amp=1)$cn
       prms_3$land_mat <- abs(min(mat)) + mat
       prms_3$land_max <- max(prms_3$land_mat)	
  	}

  	plot_hm(prms_3, i="", pop=NA) 
	mtext("A. Landscapes", adj=0.01)
 
	
   screen(2)
	par(mai=c(.2, .05, .3, .1), ps=18) #(bottom,left,top, right)
  	prms_2 <- base_prms(land_type = land_type, acl=0.01) #acl=acl[round(length(acl)/2)])
  	if(land_type=="noisy"){
  		#browser()
  	   mat <- cn_2D(acl=prms_2$acl, n=128, amp=1)$cn
       prms_2$land_mat <- abs(min(mat)) + mat
       prms_2$land_max <- max(prms_2$land_mat)	
  	}
  
  	plot_hm(prms_2, i="", pop=NA)

 screen(3)

 	par(mai=c(.2, .05, .3, .1), ps=18) #(bottom,left,top, right)
  	prms_1 <- base_prms(land_type = land_type, acl=0.1)
  	if(land_type=="noisy"){
  		#browser()
  	   mat <- cn_2D(acl=prms_1$acl, n=128, amp=1)$cn
       prms_1$land_mat <- abs(min(mat)) + mat
       prms_1$land_max <- max(prms_1$land_mat)	
  	}

  	plot_hm(prms_1, i="", pop=NA)
	  	  
  screen(4)
  par(mai=c(.3, .9, .2, .1), ps=18) #(bottom,left,top, right)
  plot_multi_noisy(summary=disp_list[[1]], ylim=c(0,1.6), ylab="Dispersal distance", header="B. Coevolution - dispersal", legend = F)
  	lines(x=c(log10(0.001),log10(0.001)), y=c(-0.5,2), lty=3)
  	lines(x=c(log10(0.01),log10(0.01)), y=c(-0.5,2), lty=3)
  	lines(x=c(log10(0.1),log10(0.1)), y=c(-0.5,2), lty=3)
  	
  screen(5)
   par(mai=c(.9, .9, .5, .1), ps=18) #(bottom,left,top, right)
  plot_multi_noisy(summary=diap_list[[1]], ylim=c(0,0.4), ylab="Dormancy probability", header="C. Coevolution - dormancy", legend = F)  	
    lines(x=c(log10(0.001),log10(0.001)), y=c(-0.5,2), lty=3)
  	lines(x=c(log10(0.01),log10(0.01)), y=c(-0.5,2), lty=3)
  	lines(x=c(log10(0.1),log10(0.1)), y=c(-0.5,2), lty=3)
  
  	
  screen(6) 	 
  	par(mai=c(.8, .9, .2, .1), ps=18) #(bottom,left,top, right)
	color = viridis(length(acl), alpha=0.4)
  	plot(NA, xlim=c(0.05,2.45), yaxt='n', ylim=c(0,0.65),
       xlab="Dispersal distance",
       ylab="Dormancy probability", las=1)
       axis(side=2, at=seq(0:0.6, by=0.2), las=1)
 	mtext("D. Coevolution tradeoff", adj=0.01)      
  	for(var in acl) {
    	name <- as.character(var)
    	norm_disp<- correlation_list[name][[1]][,'dispersal']                               	
    	points(correlation_list[name][[1]][,'diapause']~norm_disp, col=color[which(acl==var)], pch=19)
  	}	 
 	
  close.screen(all.screens=TRUE)
  dev.off()
  }

write_diap_summary_figure_four_panel_noisy(acl=c(0.001, 0.0014, 0.0019, 0.0025, 0.0035, 0.005, 0.006, 0.0075, 0.0082, 0.0091, 0.010, 0.013, 0.016, 0.0200, 0.025, 0.035, 0.05, 0.075, 0.100))












####Supplemental Figure 1
setwd("~/Dropbox/diapause_ibm/")
write_fixed_disp_figure <- function(write_to="./figures", load_from="./model_output/noisy", acl=c(0.001, 0.01, 0.1), land_type="noisy", legend_text_1, legend_text_2){
	diap_list <- make_summary_multi(variable="diap", acl=acl, load_from=load_from, stripes=NA)
	disp_list <- make_summary_multi(variable="disp", acl=acl, load_from=load_from, stripes=NA)
	#write to pdf
  	setwd(write_to)
	pdf(file="fixed_noisy_two_panel.pdf", width=5.5, height=8.5)
	split.screen(rbind(c(0,1,0.53,1), c(0,1,0,0.53)))
	   par(mai=c(.3, .9, .5, .1), ps=18) #(bottom,left,top, right)
	   screen(1)
  			plot_multi_noisy(summary=disp_list[[5]], summaries_list=disp_list[6:7], ylim=c(0,2), ylab="Dispersal distance", header="", legend=T, legend_text=legend_text_2, legend_title=expression(paste(italic("q")["f"])))
  			 mtext("A. Dispersal with fixed dormancy", adj=0.01)	
  		   	lines(x=c(log10(0.001),log10(0.001)), y=c(-0.5,2.5), lty=3)
  			lines(x=c(log10(0.01),log10(0.01)), y=c(-0.5,2.5), lty=3)
  			lines(x=c(log10(0.1),log10(0.1)), y=c(-0.5,2.5), lty=3)  			  	
  		screen(2)
	   	par(mai=c(.8, .9, .5, .1), ps=18) #(bottom,left,top, right)
			plot_multi_noisy(summary=diap_list[[2]], summaries_list=diap_list[3:4], ylim=c(0,.55), ylab="Dormancy probability", header="", legend=T, legend_text=legend_text_1, legend_title=expression(paste(sigma["m"])))
  			mtext("B. Dormancy with fixed dispersal", adj=0.01)
  			lines(x=c(log10(0.001),log10(0.001)), y=c(-0.5,2), lty=3)
  			lines(x=c(log10(0.01),log10(0.01)), y=c(-0.5,2), lty=3)
  			lines(x=c(log10(0.1),log10(0.1)), y=c(-0.5,2), lty=3)

  
	close.screen(all.screens = TRUE)
	dev.off()
}
write_fixed_disp_figure(acl=c(0.001, 0.0014, 0.0019, 0.0025, 0.0035, 0.005, 0.006, 0.0075, 0.0082, 0.0091, 0.010, 0.013, 0.016, 0.0200, 0.025, 0.035, 0.05, 0.075, 0.100), legend_text_1 = c("0.01", "0.05", "0.5"), legend_text_2=c("0.1","0.3", "0.5"))
