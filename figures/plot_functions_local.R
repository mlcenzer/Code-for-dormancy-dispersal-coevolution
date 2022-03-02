setwd("~/Documents/GitHub/Code-for-dormancy-dispersal-coevolution")
source('src/initialize.R')
library(abind)
library(scales)
library(viridis)
library(autoimage)

##Summarize runs for multi-panel diap/disp summary figures
  make_summary_multi<-function(variable="diap", stripes=stripes_all, load_from=load_from, acl=NA){
    folders_list=list.files(path=load_from)
    var_list=list()
    list_num <- 1
    for(folder in folders_list){
      all_files<-list.files(path=paste(load_from, "/", folder, sep=""))	
	if(!is.na(stripes[1]))  var_summary<-data.frame(mean=rep(NA, times=length(all_files)), CI=rep(NA, times=length(all_files)), stripes=rep(NA, times=length(all_files)))
	if(!is.na(acl[1]))  var_summary<-data.frame(mean=rep(NA, times=length(all_files)), CI=rep(NA, times=length(all_files)), acl=rep(NA, times=length(all_files)))
	#browser()
      file_num <- 1
      for(file in all_files){
        load(paste(load_from, "/", folder, "/", file, sep=""))
		if(!is.na(stripes[1]))  var_summary$stripes[file_num]<-as.numeric(strsplit(strsplit(file, "stripes_")[[1]][2], "_")[[1]][1])
		if(!is.na(acl[1]))  var_summary$acl[file_num]<-as.numeric(strsplit(strsplit(file, "acl_")[[1]][2], "_")[[1]][1])
        if(variable=="diap") var <- 1
        if(variable=="disp") var <- 2
        if(variable=="pop_size") var <- 3
        var_summary$mean[file_num]<-mean(summary_list[[var]]$mean)
        var_summary$CI[file_num]<-sd(summary_list[[var]]$mean)/sqrt(nrow(summary_list[[var]]))
        file_num <- file_num + 1
      }
                                        #browser()
      if(!is.na(stripes[1])) restricted_var_summary <- var_summary[var_summary$stripes%in%stripes,]
      if(!is.na(acl[1])) restricted_var_summary <- var_summary[var_summary$acl%in%acl,]

      var_list[[list_num]] <- restricted_var_summary
      names(var_list)[[list_num]] <- folder
      list_num <- list_num + 1
    }
    var_list
  }


##summarize runs to plot for diap/disp corr plot, others?
  make_summary_corr<-function(stripes=NA, acl=NA, land_var, load_from=load_from) {
    corr_list <- list()
    #browser()
    all_files<-list.files(path=paste(load_from, sep=""))
    
    list_names <- rep(NA, length=length(land_var))
    
    list_num <- 1
    for(file in all_files){
     if(stripes==T) var_num <- as.numeric(strsplit(strsplit(file, "stripes_")[[1]][2], "_")[[1]][1])
     if(acl==T) var_num <- as.numeric(strsplit(strsplit(file, "acl_")[[1]][2], "_")[[1]][1])
	#browser()
      if(var_num%in%land_var){
        load(paste(load_from, "/", file, sep=""))
        temp_summary <- cbind(summary_list[[1]]$mean, summary_list[[2]]$mean)
        colnames(temp_summary) <- c("diapause", "dispersal")
        corr_list[[list_num]] <- temp_summary
        list_names[list_num] <- var_num
        list_num <- list_num + 1
      }
    }
    if(acl==T) list_names<--log10(list_names)
    names(corr_list)<-list_names
    corr_list
  }

##function for adding color scale bar  
  add.gray.scale <- function() {
    vp <- baseViewports()
    pushViewport(vp$inner,vp$figure,vp$plot)
    ## adjust x, y, and height to move location around, etc
    pushViewport(viewport(x=0.9, y=0.95, width=0.06, height=.4,
                          just=c("left","top")))
    op <- par(plt=gridPLT(),new=T)
    on.exit(par(op))
    cols <- viridis(100, alpha=0.9)
    image(matrix(1:100, ncol=100),
          col=cols, xaxt="n", yaxt="n")
    popViewport(4)
  }    

## (noisy) plot output from multiple simulations in one plot
plot_multi_noisy<-function(summary, summaries_list=NA, ylim=c(0,0.15), ylab="", header="", legend=T, legend_text=NA, legend_title=NA){
  summary$acl_flip<--log10(summary$acl)
  summary <- summary[order(summary$acl_flip),]
	xlim=c(min(summary$acl_flip),max(summary$acl_flip))

  plot(NA, xlim=xlim, ylim=ylim, xlab=expression(paste("Spatial autocorrelation (", italic("-log"["10"]), italic("(acl)"), ")")), ylab=ylab, las=1)

  color<-viridis(length(summary$acl), alpha=0.5)
	#browser()
  points(summary$mean~summary$acl_flip, pch=21, bg=color[as.factor(summary$acl_flip)], col=color[as.factor(summary$acl_flip)], cex=2)

  for(row in 1:nrow(summary)){
                                        #browser()
    lines(x=rep(summary$acl_flip[row],2), y=c((summary$mean[row]+summary$CI[row]), (summary$mean[row]-summary$CI[row])), col=color[row], lwd=2)
  } 
    pch_vec <- c(22:25) 
  mtext(header, side=3, adj=0.01)
  if(!is.na(summaries_list)){

    for(summary_num in 1:length(summaries_list)){
                                        #browser()
    	summary_2 <- summaries_list[[summary_num]]
  		summary_2$acl_flip<--log10(summary_2$acl)
		summary_2 <- summary_2[order(summary_2$acl_flip),]
      
   		points(summary_2$mean~summary_2$acl_flip, pch=pch_vec[summary_num], bg=color[as.factor(summary_2$acl_flip)], col=color[as.factor(summary_2$acl_flip)], cex=2)

      for(row in 1:nrow(summary_2)){
        lines(x=rep(summary_2$acl_flip[row],2), y=c((summary_2$mean[row]+summary_2$CI[row]), (summary_2$mean[row]-summary_2$CI[row])), col=color[row], lwd=2)
      		}
    	}
 	}
	if(legend==TRUE) 	legend(xlim[1], (ylim[2]-0.005), legend=legend_text, col=1, pt.bg='grey', pch=c(21,pch_vec), y.intersp = 1.3, cex=0.8, pt.cex=1.2, title = legend_title) 
 	
}

## plot output from multiple simulations in one plot
plot_multi_stripes<-function(summary, summaries_list=NA, ylim=c(0,0.15), ylab="", header="", legend=TRUE, legend_text=NA, legend_title=NA, add_disp=FALSE, disp_summary, outlier=NA){
	temp_summary <- summary[summary$stripes!=outlier,]
	xlim=c(min(summary$stripes),max(summary$stripes+1))
	if(!is.na(outlier[1])) {
		xlim=c(min(temp_summary$stripes), (max(temp_summary$stripes)+3))
		outlier_color <- rgb(1,0,0,alpha=0.5)
		}
  plot(NA, xlim=xlim, ylim=ylim, xlab=expression(paste("Patch number (", italic("n"["p"]), ")")), ylab=ylab, las=1)

  color<-viridis(max(temp_summary), alpha=0.5)
      
  points(temp_summary$mean~temp_summary$stripes, pch=21, bg=color[temp_summary$stripes], col=color[temp_summary$stripes], cex=2)

  for(row in 1:nrow(temp_summary)){
    lines(x=rep(temp_summary$stripes[row],2), y=c((temp_summary$mean[row]+temp_summary$CI[row]), (temp_summary$mean[row]-temp_summary$CI[row])), col=color[temp_summary$stripes[row]], lwd=2)
  }

	if(!is.na(outlier[1])){
		outlier_summary <- summary[summary$stripes==outlier,]

		points(outlier_summary$mean~max(xlim), pch=21, bg=outlier_color, col=outlier_color, cex=2)
		lines(x=rep(max(xlim), 2), y=c((outlier_summary$mean+outlier_summary$CI),(outlier_summary$mean-outlier_summary$CI)), col=outlier_color, lwd=2)
	}
  
  if(!is.na(summaries_list[1])){
    pch_vec <- c(22:25)
    for(summary_num in 1:length(summaries_list)){
                                        #browser()
      summary_2 <- summaries_list[[summary_num]]
      
      points(summary_2$mean~summary_2$stripes, pch=pch_vec[summary_num], bg=color[summary_2$stripes], col=color[summary_2$stripes], cex=2)

      for(row in 1:nrow(summary_2)){
        lines(x=rep(summary_2$stripes[row],2), y=c((summary_2$mean[row]+summary_2$CI[row]), (summary_2$mean[row]-summary_2$CI[row])), col=color[summary_2$stripes[row]], lwd=2)
      		}
      	if(!is.na(outlier[1])){
      		points(summary_2$mean[summary_2$stripes==outlier] ~ max(xlim), pch=pch_vec[summary_num], bg=outlier_color, col=outlier_color, cex=2)
      		lines(x=rep(max(xlim), 2), y=c((summary_2$mean[summary_2$stripes==outlier]+summary_2$CI[summary_2$stripes==outlier]), (summary_2$mean[summary_2$stripes==outlier]-summary_2$CI[summary_2$stripes==outlier])), col=outlier_color, lwd=2)
      		}
    	}
 	}  
 
  mtext(header, side=3, adj=0.01)
  #par(ps=14)
  if(legend==TRUE) 	legend(xlim[1], (ylim[2]-0.005), legend=legend_text, col=1, pt.bg='grey', pch=c(21,pch_vec), y.intersp = 1.3, cex=0.8, pt.cex=1.2, title = legend_title)

}

#plot tau change
plot_kin_comp_tau<-function(summary_list, stripes_all=c(1,2,4,8), disp_ev=T, header="Header here", diap_length=FALSE, legend=TRUE, ylim, outlier=F, yaxt="s", ylab=T){
	if(outlier==T) {
		stripes=stripes_all[1:(length(stripes_all)-1)]
		outlier_stripe = stripes_all[length(stripes_all)]
		#browser()
	}
	if(diap_length) {
		plot(NA, xlim=c(0,50), ylim=ylim, xlab="Dormancy probability", ylab=expression(paste("Change in competition (", delta, ")")), las=1, yaxt=yaxt)
	}
	else {
	plot(NA, xlim=c(0,1), ylim=ylim, xlab="Dormancy probability", ylab=NA, las=1, yaxt=yaxt)
	if(ylab==T) mtext(expression(paste("Change in competition (", Delta, tau, ")")), side = 2, line=4)
	}
	#browser()
	color<-viridis(max(stripes))
	for(stripe in names(summary_list)){
		#browser()
		summary<-summary_list[[as.character(stripe)]]
		lines(mean_tau_change~diap_bins, data=summary, col=color[as.numeric(stripe)], lwd=3)	
	}
	if(outlier==T){
		summary<-summary_list[[as.character(outlier_stripe)]]
		lines(mean_tau_change~diap_bins, data=summary, col=rgb(1,0,0), lwd=3)
	}
	if(legend==TRUE){
	#legend(0, ylim[2], legend=c(1, (max(stripes)/2), max(stripes)), col=color[c(1, (max(stripes)/2), max(stripes))], pch=19, y.intersp = 1.3, pt.cex=1.2)

		}

	mtext(header, side=3, adj=0.01)
}





