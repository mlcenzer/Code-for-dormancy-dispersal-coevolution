##Update to your local directory here
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
      if(!is.na(stripes[1])) restricted_var_summary <- var_summary[var_summary$stripes%in%stripes,]
      if(!is.na(acl[1])) restricted_var_summary <- var_summary[var_summary$acl%in%acl,]

      var_list[[list_num]] <- restricted_var_summary
      names(var_list)[[list_num]] <- folder
      list_num <- list_num + 1
    }
    var_list
  }


##summarize runs to plot for diap/disp corr plots
  make_summary_corr<-function(stripes=NA, acl=NA, land_var, load_from=load_from) {
    corr_list <- list()
    #browser()
    all_files<-list.files(path=paste(load_from, sep=""))
    
    list_names <- rep(NA, length=length(land_var))
    
    list_num <- 1
    for(file in all_files){
     if(stripes==T) var_num <- as.numeric(strsplit(strsplit(file, "stripes_")[[1]][2], "_")[[1]][1])
     if(acl==T) var_num <- as.numeric(strsplit(strsplit(file, "acl_")[[1]][2], "_")[[1]][1])
	
      if(var_num%in%land_var){
        load(paste(load_from, "/", file, sep=""))
        temp_summary <- cbind(summary_list[[1]]$mean, summary_list[[2]]$mean)
        colnames(temp_summary) <- c("diapause", "dispersal")
        corr_list[[list_num]] <- temp_summary
        list_names[list_num] <- var_num
        list_num <- list_num + 1
      }
    }
    if(acl==T) list_names<-log10(list_names)
    names(corr_list)<-list_names
    corr_list
  }

## (noisy) plot output from multiple simulations in one plot
plot_multi_noisy<-function(summary, summaries_list=NA, ylim=c(0,0.15), ylab="", header="", legend=T, legend_text=NA, legend_title=NA){
  summary$acl_flip<-log10(summary$acl)
  #summary$acl_flip<-summary$acl
   summary <- summary[order(summary$acl_flip, decreasing=T),]
	xlim=c(min(summary$acl_flip),max(summary$acl_flip))

  plot(NA, xlim=xlim, ylim=ylim, xlab=expression(paste("Spatial autocorrelation (", "log"["10"], italic("(acl)"), ")")), ylab=ylab, las=1)

  color<-viridis(length(summary$acl), alpha=0.5)

  for(row in 1:nrow(summary)){
   points(summary$mean[row]~summary$acl_flip[row], pch=21, bg=color[row], col=color[row], cex=2)                                       #browser()
    lines(x=rep(summary$acl_flip[row],2), y=c((summary$mean[row]+summary$CI[row]), (summary$mean[row]-summary$CI[row])), col=color[row], lwd=2)
  } 
    pch_vec <- c(22:25) 
  mtext(header, side=3, adj=0.01)
  if(!is.na(summaries_list[1])){

    for(summary_num in 1:length(summaries_list)){
    	summary_2 <- summaries_list[[summary_num]]
  		summary_2$acl_flip<-log10(summary_2$acl)
		summary_2 <- summary_2[order(summary_2$acl_flip, decreasing=T),]

      for(row in 1:nrow(summary_2)){
      	points(summary_2$mean[row]~summary_2$acl_flip[row], pch=pch_vec[summary_num], bg=color[row], col=color[row], cex=2)
        lines(x=rep(summary_2$acl_flip[row],2), y=c((summary_2$mean[row]+summary_2$CI[row]), (summary_2$mean[row]-summary_2$CI[row])), col=color[row], lwd=2)
      		}
    		}
 	}
	if(legend==TRUE) 	legend(xlim[2]-0.4, (ylim[2]-0.005), legend=legend_text, col=1, pt.bg='grey', pch=c(21,pch_vec), y.intersp = 1.3, cex=0.8, pt.cex=1.2, title = legend_title) 
 	
}


## plot (stripes) output from multiple simulations in one plot
plot_multi_stripes<-function(summary, summaries_list=NA, ylim=c(0,0.15), ylab="", header="", legend=TRUE, legend_text=NA, legend_title=NA, add_disp=FALSE, disp_summary, outlier=NA, zero=NA, width=F){
	summary<-summary[order(summary$stripes, decreasing=F),]
	summary$plot_by<-summary$stripes
	summary$labels<-summary$stripes
	if(width==T) {
		plot_by<-sort(summary$stripes, decreasing=T)
		outlier_loc <- plot_by[(length(plot_by)-1)]-(3/16)*plot_by[2]
		zero_loc <- plot_by[2]+(3/16)*plot_by[2]
		summary$plot_by <- c(zero_loc, plot_by[2:(length(plot_by)-1)], outlier_loc)
		summary$labels <- round(1/(2*summary$stripes), digits=3)
		summary$labels[summary$labels==-Inf | summary$labels==Inf] <- 1
		}

	temp_summary <- summary[summary$stripes!=outlier & summary$stripes!=zero,]
	xlim=c(min(temp_summary$stripes),max(temp_summary$stripes))
	
	if(!is.na(outlier[1])) {
		xlim[2]=(max(temp_summary$plot_by)+(3/16)*max(temp_summary$plot_by))
		outlier_color <- rgb(1,0,0,alpha=0.5)
		}
		
	if(!is.na(zero)) {
		xlim[1]=(min(temp_summary$plot_by)-(3/16)*max(temp_summary$plot_by))
		zero_color <- rgb(0.5, 0.5, 0.5, alpha=0.5)
		}
		
  plot(NA, xlim=xlim, ylim=ylim, xlab=expression(paste("Patch width (1/2", italic("n")["p"], "))")), ylab=ylab, xaxt='n', las=1)
	axis(1, at=temp_summary$plot_by[c(2,5,8,11,14)], labels=temp_summary$labels[c(2,5,8,11,14)])
	
  color<-viridis(nrow(temp_summary), alpha=0.5)
      
  points(temp_summary$mean~temp_summary$plot_by, pch=21, bg=color, col=color, cex=2)

  for(row in 1:nrow(temp_summary)){
    lines(x=rep(temp_summary$plot_by[row],2), y=c((temp_summary$mean[row]+temp_summary$CI[row]), (temp_summary$mean[row]-temp_summary$CI[row])), col=color[row], lwd=2)
  }

	if(!is.na(outlier[1])){
		outlier_summary <- summary[summary$stripes==outlier,]

		points(outlier_summary$mean~outlier_summary$plot_by, pch=21, bg=outlier_color, col=outlier_color, cex=2)
		lines(x=rep(outlier_summary$plot_by, 2), y=c((outlier_summary$mean+outlier_summary$CI),(outlier_summary$mean-outlier_summary$CI)), col=outlier_color, lwd=2)
	}
	if(!is.na(zero)){
		zero_summary <- summary[summary$stripes==zero,]

		points(zero_summary$mean~zero_summary$plot_by, pch=21, bg=zero_color, col=zero_color, cex=2)
		lines(x=rep(zero_summary$plot_by, 2), y=c((zero_summary$mean+zero_summary$CI),(zero_summary$mean-zero_summary$CI)), col=zero_color, lwd=2)
	}
  
  if(!is.na(summaries_list[1])){
    pch_vec <- c(22:25)
    for(summary_num in 1:length(summaries_list)){
                                       
      	summary_2 <- summaries_list[[summary_num]]
      	summary_2<-summary_2[order(summary_2$stripes, decreasing=F),]
      	summary_2$plot_by<-summary_2$stripes
		summary_2$labels<-summary_2$stripes
		
	if(width==T) {
		plot_by<-sort(summary_2$stripes, decreasing=T)
		outlier_loc <- plot_by[(length(plot_by)-1)]-(3/16)*plot_by[2]
		zero_loc <- plot_by[2]+(3/16)*plot_by[2]
		summary_2$plot_by <- c(zero_loc, plot_by[2:(length(plot_by)-1)], outlier_loc)
		summary_2$labels <- round(1/(2*summary_2$stripes), digits=2)
		summary_2$labels[summary_2$labels==-Inf | summary_2$labels==Inf] <- 1
		} 
            
      temp_summary_2 <- summary_2[summary_2$stripes!=zero & summary_2$stripes!=outlier,]
      
      points(temp_summary_2$mean~temp_summary_2$plot_by, pch=pch_vec[summary_num], bg=color, col=color, cex=2)

      for(row in 1:nrow(temp_summary_2)){
        lines(x=rep(temp_summary_2$plot_by[row],2), y=c((temp_summary_2$mean[row]+temp_summary_2$CI[row]), (temp_summary_2$mean[row]-temp_summary_2$CI[row])), col=color[row], lwd=2)
      		}
      	if(!is.na(outlier[1])){
      		outlier_summary_2 <- summary_2[summary_2$stripes==outlier,]
      		points(outlier_summary_2$mean ~ outlier_summary_2$plot_by, pch=pch_vec[summary_num], bg=outlier_color, col=outlier_color, cex=2)
      		lines(x=rep(outlier_summary_2$plot_by, 2), y=c((outlier_summary_2$mean+outlier_summary_2$CI), (outlier_summary_2$mean-outlier_summary_2$CI)), col=outlier_color, lwd=2)
      		}
      	if(!is.na(zero)){
      		zero_summary_2 <- summary_2[summary_2$stripes==zero,]
      		points(zero_summary_2$mean ~ zero_summary_2$plot_by, pch=pch_vec[summary_num], bg=zero_color, col=zero_color, cex=2)
      		lines(x=rep(zero_summary_2$plot_by, 2), y=c((zero_summary_2$mean+zero_summary_2$CI), (zero_summary_2$mean-zero_summary_2$CI)), col=zero_color, lwd=2)
      		}
    		}
 	}  
 
  mtext(header, side=3, adj=0.01)
  #par(ps=14)
  if(legend==TRUE) 	legend((xlim[1]+12), (ylim[2]), legend=legend_text, col=1, pt.bg='grey', pch=c(21,pch_vec), y.intersp = 1.3, cex=0.8, pt.cex=1.2, title = legend_title)

}

#plot tau change
plot_kin_comp_tau<-function(summary_list, stripes_all=c(1,2,4,8), disp_ev=T, header="Header here", diap_length=FALSE, ylim, outlier=NA, zero=NA, yaxt="s", ylab=T){
	if(!is.na(outlier)) stripes_all=stripes_all[stripes_all!=outlier]
	if(!is.na(zero)) stripes_all=stripes_all[stripes_all!=zero]

	plot(NA, xlim=c(0,1), ylim=ylim, xlab="Dormancy probability", ylab=NA, las=1, yaxt=yaxt)
	if(ylab==T) mtext(expression(paste("Change in competition (", Delta, tau, ")")), side = 2, line=4)
	
	color<-viridis(max(stripes_all))
	for(stripe in names(summary_list)){
		#browser()
		summary<-summary_list[[as.character(stripe)]]
		lines(mean_tau_change~diap_bins, data=summary, col=color[as.numeric(stripe)], lwd=3)	
	}
	if(!is.na(outlier)){
		summary<-summary_list[[as.character(outlier)]]
		lines(mean_tau_change~diap_bins, data=summary, col=rgb(1,0,0), lwd=3)
	}
	if(!is.na(zero)){
		summary<-summary_list[[as.character(zero)]]
		lines(mean_tau_change~diap_bins, data=summary, col=rgb(0.5,0.5,0.5), lwd=3)
	}
	mtext(header, side=3, adj=0.01)
}
