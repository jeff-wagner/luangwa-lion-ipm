
library(RMark)
#setwd("C:/Users/pcadmin/Desktop/rtemp")

setwd("C:/Users/davealanchris/Desktop/rtemp")

#exclude dead animals but add at the end
# add sex and age to add these animals to totals at end.
deaths<-data.frame(lionid=lionids$MARKID[!is.na(lionids$Date.of.Death)],
                   deathdate=lionids$Date.of.Death[!is.na(lionids$Date.of.Death)])
deaths$deathyear<-as.POSIXlt(deaths$deathdate)$year+1900
deaths$deathmonth<-as.POSIXlt(deaths$deathdate)$mon+1
deaths$sex<-inputdata$sex[match(deaths$lionid,inputdata$lionid)]
deaths$ageatfirst<-inputdata$ageatfirst[match(deaths$lionid,inputdata$lionid)]
deaths$yearatfirst<-inputdata$yearatfirst[match(deaths$lionid,inputdata$lionid)]
deaths$monthatfirst<-inputdata$monthatfirst[match(deaths$lionid,inputdata$lionid)]
deaths$ageatdeath<-deaths$ageatfirst+(difftime(strptime(paste(deaths$deathyear,deaths$deathmonth,"28",sep="-"),format="%Y-%m-%d"),
                                     strptime(paste(deaths$yearatfirst,deaths$monthatfirst,"28",sep="-"),format="%Y-%m-%d"),
                            units="days")/365)


## Huggins full likelihood heterogeneity models
## fit using model averaged pi estimates from best models of CJS
## note that rosenblatt et al fit mixture * season models with p=c for all lions, (despite what the Methods say).
## similar parameter estimates to rosenblatt et al. are retrieved wether or not an intercept is fit and whether
##  each mixture:season combination is modelled as a single categorical variable 
## (with 4 categories, e.g. mix1season1, mix1season2, mix2season1 and mix2season2)

## 2008 still has fitting problems.  Estimation of the capture probability in the second season for mixture 2 is primary 
## problem.  I fixed this real parameter using the mean for all other years.
#pmix2season1<-mean(do.call("rbind",realall)[2+seq(1,36,by=5),][2:8,]$estimate)
	pmix2season1<-0.1203006
	

	k=0
	nall<-list()
	betaall<-list()
	realall<-list()
	for (i in 2008:2015){
		k=k+1
		chyear<-paste0("ch",i)
		popdata<-inputdata[as.numeric(inputdata[[chyear]])>=0,]
		popdata$ch<-popdata[[chyear]]
		popdata<-popdata[as.numeric(popdata$ch)>0,]
		# exclude dead lions
		if(i<2014){
		popdata<-popdata[!popdata$lionid %in% deaths$lionid[deaths$deathyear==i],]
		}
		
		## to identify cub groups,the age at first detection must be calculated from the dob and the first 1 in the ch for each year
		dobs<-lionids$DOB[match(popdata$lionid,lionids$MARKID)]
		imonthatfirst<-3+(regexpr("1",popdata$ch))
		popdata$chageatfirst<-as.numeric(difftime(strptime(paste(i, imonthatfirst,"28",sep="-"),"%Y-%m-%d"),dobs)/365)
        # lions with no dob are assigned age 3 in the first detection for CJS. To maintian consistency, these animals must
		# age
		for(j in 1:length(dobs)){
		    if(is.na(dobs[j])){
		        if(i==popdata$yearatfirst[j]){
		            popdata$chageatfirst[j]<-4.7
		        }
		        else{
		           popdata$chageatfirst[j]<-popdata$ageatfirst[j]+(difftime(strptime(paste(i,imonthatfirst[j],"28",sep="-"),"%Y-%m-%d"),strptime(paste(popdata$yearatfirst[j],popdata$monthatfirst[j],"28",sep="-"),"%Y-%m-%d"))/365)
		        }
		    } 
		}
		## define groups as needed
		
		popdata$agegroup<-"noncub0"
		#popdata$agegroup[popdata$chageatfirst>=4&popdata$sex==0]<-"fold4"
		popdata$agegroup[popdata$chageatfirst<1]<-"cub0"
		
		#hughet
		huggins.process<-process.data(popdata,model="HugFullHet",
                         mixtures=2,groups=c("agegroup"))
		
		huggins.dll<-make.design.data(huggins.process)
		
		seasonbins<-c(1,7,9)
		huggins.dll<-add.design.data(huggins.process,huggins.dll,parameter='p',type='time',bins=seasonbins,name='season',right=F)
		huggins.dll<-add.design.data(huggins.process,huggins.dll,parameter='c',type='time',bins=seasonbins,name='season',right=F)
		
		# formulas
		p.mixture=list(formula=~mixture*season,share=T)
		pi.mix=list(formula=~1)
		huggins.dll$pi$fixed<-c(0.5227588)
		if (i==2008){
		  p.mixture=list(formula=~mixture+season,share=T)
		}

		pop.huggins<-mark(huggins.process,huggins.dll,model="HugFullHet", model.parameters=list(p=p.mixture,pi=pi.mix))

		nall[[k]]<-pop.huggins$results$derived$`N Population Size`
		realall[[k]]<-pop.huggins$results$real
		betaall[[k]]<-pop.huggins$results$beta
	}
 
# no groups
	k=0
	nallng<-list()
	betaallng<-list()
	realallng<-list()
	for (i in 2008:2015){
	  k=k+1
	  chyear<-paste0("ch",i)
	  popdata<-inputdata[as.numeric(inputdata[[chyear]])>=0,]
	  popdata$ch<-popdata[[chyear]]
	  popdata<-popdata[as.numeric(popdata$ch)>0,]
	  # exclude dead lions
	  if(i<2014){
	    popdata<-popdata[!popdata$lionid %in% deaths$lionid[deaths$deathyear==i],]
	  }
	  
	   # age
	  for(j in 1:length(dobs)){
	    if(is.na(dobs[j])){
	      if(i==popdata$yearatfirst[j]){
	        popdata$chageatfirst[j]<-4.7
	      }
	      else{
	        popdata$chageatfirst[j]<-popdata$ageatfirst[j]+(difftime(strptime(paste(i,imonthatfirst[j],"28",sep="-"),"%Y-%m-%d"),strptime(paste(popdata$yearatfirst[j],popdata$monthatfirst[j],"28",sep="-"),"%Y-%m-%d"))/365)
	      }
	    } 
	  }
	  print(i)
    print(dim(popdata))
	  #hughet
	  huggins.process<-process.data(popdata,model="HugFullHet",
	                                mixtures=2)
	  
	  huggins.dll<-make.design.data(huggins.process)
	  
	  seasonbins<-c(1,7,9)
	  huggins.dll<-add.design.data(huggins.process,huggins.dll,parameter='p',type='time',bins=seasonbins,name='season',right=F)
	  huggins.dll<-add.design.data(huggins.process,huggins.dll,parameter='c',type='time',bins=seasonbins,name='season',right=F)
	  
	  # formulas
	  p.mixture=list(formula=~mixture*season,share=T)
	  pi.mix=list(formula=~1)
	  huggins.dll$pi$fixed<-c(0.5227588)
	  if (i==2008){
	    p.mixture=list(formula=~mixture+season,share=T)
	  }
	  pop.huggins<-mark(huggins.process,huggins.dll,model="HugFullHet", model.parameters=list(p=p.mixture,pi=pi.mix))
	  
	  nallng[[k]]<-pop.huggins$results$derived$`N Population Size`
	  realallng[[k]]<-pop.huggins$results$real
	  betaallng[[k]]<-pop.huggins$results$beta
	}
	

	## summarize population data and add dead animals
	deaths$noncub0<-(deaths$ageatdeath<1|is.na(deaths$ageatdeath))==F
	
	## these estimates use all data but with groups defined
	popest<-do.call("rbind",nall)
	popest$agegroup<-c("cub0","noncub0")
	popest$year<-rep(seq(2008,2015),each=2)
	popest$estimate<-popest$estimate + c(with(deaths,c(table(noncub0,deathyear))),0,0,0,0)
	popest$lcl<-popest$lcl + c(with(deaths,c(table(noncub0,deathyear))),0,0,0,0)
	popest$ucl<-popest$ucl + c(with(deaths,c(table(noncub0,deathyear))),0,0,0,0)
	write.table(popest,"clipboard",sep="\t")

	# these estimates use all data with no groups defined
	popestallng<-do.call("rbind",nallng)
	popestallng$year<-rep(seq(2008,2015),each=1)
	# add dead animals
	popestallng$estimate<-popestallng$estimate + with(deaths,c(table(deathyear),0,0))
	popestallng$lcl<-popestallng$lcl + with(deaths,c(table(deathyear),0,0))
	popestallng$ucl<-popestallng$ucl + with(deaths,c(table(deathyear),0,0))
	write.table(popestallng,"clipboard",sep="\t")
	

	

# the population estimate using all lions, no groups is more accurate than summing the estimate using all lions with
# groups.  The estimate using data that excludes a group (like first year cubs) before fitting the model also has larger variance
# than the estimate produced fitting all the data with all groups in the data and groups defined.
		
#	jpeg("D:/Dropbox/projectdave/zambia/lions/paper/hugfullhetcounts_all.jpg",width=6.5,height=3.25,units="in",res=300)
	jpeg("C:/Dropbox/projectdave/zambia/lions/paper/hugfullhetcounts_all.jpg",width=6.5,height=6.5,units="in",res=600)
	
	xoffset=0.1
	larrow=0.05
	cub0col="blue"
	noncub0col="green1"
	xtcl=0.5
	ytcl=0.5
	htitlebox=0.15
	par(oma=c(0,0,0,0),mar=c(2,4,0.5,0.5))	
	plot(c(2008,2015+xoffset),c(0,250),type="n",xlab="",ylab="N-hat",yaxt="n",xaxt="n",yaxs="i")
	grid()
	points(2008:2015-xoffset,popestallng$estimate,pch=19)
	arrows(2008:2015-xoffset,popestallng$lcl,2008:2015-xoffset,popestallng$ucl,code=3,angle=90,length=larrow)
	points(2008:2015+xoffset,popest$estimate[popest$agegroup=="cub0"],col=cub0col,pch=19)
	arrows(2008:2015+xoffset,popest$lcl[popest$agegroup=="cub0"],2008:2015+xoffset,popest$ucl[popest$agegroup=="cub0"],code=3,angle=90,length=larrow,col=cub0col)
	points(2008:2015+xoffset,popest$estimate[popest$agegroup=="noncub0"],col=noncub0col,pch=19)
	arrows(2008:2015+xoffset,popest$lcl[popest$agegroup=="noncub0"],2008:2015+xoffset,popest$ucl[popest$agegroup=="noncub0"],code=3,angle=90,length=larrow,col=noncub0col)
	
	
	axis(2,at=c(0,50,100,150,200,250),tcl=ytcl,las=2)
	axis(1,at=2007:2016,tcl=xtcl,labels=NA)
	axis(1,at=2008:2015,tcl=xtcl)
	axis(3,at=c(2007,2016),tcl=0)
	#rect(par("usr")[1],par("usr")[4],par("usr")[2],par("usr")[4]+(htitlebox*(par("usr")[4]-par("usr")[3])),xpd=T,col="grey",border="black")
	#text(2011.5,par("usr")[4]+((htitlebox/2)*(par("usr")[4]-par("usr")[3])),"All Lions",xpd=T,cex=2)
	legend(2008.25,240,pch=c(19,19,19),col=c("black",noncub0col,cub0col),legend=c("all lions","\u22651 yr","<1 yr"),box.col="light grey",bg="light grey",cex=1.25)
	lines(c(2012.5,2012.5),par("usr")[3:4],col="red",lwd=1.5)
	text(2014,par("usr")[3]+((htitlebox/2)*(par("usr")[4]-par("usr")[3])),"hunting\nmoratorium",col="red",cex=1.5)
	
	dev.off()


write.table(do.call("rbind",betaallng),"clipboard",sep="\t")
write.table(do.call("rbind",realallng),"clipboard",sep="\t")

	