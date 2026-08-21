#rm(list=ls())  
library(RMark)

setwd("C:/Users/David/Desktop/rtemp")
#setwd("C:/Users/pcadmin/Desktop/rtemp")

## the rmark_input_script must be run before running this script
inputdata$freq<-1

inputdata<-inputdata[,c("ch","freq","sex","ageatfirst")]

inputdata<-inputdata[inputdata$ageatfirst>=2,]

# of occasions and years
primary<-4
years<-8

# round age to integer after adding the length of 1 or 2 months to make certain that an individual hasn't changed age classes within an interval
# thus, individuals can be grouped based on ageclass at the end of each primary occasion (May 31, July 31, September 30, November 30)
inputdata$ageatfirst<-floor(inputdata$ageatfirst*2)
#collapse captures histories into new intervals, as in Rosenblatt et al 2014, if desired
cjschlist<-list(length=0)

for (k in 1:dim(inputdata)[1]){
    cjsch<-vector(length=0)
    for(i in 1:(primary*years)){
        cjsch<-paste0(cjsch,as.numeric(as.numeric(as.character(substr(inputdata$ch[k],1+((i-1)*(8/primary)),(8/primary)+((i-1)*(8/primary)))))>0))
    }
    cjschlist[[k]]<-cjsch
}


inputdata$ch<-as.character(as.data.frame(unlist(cjschlist))[,1])




#define time intervals 
interval1<-2/12
interval2<-6/12
time.intervals=c(rep(c(rep(interval1,primary-1),interval2),years-1),rep(interval1,primary-1))

#make 4 processed data types, with no mixture, mixture of 2 groups, mixture of 3 groups, 4 grups
# note the begin time is the end of the fist capture occasion


lions.process.cjsh2<-process.data(inputdata,
                                  model="CJSMixture",
                                  begin.time=(2008+(5/12)),
                                  time.intervals=time.intervals,
                                  groups=c("sex","ageatfirst"),
                                  age.var=2,
                                  initial.ages=c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,23,24,25,27),
                                  allgroups=T,
                                  mixtures=2,age.unit=2)


## design data
lions.ddl.cjsh2<-make.design.data(lions.process.cjsh2,common.zero=T)


# for eli's categories

agebins1246810<-c(2,4,6,8,10,22)*2
agebinseli<-c(2,4,6,8,22)*2

timebins<-c(2008,2009,2010,2011,2012,2013,2014,2015,2016)
huntbins<-c(2008,2012.8,2016)
seasonbins<-c(0,6/12,8/12,10/12,12/12)
eliseasonbins<-c(0,10/12,12/12)

# age
lions.ddl.cjsh2<-add.design.data(lions.process.cjsh2,lions.ddl.cjsh2,parameter="Phi", type='age', bins=agebins1246810,name='agebins1246810',right=F,replace=T)
lions.ddl.cjsh2<-add.design.data(lions.process.cjsh2,lions.ddl.cjsh2,parameter="Phi", type='age', bins=agebinseli,name='agebinseli',right=F,replace=T)
lions.ddl.cjsh2<-add.design.data(lions.process.cjsh2,lions.ddl.cjsh2,parameter="p", type='age', bins=agebins1246810,name='agebins1246810',right=F,replace=T)
lions.ddl.cjsh2<-add.design.data(lions.process.cjsh2,lions.ddl.cjsh2,parameter="p", type='age', bins=agebinseli,name='agebinseli',right=F,replace=T)

## hunting 
lions.ddl.cjsh2<-add.design.data(lions.process.cjsh2,lions.ddl.cjsh2,parameter="Phi", type='time', bins=huntbins,name='hunting',right=F,replace=T)
lions.ddl.cjsh2<-add.design.data(lions.process.cjsh2,lions.ddl.cjsh2,parameter="p", type='time', bins=huntbins,name='hunting',right=F,replace=T)

# seasons and times for p (effort varied by year)
lions.ddl.cjsh2$p$Season<-as.numeric(cut(as.numeric(as.character(lions.ddl.cjsh2$p$time))-floor(as.numeric(as.character(lions.ddl.cjsh2$p$time))),breaks=seasonbins,labels=as.character(seq(0,length(seasonbins)-2)),right=F,include.lowest=T))
lions.ddl.cjsh2$p$eliSeason<-as.numeric(cut(as.numeric(as.character(lions.ddl.cjsh2$p$time))-floor(as.numeric(as.character(lions.ddl.cjsh2$p$time))),breaks=eliseasonbins,labels=as.character(seq(0,length(eliseasonbins)-2)),right=F,include.lowest=T))
lions.ddl.cjsh2$p$season<-as.factor(lions.ddl.cjsh2$p$Season)

lions.ddl.cjsh2$p$eliseason<-as.factor(lions.ddl.cjsh2$p$eliSeason)

# year    
lions.ddl.cjsh2<-add.design.data(lions.process.cjsh2,lions.ddl.cjsh2,parameter="p", type='time', bins=timebins,name='studyyear',right=T,replace=T)

#######################################################

#subadults

lions.ddl.cjsh2$Phi$sub<-F
lions.ddl.cjsh2$Phi$sub[lions.ddl.cjsh2$Phi$agebinseli=="[4,8)"]<-T

lions.ddl.cjsh2$Phi$msub<-F
lions.ddl.cjsh2$Phi$msub[lions.ddl.cjsh2$Phi$agebinseli=="[4,8)"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$fsub<-F
lions.ddl.cjsh2$Phi$fsub[lions.ddl.cjsh2$Phi$agebinseli=="[4,8)"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$msubhunt<-F
lions.ddl.cjsh2$Phi$msubhunt[lions.ddl.cjsh2$Phi$agebinseli=="[4,8)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

# young adults
lions.ddl.cjsh2$Phi$yad<-F
lions.ddl.cjsh2$Phi$yad[lions.ddl.cjsh2$Phi$agebinseli=="[8,12)"]<-T

lions.ddl.cjsh2$Phi$myad<-F
lions.ddl.cjsh2$Phi$myad[lions.ddl.cjsh2$Phi$agebinseli=="[8,12)"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$fyad<-F
lions.ddl.cjsh2$Phi$fyad[lions.ddl.cjsh2$Phi$agebinseli=="[8,12)"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$myadhunt<-F
lions.ddl.cjsh2$Phi$myadhunt[lions.ddl.cjsh2$Phi$agebinseli=="[8,12)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

# adults	
lions.ddl.cjsh2$Phi$adu<-F
lions.ddl.cjsh2$Phi$adu[lions.ddl.cjsh2$Phi$agebinseli=="[12,16)"]<-T

lions.ddl.cjsh2$Phi$madu<-F
lions.ddl.cjsh2$Phi$madu[lions.ddl.cjsh2$Phi$agebinseli=="[12,16)"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$fadu<-F
lions.ddl.cjsh2$Phi$fadu[lions.ddl.cjsh2$Phi$agebinseli=="[12,16)"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$maduhunt<-F
lions.ddl.cjsh2$Phi$maduhunt[lions.ddl.cjsh2$Phi$agebinseli=="[12,16)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T


#adu4 refers to 48 month yng and adu age class
lions.ddl.cjsh2$Phi$adu4<-F
lions.ddl.cjsh2$Phi$adu4[lions.ddl.cjsh2$Phi$agebinseli %in% c("[8,12)","[12,16)")]<-T

lions.ddl.cjsh2$Phi$fadu4<-F
lions.ddl.cjsh2$Phi$fadu4[lions.ddl.cjsh2$Phi$agebinseli %in% c("[8,12)","[12,16)")&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$madu4<-F
lions.ddl.cjsh2$Phi$madu4[lions.ddl.cjsh2$Phi$agebinseli %in% c("[8,12)","[12,16)")&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$madu4hunt<-F
lions.ddl.cjsh2$Phi$madu4hunt[lions.ddl.cjsh2$Phi$agebinseli %in% c("[8,12)","[12,16)")&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

lions.ddl.cjsh2$Phi$fadu468<-F
lions.ddl.cjsh2$Phi$fadu468[lions.ddl.cjsh2$Phi$agebins1246810 %in% c("[8,12)","[12,16)","[16,20)")&lions.ddl.cjsh2$Phi$sex==0]<-T


#adu8 refers to another adult age class extend from 8-10
lions.ddl.cjsh2$Phi$adu8<-F
lions.ddl.cjsh2$Phi$adu8[lions.ddl.cjsh2$Phi$agebins1246810=="[16,20)"]<-T

lions.ddl.cjsh2$Phi$fadu8<-F
lions.ddl.cjsh2$Phi$fadu8[lions.ddl.cjsh2$Phi$agebins1246810=="[16,20)"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$madu8<-F
lions.ddl.cjsh2$Phi$madu8[lions.ddl.cjsh2$Phi$agebins1246810=="[16,20)"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$madu8hunt<-F
lions.ddl.cjsh2$Phi$madu8hunt[lions.ddl.cjsh2$Phi$agebins1246810=="[16,20)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

#NONCUBS
#mnoncub refers to all males >=2 years (all non cub age classes)
lions.ddl.cjsh2$Phi$noncub<-F
lions.ddl.cjsh2$Phi$noncub[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"]<-T

lions.ddl.cjsh2$Phi$noncubhunt<-F
lions.ddl.cjsh2$Phi$noncubhunt[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

lions.ddl.cjsh2$Phi$mnoncub<-F
lions.ddl.cjsh2$Phi$mnoncub[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$mnoncubhunt<-F
lions.ddl.cjsh2$Phi$mnoncubhunt[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

lions.ddl.cjsh2$Phi$fnoncubhunt<-F
lions.ddl.cjsh2$Phi$fnoncubhunt[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"&lions.ddl.cjsh2$Phi$sex==0&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

#mnoncub0 refers to all males >=1 years (all non 1st year cub age classes)
lions.ddl.cjsh2$Phi$mnoncub0<-F
lions.ddl.cjsh2$Phi$mnoncub0[lions.ddl.cjsh2$Phi$agebins1246810!="[0,2)"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$mnoncub0hunt<-F
lions.ddl.cjsh2$Phi$mnoncub0hunt[lions.ddl.cjsh2$Phi$agebins1246810!="[0,2)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T


# old 	
lions.ddl.cjsh2$Phi$old<-F
lions.ddl.cjsh2$Phi$old[lions.ddl.cjsh2$Phi$agebinseli=="[16,44]"]<-T

lions.ddl.cjsh2$Phi$mold<-F
lions.ddl.cjsh2$Phi$mold[lions.ddl.cjsh2$Phi$agebinseli=="[16,44]"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$fold<-F
lions.ddl.cjsh2$Phi$fold[lions.ddl.cjsh2$Phi$agebinseli=="[16,44]"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$moldhunt<-F
lions.ddl.cjsh2$Phi$moldhunt[lions.ddl.cjsh2$Phi$agebinseli=="[16,44]"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

# old 10 refers to a very old age class 10-22
lions.ddl.cjsh2$Phi$old10<-F
lions.ddl.cjsh2$Phi$old10[lions.ddl.cjsh2$Phi$agebins1246810=="[20,44]"]<-T

lions.ddl.cjsh2$Phi$mold10<-F
lions.ddl.cjsh2$Phi$mold10[lions.ddl.cjsh2$Phi$agebins1246810=="[20,44]"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$fold10<-F
lions.ddl.cjsh2$Phi$fold10[lions.ddl.cjsh2$Phi$agebins1246810=="[20,44]"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$mold10hunt<-F
lions.ddl.cjsh2$Phi$mold10hunt[lions.ddl.cjsh2$Phi$agebins1246810=="[20,44]"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

## hunting as a factor
lions.ddl.cjsh2$Phi$hunt<-F
lions.ddl.cjsh2$Phi$hunt[lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

lions.ddl.cjsh2$Phi$mhunt<-F
lions.ddl.cjsh2$Phi$mhunt[lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"&lions.ddl.cjsh2$Phi$sex==1]<-T


# old 6 refers to an upper age class 6-22 in males
lions.ddl.cjsh2$Phi$old6<-F
lions.ddl.cjsh2$Phi$old6[lions.ddl.cjsh2$Phi$agebinseli %in% c("[12,16)","[16,44]")]<-T

lions.ddl.cjsh2$Phi$mold6<-F
lions.ddl.cjsh2$Phi$mold6[lions.ddl.cjsh2$Phi$agebinseli %in% c("[12,16)","[16,44]")&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$fold6<-F
lions.ddl.cjsh2$Phi$fold6[lions.ddl.cjsh2$Phi$agebinseli %in% c("[12,16)","[16,44]")&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$mold6hunt<-F
lions.ddl.cjsh2$Phi$mold6hunt[lions.ddl.cjsh2$Phi$agebinseli %in% c("[12,16)","[16,44]")&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T


# old 4 refers to an upper age class 4-22 in males
lions.ddl.cjsh2$Phi$old4<-F
lions.ddl.cjsh2$Phi$old4[lions.ddl.cjsh2$Phi$agebinseli %in% c("[8,12)","[12,16)","[16,44]")]<-T

lions.ddl.cjsh2$Phi$mold4<-F
lions.ddl.cjsh2$Phi$mold4[lions.ddl.cjsh2$Phi$agebinseli %in% c("[8,12)","[12,16)","[16,44]")&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$fold4<-F
lions.ddl.cjsh2$Phi$fold4[lions.ddl.cjsh2$Phi$agebinseli %in% c("[8,12)","[12,16)","[16,44]")&lions.ddl.cjsh2$Phi$sex==0]<-T


lions.ddl.cjsh2$Phi$mold4hunt<-F
lions.ddl.cjsh2$Phi$mold4hunt[lions.ddl.cjsh2$Phi$agebinseli %in% c("[8,12)","[12,16)","[16,44]")&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

# adu6 refers to an upper age class 6-10
lions.ddl.cjsh2$Phi$adu6<-F
lions.ddl.cjsh2$Phi$adu6[lions.ddl.cjsh2$Phi$agebins1246810 %in% c("[12,16)","[16,20)")]<-T

lions.ddl.cjsh2$Phi$madu6<-F
lions.ddl.cjsh2$Phi$madu6[lions.ddl.cjsh2$Phi$agebins1246810 %in% c("[12,16)","[16,20)")&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$fadu6<-F
lions.ddl.cjsh2$Phi$fadu6[lions.ddl.cjsh2$Phi$agebins1246810 %in% c("[12,16)","[16,20)")&lions.ddl.cjsh2$Phi$sex==0]<-T

# codes for posthoc model
lions.ddl.cjsh2$Phi$posthoc<-"cubs"
lions.ddl.cjsh2$Phi$posthoc[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"&lions.ddl.cjsh2$Phi$sex==0&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-"fnoncubhunt"
lions.ddl.cjsh2$Phi$posthoc[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-"mnoncubhunt"
lions.ddl.cjsh2$Phi$posthoc[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"&lions.ddl.cjsh2$Phi$sex==0&lions.ddl.cjsh2$Phi$hunting!="[2008,2013)"]<-"fnoncubban"
lions.ddl.cjsh2$Phi$posthoc[lions.ddl.cjsh2$Phi$agebinseli!="[0,4)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting!="[2008,2013)"]<-"mnoncubban"



postsub<-function()
{
    Phi.apost<-list(formula=~-1+posthoc)
    p.cjs.mixbyseasonbyyear=list(formula=~mixture*eliseason*studyyear)
    pi.mix<-list(formula=~1)
    pml <- create.model.list("CJSMixture")
    return(mark.wrapper(pml, data = lions.process.cjsh2, ddl = lions.ddl.cjsh2, adjust = FALSE, output = FALSE))
}


postmodel<-postsub()

postmodel
postmodel[[1]]$results$real

postmodel[[1]]$results$beta




##95% CI's are overlapping for males.  must calculate cl where intervals do not overlap.  
## t value used for 95% is 1.96,qt(0.975,300,lower.tail=T)


betas<-postmodel[[1]]$results$beta$estimate[3:6]
betase<-postmodel[[1]]$results$beta$se[3:6]

t<-qt((1-(0.10/2)),161,lower.tail=T)
survs<-exp(betas)/(1+exp(betas))
lcls<-print(exp(betas-t*betase)/(1+exp(betas-t*betase)))
ucls<-print(exp(betas+t*betase)/(1+exp(betas+t*betase)))


jpeg("D:/Dropbox/projectdave/zambia/lions/paper/post_reals.jpg",width=3,height=3,units="in",res=600)

pcex=0.75
par(mar=c(3,5.5,0.5,0.5))
plot(c(0,6),c(0,1),type="n",ylab="Apparent Survival",yaxt="n",xaxt="n",xlab="",yaxs="i",xaxs="i",bty="n")


points(c(2,1), survs[c(2,1)],pch=c(17,15),col=c("red","blue"))

points(c(5,4), survs[c(4,3)],pch=c(17,15),col=c("red","blue"))


arrows(c(2,1),
       lcls[c(2,1)],
       c(2,1),
       ucls[c(2,1)],
       code=0,col=c("red","blue")) 

arrows(c(5,4),
	   lcls[c(4,3)],
	   c(5,4),
	   ucls[c(4,3)],
	   code=0,col=c("red","blue")) 

axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0),las=2,tcl=0.25)

axis(1,at=c(-1,1.5,4.5,6),labels=c("","female","male",""),tcl=0.25,padj=-1.25,lwd=1)
legend(1.5,0.4,pch=c(15,17),col=c("blue","red"),legend=c("moratorium","hunting"),border=NA,bty="n",ncol=1)


dev.off()










