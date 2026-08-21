#rm(list=ls())  
library(RMark)

setwd("C:/Users/David/Desktop/rtemp")
#setwd("C:/Users/pcadmin/Desktop/rtemp")

## the rmark_input_script must be run before running this script
inputdata$freq<-1

inputdata<-inputdata[,c("ch","freq","sex","ageatfirst")]


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
                                  initial.ages=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,23,24,25,27),
                                  allgroups=T,
                                  mixtures=2,age.unit=2)


## design data
lions.ddl.cjsh2<-make.design.data(lions.process.cjsh2,common.zero=T)


# for eli's categories

agebins1246810<-c(0,1,2,4,6,8,10,22)*2
agebinseli<-c(0,2,4,6,8,22)*2

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
# age structure on detection
lions.ddl.cjsh2$p$sub<-F
lions.ddl.cjsh2$p$sub[lions.ddl.cjsh2$p$agebinseli=="[4,8)"]<-T

lions.ddl.cjsh2$p$old4<-F
lions.ddl.cjsh2$p$old4[lions.ddl.cjsh2$p$agebinseli %in% c("[8,12)","[12,16)","[16,44]")]<-T

lions.ddl.cjsh2$p$fcub2<-F
lions.ddl.cjsh2$p$fcub2[lions.ddl.cjsh2$p$agebinseli=="[0,4)"&lions.ddl.cjsh2$p$sex==0]<-T

lions.ddl.cjsh2$p$fsub<-F
lions.ddl.cjsh2$p$fsub[lions.ddl.cjsh2$p$agebinseli=="[4,8)"&lions.ddl.cjsh2$p$sex==0]<-T

lions.ddl.cjsh2$p$fyad<-F
lions.ddl.cjsh2$p$fyad[lions.ddl.cjsh2$p$agebinseli=="[8,12)"&lions.ddl.cjsh2$p$sex==0]<-T

lions.ddl.cjsh2$p$fadu<-F
lions.ddl.cjsh2$p$fadu[lions.ddl.cjsh2$p$agebinseli=="[12,16)"&lions.ddl.cjsh2$p$sex==0]<-T

lions.ddl.cjsh2$p$fold<-F
lions.ddl.cjsh2$p$fold[lions.ddl.cjsh2$p$agebinseli=="[16,44]"&lions.ddl.cjsh2$p$sex==0]<-T


#Age structure on survival
#cubs
lions.ddl.cjsh2$Phi$cub0<-F
lions.ddl.cjsh2$Phi$cub0[lions.ddl.cjsh2$Phi$agebins1246810=="[0,2)"]<-T

lions.ddl.cjsh2$Phi$cub0hunt<-F
lions.ddl.cjsh2$Phi$cub0hunt[lions.ddl.cjsh2$Phi$agebins1246810=="[0,2)"&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

lions.ddl.cjsh2$Phi$cub1<-F
lions.ddl.cjsh2$Phi$cub1[lions.ddl.cjsh2$Phi$agebins1246810=="[2,4)"]<-T

lions.ddl.cjsh2$Phi$cub1hunt<-F
lions.ddl.cjsh2$Phi$cub1hunt[lions.ddl.cjsh2$Phi$agebins1246810=="[2,4)"&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

lions.ddl.cjsh2$Phi$fcub0<-F
lions.ddl.cjsh2$Phi$fcub0[lions.ddl.cjsh2$Phi$agebins1246810=="[0,2)"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$fcub1<-F
lions.ddl.cjsh2$Phi$fcub1[lions.ddl.cjsh2$Phi$agebins1246810=="[2,4)"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$mcub0<-F
lions.ddl.cjsh2$Phi$mcub0[lions.ddl.cjsh2$Phi$agebins1246810=="[0,2)"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$mcub0hunt<-F
lions.ddl.cjsh2$Phi$mcub0hunt[lions.ddl.cjsh2$Phi$agebins1246810=="[0,2)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

lions.ddl.cjsh2$Phi$mcub1<-F
lions.ddl.cjsh2$Phi$mcub1[lions.ddl.cjsh2$Phi$agebins1246810=="[2,4)"&lions.ddl.cjsh2$Phi$sex==1]<-T

lions.ddl.cjsh2$Phi$mcub1hunt<-F
lions.ddl.cjsh2$Phi$mcub1hunt[lions.ddl.cjsh2$Phi$agebins1246810=="[2,4)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

# cub2 refers to 24 month cub age class
lions.ddl.cjsh2$Phi$cub2<-F
lions.ddl.cjsh2$Phi$cub2[lions.ddl.cjsh2$Phi$agebinseli=="[0,4)"]<-T

lions.ddl.cjsh2$Phi$cub2hunt<-F
lions.ddl.cjsh2$Phi$cub2hunt[lions.ddl.cjsh2$Phi$agebinseli=="[0,4)"&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

lions.ddl.cjsh2$Phi$fcub2<-F
lions.ddl.cjsh2$Phi$fcub2[lions.ddl.cjsh2$Phi$agebinseli=="[0,4)"&lions.ddl.cjsh2$Phi$sex==0]<-T

lions.ddl.cjsh2$Phi$mcub2<-F
lions.ddl.cjsh2$Phi$mcub2[lions.ddl.cjsh2$Phi$agebinseli=="[0,4)"&lions.ddl.cjsh2$Phi$sex==1]<-T


lions.ddl.cjsh2$Phi$mcub2hunt<-F
lions.ddl.cjsh2$Phi$mcub2hunt[lions.ddl.cjsh2$Phi$agebinseli=="[0,4)"&lions.ddl.cjsh2$Phi$sex==1&lions.ddl.cjsh2$Phi$hunting=="[2008,2013)"]<-T

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



# detection model
# fit elis with a hunting term
combos <- function()
{
    # Phi
    Phi.cjs.age=list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+mnoncubhunt)
    
    # p
    #1 single effect of sex
    p.cjs.sex=list(formula=~sex)#
    p.cjs.mixsex=list(formula=~mixture+sex)
    p.cjs.mixbysex=list(formula=~mixture*sex)
    #2 single effect of season
    p.cjs.season=list(formula=~eliseason)  #
    p.cjs.mixseason=list(formula=~mixture+eliseason)#
    p.cjs.mixbyseason=list(formula=~mixture*eliseason)#
    #3 single effect of year
    p.cjs.year=list(formula=~studyyear) #
    p.cjs.mixyear=list(formula=~mixture+studyyear)#
    p.cjs.mixbyyear=list(formula=~mixture*studyyear)#
    #4 year and seaon
    p.cjs.yearseason=list(formula=~studyyear+eliseason)#
    p.cjs.mixyearseason=list(formula=~mixture+studyyear+eliseason)
    p.cjs.mixbyseasonyear=list(formula=~studyyear+mixture*eliseason)
    p.cjs.mixbyyearseason=list(formula=~mixture*studyyear+eliseason)
    p.cjs.mixbyseasonmixbyyear=list(formula=~mixture*eliseason+studyyear*mixture)
    #5 year and season and interaction
    p.cjs.yearbyseason=list(formula=~studyyear*eliseason)#
    p.cjs.mixyearbyseason=list(formula=~mixture+studyyear*eliseason)
    p.cjs.mixbyseasonseasonbyyear=list(formula=~mixture*eliseason+studyyear*eliseason)
    p.cjs.mixbyyearseasonbyyear=list(formula=~mixture*studyyear+studyyear*eliseason)
    p.cjs.mixbyseasonbyyear=list(formula=~mixture*eliseason*studyyear)
    #6  year and sex
    p.cjs.yearsex=list(formula=~studyyear+sex)
    p.cjs.mixyearsex=list(formula=~mixture+studyyear+sex)
    p.cjs.mixbyyearmixbysex=list(formula=~mixture*studyyear+mixture*sex)
    p.cjs.mixbyyearsex=list(formula=~mixture*studyyear+sex)
    p.cjs.mixbysexyear=list(formula=~studyyear+mixture*sex)
    #7 season and sex
    p.cjs.seasonsex=list(formula=~eliseason+sex)
    p.cjs.mixseasonsex=list(formula=~mixture+eliseason+sex)
    p.cjs.mixbyseasonmixbysex=list(formula=~mixture*eliseason+mixture*sex)
    p.cjs.mixbyseasonsex=list(formula=~mixture*eliseason+sex)
    p.cjs.mixbysexseason=list(formula=~eliseason+mixture*sex)
    #8 eli's 8 age/sex classes but modifed formulation to reflect early mort of males
    p.cjs.agesex=list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold)
    p.cjs.mixagesex=list(formula=~mixture+sub+old4+fcub2+fsub+fyad+fadu+fold)
    #9 eli's 8 age/sex classes (modified) and season
    p.cjs.agesexseason=list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+eliseason)
    p.cjs.mixagesexseason=list(formula=~mixture+sub+old4+fcub2+fsub+fyad+fadu+fold+eliseason) 
    p.cjs.mixbyseasonagesex=list(formula=~mixture*season+sub+old4+fcub2+fsub+fyad+fadu+fold) 
    #10 eli's 8 age/sex classes (modified) and year
    p.cjs.agesexyear=list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+studyyear)
    p.cjs.mixagesexyear=list(formula=~mixture+sub+old4+fcub2+fsub+fyad+fadu+fold+studyyear)
    p.cjs.mixbyyearagesex=list(formula=~mixture*studyyear+sub+old4+fcub2+fsub+fyad+fadu+fold)
    #11 eli's 8 age/sex classes (modified) and year and season
    p.cjs.agesexyearseason=list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+studyyear+eliseason)
    p.cjs.mixagesexyearseason=list(formula=~mixture+sub+old4+fcub2+fsub+fyad+fadu+fold+studyyear+eliseason)
    p.cjs.mixbyyearagesexseason=list(formula=~mixture*studyyear+sub+old4+fcub2+fsub+fyad+fadu+fold+eliseason)
    p.cjs.mixbyseasonagesexyear=list(formula=~mixture*eliseason+sub+old4+fcub2+fsub+fyad+fadu+fold+studyyear)
    #12 eli age classes will all sex interactions (10 age/sex classes) 
    p.cjs.eliagebysex=list(formula=~agebinseli*sex)
    p.cjs.mixeliagebysex=list(formula=~mixture+agebinseli*sex)
    p.cjs.mixbysexeliagebysex=list(formula=~mixture*sex+agebinseli*sex)
    #13 null models
    p.cjs.dot=list(formula=~1)
    p.cjs.mix=list(formula=~mixture)
    
    #pi
    pi.cjs.fix=list(formula=~1,fixed=1)
    pi.cjs.mix=list(formula=~1)
    
    
    cml <- create.model.list("CJSMixture")
    
    ## remove models where pi is fixed but mixture is specified in the detection model
    ## or if mixtures are defined but not in the detection model
    cml<-cml[!((grepl("fix",cml$pi)&grepl("mix",cml$p))|(!grepl("fix",cml$pi)&!grepl("mix",cml$p))),]
    
    
    return(mark.wrapper(cml, data = lions.process.cjsh2, ddl = lions.ddl.cjsh2, adjust = FALSE, output = FALSE))
}

mods2 <- combos()
mods2
save(mods2,file="mods2.RData")

model.table(mods2,model.name=F)

#write.csv(print(mods2),"D:/Dropbox/projectdave/zambia/lions/paper/detectionmods.csv",row.names=T)
#write.csv(print(model.table(mods2,model.name=F)),"D:/Dropbox/projectdave/zambia/lions/paper/detectionmods_modelnames.csv",row.names=T)

write.csv(print(mods2),"C:/Dropbox/projectdave/zambia/lions/paper/detectionmods.csv",row.names=T)
write.csv(print(model.table(mods2,model.name=F)),"C:/Dropbox/projectdave/zambia/lions/paper/detectionmods_modelnames.csv",row.names=T)

mods2[[24]]$results$beta
mods2[[24]]$results$real

rm(mods2)
gc()
### MODELLING
## only include models with time as efort likley varied with time and Pledger et al 2003 indicated that failure
## to model effort could lead to false detection of heterogeneity




## survival model

basemod<-function()
{
    #1 base model from rosenblatt et al. one cub age class, same gender effect on both cubs, noncub males split at 4, noncub females split at 4,6 & 8
    Phi.f02468.024       <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold)
    Phi.f02468.024.h0    <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+hunt)
    Phi.f02468.024.hm0   <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+mhunt)
    Phi.f02468.024.hm2   <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+mnoncubhunt)
    Phi.f02468.024.hm4   <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+mold4hunt)
    Phi.f02468.024.hm24  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+msubhunt+mold4hunt)
    Phi.f02468.024.hc2   <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+noncubhunt)
    Phi.f02468.024.hcm2  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt)
    Phi.f02468.024.hcm2f2<-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f02468.024.hcm24 <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+msubhunt+mold4hunt)
    Phi.f02468.024.hcm4  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+mold4hunt)
    Phi.f02468.024.h0m2  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+cub0hunt+mnoncubhunt)
    Phi.f02468.024.h0m4  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fold+cub0hunt+mold4hunt)
    
    #2 one cub age class, same gender effect on both cubs, noncub males split at 4, noncub females split at 4,6 8 & 10
    Phi.f0246810.024       <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10)
    Phi.f0246810.024.h0    <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f0246810.024.hm0   <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f0246810.024.hm2   <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f0246810.024.hm4   <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f0246810.024.hm24  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f0246810.024.hc2   <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f0246810.024.hcm2  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f0246810.024.hcm2f2<-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f0246810.024.hcm4  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f0246810.024.hcm24 <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f0246810.024.h0m2  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f0246810.024.h0m4  <-list(formula=~sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    #3 two cub age classes, same gender effect on both cubs, noncub males split at 4,and noncub females split at 4,6 & 8
    Phi.f02468.0124       <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold)
    Phi.f02468.0124.h0    <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+hunt)
    Phi.f02468.0124.hm0   <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+mhunt)
    Phi.f02468.0124.hm2   <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+mnoncubhunt)
    Phi.f02468.0124.hm4   <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+mold4hunt)
    Phi.f02468.0124.hm24  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+msubhunt+mold4hunt)
    Phi.f02468.0124.hc2   <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+noncubhunt)
    Phi.f02468.0124.hcm2  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt)
    Phi.f02468.0124.hcm2f2<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f02468.0124.hcm4  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+mold4hunt)
    Phi.f02468.0124.hcm24 <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+cub2hunt+msubhunt+mold4hunt)
    Phi.f02468.0124.h0m2  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+cub0hunt+mnoncubhunt)
    Phi.f02468.0124.h0m4  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fold+cub0hunt+mold4hunt)
    
    #4 two cub age classes, same gender effect on both cubs, noncub males split at 4,and noncub females split at 4,6,8 & 10
    Phi.f0246810.0124       <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10)
    Phi.f0246810.0124.h0    <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f0246810.0124.hm0   <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f0246810.0124.hm2   <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f0246810.0124.hm4   <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f0246810.0124.hm24  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f0246810.0124.hc2   <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f0246810.0124.hcm2  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f0246810.0124.hcm2f2<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f0246810.0124.hcm24 <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f0246810.0124.hcm4  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f0246810.0124.h0m2  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f0246810.0124.h0m4  <-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    #5 two cub age classes, gender effect on 2nd year cubs, noncub males split at 4, noncub females split 4,6, 8 & 10
    Phi.f1246810.0124       <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10)
    Phi.f1246810.0124.h0    <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f1246810.0124.hm0   <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f1246810.0124.hm2   <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f1246810.0124.hm4   <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f1246810.0124.hm24  <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f1246810.0124.hc2   <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f1246810.0124.hcm2  <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f1246810.0124.hcm2f2<-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f1246810.0124.hcm4  <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f1246810.0124.hcm24 <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f1246810.0124.h0m2  <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f1246810.0124.h0m4  <-list(formula=~cub1+sub+old4+fcub1+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    
    #6 two cub age classes, no gender effect on cubs, noncub males split at 4, noncub females split 4,6, 8 & 10
    Phi.f246810.0124       <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10)
    Phi.f246810.0124.h0    <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f246810.0124.hm0   <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f246810.0124.hm2   <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f246810.0124.hm4   <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f246810.0124.hm24  <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f246810.0124.hc2   <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f246810.0124.hcm2  <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f246810.0124.hcm2f2<-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f246810.0124.hcm4  <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f246810.0124.hcm24 <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f246810.0124.h0m2  <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f246810.0124.h0m4  <-list(formula=~cub1+sub+old4+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    #7 two cub age classes, gender effect on second years cubs, noncub males split at 4 and 6, noncub females split at 4
    Phi.0124.m1246       <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6)
    Phi.0124.m1246.h0    <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+hunt)
    Phi.0124.m1246.hm0   <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+mhunt)
    Phi.0124.m1246.hm2   <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+mnoncubhunt)
    Phi.0124.m1246.hm4   <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+mold4hunt)
    Phi.0124.m1246.hm24  <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+msubhunt+mold4hunt)
    Phi.0124.m1246.hc2   <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+cub2hunt+noncubhunt)
    Phi.0124.m1246.hcm2  <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+cub2hunt+mnoncubhunt)
    Phi.0124.m1246.hcm2f2<-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.0124.m1246.hcm4  <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+cub2hunt+mold4hunt)
    Phi.0124.m1246.hcm24 <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+cub2hunt+msubhunt+mold4hunt)
    Phi.0124.m1246.h0m2  <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+cub0hunt+mnoncubhunt)
    Phi.0124.m1246.h0m4  <-list(formula=~cub1+sub+old4+mcub1+msub+myad+mold6+cub0hunt+mold4hunt)
    
    #8 two cub age classes, no gender effect on cubs, noncub males split at 4 and 6 and noncub females split at 4, 6, 8 & 10
    Phi.f246810.01246       <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10)
    Phi.f246810.01246.h0    <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f246810.01246.hm0   <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f246810.01246.hm2   <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f246810.01246.hm4   <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f246810.01246.hm24  <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f246810.01246.hc2   <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f246810.01246.hcm2  <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f246810.01246.hcm2f2<-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f246810.01246.hcm4  <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f246810.01246.hcm24 <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f246810.01246.h0m2  <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f246810.01246.h0m4  <-list(formula=~cub1+sub+yad+old6+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    
    #9 two cub age classes, gender effect on 2nd year cubs, noncub males split at 4 and 6, noncub females split at 4, 6, 8, and 10
    Phi.f1246810.01246       <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10)
    Phi.f1246810.01246.h0    <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f1246810.01246.hm0   <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f1246810.01246.hm2   <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f1246810.01246.hm4   <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f1246810.01246.hm24  <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f1246810.01246.hc2   <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f1246810.01246.hcm2  <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f1246810.01246.hcm2f2<-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f1246810.01246.hcm4  <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f1246810.01246.hcm24 <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f1246810.01246.h0m2  <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f1246810.01246.h0m4  <-list(formula=~cub1+sub+yad+old6+fcub1+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    #10 two cub age classes, same gender effect on both cubs, noncub males split at 4 & 6, noncub females split at 4, 6, 8 & 10
    Phi.f0246810.01246			<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10)
    Phi.f0246810.01246.h0		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f0246810.01246.hm0		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f0246810.01246.hm2		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f0246810.01246.hm4		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f0246810.01246.hm24		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f0246810.01246.hc2		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f0246810.01246.hcm2		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f0246810.01246.hcm2f2	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f0246810.01246.hcm4		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f0246810.01246.hcm24	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f0246810.01246.h0m2		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f0246810.01246.h0m4		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    #11 two cub age classes, same gender effect on both cubs, noncub males split at 4 and 6, noncub femalessplit at 4, 6 & 10
    Phi.f024610.01246		<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10)
    Phi.f024610.01246.h0	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+hunt)
    Phi.f024610.01246.hm0	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+mhunt)
    Phi.f024610.01246.hm2	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+mnoncubhunt)
    Phi.f024610.01246.hm4	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+mold4hunt)
    Phi.f024610.01246.hm24	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+msubhunt+mold4hunt)
    Phi.f024610.01246.hc2	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+cub2hunt+noncubhunt)
    Phi.f024610.01246.hcm2	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+cub2hunt+mnoncubhunt)
    Phi.f024610.01246.hcm2f2<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f024610.01246.hcm4	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+cub2hunt+mold4hunt)
    Phi.f024610.01246.hcm24	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f024610.01246.h0m2	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+cub0hunt+mnoncubhunt)
    Phi.f024610.01246.h0m4	<-list(formula=~cub1+sub+yad+old6+fcub2+fsub+fyad+fadu6+fold10+cub0hunt+mold4hunt)
    
    #12 two cub age classes, same gender effect on both cubs, noncub males split at 4, and noncub females split at 4, 6, & 10
    Phi.f024610.0124 		<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10)
    Phi.f024610.0124.h0		<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+hunt)
    Phi.f024610.0124.hm0	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+mhunt)
    Phi.f024610.0124.hm2	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+mnoncubhunt)
    Phi.f024610.0124.hm4	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+mold4hunt)
    Phi.f024610.0124.hm24	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+msubhunt+mold4hunt)
    Phi.f024610.0124.hc2	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+cub2hunt+noncubhunt)
    Phi.f024610.0124.hcm2	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+cub2hunt+mnoncubhunt)
    Phi.f024610.0124.hcm2f2	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f024610.0124.hcm4	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+fold6+cub2hunt+mold4hunt)
    Phi.f024610.0124.hcm24	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+fold6+cub2hunt+msubhunt+mold4hunt)
    Phi.f024610.0124.h0m2	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+cub0hunt+mnoncubhunt)
    Phi.f024610.0124.h0m4	<-list(formula=~cub1+sub+old4+fcub2+fsub+fyad+fadu6+fold10+cub0hunt+mold4hunt)
    
    #13 two cub age classes, same gender effect on both cubs, noncub males split at 4 & 8, noncub females split at 4, 6, 8, & 10
    Phi.f0246810.01248			<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10)
    Phi.f0246810.01248.h0		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f0246810.01248.hm0		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f0246810.01248.hm2		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f0246810.01248.hm4		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f0246810.01248.hm24		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f0246810.01248.hm24  	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f0246810.01248.hc2		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f0246810.01248.hcm2		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f0246810.01248.hcm2f2	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f0246810.01248.hcm4		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f0246810.01248.hcm24	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f0246810.01248.h0m2		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f0246810.01248.h0m4		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    #14 two cub age classes, gender effect on 2nd year cubs, noncub males split at 4 & 8, noncub females split at 4 6, 8, & 10
    Phi.f1246810.01248			<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10)
    Phi.f1246810.01248.h0		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+hunt)
    Phi.f1246810.01248.hm0		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f1246810.01248.hm2		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f1246810.01248.hm4		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f1246810.01248.hm24		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f1246810.01248.hc2		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f1246810.01248.hcm2		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f1246810.01248.hcm2f2	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f1246810.01248.hcm4		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f1246810.01248.hcm24	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f1246810.01248.h0m2		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f1246810.01248.h0m4		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    #15 two cub age classes, no gender effect on cubs, noncub males split at 4 & 8, noncub females split at 4, 6, 8 & 10
    Phi.f246810.01248		<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+fold)
    Phi.f246810.01248.h0	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+fold+hunt)
    Phi.f246810.01248.hm0	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+mhunt)
    Phi.f246810.01248.hm2	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+mnoncubhunt)
    Phi.f246810.01248.hm4	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+mold4hunt)
    Phi.f246810.01248.hm24	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f246810.01248.hc2	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f246810.01248.hcm2	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f246810.01248.hcm2f2<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f246810.01248.hcm4	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f246810.01248.hcm24	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f246810.01248.h0m2	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f246810.01248.h0m4	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fadu8+fold10+cub0hunt+mold4hunt)
    
    #16 two cub age classes, no gender effect on cubs, noncub males split at 4 & 8, noncub females split at 4 & 6 & 8
    Phi.f2468.01248			<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold)
    Phi.f2468.01248.h0		<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+hunt)
    Phi.f2468.01248.hm0		<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+mhunt)
    Phi.f2468.01248.hm2		<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+mnoncubhunt)
    Phi.f2468.01248.hm4		<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+mold4hunt)
    Phi.f2468.01248.hm24	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+msubhunt+mold4hunt)
    Phi.f2468.01248.hc2		<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+cub2hunt+noncubhunt)
    Phi.f2468.01248.hcm2	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt)
    Phi.f2468.01248.hcm2f2	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f2468.01248.hcm4	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+cub2hunt+mold4hunt)
    Phi.f2468.01248.hcm24	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+cub2hunt+msubhunt+mold4hunt)
    Phi.f2468.01248.h0m2	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+cub0hunt+mnoncubhunt)
    Phi.f2468.01248.h0m4	<-list(formula=~cub1+sub+adu4+old+fsub+fyad+fadu+fold+cub0hunt+mold4hunt)
    
    #17 two cub age classes, gender effect on 2nd year cubs, noncub males split at 4 & 8, noncub females slit at 4 & 6 & 8
    Phi.f12468.01248		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold)
    Phi.f12468.01248.h0		<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+hunt)
    Phi.f12468.01248.hm0	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+mhunt)
    Phi.f12468.01248.hm2	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+mnoncubhunt)
    Phi.f12468.01248.hm4	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+mold4hunt)
    Phi.f12468.01248.hm24	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+msubhunt+mold4hunt)
    Phi.f12468.01248.hc2	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+cub2hunt+noncubhunt)
    Phi.f12468.01248.hcm2	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt)
    Phi.f12468.01248.hcm2f2	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f12468.01248.hcm4	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+cub2hunt+mold4hunt)
    Phi.f12468.0248.hcm24	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+cub2hunt+msubhunt+mold4hunt)
    Phi.f12468.01248.h0m2	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+cub0hunt+mnoncubhunt)
    Phi.f12468.01248.h0m4	<-list(formula=~cub1+sub+adu4+old+fcub1+fsub+fyad+fadu+fold+cub0hunt+mold4hunt)
    
    #18 two cub age classes, same gender effect on both cubs, males split on 4 & 8, nocub females split at 4 & 6 & 8
    Phi.f02468.01248		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold)
    Phi.f02468.01248.h0		<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+hunt)
    Phi.f02468.01248.hm0	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+mhunt)
    Phi.f02468.01248.hm2	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+mnoncubhunt)
    Phi.f02468.01248.hm4	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+mold4hunt)
    Phi.f02468.01248.hm24	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+msubhunt+mold4hunt)
    Phi.f02468.01248.hc2	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+cub2hunt+noncubhunt)
    Phi.f02468.01248.hcm2	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt)
    Phi.f02468.01248.hcm2f2	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f02468.01248.hcm4	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+cub2hunt+mold4hunt)
    Phi.f02468.01248.hcm24	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+cub2hunt+msubhunt+mold4hunt)
    Phi.f02468.01248.h0m2	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+cub0hunt+mnoncubhunt)
    Phi.f02468.01248.h0m4	<-list(formula=~cub1+sub+adu4+old+fcub2+fsub+fyad+fadu+fold+cub0hunt+mold4hunt)
    
    #19 two cub age classes, no gender effect on cubs, noncub males split at 4, noncub males split at 4 & 10
    Phi.f2410.0124			<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10)
    Phi.f2410.0124.h0		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+hunt)
    Phi.f2410.0124.hm0		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+mhunt)
    Phi.f2410.0124.hm2		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+mnoncubhunt)
    Phi.f2410.0124.hm4		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+mold4hunt)
    Phi.f2410.0124.hm24		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+msubhunt+mold4hunt)
    Phi.f2410.0124.hc2		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub2hunt+noncubhunt)
    Phi.f2410.0124.hcm2		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub2hunt+mnoncubhunt)
    Phi.f2410.0124.hcm2f2	<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f2410.0124.hcm4		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub2hunt+mold4hunt)
    Phi.f2410.0124.hcm24	<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f2410.0124.h0m2		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub0hunt+mnoncubhunt)
    Phi.f2410.0124.h0m4		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub0hunt+mold4hunt)
    
    #20 two cub age classes, gender effect on 2nd year cubs, noncub males split at 4, noncubs females split at 4 & 10
    Phi.f12410.0124			<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10)
    Phi.f12410.0124.h0		<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+hunt)
    Phi.f12410.0124.hm0		<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+mhunt)
    Phi.f12410.0124.hm2		<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+mnoncubhunt)
    Phi.f12410.0124.hm4		<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+mold4hunt)
    Phi.f12410.0124.hm24	<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+msubhunt+mold4hunt)
    Phi.f12410.0124.hc2		<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+cub2hunt+noncubhunt)
    Phi.f12410.0124.hcm2	<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+cub2hunt+mnoncubhunt)
    Phi.f12410.0124.hcm2f2	<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f12410.0124.hcm4	<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+cub2hunt+mold4hunt)
    Phi.f12410.0124.hcm24	<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f12410.0124.h0m2	<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+cub0hunt+mnoncubhunt)
    Phi.f12410.0124.h0m4	<-list(formula=~cub1+sub+old4+fcub1+fsub+fadu468+fold10+cub0hunt+mold4hunt)
    
    #21 two cub age classes, no gender effect on cubs, noncub males split at 4, noncub females split at 4 & 8 & 10
    Phi.f24810.0124			<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10)
    Phi.f24810.0124.h0		<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+hunt)
    Phi.f24810.0124.hm0		<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+mhunt)
    Phi.f24810.0124.hm2		<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+mnoncubhunt)
    Phi.f24810.0124.hm4		<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+mold4hunt)
    Phi.f24810.0124.hm24	<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+msubhunt+mold4hunt)
    Phi.f24810.0124.hc2		<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+cub2hunt+noncubhunt)
    Phi.f24810.0124.hcm2	<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+cub2hunt+mnoncubhunt)
    Phi.f24810.0124.hcm2f2	<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f24810.0124.hcm4	<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+cub2hunt+mold4hunt)
    Phi.f24810.0124.hcm24	<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+cub2hunt+msubhunt+mold4hunt)
    Phi.f24810.0124.h0m2	<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+cub0hunt+mnoncubhunt)
    Phi.f24810.0124.h0m4	<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+cub0hunt+mold4hunt)
    
    #22 two cub age classes, no gender effect on cubs, noncub males split at 4, noncub femlaes split at 4 & 8
    Phi.f248.0124			<-list(formula=~cub1+sub+old4+fsub+fadu4+fold)
    Phi.f248.0124.h0		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+hunt)
    Phi.f248.0124.hm0		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+mhunt)
    Phi.f248.0124.hm2		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+mnoncubhunt)
    Phi.f248.0124.hm4		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+mold4hunt)
    Phi.f248.0124.hm24		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+msubhunt+mold4hunt)
    Phi.f248.0124.hc2		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+cub2hunt+noncubhunt)
    Phi.f248.0124.hcm2		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+cub2hunt+mnoncubhunt)
    Phi.f248.0124.hcm2f2	<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+cub2hunt+mnoncubhunt+fnoncubhunt)
    Phi.f248.0124.hcm4		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+cub2hunt+mold4hunt)
    Phi.f248.0124.hcm24		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+cub2hunt+msubhunt+mold4hunt)
    Phi.f248.0124.h0m2		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+cub0hunt+mnoncubhunt)
    Phi.f248.0124.h0m4		<-list(formula=~cub1+sub+old4+fsub+fadu4+fold+cub0hunt+mold4hunt)
    
    # detection model
    p.cjs.mixbyseasonbyyear=list(formula=~mixture*eliseason*studyyear)
    
    # mixture model
    pi.mix<-list(formula=~1)
    
    cml <- create.model.list("CJSMixture")
    ## code to subset top models
    #cml<-cml[c(205,200,283),]
    
    return(mark.wrapper(cml, data = lions.process.cjsh2, ddl = lions.ddl.cjsh2, adjust = FALSE, output = FALSE))
    
    
}

basemods <- basemod()
basemods
save(basemods,file="basemods.RData")
load("basemods.RData")

#save(basemods,file="D:/Dropbox/projectdave/zambia/lions/paper/basemods.RData" )
#load("D:/Dropbox/projectdave/zambia/lions/paper/basemods.RData")


#write.csv(print(basemods),"D:/Dropbox/projectdave/zambia/lions/paper/basemods.csv",row.names=T)
write.csv(print(basemods),"C:/Dropbox/projectdave/zambia/lions/paper/basemods.csv",row.names=T)


basemods$model.table<-model.table(basemods,model.name=F)
#write.csv(print(basemods),"D:/Dropbox/projectdave/zambia/lions/paper/basemods_modelname.csv",row.names=T)
write.csv(print(basemods),"C:/Dropbox/projectdave/zambia/lions/paper/basemods_modelname.csv",row.names=T)


write.table(print(basemods),"clipboard",sep="\t", row.names=FALSE)
write.table(basemods[[205]]$results$beta,"clipboard",sep="\t")
## examine a model

basemods[[205]]$results$beta

write.table(print(basemods[[205]]$results$beta),"clipboard",sep="\t")
write.table(print(basemods[[200]]$results$beta),"clipboard",sep="\t")
write.table(print(basemods[[283]]$results$beta),"clipboard",sep="\t")




# est for overdispersion
release.gof(lions.process.cjsh2)


# identify and average top models
topmodids<-row.names(basemods$model.table)[basemods$model.table$DeltaAICc<=2.0]

excludemods<-seq(1,(length(basemods)-1))[!(seq(1,(length(basemods)-1)) %in% as.numeric(topmodids))]

topmods<-remove.mark(basemods,excludemods)

#topmods<-basemods
topmods$model.table<-model.table(topmods,model.name=F)
#pi.mix.Phi.f012410.m24.hm2.p.cjs.mixbyseasonbyyear
#pi.mix.Phi.f012410.m24.hm2.p.cjs.mixbyseasonbyyear

write.table(print(topmods),"clipboard",sep="\t", row.names=FALSE)

topmods<-topsubs

newpiresults<-model.average(topmods,"pi",drop=TRUE)
newpresults<-model.average(topmods,"p",drop=TRUE)


newresults<-model.average(topmods,"Phi",drop=TRUE)
newresults<-newresults[match(unique(newresults$estimate),newresults$estimate),]
newresults$age<-as.numeric(as.character(newresults$age))/2
newresults$hunt<-as.numeric(as.character(newresults$time))<2012.5
newresults$sex<-as.numeric(as.character(newresults$sex))

topmods<-topsubs

modelindex<-NA
k<-0
for(i in c(0,1)){
    for (j in c(0,1,2,4,6,8,10)*2){
        for(m in c(T,F)){
            k=k+1
            modelindex[k]<-lions.ddl.cjsh2$Phi$model.index[which(lions.ddl.cjsh2$Phi$sex==i&lions.ddl.cjsh2$Phi$Age==j&lions.ddl.cjsh2$Phi$hunt==m)[1]] 
        }
    }
}
modelindex<-data.frame(model.index=modelindex,hunt=rep(c(T,F),14),age=rep(c(0,0,1,1,2,2,4,4,6,6,8,8,10,10)*2,2),sex=c(rep(0,14),rep(1,14)))
modavg<-cbind(modelindex,model.average(topmods,indices=modelindex$model.index,vcv=T,drop=TRUE)$estimates)
modavg$age<-modavg$age/2
modavg[modavg$age==0,c("estimate","lcl","ucl")]<-modavg[modavg$age==0,c("estimate","lcl","ucl")]^(12/(12-(0.367*12)))

parorder<-order(modelindex$model.index)





top1<-get.real(topmods[[3]],"Phi",se=T)[(lions.ddl.cjsh2$Phi$model.index %in% modelindex$model.index),][order(parorder),c("estimate","se","lcl","ucl")]
top1<-cbind(modelindex,top1)
top1$age<-as.numeric(as.character(top1$age))/2
top1[top1$age==0,c("estimate","lcl","ucl")]<-top1[top1$age==0,c("estimate","lcl","ucl")]^(12/(12-(0.367*12)))

top2<-get.real(topmods[[2]],"Phi",se=T)[(lions.ddl.cjsh2$Phi$model.index %in% modelindex$model.index),][order(parorder),c("estimate","se","lcl","ucl")]
top2<-cbind(modelindex,top2)
top2$age<-as.numeric(as.character(top2$age))/2
top2[top2$age==0,c("estimate","lcl","ucl")]<-top2[top2$age==0,c("estimate","lcl","ucl")]^(12/(12-(0.367*12)))

top3<-get.real(topmods[[4]],"Phi",se=T)[(lions.ddl.cjsh2$Phi$model.index %in% modelindex$model.index),][order(parorder),c("estimate","se","lcl","ucl")]
top3<-cbind(modelindex,top3)
top3$age<-as.numeric(as.character(top3$age))/2

top4<-get.real(topmods[[1]],"Phi",se=T)[(lions.ddl.cjsh2$Phi$model.index %in% modelindex$model.index),][order(parorder),c("estimate","se","lcl","ucl")]
top4<-cbind(modelindex,top4)
top4$age<-as.numeric(as.character(top4$age))/2



dev.off()
jpeg("D:/Dropbox/projectdave/zambia/lions/paper/CJS_reals.jpg",width=6.5,height=6.5,units="in",res=600)

hoffset=0.15
moffset=0.50
moratoriumcol="blue"
pcex=0.75
par(mfrow=c(2,1),oma=c(1,0,0.5,0.5),mar=c(3,5,0,0))
    plot(c(-0.5,11),c(0,1),type="n",ylab="Apparent Survival",yaxt="n",xaxt="n",xlab="",yaxs="i",xaxs="i",bty="n")
#        with(modavg[modavg$sex==0&modavg$hunt==F,],points(age[age %in% c(0,1,2,4,8,10)],
#                                                          estimate[age %in% c(0,1,2,4,8,10)],
#                                                          pch=19,col=moratoriumcol))
#            with(modavg[modavg$sex==0&modavg$hunt==F,],arrows(age[age %in% c(0,1,2,4,8,10)],
#                                                              lcl[age %in% c(0,1,2,4,8,10)],
#                                                              age[age %in% c(0,1,2,4,8,10)],
#                                                              ucl[age %in% c(0,1,2,4,8,10)],
#                                                              code=0,col=moratoriumcol))

        with(top1[top1$sex==0&top1$hunt==F,],points(c(0.5,1.5,3,7,10.5),
                                                              estimate[age %in% c(0,1,2,4,10)],
                                                              pch=15,col=moratoriumcol,cex=1.5))
            with(top1[top1$sex==0&top1$hunt==F,],arrows(c(0.5,1.5,3,7,10.5),
                                                              lcl[age %in% c(0,1,2,4,10)],
                                                        c(0.5,1.5,3,7,10.5),
                                                              ucl[age %in% c(0,1,2,4,10)],
                                                              code=0,col=moratoriumcol)) 
        # with(top2[top2$sex==0&top2$hunt==F,],points(age[age %in% c(0,1,2,4,10)]+moffset*2,
        #                                    estimate[age %in% c(0,1,2,4,10)],
        #                                    pch=18))
        #     with(top2[top2$sex==0&top2$hunt==F,],arrows(age[age %in% c(0,1,2,4,10)]+moffset*2,
        #                                    lcl[age %in% c(0,1,2,4,10)],
        #                                    age[age %in% c(0,1,2,4,10)]+moffset*2,
        #                                    ucl[age %in% c(0,1,2,4,10)],
        #                                    code=0)) 
        # with(top3[top3$sex==0&top3$hunt==F,],points(age[age %in% c(0,1,2,4,10)]+moffset*3,
        #                                                 estimate[age %in% c(0,1,2,4,10)],
        #                                                 pch=16))
        #     with(top3[top3$sex==0&top3$hunt==F,],arrows(age[age %in% c(0,1,2,4,10)]+moffset*3,
        #                                                 lcl[age %in% c(0,1,2,4,10)],
        #                                                 age[age %in% c(0,1,2,4,10)]+moffset*3,
        #                                                 ucl[age %in% c(0,1,2,4,10)],
        #                                                 code=0)) 
        # with(top4[top4$sex==0&top4$hunt==F,],points(age[age %in% c(0,1,2,4,10)]+moffset*4,
        #                                                 estimate[age %in% c(0,1,2,4,10)],
        #                                                 pch=15))
        #     with(top4[top3$sex==0&top4$hunt==F,],arrows(age[age %in% c(0,1,2,4,10)]+moffset*4,
        #                                                 lcl[age %in% c(0,1,2,4,10)],
        #                                                 age[age %in% c(0,1,2,4,10)]+moffset*4,
        #                                                 ucl[age %in% c(0,1,2,4,10)],
        #                                                 code=0)) 
            
#        with(modavg[modavg$sex==0&modavg$hunt==T,],points(age[age %in% c(0,1,2,4,8,10)]+hoffset,
#                                                              estimate[age %in% c(0,1,2,4,8,10)],
#                                                              pch=19,col="red"))
#            with(modavg[modavg$sex==0&modavg$hunt==T,],arrows(age[age %in% c(0,1,2,4,8,10)]+hoffset,
#                                                              lcl[age %in% c(0,1,2,4,8,10)],
#                                                              age[age %in% c(0,1,2,4,8,10)]+hoffset,
#                                                              ucl[age %in% c(0,1,2,4,8,10)],
#                                                              code=0,col="red"))
        # with(top2[top2$sex==0&top2$hunt==T,],points(age[age %in% c(0,1)]+moffset*2+hoffset,
        #                                                 estimate[age %in% c(0,1)],
        #                                                 pch=18,col="red"))
        #     with(top2[top2$sex==0&top2$hunt==T,],arrows(age[age %in% c(0,1)]+moffset*2+hoffset,
        #                                                 lcl[age %in% c(0,1)],
        #                                                 age[age %in% c(0,1)]+moffset*2+hoffset,
        #                                                 ucl[age %in% c(0,1)],
        #                                                 code=0,col="red")) 
        # with(top4[top4$sex==0&top4$hunt==T,],points(age[age %in% c(0,1,2,4,10)]+moffset*4+hoffset,
        #                                                 estimate[age %in% c(0,1,2,4,10)],
        #                                                 pch=15,col="red"))
        #     with(top4[top4$sex==0&top4$hunt==T,],arrows(age[age %in% c(0,1,2,4,10)]+moffset*4+hoffset,
        #                                                 lcl[age %in% c(0,1,2,4,10)],
        #                                                 age[age %in% c(0,1,2,4,10)]+moffset*4+hoffset,
        #                                                 ucl[age %in% c(0,1,2,4,10)],
        #                                                 code=0,col="red")) 
        axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0),las=2)
        arrows(c(0,1,2,4,10),0,c(1,2,4,10,10.2)-0.05,0,code=1,xpd=T,length=0.045,angle=c(90))
        arrows(c(0,1,2,4),0,c(1,2,4,10)-0.05,0,code=2,xpd=T,length=0.075,angle=c(35))
        arrows(c(10.8),0,c(11),0,code=2,xpd=T,length=0.075,angle=c(35))
        lines(c(10.2,10.4,10.6,10.8),c(0,-0.05,0.05,0),xpd=T)
        #arrows(c(4,8),-0.2,c(8,10),-0.2,code=3,xpd=T,length=0.05,angle=c(90))
        axis(1,c(0,1,2,4,10,11),labels=as.character(c(0,1,2,4,10,14)),tcl=0,padj=-1.25,lwd=0)
        #axis(1,c(8),labels=as.character(c(8)),tcl=0,padj=0.5)
        
        text(5,0.2,"Female",cex=1.5)
        
            
            
    
    plot(c(-0.5,11),c(0,1),type="n",ylab="Apparent Survival",yaxt="n",xaxt="n",xlab="",yaxs="i",xaxs="i",bty="n")
#        with(modavg[modavg$sex==1&modavg$hunt==F,],points(age[age %in% c(0,1,2,4)],
#                                                              estimate[age %in% c(0,1,2,4)],
#                                                              pch=19,col=moratoriumcol))
#            with(modavg[modavg$sex==1&modavg$hunt==F,],arrows(age[age %in% c(0,1,2,4)],
#                                                              lcl[age %in% c(0,1,2,4)],
#                                                              age[age %in% c(0,1,2,4)],
#                                                              ucl[age %in% c(0,1,2,4)],
#                                                              code=0,col=moratoriumcol))
            
#        with(modavg[modavg$sex==1&modavg$hunt==T,],points(age[age %in% c(0,1,2,4)]+hoffset,
#                                                              estimate[age %in% c(0,1,2,4)],
#                                                              pch=19,col="red"))
#            with(modavg[modavg$sex==1&modavg$hunt==T,],arrows(age[age %in% c(0,1,2,4)]+hoffset,
#                                                              lcl[age %in% c(0,1,2,4)],
#                                                              age[age %in% c(0,1,2,4)]+hoffset,
#                                                              ucl[age %in% c(0,1,2,4)],
#                                                              code=0,col="red"))
            
        with(top1[top1$sex==1&top1$hunt==F,],points(c(0.5,1.5,3,7),
                                           estimate[age %in% c(0,1,2,4)],
                                           pch=15,col=moratoriumcol,cex=1.5))
            with(top1[top1$sex==1&top1$hunt==F,],arrows(c(0.5,1.5,3,7),
                                           lcl[age %in% c(0,1,2,4)],
                 c(0.5,1.5,3,7),
                                           ucl[age %in% c(0,1,2,4)],
                                           code=0,col=moratoriumcol)) 
        # with(top2[top2$sex==1&top2$hunt==F,],points(age[age %in% c(0,1,2,4)]+moffset*2,
        #                                                 estimate[age %in% c(0,1,2,4)],
        #                                                 pch=18))
        #     with(top2[top2$sex==1&top2$hunt==F,],arrows(age[age %in% c(0,1,2,4)]+moffset*2,
        #                                                 lcl[age %in% c(0,1,2,4)],
        #                                                 age[age %in% c(0,1,2,4)]+moffset*2,
        #                                                 ucl[age %in% c(0,1,2,4)],
        #                                                 code=0)) 
        # with(top3[top3$sex==1&top3$hunt==F,],points(age[age %in% c(0,1,2,4)]+moffset*3,
        #                                                 estimate[age %in% c(0,1,2,4)],
        #                                                 pch=16))
        #     with(top3[top3$sex==1&top3$hunt==F,],arrows(age[age %in% c(0,1,2,4)]+moffset*3,
        #                                                 lcl[age %in% c(0,1,2,4)],
        #                                                 age[age %in% c(0,1,2,4)]+moffset*3,
        #                                                 ucl[age %in% c(0,1,2,4)],
        #                                                 code=0)) 
        # with(top4[top4$sex==1&top4$hunt==F,],points(age[age %in% c(0,1,2,4)]+moffset*4,
        #                                                 estimate[age %in% c(0,1,2,4)],
        #                                                 pch=15))
        #     with(top4[top4$sex==1&top4$hunt==F,],arrows(age[age %in% c(0,1,2,4)]+moffset*4,
        #                                                 lcl[age %in% c(0,1,2,4)],
        #                                                 age[age %in% c(0,1,2,4)]+moffset*4,
        #                                                 ucl[age %in% c(0,1,2,4)],
        #                                                 code=0))  
            
            
        with(top1[top1$sex==1&top1$hunt==T,],points(c(3,7)+hoffset,
                                                        estimate[age %in% c(2,4)],
                                                        pch=17,col="red",cex=1.5))
            with(top1[top1$sex==1&top1$hunt==T,],arrows(c(3,7)+hoffset,
                                                        lcl[age %in% c(2,4)],
                                                        c(3,7)+hoffset,
                                                        ucl[age %in% c(2,4)],
                                                        code=0,col="red"))   
        # with(top2[top2$sex==1&top2$hunt==T,],points(age[age %in% c(0,1,2,4)]+moffset*2+hoffset,
        #                                                 estimate[age %in% c(0,1,2,4)],
        #                                                 pch=18,col="red"))
        #     with(top2[top2$sex==1&top2$hunt==T,],arrows(age[age %in% c(0,1,2,4)]+moffset*2+hoffset,
        #                                                 lcl[age %in% c(0,1,2,4)],
        #                                                 age[age %in% c(0,1,2,4)]+moffset*2+hoffset,
        #                                                 ucl[age %in% c(0,1,2,4)],
        #                                                 code=0,col="red"))   
        # with(top3[top3$sex==1&top3$hunt==T,],points(age[age %in% c(2,4)]+moffset*3+hoffset,
        #                                                 estimate[age %in% c(2,4)],
        #                                                 pch=16,col="red"))
        #     with(top3[top3$sex==1&top3$hunt==T,],arrows(age[age %in% c(2,4)]+moffset*3+hoffset,
        #                                                 lcl[age %in% c(2,4)],
        #                                                 age[age %in% c(2,4)]+moffset*3+hoffset,
        #                                                 ucl[age %in% c(2,4)],
        #                                                 code=0,col="red"))     
        # with(top4[top4$sex==1&top4$hunt==T,],points(age[age %in% c(2,4)]+moffset*4+hoffset,
        #                                                 estimate[age %in% c(2,4)],
        #                                                 pch=15,col="red"))
        #     with(top4[top4$sex==1&top4$hunt==T,],arrows(age[age %in% c(2,4)]+moffset*4+hoffset,
        #                                                 lcl[age %in% c(2,4)],
        #                                                 age[age %in% c(2,4)]+moffset*4+hoffset,
        #                                                 ucl[age %in% c(2,4)],
        #                                                 code=0,col="red"))  
        # 
        #     
            axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0),las=2)
            
            
            arrows(c(0,1,2,4),0,c(1,2,4,10.21)-0.05,0,code=1,xpd=T,length=0.045,angle=c(90))
            arrows(c(0,1,2),0,c(1,2,4)-0.05,0,code=2,xpd=T,length=0.075,angle=c(35))
            arrows(c(10.8),0,c(11),0,code=2,xpd=T,length=0.075,angle=c(35))
            lines(c(10.2,10.4,10.6,10.8),c(0,-0.05,0.05,0),xpd=T)
            #arrows(c(4,8),-0.2,c(8,10),-0.2,code=3,xpd=T,length=0.05,angle=c(90))
            axis(1,c(0,1,2,4,11),labels=as.character(c(0,1,2,4,12)),tcl=0,padj=-1.25,lwd=0)
            #axis(1,c(8),labels=as.character(c(8)),tcl=0,padj=0.5)
            
            mtext("age class",side=1,xpd=T,line=2.75,cex=2)
            text(5,0.2,"Male",cex=1.5)
            
            legend(7.2,0.45,pch=c(15,17),legend=c("moratorium","hunting"),cex=1.5,col=c("blue","red"),bty="n",ncol=1)
    

dev.off()













