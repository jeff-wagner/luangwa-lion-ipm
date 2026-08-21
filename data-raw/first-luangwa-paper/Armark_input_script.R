rm(list=ls())  
library(RMark)
## This script reads in monthly capture histories from .inp files created in Excel and
## pairs those capture histories with demographic, behavioral, and social data on each lion
## from the MS Access db.  Alignment of these different databases is made on the 7 digit lion ID
##

## seperate R scripts that load the necessary data from the MS Access db must  be run prior to 
## running this script.

## sightings_query
## animasatsightings_query
## slnp_animals_query

# set a working directory, if needed


## specify directory for all .inp files of raw capture histories
#inputdirectory<-"rmark_input_files"
inputdirectory<-"D:/Dropbox/projectdave/zambia/lions/paper/revision/datafiles/rmark_input_files"

#create a list of ch files
inputlist<-lapply(list.files(inputdirectory,pattern=".txt",full.names=T),convert.inp,use.comments=T)

#extract all unique lion names from the inputed list of ch files
uniquelions<-sort(unique(unlist(lapply(inputlist,rownames))))

#create and empty data.frame to hold ch's
inputdata<-data.frame(lionid=uniquelions,ch2008=NA,
					  ch2009=NA,
					  ch2010=NA,
					  ch2011=NA,
					  ch2012=NA,
					  ch2013=NA,
					  ch2014=NA,
					  ch2015=NA)

##  population inputdata with capture histories year by year 
for(i in 1:8){
	inputdata[,i+1]<-inputlist[[i]][match(inputdata$lionid,rownames(inputlist[[i]])),1]
}

## mmight need to trim 2008, 2013, 2014, 2015 which can have 12-month histories 
#inputdata[,"ch2008"]<-substr(inputdata[,"ch2008"],4,11)
#inputdata[,"ch2013"]<-substr(inputdata[,"ch2013"],4,11)
#inputdata[,"ch2014"]<-substr(inputdata[,"ch2014"],4,11)
#inputdata[,"ch2015"]<-substr(inputdata[,"ch2015"],4,11)

## replace all ch's that are NA's with a string of eight zeroes
inputdata[is.na(inputdata)] <- "00000000"

## Read in Lion ID .xlsx file for covariates
library(XLConnect)

#lionid<-loadWorkbook("original_excel_files/Lion ID 2008-2016 with DOBs.xlsx")
lionid<-loadWorkbook("D:/Dropbox/projectdave/zambia/lions/paper/revision/datafiles/original_excel_files/Lion ID 2008-2016 with DOBs.xlsx")

lionids<-readWorksheet(lionid,sheet="All Lions",header=T)

## merge capture histories into one 
chcols<-c("ch2008","ch2009","ch2010","ch2011","ch2012","ch2013","ch2014","ch2015")
inputdata$ch <- apply( inputdata[ , chcols ] , 1 , paste , collapse = "" )

##age at first detection
## this is estimated from the first detection in the capture histories and the DOB in the master Lion ID file
inputdata$yearatfirst<-2007+ceiling(8*(regexpr("1",inputdata$ch)/64))
inputdata$monthatfirst<-3+(regexpr("1",inputdata$ch)-(ceiling(8*(regexpr("1",inputdata$ch)/64))-1)*8)

# convert age in master ID file to an R date
lionids$DOB<-strptime(lionids$DOB,"%Y-%m-%d")

# estimate age at first detection using the month and year of first detection and DOB.  Assume day of month at first detection was the 15th to avoid negative ages.
inputdata$ageatfirst<-as.numeric(difftime(strptime(paste(inputdata$yearatfirst, inputdata$monthatfirst,"28",sep="-"),"%Y-%m-%d"),lionids$DOB[match(inputdata$lionid,lionids$MARKID)],units="days")/365)

# get sex of each lion from master Lion ID file, 1 for male, 0 for female, NA for unknown
inputdata$sex<-NA
inputdata$sex<-lionids$Sex[match(inputdata$lionid,lionids$MARKID)]
inputdata$sex[inputdata$sex=="Male"]<-1
inputdata$sex[inputdata$sex=="Female"]<-0
inputdata$sex<-as.numeric(as.character(inputdata$sex))
inputdata$knownsex<-!is.na(inputdata$sex)
# for lions with unknown gender, assign a gender using a consistent rule to ensure assigned gender will not change from one simulation to another.  ALso want about a 50:50 sex ratio 
# as of 11/14/2016 the rule used was based on the Lion ID number: if the ID number was even, sex was male, if odd, female.
# this resulted in 23 unknown gender lions classified as male and 22 lions classified as female.
# most lions of unknown gender were cubs
inputdata$sex[which(lionids$Sex[match(inputdata$lionid,lionids$MARKID)]=="NA")]<-1*((as.numeric(substr(inputdata$lionid[which(lionids$Sex[match(inputdata$lionid,lionids$MARKID)]=="NA")],7,7)) %%2)==0)
# turn any NA ages at first detection into 3 year-old lions (assuming these were not cubs)
inputdata$ageatfirst[is.na(inputdata$ageatfirst)&inputdata$sex==0]<-median(inputdata$ageatfirst[inputdata$sex==0&inputdata$ageatfirst>=2],na.rm=T)
inputdata$ageatfirst[is.na(inputdata$ageatfirst)&inputdata$sex==1]<-median(inputdata$ageatfirst[inputdata$sex==1&inputdata$ageatfirst>=2],na.rm=T)

mean(inputdata$ageatfirst[inputdata$ageatfirst<1],na.rm=T)

# ELI-HW2 appears to be a 'counted' cub (known litter size but not id'ed cubs) Will give the median age of cubs

inputdata$ageatfirst[inputdata$lionid=="ELI-HW2"]<-0.367
inputdata$sex[inputdata$lionid=="ELI-HW2"]<-1


write.csv(inputdata,"D:/Dropbox/projectdave/zambia/lions/paper/revision/datafiles/inputdata_for_plos.csv")


agehist2008<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2008)>0],breaks=c(0,2,4,6,8,10,20),right=F)
agehist2009<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2009)>0]+(2009-inputdata$yearatfirst[as.numeric(inputdata$ch2009)>0]),breaks=c(0,2,4,6,8,10,20),right=F)
agehist2010<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2010)>0]+(2010-inputdata$yearatfirst[as.numeric(inputdata$ch2010)>0]),breaks=c(0,2,4,6,8,10,20),right=F)
agehist2011<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2011)>0]+(2011-inputdata$yearatfirst[as.numeric(inputdata$ch2011)>0]),breaks=c(0,2,4,6,8,10,20),right=F)
agehist2012<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2012)>0]+(2012-inputdata$yearatfirst[as.numeric(inputdata$ch2012)>0]),breaks=c(0,2,4,6,8,10,20),right=F)
agehist2013<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2013)>0]+(2013-inputdata$yearatfirst[as.numeric(inputdata$ch2013)>0]),breaks=c(0,2,4,6,8,10,20),right=F)
agehist2014<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2014)>0]+(2014-inputdata$yearatfirst[as.numeric(inputdata$ch2014)>0]),breaks=c(0,2,4,6,8,10,20),right=F)
agehist2015<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2015)>0]+(2015-inputdata$yearatfirst[as.numeric(inputdata$ch2015)>0]),breaks=c(0,2,4,6,8,10,20),right=F)

magehist2008<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2008)>0&inputdata$sex==1],breaks=c(0,2,4,6,8,10,20),right=F)
magehist2009<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2009)>0&inputdata$sex==1]+(2009-inputdata$yearatfirst[as.numeric(inputdata$ch2009)>0&inputdata$sex==1]),breaks=c(0,2,4,6,8,10,20),right=F)
magehist2010<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2010)>0&inputdata$sex==1]+(2010-inputdata$yearatfirst[as.numeric(inputdata$ch2010)>0&inputdata$sex==1]),breaks=c(0,2,4,6,8,10,20),right=F)
magehist2011<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2011)>0&inputdata$sex==1]+(2011-inputdata$yearatfirst[as.numeric(inputdata$ch2011)>0&inputdata$sex==1]),breaks=c(0,2,4,6,8,10,20),right=F)
magehist2012<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2012)>0&inputdata$sex==1]+(2012-inputdata$yearatfirst[as.numeric(inputdata$ch2012)>0&inputdata$sex==1]),breaks=c(0,2,4,6,8,10,20),right=F)
magehist2013<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2013)>0&inputdata$sex==1]+(2013-inputdata$yearatfirst[as.numeric(inputdata$ch2013)>0&inputdata$sex==1]),breaks=c(0,2,4,6,8,10,20),right=F)
magehist2014<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2014)>0&inputdata$sex==1]+(2014-inputdata$yearatfirst[as.numeric(inputdata$ch2014)>0&inputdata$sex==1]),breaks=c(0,2,4,6,8,10,20),right=F)
magehist2015<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2015)>0&inputdata$sex==1]+(2015-inputdata$yearatfirst[as.numeric(inputdata$ch2015)>0&inputdata$sex==1]),breaks=c(0,2,4,6,8,10,20),right=F)

fagehist2008<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2008)>0&inputdata$sex==0],breaks=c(0,2,4,6,8,10,20),right=F)
fagehist2009<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2009)>0&inputdata$sex==0]+(2009-inputdata$yearatfirst[as.numeric(inputdata$ch2009)>0&inputdata$sex==0]),breaks=c(0,2,4,6,8,10,20),right=F)
fagehist2010<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2010)>0&inputdata$sex==0]+(2010-inputdata$yearatfirst[as.numeric(inputdata$ch2010)>0&inputdata$sex==0]),breaks=c(0,2,4,6,8,10,20),right=F)
fagehist2011<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2011)>0&inputdata$sex==0]+(2011-inputdata$yearatfirst[as.numeric(inputdata$ch2011)>0&inputdata$sex==0]),breaks=c(0,2,4,6,8,10,20),right=F)
fagehist2012<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2012)>0&inputdata$sex==0]+(2012-inputdata$yearatfirst[as.numeric(inputdata$ch2012)>0&inputdata$sex==0]),breaks=c(0,2,4,6,8,10,20),right=F)
fagehist2013<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2013)>0&inputdata$sex==0]+(2013-inputdata$yearatfirst[as.numeric(inputdata$ch2013)>0&inputdata$sex==0]),breaks=c(0,2,4,6,8,10,20),right=F)
fagehist2014<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2014)>0&inputdata$sex==0]+(2014-inputdata$yearatfirst[as.numeric(inputdata$ch2014)>0&inputdata$sex==0]),breaks=c(0,2,4,6,8,10,20),right=F)
fagehist2015<-hist(inputdata$ageatfirst[as.numeric(inputdata$ch2015)>0&inputdata$sex==0]+(2015-inputdata$yearatfirst[as.numeric(inputdata$ch2015)>0&inputdata$sex==0]),breaks=c(0,2,4,6,8,10,20),right=F)


max(inputdata$ageatfirst[as.numeric(inputdata$ch2015)>0&inputdata$sex==0]+(2015-inputdata$yearatfirst[as.numeric(inputdata$ch2015)>0&inputdata$sex==0]))
max(inputdata$ageatfirst[as.numeric(inputdata$ch2013)>0&inputdata$sex==1]+(2013-inputdata$yearatfirst[as.numeric(inputdata$ch2013)>0&inputdata$sex==1]))



allcounts<-rbind(agehist2008$counts,agehist2009$counts,agehist2010$counts,agehist2011$counts,agehist2012$counts,agehist2013$counts,agehist2014$counts,agehist2015$counts)


femalecount<-rbind(fagehist2008$counts,fagehist2009$counts,fagehist2010$counts,fagehist2011$counts,fagehist2012$counts,fagehist2013$counts,fagehist2014$counts,fagehist2015$counts)
femalecounts<-femalecount/(rowSums(allcounts))

malecount<-rbind(magehist2008$counts,magehist2009$counts,magehist2010$counts,magehist2011$counts,magehist2012$counts,magehist2013$counts,magehist2014$counts,magehist2015$counts)
malecounts<-malecount/(rowSums(allcounts))

inputdata2<-inputdata
inputdata2$agegroups<-"sub"
inputdata2$agegroups[inputdata2$ageatfirst<2]<-"cub"
inputdata2$agegroups[inputdata2$ageatfirst>=4]<-"zadu"
newmales<-with(inputdata2[inputdata2$sex==1,],table(yearatfirst,agegroups))
newfemales<-with(inputdata2[inputdata2$sex==0,],table(yearatfirst,agegroups))

newassignedmales<-with(inputdata2[inputdata2$sex==1&inputdata2$knownsex==F,],table(yearatfirst,agegroups))
newassignedfemales<-with(inputdata2[inputdata2$sex==0&inputdata2$knownsex==F,],table(yearatfirst,agegroups))


jpeg("D:/Dropbox/projectdave/zambia/lions/paper/sexagedist.jpg",width=6.5,height=6.5,units="in",res=600)
par(mfrow=c(3,1),mar=c(0,0,2,0),oma=c(3,5,0.5,0.5),mgp=c(2,0.5,0))
plot(c(0.75,8.75),c(0,0.5),type="n",xaxt="n",yaxt="n",xaxs="i",yaxs="i")
  grid(NA,NULL)  
  rect(1:8,0,1:8+0.25,femalecounts[,1],col="red",border=NA)
  rect(1:8+0.25,0,1:8+0.5,malecounts[,1],col="blue",border=NA)
  axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5),las=2,tcl=0.25)
  axis(1,at=c(-1,10))
  axis(3,at=c(-1,10))
  lines(c(5.75,5.75),c(0,1),lty=2)
  text(1,0.45,"0.00 - 1.99 years",pos=4)
  text(1:8+0.125,femalecounts[,1]+0.06,femalecount[,1])
  text(1:8+0.375,malecounts[,1]+0.06,malecount[,1])
  text(1:8+0.125,femalecounts[,1]+0.02,paste0("(",newfemales[,1],")"))
  text(1:8+0.375,malecounts[,1]+0.02,paste0("(",newmales[,1],")"))
  legend(6.5,0.5,fil=c("red","blue"),legend=c("female","male"),border=NA,bty="n",ncol=2)
  
  
plot(c(0.75,8.75),c(0,0.5),type="n",xaxt="n",yaxt="n",xaxs="i",yaxs="i")
  grid(NA,NULL) 
  rect(1:8,0,1:8+0.25,femalecounts[,2],col="red",border=NA)
  rect(1:8+0.25,0,1:8+0.5,malecounts[,2],col="blue",border=NA)
  axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5),las=2,tcl=0.25)
  mtext("Observed Propoprtion of Known Lions",side=2,line=3) 
  axis(1,at=c(-1,10))
  axis(3,at=c(-1,10))
  lines(c(5.75,5.75),c(0,1),lty=2)
  text(1,0.45,"2.00 - 3.99 years",pos=4)
  text(1:8+0.125,femalecounts[,2]+0.06,femalecount[,2])
  text(1:8+0.375,malecounts[,2]+0.06,malecount[,2])
  text(1:8+0.125,femalecounts[,2]+0.02,paste0("(",newfemales[,2],")"))
  text(1:8+0.375,malecounts[,2]+0.02,paste0("(",newmales[,2],")"))
  
plot(c(0.75,8.75),c(0,0.5),type="n",xaxt="n",yaxt="n",xaxs="i",yaxs="i")
  grid(NA,NULL) 
  rect(1:8,0,1:8+0.25,colSums(rbind(femalecounts[,3],femalecounts[,4],femalecounts[,5])),col="red",border=NA)
  rect(1:8+0.25,0,1:8+0.5,colSums(rbind(malecounts[,3],malecounts[,4],malecounts[,5])),col="blue",border=NA)
  axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5),las=2,tcl=0.25)
  axis(1,at=c(-1,10))
  axis(3,at=c(-1,10))
  axis(1,at=c(1:8+0.25),labels=as.character(2008:2015),las=1,tcl=0)
  lines(c(5.75,5.75),c(0,1),lty=2)
  text(1,0.45,"\u22654.00 years",pos=4)
  text(1:8+0.125,colSums(rbind(femalecounts[,3],femalecounts[,4],femalecounts[,5]))+0.06,colSums(rbind(femalecount[,3],femalecount[,4],femalecount[,5])))
  text(1:8+0.375,colSums(rbind(malecounts[,3],malecounts[,4],malecounts[,5]))+0.06,colSums(rbind(malecount[,3],malecount[,4],malecount[,5])))
  text(1:8+0.125,colSums(rbind(femalecounts[,3],femalecounts[,4],femalecounts[,5]))+0.02,paste0("(",newfemales[,3],")"))
  text(1:8+0.375,colSums(rbind(malecounts[,3],malecounts[,4],malecounts[,5]))+0.02,paste0("(",newmales[,3],")"))
  
dev.off()






## mean age of observed animals
t.test(c(inputdata$ageatfirst[as.numeric(inputdata$ch2008)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2008-inputdata$yearatfirst[as.numeric(inputdata$ch2008)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2009)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2009-inputdata$yearatfirst[as.numeric(inputdata$ch2009)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2010)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2010-inputdata$yearatfirst[as.numeric(inputdata$ch2010)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2011)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2011-inputdata$yearatfirst[as.numeric(inputdata$ch2011)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2012)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2012-inputdata$yearatfirst[as.numeric(inputdata$ch2012)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2013)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2013-inputdata$yearatfirst[as.numeric(inputdata$ch2013)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2014)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2014-inputdata$yearatfirst[as.numeric(inputdata$ch2014)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2015)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2015-inputdata$yearatfirst[as.numeric(inputdata$ch2015)!=0&inputdata$sex==0&inputdata$ageatfirst>=4])))

range(c(inputdata$ageatfirst[as.numeric(inputdata$ch2008)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2008-inputdata$yearatfirst[as.numeric(inputdata$ch2008)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2009)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2009-inputdata$yearatfirst[as.numeric(inputdata$ch2009)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2010)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2010-inputdata$yearatfirst[as.numeric(inputdata$ch2010)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2011)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2011-inputdata$yearatfirst[as.numeric(inputdata$ch2011)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2012)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2012-inputdata$yearatfirst[as.numeric(inputdata$ch2012)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2013)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2013-inputdata$yearatfirst[as.numeric(inputdata$ch2013)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2014)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2014-inputdata$yearatfirst[as.numeric(inputdata$ch2014)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2015)!=0&inputdata$sex==0&inputdata$ageatfirst>=4]+
           (2015-inputdata$yearatfirst[as.numeric(inputdata$ch2015)!=0&inputdata$sex==0&inputdata$ageatfirst>=4])))


t.test(c(inputdata$ageatfirst[as.numeric(inputdata$ch2008)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2008-inputdata$yearatfirst[as.numeric(inputdata$ch2008)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2009)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2009-inputdata$yearatfirst[as.numeric(inputdata$ch2009)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2010)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2010-inputdata$yearatfirst[as.numeric(inputdata$ch2010)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2011)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2011-inputdata$yearatfirst[as.numeric(inputdata$ch2011)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2012)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2012-inputdata$yearatfirst[as.numeric(inputdata$ch2012)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2013)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2013-inputdata$yearatfirst[as.numeric(inputdata$ch2013)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2014)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2014-inputdata$yearatfirst[as.numeric(inputdata$ch2014)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2015)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2015-inputdata$yearatfirst[as.numeric(inputdata$ch2015)!=0&inputdata$sex==1&inputdata$ageatfirst>=4])))

range(c(inputdata$ageatfirst[as.numeric(inputdata$ch2008)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2008-inputdata$yearatfirst[as.numeric(inputdata$ch2008)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2009)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2009-inputdata$yearatfirst[as.numeric(inputdata$ch2009)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2010)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2010-inputdata$yearatfirst[as.numeric(inputdata$ch2010)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2011)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2011-inputdata$yearatfirst[as.numeric(inputdata$ch2011)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2012)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2012-inputdata$yearatfirst[as.numeric(inputdata$ch2012)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2013)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2013-inputdata$yearatfirst[as.numeric(inputdata$ch2013)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2014)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2014-inputdata$yearatfirst[as.numeric(inputdata$ch2014)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]),
         inputdata$ageatfirst[as.numeric(inputdata$ch2015)!=0&inputdata$sex==1&inputdata$ageatfirst>=4]+
           (2015-inputdata$yearatfirst[as.numeric(inputdata$ch2015)!=0&inputdata$sex==1&inputdata$ageatfirst>=4])))


## sex ratios in cubs
newassignedfemales
newassignedmales
## these count reflect the number of known sex cubs (thats why different from fig 2 in paper)
femcubhunt<-c(15,9,2,11,7)
malecubhunt<-c(11,9,8,6,7)

femcubban<-c(21,27,13)
malecubban<-c(10,10,17)

sum(femcubhunt)/sum(malecubhunt)

sum(femcubban)/sum(malecubban)

binom.test(sum(femcubhunt),sum(femcubhunt+malecubhunt))
#lcl
1/((1/0.4066)-1)
1/((1/0.6274)-1)

binom.test(sum(femcubban),sum(femcubban+malecubban))
#lcl ucl
1/((1/0.5189)-1)
1/((1/0.7184)-1)


