#Summary of detection data
lapply(lapply(lions.process.cjsh2$data$ch,unlist),grep,"1")
do.call("max",c(lapply(lapply(lions.process.cjsh2$data$ch,unlist),str_count,"1")))
do.call("sum",c(lapply(lapply(lions.process.cjsh2$data$ch,unlist),str_count,"1")))

do.call("sum",c(lapply(lapply(lions.process.cjsh2$data$ch,unlist),str_count,"1")))

#number of intervals from time of First detection and end of study
sum(32-unlist(lapply(gregexpr(lions.process.cjsh2$data$ch,pattern="1"),min)))


## litter size tests



#Top models
library(RMark)

topsub<-function()
{

Phi.f2410.0124.hm2		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+mnoncubhunt)
Phi.f2410.0124.hc2		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub2hunt+noncubhunt)
Phi.f2410.0124.hcm2		<-list(formula=~cub1+sub+old4+fsub+fadu468+fold10+cub2hunt+mnoncubhunt)
Phi.f24810.0124.hm2		<-list(formula=~cub1+sub+old4+fsub+fadu4+fadu8+fold10+mnoncubhunt)

  # detection model
p.cjs.mixbyseasonbyyear=list(formula=~mixture*eliseason*studyyear)

# mixture model
pi.mix<-list(formula=~1)

cml <- create.model.list("CJSMixture")
## code to subset top models
#cml<-cml[c(205,200,283),]

return(mark.wrapper(cml, data = lions.process.cjsh2, ddl = lions.ddl.cjsh2, threads=1, adjust = FALSE, output = FALSE))


}

topsubs <- topsub()
topsubs
save(topsubs,file="C:/Dropbox/projectdave/zambia/lions/paper/topsubs.RData")
#save(topsubs,file="D:/Dropbox/projectdave/zambia/lions/paper/topsubs.RData" )
load("D:/Dropbox/projectdave/zambia/lions/paper/topsubs.RData")

write.csv(print(topsubs),"C:/Dropbox/projectdave/zambia/lions/paper/topsubs.csv",row.names=T)


topsubs$model.table<-model.table(topsubs,model.name=F)
write.csv(print(topsubs),"C:/Dropbox/projectdave/zambia/lions/paper/topsubs_modelname.csv",row.names=T)


write.table(print(topsubs),"clipboard",sep="\t", row.names=FALSE)






# a posterior model
postsub<-function()
{
Phi.apost<-list(formula=~cub1+noncub+mnoncub+mnoncubhunt+fnoncubhunt)
p.cjs.mixbyseasonbyyear=list(formula=~mixture*eliseason*studyyear)
pi.mix<-list(formula=~1)
pml <- create.model.list("CJSMixture")
return(mark.wrapper(pml, data = lions.process.cjsh2, ddl = lions.ddl.cjsh2, adjust = FALSE, output = FALSE))
}

postmodel<-postsub()

postmodel
postmodel[[1]]$results$real

postmodel[[1]]$results$beta

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
parorder<-order(modelindex$model.index)

post1<-get.real(postmodel[[1]],"Phi",se=T)[(lions.ddl.cjsh2$Phi$model.index %in% modelindex$model.index),][order(parorder),c("estimate","se","lcl","ucl")]
post1<-cbind(modelindex,post1)
post1$age<-as.numeric(as.character(post1$age))/2
post1[post1$age==0,c("estimate","lcl","ucl")]<-post1[post1$age==0,c("estimate","lcl","ucl")]^(12/(12-(0.367*12)))




jpeg("D:/Dropbox/projectdave/zambia/lions/paper/post_reals.jpg",width=3,height=3,units="in",res=600)

pcex=0.75
par(mar=c(3,5.5,0.5,0.5))
plot(c(0,6),c(0,1),type="n",ylab="Apparent Survival",yaxt="n",xaxt="n",xlab="",yaxs="i",xaxs="i",bty="n")


points(c(2,1), post1$estimate[post1$sex==0&post1$age==4],pch=c(17,15),col=c("red","blue"))

points(c(5,4), post1$estimate[post1$sex==1&post1$age==4],pch=c(17,15),col=c("red","blue"))


arrows(c(2,1),
       post1$lcl[post1$sex==0&post1$age==4],
       c(2,1),
       post1$ucl[post1$sex==0&post1$age==4],
       code=0,col=c("red","blue")) 

arrows(c(5,4),
       post1$lcl[post1$sex==1&post1$age==4],
       c(5,4),
       post1$ucl[post1$sex==1&post1$age==4],
       code=0,col=c("red","blue")) 

axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0),las=2,tcl=0.25)

axis(1,at=c(-1,1.5,4.5,6),labels=c("","female","male",""),tcl=0.25,padj=-1.25,lwd=1)
legend(1.5,0.4,pch=c(15,17),col=c("blue","red"),legend=c("moratorium","hunting"),border=NA,bty="n",ncol=1)


dev.off()

