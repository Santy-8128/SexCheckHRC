
initial.dir<-getwd()


# INPUT VARIABLES
                 
Args <- commandArgs(TRUE)
if (length(Args)!=3){
        print ("Usage: StudyName BamFileList PathToSamtools")
        quit()}
cohort=Args[1]
dat <- Args[2]
path <- Args[3]



dat=read.table(dat,colClasses="character")[,1]

l=length(dat)
system(paste("rm ",cohort,"_Orig_sample",sep=""))

print(" Generating Summary Statistics for Each Sample ... \n\n")

for(i in 1:l)
{

	if(i%%20==0)
	{
		print(paste(" Analyzing Sample ",i,sep=""))
	}
	system(paste(path," idxstats ",dat[i]," > ",cohort,"_sample_",i," ",sep=""))
	system(paste(path," view -H ",dat[i]," | tr '\t' '\n' | grep -m1 -w SM | cut -d: -f 2 >> ",cohort,"_Orig_sample",sep=""))

}

	
print(paste(" Generating Graph for Study : SEX_CHECK_",cohort,".pdf",sep=""))



pdf(paste("SEX_CHECK_",cohort,".pdf",sep=""))



orignam=read.table(paste(cohort,"_Orig_sample",sep=""),colClasses="character")[,1]
	
		
	
t="cat "
for(i in 1:l)
t=paste(t," ",cohort,"_sample_",i,sep="")

t=paste(t," | grep -v GL | grep -v NC | grep -v MT > DERP",sep="")
system(t)
	
	



for(j in c(1:22,"X","Y"))
{

write.table(j,"FILE",quote=FALSE,col.names=FALSE,row.names=FALSE)
 system(paste("perl GREP_CODE.pl -f FILE -t DERP -c 1 -o ",cohort,"_Summary_Chr_",j," ",sep=""))
}

FINAL=dat

dat=list()
LIST=c(1:22,"X","Y")

for(jj in 1:length(LIST))
{

j=LIST[jj]

dat[[jj]]=apply(read.table(paste(cohort,"_Summary_Chr_",j,sep=""),colClasses="character")[,c(2:4)],2,as.numeric)
}
	

Frac=list()

for(jj in 1:length(LIST))
{

Frac[[jj]]= dat[[jj]][,2]/dat[[jj]][,1]

}

FracSum=0*(1:length(Frac[[1]]))
for(jj in 1:22)
{

FracSum=FracSum+Frac[[jj]]

}

FracXNormFull=Frac[[23]]/(FracSum)
FracYNormFull=Frac[[24]]/(FracSum)
mydata=data.frame(cbind(FracXNormFull,FracYNormFull))

gg=union(which(is.na(mydata[,1])==TRUE) , which(is.na(mydata[,1])==TRUE) )



if(length(gg)>0)
mydata=mydata[-gg,]

fit <- kmeans(mydata, 2) # 5 cluster solution
# get cluster means
#aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster) 


# vall=c(">10",">5-10",">1-5",">0.5-1",">0.1-0.5",">0.07-0.1",">0.03-0.07","<=0.03")
    plot(FracXNormFull[-gg],FracYNormFull[-gg],col=mydata[,3],xlab="Average X Coverage / Autosome",ylab="Average Y Coverage / Autosome",main=paste("SEX CHECK - ",cohort,sep=""))
	
	# legend(1,9,vall,col=c("red","orangered","orange","green4","lawngreen","navyblue","royalblue1","purple3"),title="% of Variants",pch=15,cex=1.2)

	dev.off()

male=which(mydata[,3]==1)

g=aggregate(mydata,by=list(fit$cluster),FUN=mean)
if(g[2,3]>g[1,3])
	male=which(mydata[,3]==2)


print(paste(" Generating Summary Results  for Study : SEX_CHECK_",cohort,".txt",sep=""))


	TOP=NULL
if(length(gg)>0)
{
	TOP=cbind(FINAL[gg],orignam[gg],"NA")
	BOT=cbind( (FINAL[-gg])[male] ,(orignam[-gg])[male] ,"M")
	BOT2=cbind( (FINAL[-gg])[-male] ,(orignam[-gg])[-male] ,"F")
} else
{
	BOT=cbind( FINAL[male] ,orignam[male] ,"M")
	BOT2=cbind( FINAL[-male] ,orignam[-male] ,"F")

}




write.table(rbind(TOP,BOT,BOT2),paste("SEX_CHECK_",cohort,".txt",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)




