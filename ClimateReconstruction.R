# R codes below are used for climate reconstruction of past climate using Hinoki cypress tree ring width:
# Tree ring ID: '100$19M', Period: 1734-2006
# Climate Records: 
#   Site: Matsumoto

# Load R packages:
library(RMySQL)
library(dplyr)
library(tseries)
library(bootRes)
library(reshape)
library(reshape2)
library(dplR)



# clear workspace
rm(list=ls())

# Tree ring data:
con=dbConnect(RMySQL::MySQL(), host = "localhost",user = "root", password = "timberland",dbname="HinokiCypress")
dbListTables(con)
dbListFields(con,'treeringsource')
dbListFields(con,'treeringwidthdata')

rs=dbSendQuery(con,"SELECT Years,RawWidth,TreeRingName FROM treeringsource, treeringwidthdata WHERE TreeRingID=ringID AND TreeRingName='100$19M' ORDER BY Years;")
ring=dbFetch(rs,n=-1); ring=ring[,c(1,2)]; colnames(ring)=c("year","width")
dbClearResult(rs)
dbDisconnect(con)

# Climate Data:
con=dbConnect(RMySQL::MySQL(),host='localhost',user='root',password='timberland',dbname='HinokiCypress')
dbListTables(con)
dbListFields(con,'climatedata'); dbListFields(con,'climatestation'); dbListFields(con,'climatevariable')

rs=dbSendQuery(con,"SELECT Years,Months,variable,yvalue,SiteNames FROM climatedata,climatestation,climatevariable 
               WHERE sitid=SiteID AND varid=VariableID AND SiteNames='Matsumoto' ORDER BY Years;")
clim=dbFetch(rs,n=-1)
dbClearResult(rs)
dbDisconnect(con)


# Ring width:----
# Detrend a raw ring width by Ar and get residual chronology
year=ring[,"year"]
width=ring[,"width"]

series=data.frame(width,row.names=year)
names(width)=rownames(series)

#series.rwi=detrend.series(y=width,y.name="100-19M",verbose=TRUE,nyrs=20)

pdf(file="detrend_width.pdf",width=12,height=9)
series.rwi=detrend.series(y=width,y.name="100-19M",verbose=TRUE,nyrs=20)
dev.off()

gg=data.frame(year=rownames(series.rwi),series.rwi,row.names=NULL)
ring=merge(ring,gg,by="year") # Use this ring data for all analyses


# Fit AR1 for Residual Series
AR=select.list(names(ring[,2:ncol(ring)]),multiple=TRUE,title="Select Type of Rind Width Series For Residual Series:",graphics=TRUE)

x=ring[,names(ring) %in% AR]
M=arma(x,order=c(1,0))

acf(residuals(M),na.action=na.remove)
X=data.frame(ring,ar1=residuals(M)) # Use "residual chronology" derived from selected type of ring width series
colnames(X)[c(2:ncol(X))]=c("raw","spline","modnegexp","means","ar","res")



# Climate Data:----
# NDM0 = the number of days of daily min. temp. below 0 Celsius (excluding 0)　（日最低気温0度未満の日数）
# NDA0 = the number of days of average daily temp. below 0 Celsius　（日平均気温0度未満の日数）
# ADMM = Monthly average daily min. temp. (日最低気温の月平均）
# ADMX = Monthly average daily max. temp. (日最高気温の月平均）
# ADMMX = ADMX - ADMM
# ADTM = Average temp.　（平均気温）
# MTPP = Monthly total precipitation

# Delete missing values (i.e., -999.0):
Y=clim
Y$yvalue[Y$yvalue==-999.0]=NA

# Choose Climate Variables:----
colnames(Y)[1:2]=c("year","month");Y=Y[,c(1:4)]
Y=dcast(Y,year+month~variable)

# BootRes:
# Choose only up to 2 cliamte variables:
Clim=select.list(names(Y[,3:ncol(Y)]),multiple=TRUE,title="Select Climate Variables of Your interest:",graphics=TRUE)
Y1=Y[,names(Y) %in% c("year","month",Clim)]

# Select Type of Ring Width (Choose only one)
Ring.list=select.list(names(X[,2:ncol(X)]),multiple=TRUE,title="Select Rind Width Series of Xour Interest:",graphics=TRUE)
X1=X[,names(X) %in% c("year",Ring.list)]

# Identify Strength of Monthly Signals with Tree Ring----
X1=data.frame(X1[,1:2],row.names="year")


# View Important Climate Variables
op=par(mar=c(5,5,6,3))
dc.corr <- dcc(X1,Y1,method = "corr")
dcplot(dc.corr)
par(op)


## Correlation and Linear Regression
# Reshape Climte dataset
## Use only one climate variable:
Clim=select.list(names(Y[,3:ncol(Y)]),multiple=TRUE,title="Select only one climate variable:",graphics=TRUE)
Y1=Y[,names(Y) %in% c("year","month",Clim)]
Y1=recast(Y1,year~variable+month,id.var=c("year","month"),na.rm=TRUE)

# Convert to time series object
Y1=ts(Y1,frequency=1,start=min(Y1$year),end=max(Y1$year))

Y1=cbind(p=Y1,c=lag(Y1))
Y1=as.data.frame(Y1)

# Rename variables for ease of interpretation
colnames(Y1)=c("year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec",
               "p.year","p.Jan","p.Feb","p.Mar","p.Apr","p.May","p.Jun","p.Jul","p.Aug","p.Sep","p.Oct","p.Nov","p.Dec")
ind=which(colnames(Y1)=="p.year")

Y1=Y1[-1,c(1:(ind-1),(ind+1):ncol(Y1))]


# Select Months to be Averaged
Month.list=select.list(names(Y1[,c(2:ncol(Y1))]),multiple=TRUE,title="Choose Months to be Averaged:",graphics=TRUE)
if(length(Month.list)==1) {Y2=transform(Y1,means=Y1[,c(Month.list)])} else{Y2=transform(Y1,means=rowMeans(Y1[,c(Month.list)]))}

# ID time period
time=paste(min(Y2$year),max(Y2$year),sep="-")
timePeriod=matrix(unlist(strsplit(time,"-")),length(time),2,byrow=TRUE)

# Ring width do not seem reliable after 1998 based on large residuals after linear fit. so remove 
# ring width after 1998
YearCut=1998

# Delete rows where missing values are observed
Y2=na.omit(Y2)
Y2=subset(Y2,year>=as.numeric(timePeriod[1,1]) & year<YearCut)

# Select Ring Width during the specified Period
#X1=subset(X,year>=as.numeric(timePeriod[1,1]) & year<=as.numeric(timePeriod[1,2]))
X1=subset(X,year>=as.numeric(timePeriod[1,1]) & year<YearCut)
Ring.list=select.list(names(X1[,2:ncol(X1)]),multiple=TRUE,title="Select Rind Width Series of Xour Interest:",graphics=TRUE)
X1=X1[,names(X1) %in% c("year",Ring.list)]

# Plot
RingName=names(X1)[2]

op=par(mar=c(5,5,4,5))

c=paste(Month.list,collapse=", ")
plot(Y2$year,Y2$means,type="l",xlab="Year",ylab=Clim,main=paste("Months: ",c,sep=""))
par(new=TRUE)
plot(X$year,X[,2],col="blue",type="l",axes = FALSE, bty = "n", xlab = "", ylab = "") # Tree ring
axis(side=4, at = pretty(range(X[,2])))
mtext(RingName, side=4, line=3)


par(op)


# Run correlation and linear regression 
XX=merge(X1,Y2,by="year")
cor(XX[,2],XX$means,method="pearson")

M0=lm(XX$means~XX[,2])
summary(M0) 

op=par(mfrow=c(2,2))
plot(M0)
par(op)
####################################################











## 1B-1: When You Chose "Daily Min. Temp. Matsumoto 1898-2006.csv-----
# Note: Minimum temperatures above 0 are already omitted from the dataset.
# Convert date to Date class
head(Climate)
Climate$date=as.Date(Climate$date)
Climate$time=as.Date(Climate$date,"%Y-%m-%d")
Climate$week=1
Climate$week[Climate$day>=15 & Climate$day<=31]=2


# Convert min. temp to numeric variable
Climate$min.temp=as.numeric(Climate$min.temp)



###########################################
## 2A: Daily Minimum Temperature Matsumoto 1898-2006:----
## 2A-1: Sum of Daily Minimum Temperature:----
library(reshape2)

X1=Climate[,c("min.temp","year","month")]
X1$min.temp[X1$min.temp>=0]=0
X1=recast(X1,year~variable+month,sum,id.var=c("year","month"),na.rm=FALSE) # Sum daily min temp
colnames(X1)=c("year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Select Type of Ring Width
Y=ring
Ring.list=select.list(names(Y[,2:ncol(Y)]),multiple=TRUE,title="Select Rind Width Series of Your Interest:",graphics=TRUE)
List="year"
Y=Y[,names(Y) %in% c(List,Ring.list)]
Y=data.frame(Y[,1:2],row.names="year")


## Prepare a table descring monthly sum of daily min. temp.
### convert to time series object
X=X1
X=ts(X,frequency=1,start=min(X$year),end=max(X$year))
X=cbind(p=X,c=lag(X))

# Convert back to data.frame
X=as.data.frame(X)

# rename variables for ease of interpretation
colnames(X)=c("p.year","p.Jan","p.Feb","p.Mar","p.Apr","p.May","p.Jun","p.Jul","p.Aug","p.Sep","p.Oct","p.Nov","p.Dec",
              "year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
ind=which(colnames(X)=="year")
X=X[,c(ind,2:(ind-1),(ind+1):ncol(X))]

# Select Months to be Summed or averaged
ind=which(colnames(X)=="year")
Month.list=select.list(names(X[,c(2:ncol(X))]),multiple=TRUE,title="Choose Months to be Summed:",graphics=TRUE)
Month.list

X=transform(X,y=rowSums(X[,c(Month.list)])) # Sum
X=transform(X,y=abs(rowSums(X[,c(Month.list)]))) # Absolute Sum
X=transform(X,y=rowMeans(X[,c(Month.list)])) # Average


# Delete rows where missing values are observed
X=na.omit(X)

# ID time period
time=paste(min(X$year),max(X$year),sep="-")
timePeriod=matrix(unlist(strsplit(time,"-")),length(time),2,byrow=TRUE)
timePeriod

# Response to be used for correlation and linear regression
# this is wrong. how to exlude "0" from calcuation? mean(sum) shoud not include "0"
X=transform(X,y1=y) # Raw: 1
X=transform(X,y1=mean(y)/y) # 2
X=transform(X,y1=y/mean(y)) # 3
X=transform(X,y1=scale(y)) # standardize: 4


# Select Type of Ring Width
Y=ring
Y=subset(ring,year>=as.numeric(timePeriod[1,1]) & year<=as.numeric(timePeriod[1,2]))
Ring.list=select.list(names(Y[,2:ncol(Y)]),multiple=TRUE,title="Select Rind Width Series of Your Interest:",graphics=TRUE)
List="year"
Y=Y[,names(Y) %in% c(List,Ring.list)]
names(Y)

# Plot
op=par(mar=c(5,5,4,5))
plot(X$year,X$y1,type="l")
par(new=TRUE)
plot(Y$year,Y[,2],col="blue",type="l",axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(Y[,2])))
mtext("Y$width", side=4, line=3)

# Run correlation and linear regression
cor(X$y1,Y[,2],method="spearman") # Spearman should be used when variables are not normally distributed
M0=lm(X$y1~Y[,2])
summary(M0) 

############################

## PCA Analysis
X=X1
X=ts(X,frequency=1,start=min(X$year),end=max(X$year))
X=cbind(p=X,c=lag(X))

# Convert back to data.frame
X=as.data.frame(X)

# rename variables for ease of interpretation
colnames(X)=c("p.year","p.Jan","p.Feb","p.Mar","p.Apr","p.May","p.Jun","p.Jul","p.Aug","p.Sep","p.Oct","p.Nov","p.Dec",
              "year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
ind=which(colnames(X)=="year")
X=X[,c(ind,2:(ind-1),(ind+1):ncol(X))]
X=na.omit(X)

# Select Months for PCA
Month.list=select.list(names(X[,c(2:ncol(X))]),multiple=TRUE,title="Choose Months Used for PCA:",graphics=TRUE)
List="year"
X=X[,names(X) %in% c(List,Month.list)]
Month.list
head(X)

# Select Type of Ring Width
Y=ring
Ring.list=select.list(names(Y[,2:ncol(Y)]),multiple=TRUE,title="Select Rind Width Series of Your Interest:",graphics=TRUE)
List="year"
Y=Y[,names(Y) %in% c(List,Ring.list)]
#Y=data.frame(Y[,1:2],row.names="year")

head(X);head(Y)
X2=merge(Y,X,by="year")

# fit PCA
library(FactoMineR)

RES=PCA(X2[,3:ncol(X2)],quanti.sup=2,scale.unit=TRUE,graph=FALSE,ncp=4)
RES$eig
RES$var$cor

Score=RES$ind$coord
X3=cbind(X2[,1:2],Score)
X3=X3[,-1]

M=lm(ModNegExp~., data=X3)
summary(M)
M1=update(M,.~.-Dim.3)
summary(M1)
M2=update(M1,.~.-Dim.4)
summary(M2)








### 2A-2-II. Analyze Bi-weekly----

# Sum daily min. temp. by year and month
library(reshape2)

X=Climate
X=X[,-c(1,5,6)]
X$min.temp[X$min.temp>=0]=0
X=recast(X,year~variable+month+week,sum,id.var=c("year","month","week"),na.rm=FALSE)
head(X)

## Prepare a table descring monthly sum of daily min. temp.
### convert to time series object
X=ts(X,frequency=1,start=min(X$year),end=max(X$year))
X=cbind(X,lag(X))

### Convert back to data.frame
X=as.data.frame(X)
names(X)

### rename variables for ease of interpretation
colnames(X)=c("lag.year","l.Jan.1","l.Jan.2","l.Feb.1","l.Feb.2","l.Mar.1","l.Mar.2","l.Apr.1","l.Apr.2",
              "l.May.1","l.May.2","l.Jun.1","l.Jun.2","l.Jul.1","l.Jul.2","l.Aug.1","l.Aug.2","l.Sep.1","l.Sep.2",
              "l.Oct.1","l.Oct.2","l.Nov.1","l.Nov.2","l.Dec.1","l.Dec.2","year",
              "Jan.1","Jan.2","Feb.1","Feb.2","Mar.1","Mar.2","Apr.1","Apr.2","May.1","May.2",
              "Jun.1","Jun.2","Ju1","Ju2","Aug.1","Aug.2","Sep.1","Sep.2","Oct.1","Oct.2","Nov.1","Nov.2","Dec.1","Dec.2")
ind=which(colnames(X)=="year")
X=X[,c(ind,2:(ind-1),(ind+1):ncol(X))]

# Select Months to be Summed
ind=which(colnames(X)=="year")
Month.list=select.list(names(X[,c(2:ncol(X))]),multiple=TRUE,title="Choose Months to be Summed:",graphics=TRUE)

if(length(Month.list)==1) {
  X=transform(X,sum=abs(X[,c(Month.list)]))
} else{
  X=transform(X,sum=abs(rowSums(X[,c(Month.list)])))
}


# Delete rows where missing values are observed
X=na.omit(X)

# ID time period
time=paste(min(X$year),max(X$year),sep="-")
timePeriod=matrix(unlist(strsplit(time,"-")),length(time),2,byrow=TRUE)

# Ceter sum with a column mean
X=transform(X,c.sum=mean(sum)/sum)
X=transform(X,c.sum=sum/mean(sum))

# Select Type of Ring Width
Y=subset(ring,year>=as.numeric(timePeriod[1,1]) & year<=as.numeric(timePeriod[1,2]))
Ring.list=select.list(names(Y[,2:ncol(Y)]),multiple=TRUE,title="Select Rind Width Series of Your Interest:",graphics=TRUE)
List="year"
Y=Y[,names(Y) %in% c(List,Ring.list)]

# Plot
head(X);head(Y)
X
op=par(mar=c(5,5,4,5))
plot(X$year,X$c.sum,type="l")
par(new=TRUE)
plot(Y$year,Y[,2],col="blue",type="l",axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(Y[,2])))
mtext("Y$width", side=4, line=3)

# Run correlation and linear regression
cor(X$c.sum,Y[,2])
M0=lm(X$c.sum~Y[,2])
summary(M0) 







# # # # # # # # # # # # # # # # # # # # # # # # 
## PCA and Random Forest
# ID important monthly climate variables
## PCA Analysis
X=Climate[,-1] # drop date

# Choose only up to 2 cliamte variables:
Clim=select.list(names(X[,3:ncol(X)]),multiple=TRUE,title="Select Climate Variables of Your interest:",graphics=TRUE)
List=c("year","month")
X=X[,names(X) %in% c(List,Clim)]

X1=recast(X,year~variable+month,id.var=c("year","month"),na.rm=TRUE)

# Convert to time series object
X1=ts(X1,frequency=1,start=min(X1$year),end=max(X1$year))
X1=cbind(p=X1,c=lag(X1))

# Convert back to data.frame
X1=as.data.frame(X1)
names(X1)

ind=which(colnames(X1)=="c.year")
X1=X1[,c(ind,2:(ind-1),(ind+1):ncol(X1))]
X1=na.omit(X1)
colnames(X1)[1]="year"


# Fit Random Forest to identify and reduce monthly variables prior to PCA
library(randomForest)

## Select Type of Ring Width
Y=ring
Ring.list=select.list(names(Y[,2:ncol(Y)]),multiple=TRUE,title="Select Rind Width Series of Your Interest:",graphics=TRUE)
List="year"
Y=Y[,names(Y) %in% c(List,Ring.list)]
#Y=data.frame(Y[,1:2],row.names="year")

## Merge ring width & Climate variables
head(X1);head(Y)
X2=merge(Y,X1,by="year")

X3=X2[,-1]

## Run RF
fit=randomForest(ModNegExp~.,data=X3,confusion=TRUE,ntree=5000,proximity=TRUE,importance=TRUE,na.action=na.omit)
print(fit)

# Variable importance
varImpPlot(fit,cex=0.7,main="Variable Importance")
box(which = "outer", lty = "solid")

# Choose top 10 important variable
rf=data.frame(importance(fit))
rf=data.frame(Var=rownames(rf),rf);rownames(rf)=NULL
rf1=rf[order(-rf$X.IncMSE),,drop=FALSE]
a=rf1[1:30,]
VarNames=as.character(a$Var)


# If you want to Select monthly Climate variables based on RandomForest:
Ring=colnames(X2[2])
Year="year"
X2=X2[,names(X2) %in% c(Year,Ring,VarNames)]

# fit PCA
library(FactoMineR)

RES=PCA(X2[,3:ncol(X2)],quanti.sup=2,scale.unit=TRUE,graph=FALSE,ncp=33)
RES$eig
RES$var$cor

Score=RES$ind$coord
X3=cbind(X2[,1:2],Score)
X3=X3[,-1]

library(RcmdrMisc)
M=lm(ModNegExp~., data=X3)
step=stepwise(M,criterion="AIC")
summary(step)










