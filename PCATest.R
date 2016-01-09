## Calculate extreme values from daily climate variables (5% and 95% percentiles) for each climate variable
## Run PCA and regress ring (or climate variables) against climate variables (or ring)

library(FactoMineR)
library(vegan)
library(dplR)
library(tseries)

# choose and import ring width data----
b=file.choose() # 100-19M Hinoki ring width 1734-2006.csv
ring=read.csv(b,stringsAsFactors=FALSE,na.strings=".")

## Detrend a raw ring width by Ar and get residual chronology
### note that the detrended series is the residuals of an Ar model divided by the mean of those residuals to yield a series with white noise and a mean of one
library(dplR)

head(ring,n=10)
year=ring[,"year"]
width=ring[,"width"]

series=data.frame(width,row.names=year)
names(width)=rownames(series)

series.rwi=detrend.series(y=width,y.name="100-19M",verbose=TRUE)
gg=data.frame(year=year,width=series.rwi,row.names=NULL)
head(gg)

ring=merge(ring,gg,by="year") # Use this ring data for all analyses
colnames(ring)[3:ncol(ring)]=c("spline","negexp","mean","ar")

# Choose and import climate data----
a=file.choose() # Monthly climate Matsumoto 1898-2006.csv
Climate=read.csv(a,stringsAsFactors=FALSE,na.strings=".")

ClimVar=select.list(names(Climate[,4:ncol(Climate)]),multiple=TRUE,title="Choose Climate Variables of Your Interest:",graphics=TRUE)
List=c("year","month")
Climate=Climate[,names(Climate) %in% c(List,ClimVar)]

library(reshape2)
Climate=recast(Climate,year~variable+month,id.var=c("year","month"))

# Merge Ring width + Climate variables
X=merge(ring,Climate,by="year")

X1=ts(X,frequency=1,start=min(X$year),end=max(X$year))
X1=cbind(p=X1,c=lag(X1))

### Convert back to data.frame
X1=as.data.frame(X1)

### Onl retain months of our interests: Nov, Dec of last years and Jan, Feb, Mar,and Apr of this year's growth
ind=which(colnames(X1)=="p.year")
X1=X1[,-c(ind:(ind+5))]

ind=which(colnames(X1)=="c.year")
X1=X1[,c(ind:(ind+5),1:(ind-1),(ind+6):ncol(X1))]
X1=X1[-1,]

head(X1)

RES=PCA(X1[,c(7:ncol(X1))],scale.unit=TRUE,quanti.sup=4,ncp=6)
RES$eig

RES$var
