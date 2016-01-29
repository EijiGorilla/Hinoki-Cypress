# R codes below are used to sum daily climate variables and calculate monthly data:

a=file.choose()
Climate=read.csv(a,stringsAsFactors=FALSE,na.strings=".")

Climate$date=as.Date(Climate$date)
Climate$time=as.Date(Climate$date,"%Y-%m-%d")
Climate$week=1
Climate$week[Climate$day>=15 & Climate$day<=31]=2

# Sum daily climate variable:----
library(reshape2)

head(Climate)
var=select.list(names(Climate),multiple=TRUE,title="Select Climate Variables to be Summed:",graphics=TRUE)

X1=Climate[,c(var,"year","month")]
X1$min.temp[X1$min.temp>=0]=0

# Sum by year and month:
library(plyr)
X2=ddply(X1,.(year,month),summarize,value=sum(min.temp))


# IF you want to convert table to wide format:
X1=recast(X1,year~variable+month,sum,id.var=c("year","month"),na.rm=FALSE) # Sum daily min temp
colnames(X1)=c("year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
head(X1)
