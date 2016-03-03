#open File
numComp<-5;

directory<-"";

filename<-commandArgs()[length(commandArgs())];
print(paste("Reading file '", filename, "'", sep=""));

#read file
a<-read.table(paste(directory, filename, sep=""), header=T, nrows=2000, skip=0);

firstStat<-grep("H_1", names(a))[1];
stats<-a[,firstStat:length(a[1,])];

p<-grep("LOG10", names(a));
if(length(p)==1){
  params <- data.frame(x=a[,p]); names(params)<-names(a)[p];
} else {
  params<-a[,p]; 
}
print(names(params));

 
#standardize the params
for(i in 1:length(params)){params[,i]<-(params[,i]-mean(params[,i]))/sd(params[,i]);}

#force stats in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stats)){
	myMax<-c(myMax, max(stats[,i]));
	myMin<-c(myMin, min(stats[,i]));
	stats[,i]<-1+(stats[,i]-myMin[i])/(myMax[i]-myMin[i]);
}

#transform statistics via boxcox  
library("MASS");	
for(i in 1:length(stats)){		  
	d<-cbind(stats[,i], params);
	mylm<-lm(as.formula(d), data=d)			
	myboxcox<-boxcox(mylm, lambda=seq(-50, 80, 1/10), plotit=T, interp=T, eps=1/50);	
	lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);			
	print(paste(names(stats)[i], myboxcox$x[myboxcox$y==max(myboxcox$y)]));
	myGM<-c(myGM, exp(mean(log(stats[,i]))));			
}

#standardize the BC-stats
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:length(stats)){
	stats[,i]<-(stats[,i]^lambda[i] - 1)/(lambda[i]*myGM[i]^(lambda[i]-1));	
	myBCSDs<-c(myBCSDs, sd(stats[,i]));
	myBCMeans<-c(myBCMeans, mean(stats[,i]));		
	stats[,i]<-(stats[,i]-myBCMeans[i])/myBCSDs[i];
}

#perform pls
library("pls");
#myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp, validation="LOO");
myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp);

#write pls to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:numComp) { myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]); } 
write.table(cbind(names(stats), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs, myPlsrDataFrame), file=paste(directory, "PLSdef_", filename, sep=""), col.names=F, row.names=F, sep="\t", quote=F);

#make RMSE plot
pdf(paste(directory, "RMSE_", filename, ".pdf", sep=""));
plot(RMSEP(myPlsr), ylim=c(0,1));
dev.off();



#obsa<-read.table("/mnt/uni/ABC/arvalis/arvalis_both.obs", header=T);
#n<-data.frame(a=1:length(names(obsa)), n=names(obsa));
#pdf(paste(directory, "stats_", filename, ".pdf", sep=""), width=9, height=12);
#par(mfrow=c(5,4), cex=0.5)		
#	for(i in c(1:13,25,26,49:51,63,64,76:80,183:227)){
#	plot(density(stats[,i]), xlim=c(min(stats[,i])-max(stats[,i])+min(stats[,i]),max(stats[,i])+max(stats[,i])-min(stats[,i])), main=names(stats)[i]);	
#	print(paste(n[n[,2]==names(stats)[i],1], obsa[n[n[,2]==names(stats)[i],1]]));
#	lines(c(obsa[,n[n[,2]==names(stats)[i],1]], obsa[,n[n[,2]==names(stats)[i],1]]), c(0,1000), col="red")
#}

#dev.off();

