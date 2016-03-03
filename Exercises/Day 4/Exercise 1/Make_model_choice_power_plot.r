library(MASS);
library(lattice);

#--------------------------------------------
#Read validation file and open pdf
#--------------------------------------------
a<-read.table("ABC_modelchoice_modelChoiceValidation.txt", header=T)

pdf(file="model_choice_power_fig.pdf", width=8, height=8);


#--------------------------------------------
#make plot
#--------------------------------------------
#functions
transformLog<-function(x, doTransformation=T){  
  if(!doTransformation){ return(x);}
  x[x>=0.99]<-1-10^(-2);
  x[x<=0.01]<-10^(-2);
  temp<-log10(2*(0.5-abs(0.5-x)));
  temp[x>0.5]<-(-temp[x>0.5]);
  return(temp);
}
getLogit<-function(x, glm){ 
  return(exp(glm$coefficients[1]+x*glm$coefficients[2])/(1+exp(glm$coefficients[1]+x*glm$coefficients[2])));
}

plotPPQQ<-function(abcPP, trueModel, maxres=2, nBins=50, windowsize=1, add=F, col='red', doTransformation=T, main="", xlab="ABC posterior probability", ylab="empirical posterior probabilities",lty=1,pch=19){
  par(cex=1.7)
  if(!add){ 
    if(doTransformation){
      plot(c(-maxres, maxres), c(-maxres, maxres),xlim=c(-maxres, maxres),ylim=c(-maxres, maxres),col='grey', lty=2, type='l', lwd=2, xaxt='n', yaxt='n', xlab=xlab, ylab=ylab, main=main); 
      l<-10^(seq(-maxres,-1));
      l<-sort(c(l, 0.5, 1-l));
      ll<-length(l)
      axis(side=1, at=transformLog(l, doTransformation), labels=c(expression(""<=0.01),l[-c(1,ll)],expression("">=0.99)));
      axis(side=2, at=transformLog(l, doTransformation), labels=c(expression(""<=0.01),l[-c(1,ll)],expression("">=0.99)));
    } else {
      plot(c(0,1), c(0,1), col=col, lty=lty,type='l', lwd=1, xaxt='n', yaxt='n', xlab=xlab, ylab=ylab, main=main); 
      l<-seq(0,1,length.out=11);
      axis(side=1, at=l, labels=l);
      axis(side=2, at=l, labels=l);      
    }
  }
  
  abcPPTrans<-transformLog(abcPP, doTransformation);  
  start<-transformLog(10^(-maxres-1), doTransformation);
  binLength<-(transformLog(1-10^(-maxres-1), doTransformation)-start)/nBins;  
  bins<-c();
  mids<-c();  
  cc<-c();
  for(i in 1:nBins){    
    x<-abcPPTrans<(start+(i+windowsize)*binLength) & abcPPTrans>=(start+(i-windowsize-1)*binLength);    
    if(sum(x)>0){
      temp<-sum(trueModel[x]==T) / sum(x);
      if(temp>0 & temp <1){
        bins<-c(bins, temp);  
        mids<-c(mids, i);
        cc<-c(cc, sum(x));
      }
    }    
  } 
  mids<-start + (mids-0.5)*binLength;    
  lines(mids, transformLog(bins, doTransformation), col=col,  type='o', pch=pch,lty=lty,lwd=1.5)  
}


par(pty="s")
add <- F;
for(modelToCheckFor in c(0,1)){
  trueModel<-rep(F, length(a[,1]));
  trueModel[a[,1]==modelToCheckFor]<-T;
  abcPP<-a[,7+modelToCheckFor]; 
  plotPPQQ(abcPP, trueModel, doTransformation=F, col="black", pch=19+2*modelToCheckFor, lty=modelToCheckFor+1,add=add, main="", xlab=expression(italic(p[ABC])), ylab=expression(italic(p[empirical])), nBins=15);
add <- T;
}

legend('topleft', bty='n', legend=c("Model 0", "Model 1"), lty=c(1,2), pch=c(19,21));

dev.off();

