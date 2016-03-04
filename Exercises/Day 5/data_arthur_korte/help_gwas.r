###how to use GWAS R script.

source('emma.r')

source('gwas.r')
source('plots_gwas.r')

load('MTMM_SAMPLE_DATA.Rdata')



a<-Sys.time()
amm_gwas(Y,X,K,include.lm=T,calculate.effect.size=FALSE)
Sys.time()-a


##  some plots 
plot_gwas(output)

local_plot(output)


op <- par(mfrow=c(2,1),mar=c(1, 4, 0.5, 4))
plot_gwas(output)
plot_gwas(output,h=9)
par(op)

qq_plot(output,farbe='blue',max.y=20)
par(new=T)
qq_plot(output,h=9,farbe='red',max.y=20)


