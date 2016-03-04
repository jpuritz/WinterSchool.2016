### help for GWAS

## source scripts we need 
source('emma.r')
source('gwas.r')
source('plots_gwas.r')

#### GWAS on dots

load("Y_dots.rda")	# phenotype data
load("X_dots.rda")	# genotype data	
load("K_dots.rda")	#kinship matrix

output <- amm_gwas(Y=Y_dots,X=X_dots,K=K_dots, n=3) 
## the script requires Y,X,K n=3 tells the script to perform the analysis on the 3. column (where our phenotype is)
output2 <- amm_gwas(Y_dots, X_dots, K_dots, n=3, include.lm=T)

par(mfrow=c(1,2))
plot_gwas(output, lower.limit=0.1)

#the linear model does not include a termed for kinship. If the linear model is different from the mixed model, it suggest that there is a structure in the data.
qq_plot(output2)



#triangle varying in size and colors
load("Y_tri.rda")
#load("X_tri.rda")   contains missing data 
load("K_tri.rda")
#missing data we have "impute" (infering missing data from linked markers)
load("X_imp.rda")

#position of the snps on chromosome 
load("SNP_tri.rda")


output_tri <- amm_gwas(Y_tri, X_imp, K_tri, n=2, include.lm=T,  gen.data="heterozygot", use.SNP_INFO=T, SNP_INFO=SNP_tri)
#### same as with dots, but two differnences : need to set gen.data='heterozygot' and need to provide SNP_INFO (otherwise NAs will be produced)
output_tri_col <- amm_gwas(Y_tri, X_imp, K_tri, n=7, include.lm=T,  gen.data="heterozygot", use.SNP_INFO=T, SNP_INFO=SNP_tri)

par(mfrow=c(2,2))
plot_gwas(output_tri, lower.limit=0.01,h=9)
plot_gwas(output_tri, lower.limit=0.01,h=8)

qq_plot(output_tri,h=9)
qq_plot(output_tri,farbe='blue')


plot_gwas(output_tri_col, lower.limit=0.1)
