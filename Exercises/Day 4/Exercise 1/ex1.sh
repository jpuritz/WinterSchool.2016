#!/bin/bash
# 
# -----------------------------------------------------------
# Exercise 1: Approximate Bayesian Computation with ABCtoolbox
# -----------------------------------------------------------
# Steps:
#   1) Compile the latest version of ABCtoolbox and check executables
#   2) Generate simulations with known parameters as truth set
#   3) Calculate summary statistics from the observed data
#   4) Generate simulations with ABCtoolbox
#   5) Estimate posterior distributions with ABCtoolbox
#   6) Plot posterior distributions in R
#   7) Generate simulations with ABCtoolbox for theta
#   8) Estimate posterior distributions with ABCtoolbox
#   9) Finding appropriate summary statistics using PLS
#   10) Transforming statistics via PLS definitions
#   11) Validation of parameter inference: are the stats reproduced?
#   12) Validation of parameter inference: power and bias
#   13) Setting up a competing model with a bottleneck
#   14) Finding statistics for model choice
#   15) Performing model choice
#   16) Model choice validation
# 
# To see the individual steps, type './ex1.sh STEP', where STEP is the number of the desired step.
# 

#1 
#1 Step 1: Compile the latest version of ABCtoolbox and check executables
#1 ----------------------------------------------------------------------
#1 - First, let us compile the latest version of ABCtoolbox. I have provided the latest code in /bin/ABCtoolbox_src/. So change into that directory. Be aware that the specific location depends on where you copied the provided folder to!
#1 - Now lets compile ABCtoolbox! This is done easily using the provided make file by simply typing
#1 
#1 > make
#1 
#1 - Next, move the resulting binary into the bin folder:
#1 
#1 > mv ABCtoolbox ../
#1 
#1 - Next make sure all other provided binaries are executable
#1 
#1 > chmod +x *
#1 
#1 - Finally, and to have easy access to these executables, add the bin folder to your PATH. Note that you need to adjust the command below to take care of the location o which you copied the files of thus exercise (replace the placeholder LOCATION in the command below).
#1 
#1 > PATH="LOCATION/bin:${PATH}"
#1 
#1 - You can easily test if that worked. Simply try to launch one of the programs from another directory_:
#1 
#1 > cd ..
#1 > ABCtoolbox
#1

#2 
#2 Step 2: Generate simulations with known parameters as truth set
#2 ---------------------------------------------------------------
#2 - To test how well the estimation is doing, we will generate some simulations under a simple model (constant population size).
#2 - To generate simulations, we will use fastsimcoal (version 2.5.2.2.1 = executable fsc25221).This program requires an input file specifying the model, and we will use here the file dna_singlepop_constsize_obs.par. Have a look at that file and locate the two important parameters i) population size and ii) mutation rate.
#2 - The manual for fastsimcoal is available here: http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal25.pdf
#2 - Use the following command to generate the simulations (run fsc25221 without arguments to get some help):
#2 
#2 > fsc25221 -i constsize_obs.par -n 1
#2 
#2 - Note: -n 1 implies that a single simulation is conducted.
#2 - fastsimcoal will generate a directory with the same name as the input will. In that directory you will find a file with the ending .arp, which contains the simulated data.
#2 

#3 
#3 Step 3: Calculate summary statistics from the observed data
#3 -----------------------------------------------------------
#3 - Here we will use the command line version of Arlequin (called arlsumstat3522_64bit) to calculate summary statistics from the output of fastsimcoal.
#3 - arlsumstat requires two setting files to run correctly: i) the file arl_run.ars that tells arlsumstats what statistics are to be calculated, ii) ssdefs.txt, which specified what summary statistics to be printed. The files provided here are set to calculate essentially everything possible. Have a look at those file, but note that both files can be generated and modified with the graphical version of Arlequin.
#3 - To run arlsumstat on the generated data, use the following command line (run arlsumstat without arguments to get some help):
#3 
#3 > arlsumstat3522_64bit constsize_obs/constsize_obs_1_1.arp constsize.obs 0 1
#3 
#3 - The results are written to the file constsize.obs. Have a look at it!
#3 


#4 
#4 Step 4: Generate simulations with ABCtoolbox
#4 --------------------------------------------
#4 - ABCtoolbox allows you to use external programs to conduct simulations with parameters taken from prior distributions and to calculate summary statistics from the resulting data.
#4 - ABCtoolbox is configured by means of an input file and a file defining the model parameters and their prior distributions. Here we will use the files constsize.input and constsize.est, respectively. Have a look at these files and try to figure out how the use of external programs is specified. Note that fastimcoal is to be run with an input file on its own. Here we will use the file constsize.par, which contains tags for the parameters that are to be sampled from the prior distributions.
#4 - Also, have a look at how the prior are specified. Why are we not using the logunif prior distribution provided by ABCtoolbox directly on N_NOW?
#4 - To run ABCtoolbox with these input files, simply run
#4 
#4 > ABCtoolbox constsize.input
#4 
#4 - The resulting simulations are found in the file sims_constsize_sampling1.txt. Also have a look at the file childOutput.txt, which contains the output of all the fastsimcoal and arlsumstat calls.
#4 - Since generating these simulations may take some time, I already prepared a file containing about 50,000 of them. You can now add those to the 100 simulations you generated with the following command:
#4 
#4 > tail -n+2 sims_constsize_50K.txt >> sims_constsize_sampling1.txt
#4 


#5 
#5 Step 5: Estimate posterior distributions with ABCtoolbox
#5 --------------------------------------------------------
#5 - Again, the settings telling ABCtoolbox to perform parameter estimation are provided in an input file. Here we will use the file estimate.input - have a look at it and try to understand what ABCtoolbox will do when using it.
#5 - Launch ABCtoolbox with this input file as follows:
#5 
#5 > ABCtoolbox estimate.input
#5 
#5 - ABCtoolbox will now complain that there are highly correlated statistics and not perform an estimation. In order to automatically prune one of any pair of statistics that are highly correlated, add the option "pruneCorrelatedStats" to the input file and rerun the estimation.
#5 - The estimation should now have finished successfully. If that was the case, the output was written to a series of files with tag "ABC_estimation_constsize_", as was specified in the input file.
#5 

#6 
#6 Step 6: Plot posterior distributions in R
#6 -----------------------------------------
#6 - ABCtoolbox ships with a bunch of R scripts to plot posterior distributions and other metrics. However, I recommend to plot them yourself so that you get a feel for the structure of the output. For this, simply open an R terminal in your working directory and follow the steps below.
#6 - Should your computer not allow X forwarding from the server, you can also output all plots to a pdf file. To open pdf output in R, simply type 'pdf("filename", width=5, height=5)', where filename is the name of the outputfile (e.g. posteriors.pdf). After plotting, you need to close the pdf output using 'dev.off()'.
#6 - Begin by plotting the marginal posterior distributions of the population size (N_NOW) and the mutation rate (MUTATION_RATE) as follows.
#6 
#6 R> post <- read.table("ABC_estimation_constsize_model0_MarginalPosteriorDensities_Obs0.txt", header=T)
#6 R> par(mfrow=c(1,2))
#6 R> plot(post$LOG10_N_NOW, post$LOG10_N_NOW.density, type='l')
#6 R> plot(post$LOG10_MUTATION_RATE, post$LOG10_MUTATION_RATE.density, type='l')
#6 
#6 - Now let's compare these marginal estimates to the 2D joint posterior we also estimated
#6 
#6 R> twoD <- read.table("ABC_estimation_constsize_model0_jointPosterior_1_2_Obs0.txt", header=T)
#6 R> x <- unique(twoD$LOG10_N_NOW)
#6 R> y <- unique(twoD$LOG10_MUTATION_RATE)
#6 R> z <- matrix(twoD$density, nrow=length(x), byrow=T)
#6 R> contour(x,y,z, xlab="Log10(N_Now)", ylab="Log10(MutRate)")
#6 
#6 - Do you notice something? Yes - we can not estimate N and mu together! So let's try to estimate theta...
#6 

#7
#7 Step 7: Generate simulations with ABCtoolbox for theta
#7 ------------------------------------------------------
#7 - The model for theta is defined in the est file constsize_theta.est. Note that the same par file can be used!
#7 - Can you generate an input file constsize_theta.input for ABCtoolbox to generate simulations under this model? Hint: copy the file from the N and Mu model and modify it. Make sure you change the output tag to sims_constsize_theta_ 
#7 - Then, you can generate simulations as follows:
#7 
#7 > ABCtoolbox constsize_theta.input
#7 
#7 - Since generating these simulations may take some time, I already prepared a file containing about 50,000 of them. You can now add those to the 100 simulations you generated with the following command:
#7 
#7 > tail -n+2 sims_constsize_theta_50K.txt >> sims_constsize_theta_sampling1.txt
#7 


#8 
#8 Step 8: Estimate posterior distributions with ABCtoolbox
#8 --------------------------------------------------------
#8 - Now estimate theta using ABCtoolbox. For this, create an input file estimate_theta.input by copying and then modifying the input file we used for the previous model.
#8 - Make sure to use the option "pruneCorrelatedStats".
#8 - Also, change the output prefix to "ABC_estimation_constsize_theta_" to avoid overwriting your results.
#8 - Launch ABCtoolbox with this input file as follows:
#8 
#8 > ABCtoolbox estimate_theta.input
#8 
#8 - You can now plot the marginal posterior for theta in R. How well does the estimate fit the parameters we used to generate the observed data? Check the file dna_singlepop_constsize_obs.par to see what parameters we used and note that theta=2*N*mu (for haploid species).
#8 
#8 - Change the parameters in constsize_obs.par and regenerate observed data, calculate statistics from them and rerun the estimation (no need to generate new simulations!). Can you accurately estimate theta?
#8 

#9 
#9 Step 9: Finding appropriate summary statistics using PLS
#9 --------------------------------------------------------
#9 - Some of the summary statistics are highly correlated, suggesting that we may reduce the summary statistics space while retain all the information.
#9 - One method to do so is via a Partial Least Squares (PLS) transformation.
#9 - For such a transformation you may use the R script find_pls.r, which is provided in the exercise directory.
#9 - However, to do so, you first need to make sure the package "pls" is installed. If it is not installed, you can install it within R as follows:
#9 
#9 R> install.packages("pls")
#9 
#9 - Once the package is installed, simply call the R script to find PLS components as follows:
#9 
#9 > Rscript find_pls.r sims_constsize_theta_sampling1.txt
#9 
#9 - This script will create two files: i) one called PLSdef_sims_dna_singlepop_constsize_theta_sampling1.txt containing the PLS definitions, and ii) one called RMSE_sims_dna_singlepop_constsize_theta_sampling1.txt.pdf, which contains the RMSE plot to determine the lowest number of PLS components that are required to retain most of the information. Have a look at it. You will probably agree that a single PLS component is capturing all the info (in fact, it is known that the number of segregating sites S is a sufficient statistics for this model).
#9 

#10 
#10 Step 10: Transforming statistics via PLS definitions
#10 ----------------------------------------------------
#10 - ABCtoolbox offers the possibility to transform the statistics into PLS components.
#10 - While you may again write an input file for that, it is time to learn how to use ABCtoolbox via the command line. To transform the statistics using the PLS definitions you generate before, use ABCtoolbox as follows:
#10 
#10 > ABCtoolbox task=transform linearComb=PLSdef_sims_constsize_theta_sampling1.txt numLinearComb=1 input=sims_constsize_theta_sampling1.txt output=sims_constsize_theta_sampling1.txt.pls boxcox verbose
#10 
#10 - This will generate a new simulation file containing all parameters and, instead of the summary statistics, PLS components. To now rerun the estimation, you also need to transform the observed data the same way. I'm sure you'll manage to do that!
#10 - Once the observed data has also been transformed, rerun the estimation. Note that you need to modify the input file to use the transformed files now!
#10 - How does the posterior look like? Any visible changes?
#10 

#11 
#11 Step 11: Validation of parameter inference: are the stats reproduced?
#11 ---------------------------------------------------------------------
#11 - Validation is an important part of any ABC application. ABCtoolbox offers plenty of validation tools - let's look at some of them!
#11 - First, let's check if the model is capable of reproducing the observed summary statistics. For this purpose, ABCtoolbox offers two statistics to compare the observed statistics with those simulated: the marginal density and the Tukey depth. To perform tests using these two statistics, you need to add the following two arguments to the input file (note: the 1000 refer to the number of retained simulations to be used when calculating p-values):
#11 
#11 marDensPValue 1000
#11 tukeyPValue 1000
#11 
#11 - Then, simply rerun the estimation and check either the output written to screen or the file with tag "modelFit.txt". Are these tests passed? Did using PLS make a difference?
#11 

#12 
#12 Step 12: Validation of parameter inference: power and bias
#12 ----------------------------------------------------------
#12 - Next, let's ask ABCtoolbox to perform some estimation on pseudo-observed data sets. This is easily done: simply add the following argument to the input file:
#12 
#12 randomValidation 1000
#12 
#12 - When rerunning the estimation with this argument, ABCtoolbox will select 1000 random simulations as pseudo-observed data sets and then perform estimation on these. The results are written to the file with tag "RandomValidation.txt". Does the content of this file make sense to you?
#12 - You can now use the content of this file to contrast the true value of the simulation with the estimated one and compute correlations between them:
#12
#12 R> d <- read.table("ABC_estimation_constsize_theta_model0_RandomValidation.txt", header=T)
#12 R> plot(d$LOG10_THETA, d$LOG10_THETA_mode)
#12 R> cor(d$LOG10_THETA, d$LOG10_THETA_mode)
#12 
#12 - Additionally, you can check for biases in the posterior distributions by using either the quantile to HDI. As discussed in class, these measures should be distributed uniformly.
#12 
#12 R> hist(d$LOG10_THETA_quantile)
#12 R> hist(d$LOG10_THETA_HDI)
#12 
#12 - Are they? What may affect these distributions and in what way?

#13 
#13 Step 13: Setting up a competing model with a bottleneck
#13 -------------------------------------------------------
#13 - Let's next generate some simulations for a model with a bottleneck. The goal is to have three parameters: the current size, the magnitude of the population size change and when that change happened.
#13 - As for the model of constant size, it does make sense to define the population size in units of Theta and to use priors on the log10 scale.
#13 - Specifically, let us have a uniform priors on i) the current theta log10(THETA_CUR) ~ U[-4, -1], ii) the old population size before the bottleneck relative to the current size log10(OLD_SIZE_RELATIVE) ~ U[-2,2] and iii) the time when the bottleneck occurred log10(T_BOT) ~ U[1,3].
#13 - Can you generate an est and par file (named bottleneck.est and bottleneck.par, respectively) for this model? Hint: copy the files of the theta example and modify them.
#13 - Then, prepare an input file bottleneck.input to generate 100 simulations under this model and run it using ABCtoolbox. Choose "sims_bottleneck" as the output file name.
#13 
#13 > ABCtoolbox bottleneck.input
#13 
#13 - Since generating these simulations may take some time, I already prepared a file containing about 50,000 of them. You can now add those to the 100 simulations you generated with the following command:
#13 
#13 > tail -n+2 sims_bottleneck_50K.txt >> sims_bottleneck_sampling1.txt
#13 

#14 
#14 Step 14: Finding statistics for model choice
#14 --------------------------------------------
#14 - Let us now check if ABC is capable of distinguishing between these two model.
#14 - To do so, we first need to find summary statistics that are informative about these models. Unfortunately, PLS does not work for this. However, ABCtoolbox implements a greedy search to find such summary statistics.
#14 - To invoke this greedy search, you can use the prepared input file "findStatsModelChoice.input" by launching ABCtoolbox with it:
#14 
#14 > ABCtoolbox findStatsModelChoice.input
#14 
#14 - This run will take quite some time. But once it finished, you should find a file called ABC_findStats_greedySearchForBestStatisticsForModelChoice.txt. This file lists the power of all tested statistic combinations. Among the combinations with highest power, choose a combinations with few statistics and low correlation among them. Which is your combination?
#14 - Prepare an obs file containing only these statistics for the observed data.

#15 
#15 Step 15: Performing model choice
#15 --------------------------------
#15 - Running model choice with ABCtoolbox is straight forward: simply add additional models to the arguments "simName" and "params".
#15 - In our case, first copy the file "estimate_theta.input"
#15 
#15 > cp estimate_theta.input modelchoice.input
#15 
#15 - Then, add the bottleneck model to simName and to params and change the output prefix:
#15 
#15 simName sims_constsize_sampling1.txt;sims_bottleneck_sampling1.txt.old
#15 params 2;2-4 
#15 outputPrefix ABC_modelchoice_
#15 
#15 - Finally, rerun ABCtoolbox on this input file and check both the output written to screen as well as the model fit file. Which model is preferred?
#15 

#16 
#16 Step 16: Model choice validation
#16 --------------------------------
#16 - Model choice validation is invoked with the argument "modelChoiceValidation", followed by the number of pseudo-observed data sets to be used.
#16 - Add the following to your input file:
#16 
#16 modelChoiceValidation 1000
#16 
#16 - And then rerun ABCtoolbox. You should now find the two files "ABC_modelchoice_confusionMatrix.txt" and "ABC_modelchoice_modelChoiceValidation.txt" containing the model choice results.
#16 - First have a look at the confusion matrix: is there power to distinguish between these models?
#16 - The file "ABC_modelchoice_modelChoiceValidation.txt" can be used to infer biases in the Bayes factor calculation. This is a rather complex analysis for which an R script is provided. Run it as follows:
#16 
#16 > Rscript Make_model_choice_power_plot.r
#16 
#16 - This scripts produces a plot called "model_choice_power_fig.pdf" comparing the posterior probability as estimated via ABC to those obtained empirically from the validation results. Can you trust the ABC posterior probabilities?

LOC=$(which ./ex1.sh)

PATTERN=#$1[[:blank:]]
grep "$PATTERN" $LOC | awk '{$1=""; print $0}'
