#!/usr/bin/env Rscript
####HiSSE - Samp Freq = 0.75_0.36####

library(hisse)

library(diversitree)

library(geiger)

library(optparse)

opts <- list(make_option('--tree'),
             make_option('--states'),
             make_option('--char-states'),
             make_option('--output'))
args <- parse_args(OptionParser(option_list=opts),
                   positional_arguments=TRUE)

tree_file       <- args$options[['tree']]
state_file      <- args$options[['states']]
hisse_dat_file  <- args$options[['char-states']]
out_prefix      <- args$options[['output']]

#Read in tree - Brendan substituted tree file variable here
# instead of hard-coded file name
tree <- read.nexus(tree_file)

#Read in states - again, state_file instead of hard-coded
# states
#states <- scan(state_file)

# Make dataframe with spp names and char states
# Variable instead of hard-coded file name
dat <- read.csv(hisse_dat_file)


#Transition rates used here are those estimated for the 'hybridizability' trait along out phylogeny. 
#We then simulate a neutral trait along the phylogeny. 
q <- list(rbind(c(-0.002692059,0.002692059), c(0.037110456,-0.037110456)))
sim.traits <- sim.char(tree, q, model='discrete', nsim = 1, root = 1)
sim.traits <- as.data.frame(sim.traits)
dat$State <- sim.traits[,1]-1

####################################
# Null 2 All Transitions All Equal #
####################################
trans.rate <- TransMatMaker(hidden.states=TRUE)
trans.rate[!is.na(trans.rate) & !trans.rate == 0] = 1

# Fixed turnover rates between 0A & 1A, and between 0B & 1B
# This runs the analysis
Null2_AllTrans_Eq <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,1,2,2), 
                        eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw", bounded.search = TRUE)

# Store likelihood, AIC, and AICc
Null2_AllTrans_Eq.logL <- Null2_AllTrans_Eq$loglik
Null2_AllTrans_Eq.AIC <- Null2_AllTrans_Eq$AIC
Null2_AllTrans_Eq.AICc <- Null2_AllTrans_Eq$AICc


# Write output to file
capture.output(Null2_AllTrans_Eq, file=paste0(out_prefix, '_Null2_AllTrans_EqRate.txt'))


##########################################
# Null 2 No Double Transitions All Equal #
##########################################
trans.rate <- TransMatMaker(hidden.states=TRUE)
trans.rate <- ParDrop(trans.rate, c(3,5,8,10))
trans.rate[!is.na(trans.rate) & !trans.rate == 0] = 1

# Fixed turnover rates between 0A & 1A, and between 0B & 1B
# This runs the analysis
Null2_NoDub_Eq <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,1,2,2), 
                        eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw", bounded.search = TRUE)

# Store likelihood, AIC, and AICc
Null2_NoDub_Eq.logL <- Null2_NoDub_Eq$loglik
Null2_NoDub_Eq.AIC <- Null2_NoDub_Eq$AIC
Null2_NoDub_Eq.AICc <- Null2_NoDub_Eq$AICc


# Write output to file
capture.output(Null2_NoDub_Eq, file=paste0(out_prefix, '_Null2_NoDub_EqRate.txt'))


#######################################
# Null 2 Three Trans Rates, No Double #
#######################################

# Define transition matrix
trans.rate <- TransMatMaker(hidden.states=TRUE)
trans.rate <- ParDrop(trans.rate, c(3,5,8,10))
to.change <- cbind(c(1,3), c(2,4))
trans.rate[to.change] = 1
# Now set all transitions from 1->0 to be governed by a single rate:
to.change <- cbind(c(2,4), c(1,3))
trans.rate[to.change] = 2
# Finally, set all transitions between the hidden state to be a single rate (essentially giving 
# you an estimate of the rate by which shifts in diversification occur:
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rate[to.change] = 3
trans.rate

# Fixed turnover rates between 0A & 1A, and between 0B & 1B
# This runs the analysis
Null2_ThreeRate_NoDub <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,1,2,2), 
                        eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw", bounded.search = TRUE)

# Store likelihood, AIC, and AICc
Null2_ThreeRate_NoDub.logL <- Null2_ThreeRate_NoDub$loglik
Null2_ThreeRate_NoDub.AIC <- Null2_ThreeRate_NoDub$AIC
Null2_ThreeRate_NoDub.AICc <- Null2_ThreeRate_NoDub$AICc


# Write output to file
capture.output(Null2_ThreeRate_NoDub, file=paste0(out_prefix, '_Null2_NoDub_ThreeRate.txt'))



###############
# BiSSE Model #
###############

# Now make the bisse model where diversification changes with hybridizability without
# the presence of hidden states. This will be the BiSSE Null model.

trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

# The transition matrix thus looks like the following
#     (0) (1)
# (0)  NA   2
# (1)   1  NA

# Given that the order of arguments in the turnover and extinction .anc commands that 
# control the number of rate classes follows this order... (0A, 1A, 0B, 1B), the
# following will set the model to only transition between states that do not include
# hidden states.
bisse <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = F, turnover.anc = c(1,2,0,0), 
               eps.anc = c(1,2,0,0), trans.rate = trans.rates.bisse, output.type = "raw", bounded.search = TRUE)

# Store likelihood, AIC, and AICc
bisse.logL <- bisse$loglik
bisse.AIC <- bisse$AIC
bisse.AICc <- bisse$AICc

# Write to csv
capture.output(bisse, file=paste0(out_prefix, '_BiSSE.txt'))


####################
# BiSSE Null Model #
####################

# Make a constrained bisse model where diversification rates are trait independent
bisse.null <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = F, turnover.anc = c(1,1,0,0), 
                    eps.anc = c(1,1,0,0), trans.rate = trans.rates.bisse, output.type = "raw", bounded.search = TRUE)

# Store likelihood, AIC, and AICc
bisse.null.logL <- bisse.null$loglik
bisse.null.AIC <- bisse.null$AIC
bisse.null.AICc <- bisse.null$AICc

# Write to csv
capture.output(bisse.null, file=paste0(out_prefix, '_BiSSE_Null.txt'))


###################################
# HiSSE Null 4 Model, Equal Rates #
###################################

# Conduct the HiSSE null-4 model that contains the same complexity as a full HiSSE 
# model
hisse.null4.equal <- hisse.null4(phy = tree, data = dat, f = c(0.75, 0.36),  turnover.anc = rep(c(1,2,3,4),2), 
                               eps.anc = rep(c(1,2,3,4),2), trans.type = "equal", output.type = "raw", bounded.search = TRUE)

#Store likelihood, AIC, and AICc
hisse.null4.equal.logL <- hisse.null4.equal$loglik
hisse.null4.equal.AIC <- hisse.null4.equal$AIC
hisse.null4.equal.AICc <- hisse.null4.equal$AICc

# Write to csv
capture.output(hisse.null4.equal, file=paste0(out_prefix, '_Null4_Equal.txt'))



###################################
# HiSSE Null 4 Model, Three Rates #
###################################

# Conduct the HiSSE null-4 model that contains the same complexity as a full HiSSE 
# model
hisse.null4.three <- hisse.null4(phy = tree, data = dat, f = c(0.75, 0.36),  turnover.anc = rep(c(1,2,3,4),2), 
                                 eps.anc = rep(c(1,2,3,4),2), trans.type = "three.rate", output.type = "raw", bounded.search = TRUE)

#Store likelihood, AIC, and AICc
hisse.null4.three.logL <- hisse.null4.three$loglik
hisse.null4.three.AIC <- hisse.null4.three$AIC
hisse.null4.three.AICc <- hisse.null4.three$AICc

# Write to csv
capture.output(hisse.null4.three, file=paste0(out_prefix, '_Null4_Three.txt'))



##############################
# HiSSE No0B All Transitions #
##############################

# Make a model that has hidden states for only state 1 - we thus cannot have transitions to 0B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(2,5,12,7,8,9))

no0B_AllTrans <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,0,3), 
                        eps.anc = c(1,2,0,3), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE)

#Store likelihood, AIC, and AICc
no0B_AllTrans.logL <- no0B_AllTrans$loglik
no0B_AllTrans.AIC <- no0B_AllTrans$AIC
no0B_AllTrans.AICc <- no0B_AllTrans$AICc

# Write to csv
capture.output(no0B_AllTrans, file=paste0(out_prefix, '_No0B_AllTrans.txt'))


####################################
# HiSSE No0B No Double Transitions #
####################################

# Make a model that has hidden states for only state 1 - we thus cannot have transitions to 0B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(2,3,5,7,8,9,10,12))

no0B_NoDub <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,0,3), 
                       eps.anc = c(1,2,0,3), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE)

#Store likelihood, AIC, and AICc
no0B_NoDub.logL <- no0B_NoDub$loglik
no0B_NoDub.AIC <- no0B_NoDub$AIC
no0B_NoDub.AICc <- no0B_NoDub$AICc

# Write to csv
capture.output(no0B_NoDub, file=paste0(out_prefix, '_No0B_NoDub.txt'))


##############################
# HiSSE No1B All Transitions #
##############################

# Make a model that has hidden states for only state 1 - we thus cannot have transitions to 0B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(3,6,9,10,11,12))

no1B_AllTrans <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,0), 
                       eps.anc = c(1,2,3,0), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE)

#Store likelihood, AIC, and AICc
no1B_AllTrans.logL <- no1B_AllTrans$loglik
no1B_AllTrans.AIC <- no1B_AllTrans$AIC
no1B_AllTrans.AICc <- no1B_AllTrans$AICc

# Write to csv
capture.output(no1B_AllTrans, file=paste0(out_prefix, '_No1B_AllTrans.txt'))


####################################
# HiSSE No0B No Double Transitions #
####################################

# Make a model that has hidden states for only state 1 - we thus cannot have transitions to 0B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(3,5,6,8,9,10,11,12))

no1B_NoDub <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,0), 
                    eps.anc = c(1,2,3,0), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE)

#Store likelihood, AIC, and AICc
no1B_NoDub.logL <- no1B_NoDub$loglik
no1B_NoDub.AIC <- no1B_NoDub$AIC
no1B_NoDub.AICc <- no1B_NoDub$AICc

# Write to csv
capture.output(no1B_NoDub, file=paste0(out_prefix, '_No1B_NoDub.txt'))


##############################
# Full HiSSE All Transitions #
##############################

# Make a model that has hidden states for all states but does not allow for a transition between 
# state 0A and 1B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)

hisse.full_AllTrans <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,4), 
                        eps.anc = c(1,2,3,4), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE)

#Store likelihood, AIC, and AICc
hisse.full_AllTrans.logL <- hisse.full_AllTrans$loglik
hisse.full_AllTrans.AIC <- hisse.full_AllTrans$AIC
hisse.full_AllTrans.AICc <- hisse.full_AllTrans$AICc

# Write to csv
# Added by Brendan
# Tried to make this work as a variable too
capture.output(hisse.full_AllTrans, file=paste0(out_prefix, '_hisse_full_AllTrans.txt'))


####################################
# Full HiSSE No Double Transitions #
####################################

# Make a model that has hidden states for all states but does not allow for a transition between 
# state 0A and 1B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))

hisse.full_NoDub <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,4), 
                             eps.anc = c(1,2,3,4), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE)

#Store likelihood, AIC, and AICc
hisse.full_NoDub.logL <- hisse.full_NoDub$loglik
hisse.full_NoDub.AIC <- hisse.full_NoDub$AIC
hisse.full_NoDub.AICc <- hisse.full_NoDub$AICc

# Write to csv
# Added by Brendan
# Tried to make this work as a variable too
capture.output(hisse.full_NoDub, file=paste0(out_prefix, '_hisse_full_NoDub.txt'))


######################
# Compile Model Fits #
######################

# Write relevent outcome for model comparison to new vectors
logL <- c(Null2_AllTrans_Eq.logL, Null2_NoDub_Eq.logL, Null2_ThreeRate_NoDub.logL, bisse.logL, bisse.null.logL, hisse.null4.equal.logL, 
          hisse.null4.three.logL, no0B_AllTrans.logL, no0B_NoDub.logL, no1B_AllTrans.logL, no1B_NoDub.logL, hisse.full_AllTrans.logL, 
          hisse.full_NoDub.logL)
AIC <- c(Null2_AllTrans_Eq.AIC, Null2_NoDub_Eq.AIC, Null2_ThreeRate_NoDub.AIC, bisse.AIC, bisse.null.AIC, hisse.null4.equal.AIC, 
         hisse.null4.three.AIC, no0B_AllTrans.AIC, no0B_NoDub.AIC, no1B_AllTrans.AIC, no1B_NoDub.AIC, hisse.full_AllTrans.AIC, 
         hisse.full_NoDub.AIC)
AICc <- c(Null2_AllTrans_Eq.AICc, Null2_NoDub_Eq.AICc, Null2_ThreeRate_NoDub.AICc, bisse.AICc, bisse.null.AICc, hisse.null4.equal.AICc, 
          hisse.null4.three.AICc, no0B_AllTrans.AICc, no0B_NoDub.AICc, no1B_AllTrans.AICc, no1B_NoDub.AICc, hisse.full_AllTrans.AICc, 
          hisse.full_NoDub.AICc)

# Make a data frame to store model comparison/results to
results <- as.data.frame(logL, row.names=c('Null2 AllTrans EqRates', 'Null2 NoDub EqRates', 'Null2 ThreeRate NoDub', 'BiSSE', 'BiSSE Null', 
                                           'Null4 EqRates', 'Null4 ThreeRate', 'No0B AllTrans', 'No0B NoDub',  'No1B AllTrans', 'No1B NoDub', 
                                           'HiSSE Full AllTrans', 'HiSSE Full NoDub'))
results$AIC <- AIC
results$AICc <- AICc

# Variable outfile name instead of hard-coded:
write.csv(results, file=paste0(out_prefix, 'ModelSupport.csv'))
