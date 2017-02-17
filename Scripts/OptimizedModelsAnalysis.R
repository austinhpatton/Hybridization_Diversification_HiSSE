### Basic script to run all models and reparameterizations and save resultant objects and results to respecive folders ###
library(hisse)
library(ggplot2)
library(reshape)

# Read in tree
tree <- read.nexus('./pw.caud.ultra.nexus')

# Read in states
dat <- read.csv('./Caudata.PyronWiens.SppStates.csv')

####################################
# Null 2 All Transitions All Equal #
####################################
trans.rate <- TransMatMaker(hidden.states=TRUE)
trans.rate[!is.na(trans.rate) & !trans.rate == 0] = 1

start.vals <- c(0.0214445204,0.9218754317,0.7819032506)
# Fixed turnover rates between 0A & 1A, and between 0B & 1B
# This runs the analysis
Null2_AllTrans_Eq <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,1,2,2),
                           eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                           starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_Null2_AllTrans_Eq <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=Null2_AllTrans_Eq$solution, hidden.states = T, aic=Null2_AllTrans_Eq$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(Null2_AllTrans_Eq, file = 'Null2_AllTrans_Eq_Object.Rsave')
save(recon_Null2_AllTrans_Eq, file='recon_Null2_AllTrans_Eq.Rsave')


##########################################
# Null 2 No Double Transitions All Equal #
##########################################
trans.rate <- TransMatMaker(hidden.states=TRUE)
trans.rate <- ParDrop(trans.rate, c(3,5,8,10))
trans.rate[!is.na(trans.rate) & !trans.rate == 0] = 1

start.vals <- c(0.2661567331,0.9992454882,0.3141280462)

# Fixed turnover rates between 0A & 1A, and between 0B & 1B
# This runs the analysis
Null2_NoDub_Eq <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,1,2,2),
                        eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                        starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_Null2_NoDub_Eq <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=Null2_NoDub_Eq$solution, hidden.states = T, aic=Null2_NoDub_Eq$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(Null2_NoDub_Eq, file = 'Null2_NoDub_Eq_Object.Rsave')
save(recon_Null2_NoDub_Eq, file='recon_Null2_NoDub_Eq.Rsave')

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

start.vals <- c(0.9240797866,0.8350535093,0.296512568)

# Fixed turnover rates between 0A & 1A, and between 0B & 1B
# This runs the analysis
Null2_ThreeRate_NoDub <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,1,2,2),
                               eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                               starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_Null2_ThreeRate_NoDub <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=Null2_ThreeRate_NoDub$solution, hidden.states = T, aic=Null2_ThreeRate_NoDub$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(Null2_ThreeRate_NoDub, file = 'Null2_ThreeRate_NoDub_Object.Rsave')
save(recon_Null2_ThreeRate_NoDub, file='recon_Null2_ThreeRate_NoDub.Rsave')

###############
# BiSSE Model #
###############

# Now make the bisse model where diversification changes with hybridizability without
# the presence of hidden states. This will be the BiSSE Null model.

trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

start.vals <- c(0.1317859475,0.038159527,0.4435103945)

# The transition matrix thus looks like the following
#     (0) (1)
# (0)  NA   2
# (1)   1  NA

# Given that the order of arguments in the turnover and extinction .anc commands that 
# control the number of rate classes follows this order... (0A, 1A, 0B, 1B), the
# following will set the model to only transition between states that do not include
# hidden states.
bisse <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = F, turnover.anc = c(1,2,0,0),
               eps.anc = c(1,2,0,0), trans.rate = trans.rates.bisse, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
               starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_bisse <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=bisse$solution, hidden.states = F, aic=bisse$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(bisse, file = 'bisse_Object.Rsave')
save(recon_bisse, file='recon_bisse.Rsave')

#################################################################
#                   Bisse Null Model - Unequal Trans            #
#################################################################

start.vals <- c(0.5700023509,0.1727187972,0.0927389375)

# Make a constrained bisse model where diversification rates are trait independent
bisse.null.uneq <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = F, turnover.anc = c(1,1,0,0),
                         eps.anc = c(1,1,0,0), trans.rate = trans.rates.bisse, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                         starting.vals=start.vals)


# Get support region by sampling likelihood surface

recon_bisse.null.uneq <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=bisse.null.uneq$solution, hidden.states = F, aic=bisse.null.uneq$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(bisse.null.uneq, file = 'bisse.null.uneq_Object.Rsave')
save(recon_bisse.null.uneq, file='recon_bisse.null.uneq.Rsave')
#################################################################
#                   Bisse Null Model - Equal Trans              #
#################################################################

trans.rate <- TransMatMaker(hidden.states=FALSE)
trans.rate <- ParEqual(trans.rate, c(1,2))

start.vals <- c(0.1110613327,0.6205475213,0.0197306661)

# Make a constrained bisse model where diversification rates are trait independent
bisse.null.eq <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = F, turnover.anc = c(1,1,0,0),
                       eps.anc = c(1,1,0,0), trans.rate = trans.rate, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                       starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_bisse.null.eq <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=bisse.null.eq$solution, hidden.states = F, aic=bisse.null.eq$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(bisse.null.eq, file = 'bisse.null.eq_Object.Rsave')
save(recon_bisse.null.eq, file='recon_bisse.null.eq.Rsave')

###################################
# HiSSE Null 4 Model, Equal Rates #
###################################

start.vals <- c(0.6365689601,0.4693103231,0.2083357648)

# Conduct the HiSSE null-4 model that contains the same complexity as a full HiSSE 
# model
hisse.null4.equal <- hisse.null4(phy = tree, data = dat, f = c(0.75, 0.36),  turnover.anc = rep(c(1,2,3,4),2),
                                 eps.anc = rep(c(1,2,3,4),2), trans.type = "equal", output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                                 starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_hisse.null4.equal <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=hisse.null4.equal$solution, hidden.states = T, aic=hisse.null4.equal$AICc, n.cores = 14, four.state.null = T, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(hisse.null4.equal, file = 'hisse.null4.equal_Object.Rsave')
save(recon_hisse.null4.equal, file='recon_hisse.null4.equal.Rsave')

###################################
# HiSSE Null 4 Model, Three Rates #
###################################

start.vals <- c(0.7951909704,0.3709173504,0.9817792378)

# Conduct the HiSSE null-4 model that contains the same complexity as a full HiSSE 
# model
hisse.null4.three <- hisse.null4(phy = tree, data = dat, f = c(0.75, 0.36),  turnover.anc = rep(c(1,2,3,4),2),
                                 eps.anc = rep(c(1,2,3,4),2), trans.type = "three.rate", output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                                 starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_hisse.null4.three <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=hisse.null4.three$solution, hidden.states = T, aic=hisse.null4.three$AICc, n.cores = 14, four.state.null = T, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(hisse.null4.three, file = 'hisse.null4.three_Object.Rsave')
save(recon_hisse.null4.three, file='recon_hisse.null4.three.Rsave')

##############################
# HiSSE No0B All Transitions #
##############################

# Make a model that has hidden states for only state 1 - we thus cannot have transitions to 0B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(2,5,12,7,8,9))

start.vals <- c(0.7027920966,0.6516460436,0.4410054008)

no0B_AllTrans <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,0,3),
                       eps.anc = c(1,2,0,3), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                       starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_no0B_AllTrans <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=no0B_AllTrans$solution, hidden.states = T, aic=no0B_AllTrans$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(no0B_AllTrans, file = 'no0B_AllTrans_Object.Rsave')
save(recon_no0B_AllTrans, file='recon_no0B_AllTrans.Rsave')

####################################
# HiSSE No0B No Double Transitions #
####################################

# Make a model that has hidden states for only state 1 - we thus cannot have transitions to 0B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(2,3,5,7,8,9,10,12))

start.vals <- c(0.199301583,0.147152632,0.7621309992)

no0B_NoDub <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,0,3),
                    eps.anc = c(1,2,0,3), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                    starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_no0B_NoDub <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=no0B_NoDub$solution, hidden.states = T, aic=no0B_NoDub$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(no0B_NoDub, file = 'no0B_NoDub_Object.Rsave')
save(recon_no0B_NoDub, file='recon_no0B_NoDub.Rsave')

##############################
# HiSSE No1B All Transitions #
##############################

# Make a model that has hidden states for only state 1 - we thus cannot have transitions to 0B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(3,6,9,10,11,12))

start.vals <- c(0.8351176892,0.4745255933,0.0574330162)

no1B_AllTrans <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,0),
                       eps.anc = c(1,2,3,0), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                       starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_no1B_AllTrans <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=no1B_AllTrans$solution, hidden.states = T, aic=no1B_AllTrans$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(no1B_AllTrans, file = 'no1B_AllTrans_Object.Rsave')
save(recon_no1B_AllTrans, file='recon_no1B_AllTrans.Rsave')

####################################
# HiSSE No1B No Double Transitions #
####################################

# Make a model that has hidden states for only state 1 - we thus cannot have transitions to 0B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(3,5,6,8,9,10,11,12))

start.vals <- c(0.9821846831,0.530902172,0.0988427541)

no1B_NoDub <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,0),
                    eps.anc = c(1,2,3,0), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                    starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_no1B_NoDub <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=no1B_NoDub$solution, hidden.states = T, aic=no1B_NoDub$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(no1B_NoDub, file = 'no1B_NoDub_Object.Rsave')
save(recon_no1B_NoDub, file='recon_no1B_NoDub.Rsave')

##############################
# Full HiSSE All Transitions #
##############################

# Make a model that has hidden states for all states but does not allow for a transition between 
# state 0A and 1B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)

start.vals <- c(0.0899926809,0.0590837358,0.1105176878)

hisse.full_AllTrans <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,4),
                             eps.anc = c(1,2,3,4), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                             starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_hisse.full_AllTrans <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=hisse.full_AllTrans$solution, hidden.states = T, aic=hisse.full_AllTrans$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(hisse.full_AllTrans, file = 'hisse.full_AllTrans_Object.Rsave')
save(recon_hisse.full_AllTrans, file='recon_hisse.full_AllTrans.Rsave')

####################################
# Full HiSSE No Double Transitions #
####################################

# Make a model that has hidden states for all states but does not allow for a transition between 
# state 0A and 1B
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))

start.vals <- c(0.1205516774,0.5354679876,0.0051798163)

hisse.full_NoDub_eps.200 <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,4),
                          eps.anc = c(1,2,3,4), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE, eps.upper = 200,
                          starting.vals=start.vals)

hisse.full_NoDub <- hisse(phy = tree, data = dat, f = c(0.75, 0.36), hidden.states = T, turnover.anc = c(1,2,3,4),
                          eps.anc = c(1,2,3,4), trans.rate = trans.rates.hisse, output.type = "raw", bounded.search = TRUE, eps.upper = 50,
                          starting.vals=start.vals)

# Get support region by sampling likelihood surface

recon_hisse.full_NoDub <- MarginRecon(tree, dat, f=c(0.75, 0.36), pars=hisse.full_NoDub$solution, hidden.states = T, aic=hisse.full_NoDub$AICc, n.cores = 14, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(hisse.full_NoDub, file = 'hisse.full_NoDub_Object.Rsave')
save(recon_hisse.full_NoDub, file='recon_hisse.full_NoDub.Rsave')
hisse.full_NoDub_SupReg <- SupportRegion(hisse.obj = hisse.full_NoDub, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
save(hisse.full_NoDub_SupReg, file='hisse.full_NoDub_SupReg.Rsave')


# Now adaptively sample for each - these have been moved later to save time. 
Null2_AllTrans_Eq_SupReg <- SupportRegion(hisse.obj = Null2_AllTrans_Eq, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
Null2_NoDub_Eq_SupReg <- SupportRegion(hisse.obj = Null2_NoDub_Eq, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
Null2_ThreeRate_NoDub_SupReg <- SupportRegion(hisse.obj = Null2_ThreeRate_NoDub, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
bisse_SupReg <- SupportRegion(hisse.obj = bisse, n.points = 3000, output.type = "raw", hidden.states = F, desired.delta = 1.92)
bisse.null.uneq_SupReg <- SupportRegion(hisse.obj = bisse.null.uneq, n.points = 3000, output.type = "raw", hidden.states = F, desired.delta = 1.92)
bisse.null.eq_SupReg <- SupportRegion(hisse.obj = bisse.null.eq, n.points = 3000, output.type = "raw", hidden.states = F, desired.delta = 1.92)
hisse.null4.equal_SupReg <- SupportRegion(hisse.obj = hisse.null4.equal, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
hisse.null4.three_SupReg <- SupportRegion(hisse.obj = hisse.null4.three, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
no0B_AllTrans_SupReg <- SupportRegion(hisse.obj = no0B_AllTrans, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
no0B_NoDub_SupReg <- SupportRegion(hisse.obj = no0B_NoDub, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
no1B_AllTrans_SupReg <- SupportRegion(hisse.obj = no1B_AllTrans, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
no1B_NoDub_SupReg <- SupportRegion(hisse.obj = no1B_NoDub, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
hisse.full_AllTrans_SupReg <- SupportRegion(hisse.obj = hisse.full_AllTrans, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)
hisse.full_NoDub_SupReg <- SupportRegion(hisse.obj = hisse.full_NoDub, n.points = 3000, output.type = "raw", hidden.states = TRUE, desired.delta = 1.92)

# Save them
save(Null2_AllTrans_Eq_SupReg, file='Null2_AllTrans_Eq_SupReg.Rsave')
save(Null2_NoDub_Eq_SupReg, file='Null2_NoDub_Eq_SupReg.Rsave')
save(Null2_ThreeRate_NoDub_SupReg, file='Null2_ThreeRate_NoDub_SupReg.Rsave')
save(bisse_SupReg, file='bisse_SupReg.Rsave')
save(bisse.null.uneq_SupReg, file='bisse.null.uneq_SupReg.Rsave')
save(bisse.null.eq_SupReg, file='bisse.null.eq_SupReg.Rsave')
save(hisse.null4.equal_SupReg, file='hisse.null4.equal_SupReg.Rsave')
save(hisse.null4.three_SupReg, file='hisse.null4.three_SupReg.Rsave')
save(no0B_AllTrans_SupReg, file='no0B_AllTrans_SupReg.Rsave')
save(no0B_NoDub_SupReg, file='no0B_NoDub_SupReg.Rsave')
save(no1B_AllTrans_SupReg, file='no1B_AllTrans_SupReg.Rsave')
save(no1B_NoDub_SupReg, file='no1B_NoDub_SupReg.Rsave')
save(hisse.full_AllTrans_SupReg, file='hisse.full_AllTrans_SupReg.Rsave')
save(hisse.full_NoDub_SupReg, file='hisse.full_NoDub_SupReg.Rsave')
