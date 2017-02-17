files <- dir(path='./', pattern='ModelSupport.csv')


lnLik.1 <- data.frame(matrix(ncol=4))
colnames(lnLik.1) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.2 <- data.frame(matrix(ncol=4))
colnames(lnLik.2) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.3 <- data.frame(matrix(ncol=4))
colnames(lnLik.3) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.4 <- data.frame(matrix(ncol=4))
colnames(lnLik.4) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.5 <- data.frame(matrix(ncol=4))
colnames(lnLik.5) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.6 <- data.frame(matrix(ncol=4))
colnames(lnLik.6) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.7 <- data.frame(matrix(ncol=4))
colnames(lnLik.7) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.8 <- data.frame(matrix(ncol=4))
colnames(lnLik.8) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.9 <- data.frame(matrix(ncol=4))
colnames(lnLik.9) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.10 <- data.frame(matrix(ncol=4))
colnames(lnLik.10) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.11 <- data.frame(matrix(ncol=4))
colnames(lnLik.11) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.12 <- data.frame(matrix(ncol=4))
colnames(lnLik.12) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.13 <- data.frame(matrix(ncol=4))
colnames(lnLik.13) <- c('Iteration', 'lnLik' , 'AIC','AICc')
lnLik.14 <- data.frame(matrix(ncol=4))
colnames(lnLik.14) <- c('Iteration', 'lnLik' , 'AIC','AICc')

# Null2 All Trans Equal Rates
for(i in 1:500){
    lnLik.1[i,] <- read.csv(file=files[i])[1, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.1 <- lnLik.1[order(abs(lnLik.1$lnLik)),]
best.1 <- print(lnLik.1[1,1])


# Null2 NoDub EqRates
for(i in 1:500){
    lnLik.2[i,] <- read.csv(file=files[i])[2, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.2 <- lnLik.2[order(abs(lnLik.2$lnLik)),]
best.2 <- print(lnLik.2[1,1])


# Null2 ThreeRate No Dub
for(i in 1:500){
    lnLik.3[i,] <- read.csv(file=files[i])[3, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.3 <- lnLik.3[order(abs(lnLik.3$lnLik)),]
best.3 <- print(lnLik.3[1,1])


# BiSSE
for(i in 1:500){
    lnLik.4[i,] <- read.csv(file=files[i])[4, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.4 <- lnLik.4[order(abs(lnLik.4$lnLik)),]
best.4 <- print(lnLik.4[1,1])


# BiSSE Null EqRates
for(i in 1:500){
    lnLik.5[i,] <- read.csv(file=files[i])[5, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.5 <- lnLik.5[order(abs(lnLik.5$lnLik)),]
best.5 <- print(lnLik.5[1,1])


#BiSSE Null UneqRates
for(i in 1:500){
    lnLik.6[i,] <- read.csv(file=files[i])[6, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.6 <- lnLik.6[order(abs(lnLik.6$lnLik)),]
best.6 <- print(lnLik.6[1,1])


# Null4 EqRates
for(i in 1:500){
    lnLik.7[i,] <- read.csv(file=files[i])[7, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.7 <- lnLik.7[order(abs(lnLik.7$lnLik)),]
best.7 <- print(lnLik.7[1,1])


#Null4 ThreeRate
for(i in 1:500){
    lnLik.8[i,] <- read.csv(file=files[i])[8, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.8 <- lnLik.8[order(abs(lnLik.8$lnLik)),]
best.8 <- print(lnLik.8[1,1])


# No 0B All Trans
for(i in 1:500){
    lnLik.9[i,] <- read.csv(file=files[i])[9, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.9 <- lnLik.9[order(abs(lnLik.9$lnLik)),]
best.9 <- print(lnLik.9[1,1])


# No 0B No Double
for(i in 1:500){
    lnLik.10[i,] <- read.csv(file=files[i])[10, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.10 <- lnLik.10[order(abs(lnLik.10$lnLik)),]
best.10 <- print(lnLik.10[1,1])


# No 1B All Trans
for(i in 1:500){
    lnLik.11[i,] <- read.csv(file=files[i])[11, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.11 <- lnLik.11[order(abs(lnLik.11$lnLik)),]
best.11 <- print(lnLik.11[1,1])


# No 1B No Double
for(i in 1:500){
    lnLik.12[i,] <- read.csv(file=files[i])[12, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.12 <- lnLik.12[order(abs(lnLik.12$lnLik)),]
best.12 <- print(lnLik.12[1,1])


#HiSSE Full All Trans
for(i in 1:500){
    lnLik.13[i,] <- read.csv(file=files[i])[13, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.13 <- lnLik.13[order(abs(lnLik.13$lnLik)),]
best.13 <- print(lnLik.13[1,1])


#HiSSE Full No Double
for(i in 1:500){
    lnLik.14[i,] <- read.csv(file=files[i])[14, c(5,2,3,4)]
}
#Sort resultant dataframe by lnLik, such that the minimum lnLik is the first row and print the iteration 
#number resulting in the best fit while saving this to a vector
lnLik.14 <- lnLik.14[order(abs(lnLik.14$lnLik)),]
best.14 <- print(lnLik.14[1,1])







#Save best fit iteration as new file in folder containing the top fits for each of the 14 models
setwd('/srv/storage/Dropbox/Dissertation_copy/PW.Caudate.HiSSE/Official Analysis/Final_Supported_Models/LHS.Optim/Frequency_75.36_Bounded_Eps100/Models/')
dest <- "/srv/storage/Dropbox/Dissertation_copy/PW.Caudate.HiSSE/Official Analysis/Final_Supported_Models/LHS.Optim/OptimalFit/freq(75.35)_Bounded_Eps100"

file.rename(from = paste0("./output.", best.1, "_Null2_AllTrans_EqRate.txt"), to = paste0(dest, best.1, "_Null2_AllTrans_EqRate.txt"))
file.rename(from = paste0("./output.", best.2, "_Null2_NoDub_EqRate.txt"), to = paste0(dest, best.2, "_Null2_NoDub_EqRate.txt"))
file.rename(from = paste0("./output.", best.3, "_Null2_NoDub_ThreeRate.txt"), to = paste0(dest, best.3, "_Null2_NoDub_ThreeRate.txt"))
file.rename(from = paste0("./output.", best.4, "_BiSSE.txt"), to = paste0(dest, best.4, "_BiSSE.txt"))
file.rename(from = paste0("./output.", best.5, "_BiSSE_Null_UnequalTrans.txt"), to = paste0(dest, best.5, "_BiSSE_Null_UnequalTrans.txt"))
file.rename(from = paste0("./output.", best.6, "_BiSSE_Null_EqualTrans.txt"), to = paste0(dest, best.6, "_BiSSE_Null_EqualTrans.txt"))
file.rename(from = paste0("./output.", best.7, "_Null4_Equal.txt"), to = paste0(dest, best.7, "_Null4_Equal.txt"))
file.rename(from = paste0("./output.", best.8, "_Null4_ThreeRate.txt"), to = paste0(dest, best.8, "_Null4_ThreeRate.txt"))
file.rename(from = paste0("./output.", best.9, "_No0B_AllTrans.txt"), to = paste0(dest, best.9, "_No0B_AllTrans.txt"))
file.rename(from = paste0("./output.", best.10, "_No0B_NoDub.txt"), to = paste0(dest, best.10, "_No0B_NoDub.txt"))
file.rename(from = paste0("./output.", best.11, "_No1B_AllTrans.txt"), to = paste0(dest, best.11, "_No1B_AllTrans.txt"))
file.rename(from = paste0("./output.", best.12, "_No1B_NoDub.txt"), to = paste0(dest, best.12, "_No1B_NoDub.txt"))
file.rename(from = paste0("./output.", best.13, "_hisse_full_AllTrans.txt"), to = paste0(dest, best.13, "_hisse_full_AllTrans.txt"))
file.rename(from = paste0("./output.", best.14, "_hisse_full_NoDub.txt"), to = paste0(dest, best.14, "_hisse_full_NoDub.txt"))
