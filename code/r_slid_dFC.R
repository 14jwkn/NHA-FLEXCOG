#Rscript --vanilla <script name> <subject> <width>

#Catches arguments.
args = commandArgs(trailingOnly=T)

#Get packages.
#install.packages('tidyverse')
#install.packages('glasso')
library(tidyverse) #Various useful functions.
library(glasso) #Graphical LASSO L1 regularization tool.

#Masked functions: cars, filter, lag, map, filter, poly, kaiser, %+%, alpha, tr, sim)

#Set code as working directory.

#Set parameters.
subject <- args[1] #Subject ID.
width <- as.numeric(args[2]) #Number of TRs in window.
nwin <- 1200 - width + 1 #Number of windows.

#Number of iterations for cross validation to find best lambda.
cviterations <- 10 

#Number of lambdas to test based on where the minimum negative likelihood is
#likely to be located for a small sample of subjects.
cvtestvals <- seq(from=0.004,to=0.07,length.out=12) 

#Define the file paths and read in as tibbles for L1 and R1.
L1path <- paste('../outputs/meants/',subject, 
                '/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv',
                sep='')
R1path <- paste('../outputs/meants/',subject, 
                '/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv',
                sep='')
L1 <- read_csv(L1path,col_names=F)
R1 <- read_csv(R1path,col_names=F)

#Sets Gaussian taper with sigma = 3.
sig <- 3
val <- c(0:(width-1)) - ((width-1)/2)
sig2 <- 2*sig*sig
taper <- exp(-(val**2)/sig2)

# Get the Windows ---------------------------------------------------------

#For each of the phase directions, get all the windows. Iterate through each 
#phase direction.
phases <- list(L1,R1)
windows <- list()
for (i in seq_along(phases)) {
  
  #Iterate through all the windows.
  win <- vector('list',nwin)
  for (j in c(1:nwin)) {
    
    #Get the window 360 region x width TR window and transpose it to width TR x
    #360 region matrix for next functions.
    start <- j
    stop <- j + width - 1
    basewin <- t(select(phases[[i]],all_of(start):all_of(stop))) 
    
    #Subtract each column timeseries by its mean and divide by sample standard
    #deviation in order to standardize into a Z-score.
    zwin <- basewin %>% as_tibble(.name_repair='universal') %>%
      purrr::map_dfc(~ (. - mean(.))/sd(.)) 
    
    #Apply the Gaussian taper to each column timeseries and add the processed
    #window into a list.
    tapwin <- taper * zwin 
    win[[j]] <- tapwin
  }
  
  #Append the windows for the phase direction to the rest of the windows.
  windows <- append(windows,win)
}


# Cross Validation for GLASSO Lambda --------------------------------------

#Set list to contain minimum lambdas.
minlams <- list()

#Make a copy of the window list to facilitate sampling of training windows 
#without replacement.
traincopy <- windows

#Repeat this the set number of times.
for (i in c(1:cviterations)) {
  
  #Select a random window as the training window and delete it from the copied list.
  indtrain <- sample(c(1:length(traincopy)),size=1)
  trainwin <- traincopy[[indtrain]]
  traincopy <- traincopy[-indtrain]
  
  #Get the covariance matrix from the training window.
  traincov <- cov(trainwin)
  
  #For each cross-validation test value, create an inverse covariance
  #matrix using GLASSO and save it with its corresponding lambda as a name.
  trainmats <- list()
  for (lam in cvtestvals) {
    trained <- glasso(traincov,rho=lam,nobs=width) #GLASSO.
    crosslam <- toString(lam)
    trainmats[[crosslam]] <- trained$wi #Inverse covariance matrix.
  }
  
  #Remove the training window from the windows to make a set of test windows
  #and concatenate.
  testwins <- windows[-indtrain]
  testcat <- data.frame(matrix(nrow=0,ncol=360))
  for (win in testwins) {
    testcat <- rbind(testcat,win)
  }
  
  #Get the covariance matrix from the concatenated test windows.
  testcov <- cov(testcat) 
  
  #Find the negative log-likelihood between training inverse covariance matrix
  #and the test covariance matrix.
  likvals <- list()
  for (lam in cvtestvals) {
    neglam <- toString(lam)
    traininv <- trainmats[[neglam]]
    likvals[[neglam]] <- -log(det(traininv)) + sum(diag(traininv%*%testcov))
  }
  
  #Find the minimum likelihood value, get its corresponding lambda,
  #convert the lambda into an integer for math, and append it for the iteration.
  minval <- as.numeric(names(which.min(likvals)))
  minlams <- append(minlams,minval)
}

#Average the minimum lambdas for each iteration to get the best lambda.
bestlam <- mean(unlist(minlams))

#Processes the full set of windows with GLASSO to produce a covariance matrix,
#r-to-z transforms it for statistics later, and adds it back to a list. Then, 
#concatenates and saves correlation matrices. Does not concatenate iteratively
#for speed.
winlist <- list()
wincat <- data.frame(matrix(nrow=360,ncol=0))

for (i in seq_along(windows)) { #Iterate through all windows.
  wincov <- cov(windows[[i]]) #Covariance matrix.
  lassocov <- glasso(wincov,rho=bestlam,nobs=width) #GLASSO.
  lassocor <- cov2cor(lassocov$w) #Correlation matrix from GLASSO covariance matrix.
  lassoz <- structure(sapply(lassocor,atanh),dim=dim(lassocor)) #r-to-z transform.
  diag(lassoz) <- 0 #Replace infinities in diagonal due to arctanh(1) with 0.
  winlist[[i]] <- data.frame(lassoz) #Saves it to the list.
}
wincat <- bind_cols(winlist) #Concatenate all by column.

#Defines outpath, creates intervening folders if it doesn't exist, and saves the dFC windows.
outpath <- paste('../outputs/r_slid_dFC/',toString(width),'/',subject,sep='')
dir.create(outpath,recursive=T)
write_csv(wincat,paste(outpath,'/r_dFC100.csv',sep=''),col_names=F)
