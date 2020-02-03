#' BESIDE-MR fitting function
#'
#' Fits BESIDE-MR model. NOTE: beta is the estimated causal effect and tau is the heterogeneity variance modelled as precision (i.e. 1/tau)
#'
#' @param tau_estimate Use DL estimate (="DL_approx") or Full Bayesian (="Full_Bayes") approach to analyse data.
#' @param N_Beta 1 parameter (=1) or 2 parameter (=2) model.
#' @param BetaXG Effect size for X-G association
#' @param BetaYG Effect size for Y-G association
#' @param seBetaXG Standard error for X-G association
#' @param seBetaYG Standard error for Y-G association
#' @param N_Ins Number of genetic instruments
#' @param N_Iter Number of iterations
#' @param Prior Hyper parameter for the prior distribution;
#'                 - For 1 parameter model
#'                    1. $hyper_Beta_mean and $hyper_Beta_sd to specify mean and standard deviation for
#'                         the normally distributed beta
#'                    2. $hyper_Prec_shape and $hyper_Prec_rate to specify shape and rate for
#'                         the gamma distribution of precision ("Full_Bayes" ONLY)
#'                    3. $Ins_prob assign prior inclusion probability for each instrument
#'                 - For 2 parameter model
#'                    1. $hyper_Beta1_mean and $hyper_Beta1_sd to specify mean and standard deviation for the normally
#'                         distributed beta1
#'                    2. $hyper_Beta2_mean and $hyper_Beta2_sd to specify mean and standard deviation for the normally
#'                         distributed beta2
#'                    3. $hyper_Prec1_shape and $hyper_Prec1_rate to specify shape and rate for the gamma distribution
#'                         of precision ("Full_Bayes" ONLY)
#'                    4. $hyper_Prec2_shape and $hyper_Prec2_rate to specify shape and rate for the gamma distribution
#'                         of precision ("Full_Bayes" ONLY)
#'                    5. $Ins1_prob assign prior inclusion probability for each instrument in set 1
#'                    6. $Ins2_prob assign prior inclusion probability for each instrument in set 2
#' @param tuning_para Tuning parameter to ensure sufficient acceptance rate (recommended between 0.25 - 0.45)
#'                 - For 1 parameter model
#'                    1. $Beta for beta
#'                    2. $Prec_LL, $Prec_UL and $Prec_gap for the upper, lower and gap for the precision respectively,
#'                         this ensures symmetry between the new and the old value ("Full_Bayes" ONLY)
#'                 - For 2 parameter model
#'                    1. $Beta1 for beta1
#'                    2. $Beta2 for beta2
#'                    3. $Prec1_LL, $Prec1_UL and $Prec1_gap for the upper, lower and gap for the precision1 respectively,
#'                         ("Full_Bayes" ONLY)
#'                    4. $Prec2_LL, $Prec2_UL and $Prec2_gap for the upper, lower and gap for the precision2 respectively,
#'                         ("Full_Bayes" ONLY)
#' @param gen_inits Initial values to start the iterations
#                 - For 1 parameter model
#                    1. $Beta for beta
#                    2. $UBPrec and $LBPrec for upper and lower limit of initial value of precision respectively,
#                         ("Full_Bayes" ONLY)
#                    3. $Ins_L, use randomS.initial.LI() generate random initial model space
#                 - For 2 parameter model
#                    1. $Beta1 for beta1
#                    2. $Beta2 for beta2
#                    3. $UBPrec1 and $LBPrec1 for upper and lower limit of initial value of precision1 respectively,
#                         ("Full_Bayes" ONLY)
#                    4. $UBPrec2 and $LBPrec2 for upper and lower limit of initial value of precision2 respectively,
#                         ("Full_Bayes" ONLY)
#                    5. $Ins1_L, use randomS.initial.LI() generate random initial model space
#                    6. $Ins2_L, use randomS.initial.LI() generate random initial model space
#' @return An object of class \code{"beside"} containing the following components:\describe{
#' \item{\code{S}}{A matrix giving the results.}
#' \item{\code{accept_rate}}{acceptance rate.}
#'}
#'@author Chin Yang Shapland; Jack Bowden.
#'@references Shapland, C.Y., et al., Profile-likelihood Bayesian model averaging for two-sample summary data Mendelian randomization in the presence of horizontal pleiotropy.
#'@export
#'@examples
#'
#' Prior choice for beta, tau and inclusion of instruments
#' Ins_prior<-rep(0.5, L)
#' Prior_DL<-list(hyper_Beta_mean=0, hyper_Beta_sd=1, Ins_prob=Ins_prior)
#' Prior_gamma<-list(hyper_Beta_mean=0, hyper_Beta_sd=1, hyper_Prec_shape=2, hyper_Prec_rate=0.00005, Ins_prob=Ins_prior)
#'
#' Tuning parameter for beta
#' H_DL<-list(Beta=0.05)
#' H_gamma<-list(Beta=0.05, Prec_LL=0, Prec_UL=1000000, Prec_gap=150000)
#'
#' Generate initial values
#' gen_inits_DL<-list(Beta=rnorm(1,0,10), Ins_L=randomS.initial.LI(rep(0,L), L,Ins_prior))
#' gen_inits_gamma<-list(Beta=rnorm(1,0,10), UBPrec=1000000, LBPrec=0, Ins_L=randomS.initial.LI(rep(0,L), L,Ins_prior))
#'
#' M-H algorithm
#' res_DL<-BMA_MRanalysis("DL_approx", data$BetaXG,data$BetaYG,data$seBetaXG,data$seBetaYG, L, nIter, Prior_DL, H_DL, gen_inits_DL)
#' res_gamma<-BMA_MRanalysis("Full_Bayes", data$BetaXG,data$BetaYG,data$seBetaXG,data$seBetaYG, L, nIter, Prior_gamma, H_gamma, gen_inits_gamma)

BMA_MRanalysis<-function(tau_estimate, N_Beta, BetaXG,BetaYG,seBetaXG,seBetaYG, N_Ins, N_Iter, Prior,
                         tuning_para, gen_inits){

  if (tau_estimate=="DL_approx" & N_Beta==1){

    #Generate initial values
    Beta  <- gen_inits$Beta
    Ins_L <- gen_inits$Ins_L

    # Saved space
    S      <- matrix(0,N_Iter,N_Ins+2)
    accept <- matrix(1,N_Iter,2)       #acceptance rate

    for(i in 1:N_Iter) {
      oldBeta  <- Beta
      oldIns_L <- Ins_L

      #1st step: Compare likelihood for old and new Beta
      Beta       <- Beta + rnorm(1,0,tuning_para$Beta)
      oldlogpost <- logpost_DL(BetaXG,BetaYG,seBetaXG,seBetaYG, oldBeta, oldIns_L, Prior$hyper_Beta_mean, Prior$hyper_Beta_sd)
      newlogpost <- logpost_DL(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, oldIns_L, Prior$hyper_Beta_mean, Prior$hyper_Beta_sd)
      oldlp      <- oldlogpost$logpost
      newlp      <- newlogpost$logpost
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Beta         <- oldBeta          #fails so return to oldb
        accept[i, 1] <- 0                #0 to indicate rejection of new value
      }

      #3nd step: Compare likelihood for old and new model space
      Ins_L      <- randomS.LI(Ins_L, N_Ins, Prior$Ins_prob)
      oldlogpost <- logpost_DL(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, oldIns_L, Prior$hyper_Beta_mean, Prior$hyper_Beta_sd)
      newlogpost <- logpost_DL(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, Ins_L, Prior$hyper_Beta_mean, Prior$hyper_Beta_sd)

      tausq_hat  <- newlogpost$tausq

      oldlp <- oldlogpost$logpost
      newlp <- newlogpost$logpost
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Ins_L        <- oldIns_L          #fails so return to oldIL
        tausq_hat    <- oldlogpost$tausq
        accept[i, 2] <- 0                 #0 to indicate rejection of new value
      }
      S[i,] <- c(Beta, tausq_hat ,Ins_L)
    }
    returnlist<-list(S=S, accept_rate=apply(accept, 2, mean))
    return(returnlist)
  }

  if (tau_estimate=="DL_approx" & N_Beta==2){

    #Generate initial values
    Beta1  <- gen_inits$Beta1
    Beta2  <- gen_inits$Beta2
    Ins_L1 <- gen_inits$Ins_L1
    Ins_L2 <- 1-Ins_L1

    # Saved space
    S      <- matrix(0,N_Iter,(N_Ins*2)+4)
    accept <- matrix(1,N_Iter,4)

    for(i in 1:N_Iter) {
      oldBeta1  <- Beta1
      oldBeta2  <- Beta2
      oldIns_L1 <- Ins_L1
      oldIns_L2 <- Ins_L2

      #1st step: Compare likelihood for old and new Beta1
      Beta1      <- Beta1 + rnorm(1,0,tuning_para$Beta1_h)
      oldlogpost <- logpost_2Beta_DL(BetaXG, BetaYG, seBetaXG, seBetaYG, oldBeta1, oldBeta2, oldIns_L1,
                                    Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd, oldIns_L2,
                                    Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd)
      newlogpost <- logpost_2Beta_DL(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, oldBeta2, oldIns_L1,
                                    Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd, oldIns_L2,
                                    Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd)

      oldlp <- oldlogpost$logpost
      newlp <- newlogpost$logpost
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Beta1        <- oldBeta1          # fails so return to oldb
        accept[i, 1] <- 0
      }

      #2nd step: Compare likelihood for old and new Beta2
      Beta2      <- Beta2 + rnorm(1,0,tuning_para$Beta2_h)
      oldlogpost <- logpost_2Beta_DL(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, oldBeta2, oldIns_L1,
                                    Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd, oldIns_L2,
                                    Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd)
      newlogpost <- logpost_2Beta_DL(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, oldIns_L1,
                                    Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd, oldIns_L2,
                                    Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd)

      oldlp <- oldlogpost$logpost
      newlp <- newlogpost$logpost
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Beta2        <- oldBeta2          # fails so return to oldb
        accept[i, 2] <- 0
      }

      #3rd step: Compare likelihood for old and new model space L1
      Ins_L1     <- randomS_cond.LI(Ins_L1, Ins_L2, N_Ins, Prior$Ins1_prob)
      oldlogpost <- logpost_2Beta_DL(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, oldIns_L1,
                                    Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd, oldIns_L2,
                                    Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd)
      newlogpost <- logpost_2Beta_DL(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Ins_L1,
                                    Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd, oldIns_L2,
                                    Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd)

      tausq1_hat <- newlogpost$tausq1
      tausq2_hat <- newlogpost$tausq2

      oldlp <- oldlogpost$logpost
      newlp <- newlogpost$logpost
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Ins_L1       <- oldIns_L1        # fails so return to oldIL
        tausq1_hat   <- oldlogpost$tausq1
        accept[i, 3] <- 0
      }

      #4th step: Compare likelihood for old and new model space L2
      Ins_L2     <- randomS_cond.LI(Ins_L2, Ins_L1, N_Ins, Prior$Ins2_prob)
      oldlogpost <- logpost_2Beta_DL(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Ins_L1,
                                    Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd, oldIns_L2,
                                    Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd)
      newlogpost <- logpost_2Beta_DL(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Ins_L1,
                                    Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd, Ins_L2,
                                    Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd)

      tausq2_hat <- newlogpost$tausq2

      oldlp <- oldlogpost$logpost
      newlp <- newlogpost$logpost
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Ins_L2       <- oldIns_L2        # fails so return to oldIL
        tausq2_hat   <-oldlogpost$tausq2
        accept[i, 4] <- 0
      }

      S[i,] <- c(Beta1, Beta2, tausq1_hat, tausq2_hat, Ins_L1, Ins_L2)
    }
    returnlist<-list(S=S, accept_rate=apply(accept, 2, mean))
    return(returnlist)
  }

  if (tau_estimate=="Full_Bayes" & N_Beta==1){

    #Generate initial values
    Beta   <- gen_inits$Beta
    UBPrec <- gen_inits$UBPrec
    LBPrec <- gen_inits$LBPrec
    Prec   <- runif(1,LBPrec,UBPrec)
    Ins_L  <- gen_inits$Ins_L

    # Saved space
    S      <- matrix(0,N_Iter,N_Ins+2)
    accept <- matrix(1,N_Iter,3)

    for(i in 1:N_Iter) {
      oldBeta   <- Beta
      oldPrec   <- Prec
      oldUBPrec <- UBPrec
      oldLBPrec <- LBPrec
      oldIns_L  <- Ins_L

      #1st step: Compare likelihood for old and new Beta
      Beta       <- Beta + rnorm(1,0,tuning_para$Beta)
      oldlogpost <- logpost_gammaTau(BetaXG,BetaYG,seBetaXG,seBetaYG, oldBeta, oldPrec, oldIns_L,
                                 Prior$hyper_Beta_mean, Prior$hyper_Beta_sd, Prior$hyper_Prec_shape, Prior$hyper_Prec_rate)
      newlogpost <- logpost_gammaTau(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, oldPrec, oldIns_L,
                                 Prior$hyper_Beta_mean, Prior$hyper_Beta_sd, Prior$hyper_Prec_shape, Prior$hyper_Prec_rate)
      oldlp <- oldlogpost+log(oldUBPrec-oldLBPrec)
      newlp <- newlogpost+log(oldUBPrec-oldLBPrec)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Beta         <- oldBeta          # fails so return to oldb
        accept[i, 1] <- 0
      }

      #2nd step: Compare likelihood for old and new tausq
      UBPrec     <- min(tuning_para$Prec_UL,Prec+tuning_para$Prec_gap)
      LBPrec     <- max(tuning_para$Prec_LL,Prec-tuning_para$Prec_gap)
      Prec       <- runif(1,LBPrec,UBPrec)
      oldlogpost <- logpost_gammaTau(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, oldPrec, oldIns_L,
                                 Prior$hyper_Beta_mean, Prior$hyper_Beta_sd, Prior$hyper_Prec_shape, Prior$hyper_Prec_rate)
      newlogpost <- logpost_gammaTau(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, Prec, oldIns_L,
                                 Prior$hyper_Beta_mean, Prior$hyper_Beta_sd, Prior$hyper_Prec_shape, Prior$hyper_Prec_rate)

      oldlp <- oldlogpost+log(UBPrec-LBPrec)
      newlp <- newlogpost+log(oldUBPrec-oldLBPrec)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Prec         <- oldPrec          # fails so return to oldPrecsq
        UBPrec       <- oldUBPrec
        LBPrec       <- oldLBPrec
        accept[i, 2] <- 0
      }

      #3nd step: Compare likelihood for old and new model space
      Ins_L      <- randomS.LI(Ins_L, N_Ins, Prior$Ins_prob)
      oldlogpost <- logpost_gammaTau(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, Prec, oldIns_L,
                                 Prior$hyper_Beta_mean, Prior$hyper_Beta_sd, Prior$hyper_Prec_shape, Prior$hyper_Prec_rate)
      newlogpost <- logpost_gammaTau(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, Prec, Ins_L,
                                 Prior$hyper_Beta_mean, Prior$hyper_Beta_sd, Prior$hyper_Prec_shape, Prior$hyper_Prec_rate)

      oldlp <- oldlogpost+log(UBPrec-LBPrec)
      newlp <- newlogpost+log(UBPrec-LBPrec)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Ins_L        <- oldIns_L        # fails so return to oldIL
        accept[i, 3] <- 0
      }
      S[i,] <- c(Beta, 1/Prec ,Ins_L)
    }
    returnlist<-list(S=S, accept_rate=apply(accept, 2, mean))
    return(returnlist)
  }

  if (tau_estimate=="Full_Bayes" & N_Beta==2){

    #Generate initial values
    Beta1   <- gen_inits$Beta1
    Beta2   <- gen_inits$Beta2

    UBPrec1 <- gen_inits$UBPrec1
    LBPrec1 <- gen_inits$LBPrec1
    Prec1   <- runif(1,LBPrec1,UBPrec1)

    UBPrec2 <- gen_inits$UBPrec2
    LBPrec2 <- gen_inits$LBPrec2
    Prec2   <- runif(1,LBPrec2,UBPrec2)

    Ins_L1  <- gen_inits$Ins_L1
    Ins_L2  <- 1-Ins_L1

    # Saved space
    S      <- matrix(0,N_Iter,(N_Ins*2)+4)
    accept <- matrix(1,N_Iter,6)

    for(i in 1:N_Iter) {
      oldBeta1   <- Beta1
      oldBeta2   <- Beta2

      oldPrec1   <- Prec1
      oldUBPrec1 <- UBPrec1
      oldLBPrec1 <- LBPrec1

      oldPrec2   <- Prec2
      oldUBPrec2 <- UBPrec2
      oldLBPrec2 <- LBPrec2

      oldIns_L1  <- Ins_L1
      oldIns_L2  <- Ins_L2

      #1st step: Compare likelihood for old and new Beta1
      Beta1      <- Beta1 + rnorm(1,0,tuning_para$Beta1)
      oldlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, oldBeta1, oldBeta2, oldPrec1, oldPrec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)

      newlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, oldBeta2, oldPrec1, oldPrec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)

      oldlp <- oldlogpost + log(oldUBPrec1-oldLBPrec1) + log(oldUBPrec2-oldLBPrec2)
      newlp <- newlogpost + log(oldUBPrec1-oldLBPrec1) + log(oldUBPrec2-oldLBPrec2)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Beta1        <- oldBeta1          # fails so return to oldb
        accept[i, 1] <- 0
      }

      #2nd step: Compare likelihood for old and new Beta2
      Beta2      <- Beta2 + rnorm(1,0,tuning_para$Beta2)
      oldlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, oldBeta2, oldPrec1, oldPrec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)

      newlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, oldPrec1, oldPrec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)

      oldlp <- oldlogpost + log(oldUBPrec1-oldLBPrec1) + log(oldUBPrec2-oldLBPrec2)
      newlp <- newlogpost + log(oldUBPrec1-oldLBPrec1) + log(oldUBPrec2-oldLBPrec2)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Beta2       <- oldBeta2          # fails so return to oldb
        accept[i, 2]<- 0
      }

      #3rd step: Compare likelihood for old and new tausq1
      UBPrec1    <- min(tuning_para$Prec1_UL,Prec1+tuning_para$Prec1_gap)
      LBPrec1    <- max(tuning_para$Prec1_LL,Prec1-tuning_para$Prec1_gap)
      Prec1      <- runif(1,LBPrec1,UBPrec1)
      oldlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, oldPrec1, oldPrec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)
      newlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Prec1, oldPrec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)

      oldlp <- oldlogpost + log(UBPrec1-LBPrec1) + log(oldUBPrec2-oldLBPrec2)
      newlp <- newlogpost + log(oldUBPrec1-oldLBPrec1) + log(oldUBPrec2-oldLBPrec2)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Prec1        <- oldPrec1          # fails so return to oldPrecsq
        UBPrec1      <- oldUBPrec1
        LBPrec1      <- oldLBPrec1
        accept[i, 3] <- 0
      }

      #4th step: Compare likelihood for old and new tausq2
      UBPrec2    <- min(tuning_para$Prec2_UL,Prec2+tuning_para$Prec2_gap)
      LBPrec2    <- max(tuning_para$Prec2_LL,Prec2-tuning_para$Prec2_gap)
      Prec2      <- runif(1,LBPrec2,UBPrec2)
      oldlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Prec1, oldPrec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)
      newlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Prec1, Prec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)

      oldlp <- oldlogpost + log(UBPrec1-LBPrec1) + log(UBPrec2-LBPrec2)
      newlp <- newlogpost + log(UBPrec1-LBPrec1) + log(oldUBPrec2-oldLBPrec2)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Prec2       <- oldPrec2          # fails so return to oldPrecsq
        UBPrec2      <- oldUBPrec2
        LBPrec2      <- oldLBPrec2
        accept[i, 4] <- 0
      }

      #5th step: Compare likelihood for old and new model space L1
      Ins_L1     <- randomS_cond.LI(Ins_L1, Ins_L2, N_Ins, Prior$Ins1_prob)
      oldlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Prec1, Prec2,
                                         oldIns_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)
      newlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Prec1, Prec2,
                                         Ins_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)

      oldlp <- oldlogpost + log(UBPrec1-LBPrec1) + log(UBPrec2-LBPrec2)
      newlp <- newlogpost + log(UBPrec1-LBPrec1) + log(UBPrec2-LBPrec2)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Ins_L1       <- oldIns_L1       # fails so return to oldIL
        accept[i, 5] <- 0
      }

      #6th step: Compare likelihood for old and new model space L2
      Ins_L2     <- randomS_cond.LI(Ins_L2, Ins_L1, N_Ins, Prior$Ins2_prob)
      oldlogpost <-logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Prec1, Prec2,
                                         Ins_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         oldIns_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)
      newlogpost <- logpost_2Beta_gammaTau(BetaXG, BetaYG, seBetaXG, seBetaYG, Beta1, Beta2, Prec1, Prec2,
                                         Ins_L1, Prior$hyper_Beta1_mean, Prior$hyper_Beta1_sd,
                                         Prior$hyper_Prec1_shape, Prior$hyper_Prec1_rate,
                                         Ins_L2, Prior$hyper_Beta2_mean, Prior$hyper_Beta2_sd,
                                         Prior$hyper_Prec2_shape, Prior$hyper_Prec2_rate)
      oldlp <- oldlogpost + log(UBPrec1-LBPrec1) + log(UBPrec2-LBPrec2)
      newlp <- newlogpost + log(UBPrec1-LBPrec1) + log(UBPrec2-LBPrec2)
      if( log(runif(1,0,1)) > (newlp-oldlp) ) {
        Ins_L2       <- oldIns_L2       # fails so return to oldIL
        accept[i, 6] <- 0
      }

      S[i,] <- c(Beta1, Beta2, 1/Prec1, 1/Prec2 ,Ins_L1, Ins_L2)
    }

    returnlist<-list(S=S, accept_rate=apply(accept, 2, mean))
    class(returnlist)<-"beside"

    return(returnlist)
    }

}

### randomS.initial.LI() ###
#Generate random initial model space that doesn't produce empty and 1 variable model space

#User specified options
#"L"         Number of instruments
#"ins_prior" prior for inclusion probability for each instrument

randomS.initial.LI <- function(L, ins_prior) {
  Ind_L <- rep(0,L)
  repeat {
    # do something
    Ind_L <- rbinom(L, 1, ins_prior)
    # exit if the condition is met
    if (sum(Ind_L) >= 5) break
  }
  return(Ind_L)
}

######################### Subfunctions that are used by BMA_MRanalysis() ###########################################

### randomS.LI() ###
#Generate random model space that doesn't produce empty and 1 variable model space

randomS.LI <- function(Ind_L, L, ins_prior) {
  repeat {
    # do something
    a<-sample(1:L, 1, prob = ins_prior)
    Ind_L[a]<-(Ind_L[a]-1)^2
    # exit if the condition is met
    if (sum(Ind_L) >= 5) break
  }
  return(Ind_L)
}

### randomS_cond.LI() ###
#Generate random model space for 2 parameter model, generation of new model space is condition on the other model
#so that an 1 instruments cannot appear in both model space

randomS_cond.LI <- function(Ind_L1, Ind_L2, L, ins_prior) {
  repeat {
    # do something
    ins_prior_sub<-ins_prior[which(Ind_L2==0)]
    a<-sample(which(Ind_L2==0), 1, prob = ins_prior_sub)
    Ind_L1[a]<-(Ind_L1[a]-1)^2
    # exit if the condition is met
    if (sum(Ind_L1) >= 5 & all((Ind_L1+Ind_L2)<2)) break
  }
  return(Ind_L1)
}

### DL() ###
#DerSimonian and Laird estimate
#From Bowden et al. 2019, IJE

DL = function(y=BIV,s=se.BIV){

  k          = length(y)
  w          = 1/s^2           ; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  ; W      = sum.w - (sum(w^2)/sum.w)
  Q          = sum(w*(y-mu.hat)^2)
  tausq.hat  = max(0,(Q-(k-1))/W)

  Isq = max(0,(Q-(k-1))/Q)

  return(list(tausq.hat=tausq.hat,Isq=Isq))

}

### logpost_DL() ###
#Calculate the log posterior with hyper parameters option and where tau is approximate from DL()

logpost_DL <- function(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, indicator, hyper_Beta_mean, hyper_Beta_sd) {

  #Subsetting instruments
  BetaXG_sub   <- BetaXG[which(indicator==1)]
  BetaYG_sub   <- BetaYG[which(indicator==1)]
  seBetaXG_sub <- seBetaXG[which(indicator==1)]
  seBetaYG_sub <- seBetaYG[which(indicator==1)]

  #Estimating tau from DL()
  BIV    <- BetaYG_sub/BetaXG_sub
  se.BIV <- sqrt((seBetaYG_sub^2)/BetaXG_sub^2)
  tausq  <- DL(BIV,se.BIV)$tausq.hat

  logpost = -0.5*sum(indicator)*log(2*pi) + (-0.5*sum(log(seBetaYG_sub^2+tausq) +
            ((BetaYG_sub-(Beta*BetaXG_sub))^2/(Beta^2*seBetaXG_sub^2+seBetaYG_sub^2+tausq)))) +  # log likelihood
            dnorm(Beta, hyper_Beta_mean, hyper_Beta_sd, log=T)                                   # priors in beta

  return(list(tausq=tausq, logpost=logpost))
}

### logpost_2Beta_DL() ###
#Calculate the log posterior for the 2 parameter model with hyper parameters option and
#where tau is approximate from DL()

logpost_2Beta_DL <- function(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta1, Beta2, indicator1, hyper_Beta1_mean,
                             hyper_Beta1_sd, indicator2, hyper_Beta2_mean, hyper_Beta2_sd) {

  #Estimating tau from DL()
  BIV    <- BetaYG/BetaXG
  se.BIV <- sqrt((seBetaYG^2)/BetaXG^2)
  tausq1 <- DL(BIV[which(indicator1==1)],se.BIV[which(indicator1==1)])$tausq.hat
  tausq2 <- DL(BIV[which(indicator2==1)],se.BIV[which(indicator2==1)])$tausq.hat

  #Weighting for the likelihood for each instrument. this allows independent probability of instrument inclusion
  w1<-ifelse(indicator1==0 & indicator2==0, 0, indicator1/(indicator1+indicator2))
  w2<-ifelse(indicator1==0 & indicator2==0, 0, indicator2/(indicator1+indicator2))

  logpost = (-0.5*sum(indicator1*w1)*log(2*pi)) + (-0.5*sum(indicator1*w1*(log(seBetaYG^2+tausq1) +
            ((BetaYG-(Beta1*BetaXG))^2/(Beta1^2*seBetaXG^2+seBetaYG^2+tausq1))))) +         # log likelihood for beta1
            dnorm(Beta1, hyper_Beta1_mean, hyper_Beta1_sd, log=T) +                         # priors in beta1
            (-0.5*sum(indicator2*w2)*log(2*pi)) + (-0.5*sum(indicator2*w2*(log(seBetaYG^2+tausq2) +
            ((BetaYG-(Beta2*BetaXG))^2/(Beta2^2*seBetaXG^2+seBetaYG^2+tausq2))))) +         # log likelihood for beta2
            dnorm(Beta2, hyper_Beta2_mean, hyper_Beta2_sd, log=T)                           # priors in beta2

  return(list(tausq1=tausq1, tausq2=tausq2, logpost=logpost))
}


### logpost_gammaTau() ###
#Calculate the log posterior with hyper parameters option

logpost_gammaTau <- function(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta, Prec, indicator, hyper_Beta_mean, hyper_Beta_sd,
                             hyper_Prec_shape, hyper_Prec_rate) {

  #Subsetting instruments
  BetaXG_sub   <-BetaXG[which(indicator==1)]
  BetaYG_sub   <-BetaYG[which(indicator==1)]
  seBetaXG_sub <-seBetaXG[which(indicator==1)]
  seBetaYG_sub <-seBetaYG[which(indicator==1)]

  #Calculate tau in terms of precision
  Tausq<-1/Prec

  logpost = -0.5*sum(indicator)*log(2*pi) + (-0.5*sum(log(seBetaYG_sub^2+Tausq) +
            ((BetaYG_sub-(Beta*BetaXG_sub))^2/(Beta^2*seBetaXG_sub^2+seBetaYG_sub^2+Tausq)))) + #log likelihood
            dnorm(Beta, hyper_Beta_mean, hyper_Beta_sd, log=T) +                                #priors in beta
            dgamma(Prec, shape=hyper_Prec_shape, rate=hyper_Prec_rate, log=T)                   #priors in tau

  return(logpost=logpost)
}

### logpost_2Beta_gammaTau() ###
#Calculate the log posterior for the 2 parameter model with hyper parameters option

logpost_2Beta_gammaTau<- function(BetaXG,BetaYG,seBetaXG,seBetaYG, Beta1, Beta2, Prec1, Prec2,
                                  indicator1, hyper_Beta1_mean, hyper_Beta1_sd, hyper_Prec1_shape, hyper_Prec1_rate,
                                  indicator2, hyper_Beta2_mean, hyper_Beta2_sd, hyper_Prec2_shape, hyper_Prec2_rate) {

  #Calculate tau in terms of precision
  Tausq1<-1/Prec1
  Tausq2<-1/Prec2

  #Weighting for the likelihood for each instrument. this allows independent probability of instrument inclusion
  w1<-ifelse(indicator1==0 & indicator2==0, 0, indicator1/(indicator1+indicator2))
  w2<-ifelse(indicator1==0 & indicator2==0, 0, indicator2/(indicator1+indicator2))

  logpost = (-0.5*sum(indicator1*w1)*log(2*pi)) + (-0.5*sum(indicator1*w1*(log(seBetaYG^2+Tausq1) +
            ((BetaYG-(Beta1*BetaXG))^2/(Beta1^2*seBetaXG^2+seBetaYG^2+Tausq1))))) +                 #log likelihood
            dnorm(Beta1, hyper_Beta1_mean, hyper_Beta1_sd, log=T) +                                 #priors in beta1
            dgamma(Prec1, shape=hyper_Prec1_shape, rate=hyper_Prec1_rate, log=T) +                  #priors in tau1
            (-0.5*sum(indicator2*w2)*log(2*pi)) + (-0.5*sum(indicator2*w2*(log(seBetaYG^2+Tausq2) +
            ((BetaYG-(Beta2*BetaXG))^2/(Beta2^2*seBetaXG^2+seBetaYG^2+Tausq2))))) +                 #log likelihood
            dnorm(Beta2, hyper_Beta2_mean, hyper_Beta2_sd, log=T) +                                 #priors in beta1
            dgamma(Prec2, shape=hyper_Prec2_shape, rate=hyper_Prec2_rate, log=T)                    #priors in tau2

  return(logpost=logpost)
}
#######################################################################################################################
