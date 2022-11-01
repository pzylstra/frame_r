#' Finds the RSE for the mean of a vector
#'
#' @param vec A vector of the values being predicted
#' @return value
#' @export

mRSE <- function(dat){
  dat<-dat[!is.na(dat)]
  m <- mean(dat)
  ssq <- 0
  for (val in dat) {
    sq <- (val - m)^2
    ssq <- ssq + sq
  }
  rse <- sqrt(ssq)
  return(rse)
}


#' Builds models for cover dynamics of surveyed Species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @param bTest Multiples of mean + mRSE for which Burr & quadratic models can predict 
#' beyond the observed mean + standard deviation
#' @param maxiter The maximum number of iterations for model fitting
#' @return dataframe
#' @export

coverDyn <- function(dat, thres = 5, pnts = 10, p = 0.05, bTest = 10, maxiter = 1000) {
  
  dat <- datClean(veg = dat,  base = "base", top = "top", he = "he", ht = "ht")
  spCov <- specCover(dat = dat, thres = thres, pnts = pnts)
  priorList <- unique(spCov$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitCov <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_Rsq' = numeric(0), 'lin_p' = numeric(0),
                       'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_Rsq' = numeric(0), 'NE_p' = numeric(0),
                       'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_Rsq' = numeric(0), 'B_p' = numeric(0),
                       'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_Rsq' = numeric(0), 'Bin_p' = numeric(0),
                       'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'q_Rsq' = numeric(0), 'Q_p' = numeric(0), 
                       'Mean' = character(0), 'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=maxiter, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- spCov %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$Cover)
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
      LM<-lm(y ~ x)
      LMSum <- base::summary(LM)
      LMa <- LMSum$coefficients[2]
      LMb <- LMSum$coefficients[1]
      LMRSE <- LMSum$sigma
      LMRsq <- LMSum$r.squared
      LMp <- if (LMRSE == 0) {
        0
      } else {
        round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
      }
      rm(LM)
    } else {
      LMa <- NA
      LMb <- NA
      LMRSE <- 100
      LMRsq <- 0
      LMp <- 1
    }
    
    #Negative exponential
    init1<-c(k=50,r=0.5)
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
      NESum <- base::summary(NE)
      k <- NESum$coefficients[1]
      r <- NESum$coefficients[2]
      NERSE <- NESum$sigma
      NERsq <- cor(predict(NE, newdata=x), y)**2
      NEp <- round(max(NESum$coefficients[7],NESum$coefficients[8]),5)
      rm(NE)
    } else {
      k <- NA
      r <- NA
      NERSE <- 100
      NERsq <- 0
      NEp <- 1
    }
    
    #Burr
    init2<-c(a=3,b=2)
    if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
      Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
      BSum <- base::summary(Burr)
      Ba <- BSum$coefficients[1]
      Bb <- BSum$coefficients[2]
      
      #Added control for Burr
      f <- function(x){Ba*Bb*((0.1*x^(Ba-1))/((1+(0.1*x)^Ba)^Bb+1))}
      if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        BRSE <- 100
        BRsq <- 0
        Bp <- 1
      } else if (bTest * (mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE)) > (optimize(f = f, interval=c(0, 150), maximum=FALSE))$objective) {
        BRSE <- 100
        BRsq <- 0
        Bp <- 1
      } else {
        BRSE <- BSum$sigma
        BRsq <- cor(predict(Burr, newdata=x), y)**2
        Bp <- round(max(BSum$coefficients[7],BSum$coefficients[8]),5)
      }
      rm(Burr)
    } else {
      Ba <- NA
      Bb <- NA
      BRSE <- 100
      BRsq <- 0
      Bp <- 1
    }
    
    #Binomial
    init3<-c(s=1,sd=3, m=20)
    if (!berryFunctions::is.error(nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control))) {
      Bin<-nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control)
      BinSum <- base::summary(Bin)
      Bs <- BinSum$coefficients[1]
      Bsd <- BinSum$coefficients[2]
      Bm <- BinSum$coefficients[3]
      BinRSE <- BinSum$sigma
      BinRsq <- cor(predict(Bin, newdata=x), y)**2
      Binp <- round(max(BinSum$coefficients[10],BinSum$coefficients[11],BinSum$coefficients[12]),5)
      rm(Bin)
    } else {
      Bs <- NA
      Bsd <- NA
      Bm <- NA
      BinRSE <- 100
      BinRsq <- 0
      Binp <- 1
    }
    
    #Quadratic
    init4 <- c(a = -1, b = 2, c = 0)
    if (!berryFunctions::is.error(nls(y ~ a*x^2 + b*x + c, data = studySpecies, 
                                      start = init4, trace = T, control = control))) {
      q <- nls(y ~ a*x^2 + b*x + c, data = studySpecies, start = init4, 
               trace = T, control = control)
      qSum <- base::summary(q)
      qa <- qSum$coefficients[1]
      qb <- qSum$coefficients[2]
      qc <- qSum$coefficients[3]
      
      # Added control for Quadratic
      # Also prevents quadratic increases
      f <- function(x){qa*x^2 + qb*x + qc}
      if ((bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective)||qa > 0) {
        qRSE <- 100
        qRsq <- 0
        qp <- 1
      } else {
        qRSE <- qSum$sigma
        qRsq <- cor(predict(q, newdata=x), y)**2
        qp <- round(max(qSum$coefficients[10], qSum$coefficients[11], 
                        qSum$coefficients[12]), 5)
      }
      rm(q)
    } else {
      qa <- NA
      qb <- NA
      qc <- NA
      qRSE <- 100
      qRsq <- 0
      qp <- 1
    }
    
    #Summary stats
    meanCov <- round(mean(y, na.rm = TRUE),1)
    m_sig <- round(mRSE(dat = y),3)
    
    #Choose the best model
    listRsq <- c(LMRsq, NERsq, BRsq, BinRsq, qRsq)
    listRSE <- c(LMRSE, NERSE, BRSE, BinRSE, qRSE)
    listMod <- c("Linear", "NegExp", "Burr", "Binomial", "Quadratic")
    
    if (listRSE[which(listRsq == max(listRsq))]<=m_sig) {
      model <- listMod[which(listRsq == max(listRsq))]
    } else {
      model <= "Mean"
    }
    
    #Record values
    fitCov[nrow(fitCov)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMRsq, LMp, k, r, NERSE, NERsq, NEp,
                                 Ba, Bb, BRSE, BRsq, Bp, Bs, Bsd, Bm, BinRSE, BinRsq, Binp, qa, qb, qc, qRSE, qRsq, qp, 
                                 meanCov, m_sig, model)
  }
  
  return(fitCov)
}



#' Builds models for top height dynamics of surveyed Species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @param bTest Multiples of mean + mRSE for which Burr & quadratic models can predict 
#' beyond the observed mean + standard deviation
#' @param maxiter The maximum number of iterations for model fitting
#' @return dataframe
#' @export

topDyn <- function(dat, base = "base", top = "top", he = "he", ht = "ht", 
                   thres = 5, pnts = 10, p = 0.05, bTest = 10, maxiter = 1000) {
  
  dat <- datClean(veg = dat,  base = "base", top = "top", he = "he", ht = "ht")
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"))
  
  priorList <- unique(dat$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitTop <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_Rsq' = numeric(0), 'lin_p' = numeric(0),
                       'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_Rsq' = numeric(0), 'NE_p' = numeric(0),
                       'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_Rsq' = numeric(0), 'B_p' = numeric(0),
                       'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_Rsq' = numeric(0), 'Bin_p' = numeric(0),
                       'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'q_Rsq' = numeric(0), 'Q_p' = numeric(0), 
                       'Mean' = character(0), 'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=maxiter, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- dat %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$top)
    
    if (length(unique(x, incomparables = FALSE))>2 & length(unique(y, incomparables = FALSE))>1) {
    
      #Linear
      if (!berryFunctions::is.error(lm(y ~ x))) {
        LM<-lm(y ~ x)
        LMSum <- base::summary(LM)
        LMa <- LMSum$coefficients[2]
        LMb <- LMSum$coefficients[1]
        LMRSE <- LMSum$sigma
        LMRsq <- LMSum$r.squared
        LMp <- if (LMRSE == 0) {
          0
        } else {
          round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
        }
        rm(LM)
      } else {
        LMa <- NA
        LMb <- NA
        LMRSE <- 100
        LMRsq <- 0
        LMp <- 1
      }
      
      #Negative exponential
      init1<-c(k=50,r=0.5)
      if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T))) {
        NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
        NESum <- base::summary(NE)
        k <- NESum$coefficients[1]
        r <- NESum$coefficients[2]
        NERSE <- NESum$sigma
        NERsq <- cor(predict(NE, newdata=x), y)**2
        NEp <- round(max(NESum$coefficients[7],NESum$coefficients[8]),5)
        rm(NE)
      } else {
        k <- NA
        r <- NA
        NERSE <- 100
        NERsq <- 0
        NEp <- 1
      }
      
      #Burr
      init2<-c(a=3,b=2)
      if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
        Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
        BSum <- base::summary(Burr)
        Ba <- BSum$coefficients[1]
        Bb <- BSum$coefficients[2]
        
        #Added control for Burr
        f <- function(x){Ba*Bb*((0.1*x^(Ba-1))/((1+(0.1*x)^Ba)^Bb+1))}
        if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
          BRSE <- 100
          BRsq <- 0
          Bp <- 1
        } else if (bTest * (mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE)) > (optimize(f = f, interval=c(0, 150), maximum=FALSE))$objective) {
          BRSE <- 100
          BRsq <- 0
          Bp <- 1
        } else {
          BRSE <- BSum$sigma
          BRsq <- cor(predict(Burr, newdata=x), y)**2
          Bp <- round(max(BSum$coefficients[7],BSum$coefficients[8]),5)
        }
        rm(Burr)
      } else {
        Ba <- NA
        Bb <- NA
        BRSE <- 100
        BRsq <- 0
        Bp <- 1
      }
      
      #Binomial
      init3<-c(s=1,sd=3, m=20)
      if (!berryFunctions::is.error(nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control))) {
        Bin<-nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control)
        BinSum <- base::summary(Bin)
        Bs <- BinSum$coefficients[1]
        Bsd <- BinSum$coefficients[2]
        Bm <- BinSum$coefficients[3]
        BinRSE <- BinSum$sigma
        BinRsq <- cor(predict(Bin, newdata=x), y)**2
        Binp <- round(max(BinSum$coefficients[10],BinSum$coefficients[11],BinSum$coefficients[12]),5)
        rm(Bin)
      } else {
        Bs <- NA
        Bsd <- NA
        Bm <- NA
        BinRSE <- 100
        BinRsq <- 0
        Binp <- 1
      }
      
      #Quadratic
      init4 <- c(a = 1, b = 2, c = 0)
      if (!berryFunctions::is.error(nls(y ~ a*x^2 + b*x + c, data = studySpecies, 
                                        start = init4, trace = T, control = control))) {
        q <- nls(y ~ a*x^2 + b*x + c, data = studySpecies, start = init4, 
                 trace = T, control = control)
        qSum <- base::summary(q)
        qa <- qSum$coefficients[1]
        qb <- qSum$coefficients[2]
        qc <- qSum$coefficients[3]
        
        #Added control for Quadratic
        f <- function(x){qa*x^2 + qb*x + qc}
        if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
          qRSE <- 100
          qRsq <- 0
          qp <- 1
        } else {
          qRSE <- qSum$sigma
          qRsq <- cor(predict(q, newdata=x), y)**2
          qp <- round(max(qSum$coefficients[10], qSum$coefficients[11], 
                          qSum$coefficients[12]), 5)
        }
        rm(q)
      } else {
        qa <- NA
        qb <- NA
        qc <- NA
        qRSE <- 100
        qRsq <- 0
        qp <- 1
      }
      
      #Summary stats
      meanTop <- round(mean(y, na.rm = TRUE),1)
      m_sig <- round(mRSE(dat = y),3)
      
      #Choose the best model
      listRsq <- c(LMRsq, NERsq, BRsq, BinRsq, qRsq)
      listRSE <- c(LMRSE, NERSE, BRSE, BinRSE, qRSE)
      listMod <- c("Linear", "NegExp", "Burr", "Binomial", "Quadratic")
      
      if (listRSE[which(listRsq == max(listRsq))]<=m_sig) {
        model <- listMod[which(listRsq == max(listRsq))]
      } else {
        model <- "Mean"
      }
      
    } else {
      LMa <- NA
      LMb <- NA
      LMRSE <- 100
      LMRsq <- 0
      LMp <- 1
      k <- NA
      r <- NA
      NERSE <- 100
      NERsq <- 0
      NEp <- 1
      Ba <- NA
      Bb <- NA
      BRSE <- 100
      BRsq <- 0
      Bp <- 1
      Bs <- NA
      Bsd <- NA
      Bm <- NA
      BinRSE <- 100
      BinRsq <- 0
      Binp <- 1
      qa <- NA
      qb <- NA
      qc <- NA
      qRSE <- 100
      qRsq <- 0
      qp <- 1
      model <- "Mean"      
      #Summary stats
      meanTop <- round(mean(y, na.rm = TRUE),1)
      m_sig <- round(mRSE(dat = y),3)
    }
    
    #Record values
    fitTop[nrow(fitTop)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMRsq, LMp, k, r, NERSE, NERsq, NEp,
                                 Ba, Bb, BRSE, BRsq, Bp, Bs, Bsd, Bm, BinRSE, BinRsq, Binp, qa, qb, qc, qRSE, qRsq, qp, 
                                 meanTop, m_sig, model)
    }
  
  return(fitTop)
}


#' Builds models for top-base height allometry dynamics of surveyed Species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @param bTest Multiples of mean + mRSE for which Burr & quadratic models can predict 
#' beyond the observed mean + standard deviation
#' @param maxiter The maximum number of iterations for model fitting
#' @return dataframe
#' @export

baseDyn <- function(dat, base = "base", top = "top", he = "he", ht = "ht", 
                    thres = 5, pnts = 10, p = 0.05, bTest = 10, maxiter = 1000) {
  
  dat <- datClean(veg = dat,  base = "base", top = "top", he = "he", ht = "ht")
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"),
           bRat = base/top)
  
  priorList <- unique(dat$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitBase <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_Rsq' = numeric(0), 'lin_p' = numeric(0),
                        'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_Rsq' = numeric(0), 'NE_p' = numeric(0),
                        'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_Rsq' = numeric(0), 'B_p' = numeric(0),
                        'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_Rsq' = numeric(0), 'Bin_p' = numeric(0),
                        'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'q_Rsq' = numeric(0), 'Q_p' = numeric(0), 
                        'Mean' = character(0), 'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=maxiter, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- dat %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$bRat)
    
    if (length(unique(x, incomparables = FALSE))>2 & length(unique(y, incomparables = FALSE))>1) {
      
      #Linear
      if (!berryFunctions::is.error(lm(y ~ x))) {
        LM<-lm(y ~ x)
        LMSum <- base::summary(LM)
        LMa <- LMSum$coefficients[2]
        LMb <- LMSum$coefficients[1]
        LMRSE <- LMSum$sigma
        LMRsq <- LMSum$r.squared
        LMp <- if (LMRSE == 0) {
          0
        } else {
          round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
        }
        rm(LM)
      } else {
        LMa <- NA
        LMb <- NA
        LMRSE <- 100
        LMRsq <- 0
        LMp <- 1
      }
      
      #Negative exponential
      init1<-c(k=50,r=0.5)
      if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T))) {
        NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
        NESum <- base::summary(NE)
        k <- NESum$coefficients[1]
        r <- NESum$coefficients[2]
        NERSE <- NESum$sigma
        NERsq <- cor(predict(NE, newdata=x), y)**2
        NEp <- round(max(NESum$coefficients[7],NESum$coefficients[8]),5)
        rm(NE)
      } else {
        k <- NA
        r <- NA
        NERSE <- 100
        NERsq <- 0
        NEp <- 1
      }
      
      #Burr
      init2<-c(a=3,b=2)
      if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
        Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
        BSum <- base::summary(Burr)
        Ba <- BSum$coefficients[1]
        Bb <- BSum$coefficients[2]
        
        #Added control for Burr
        f <- function(x){Ba*Bb*((0.1*x^(Ba-1))/((1+(0.1*x)^Ba)^Bb+1))}
        if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
          BRSE <- 100
          BRsq <- 0
          Bp <- 1
        } else if (bTest * (mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE)) > (optimize(f = f, interval=c(0, 150), maximum=FALSE))$objective) {
          BRSE <- 100
          BRsq <- 0
          Bp <- 1
        } else {
          BRSE <- BSum$sigma
          BRsq <- cor(predict(Burr, newdata=x), y)**2
          Bp <- round(max(BSum$coefficients[7],BSum$coefficients[8]),5)
        }
        rm(Burr)
      } else {
        Ba <- NA
        Bb <- NA
        BRSE <- 100
        BRsq <- 0
        Bp <- 1
      }
      
      #Binomial
      init3<-c(s=1,sd=3, m=20)
      if (!berryFunctions::is.error(nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control))) {
        Bin<-nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control)
        BinSum <- base::summary(Bin)
        Bs <- BinSum$coefficients[1]
        Bsd <- BinSum$coefficients[2]
        Bm <- BinSum$coefficients[3]
        BinRSE <- BinSum$sigma
        BinRsq <- cor(predict(Bin, newdata=x), y)**2
        Binp <- round(max(BinSum$coefficients[10],BinSum$coefficients[11],BinSum$coefficients[12]),5)
        rm(Bin)
      } else {
        Bs <- NA
        Bsd <- NA
        Bm <- NA
        BinRSE <- 100
        BinRsq <- 0
        Binp <- 1
      }
      
      #Quadratic
      init4 <- c(a = 1, b = 2, c = 0)
      if (!berryFunctions::is.error(nls(y ~ a*x^2 + b*x + c, data = studySpecies, 
                                        start = init4, trace = T, control = control))) {
        q <- nls(y ~ a*x^2 + b*x + c, data = studySpecies, start = init4, 
                 trace = T, control = control)
        qSum <- base::summary(q)
        qa <- qSum$coefficients[1]
        qb <- qSum$coefficients[2]
        qc <- qSum$coefficients[3]
        
        #Added control for Quadratic
        f <- function(x){qa*x^2 + qb*x + qc}
        if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
          qRSE <- 100
          qRsq <- 0
          qp <- 1
        } else {
          qRSE <- qSum$sigma
          qRsq <- cor(predict(q, newdata=x), y)**2
          qp <- round(max(qSum$coefficients[10], qSum$coefficients[11], 
                          qSum$coefficients[12]), 5)
        }
        rm(q)
      } else {
        qa <- NA
        qb <- NA
        qc <- NA
        qRSE <- 100
        qRsq <- 0
        qp <- 1
      }
      
      #Summary stats
      meanBase <- round(mean(y, na.rm = TRUE),1)
      m_sig <- round(mRSE(dat = y),3)
      
      #Choose the best model
      listRsq <- c(LMRsq, NERsq, BRsq, BinRsq, qRsq)
      listRSE <- c(LMRSE, NERSE, BRSE, BinRSE, qRSE)
      listMod <- c("Linear", "NegExp", "Burr", "Binomial", "Quadratic")
      
      if (listRSE[which(listRsq == max(listRsq))]<=m_sig) {
        model <- listMod[which(listRsq == max(listRsq))]
      } else {
        model <- "Mean"
      }
      
    } else {
      LMa <- NA
      LMb <- NA
      LMRSE <- 100
      LMRsq <- 0
      LMp <- 1
      k <- NA
      r <- NA
      NERSE <- 100
      NERsq <- 0
      NEp <- 1
      Ba <- NA
      Bb <- NA
      BRSE <- 100
      BRsq <- 0
      Bp <- 1
      Bs <- NA
      Bsd <- NA
      Bm <- NA
      BinRSE <- 100
      BinRsq <- 0
      Binp <- 1
      qa <- NA
      qb <- NA
      qc <- NA
      qRSE <- 100
      qRsq <- 0
      qp <- 1
      model <- "Mean"      
      #Summary stats
      meanBase <- round(mean(y, na.rm = TRUE),1)
      m_sig <- round(mRSE(dat = y),3)
    }
    
    #Record values
    fitBase[nrow(fitBase)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMRsq, LMp, k, r, NERSE, NERsq, NEp,
                                   Ba, Bb, BRSE, BRsq, Bp, Bs, Bsd, Bm, BinRSE, BinRsq, Binp, qa, qb, qc, qRSE, qRsq, qp, 
                                   meanBase, m_sig, model)
  }
  
  return(fitBase)
}


#' Builds models for top-he height allometry dynamics of surveyed Species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @return dataframe
#' @export

heDyn <- function(dat, thres = 5, pnts = 10, p = 0.05, 
                  base = "base", top = "top", he = "he", ht = "ht") {
  
  dat <- datClean(veg = dat,  base = "base", top = "top", he = "he", ht = "ht")
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"),
           bRat = he/top)
  
  priorList <- unique(dat$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fithe <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                      'linear' = character(0), 'Mean' = character(0), 'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    studySpecies <- dat %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$bRat)
    
    if (length(unique(x, incomparables = FALSE))>2 & length(unique(y, incomparables = FALSE))>1) {
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x)) && length(unique(x, incomparables = FALSE))>1) {
      LM<-lm(y ~ x)
      LMSum <- base::summary(LM)
      LMa <- LMSum$coefficients[2]
      LMb <- LMSum$coefficients[1]
      LMRSE <- LMSum$sigma
      LMp <- if (LMRSE == 0) {
        0
      } else {
        round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
      }
      LM_sig <- if (berryFunctions::is.error(LMp)) {
        ""
      }
      else if (LMp < 0.001) {
        "***"    
      } else if (LMp < 0.01) {
        "**"    
      } else if (LMp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(LM)
    } else {
      LMa <- NA
      LMb <- NA
      LRSE <- NA
      LMp <- 1
      LM_sig <- ""
    }
    }
    #Summary stats
    meanhe <- round(mean(y, na.rm = TRUE),1)
    m_sig <- round(mRSE(dat = y),3)
    model <- if (LMp < p) {
      "Linear"
    } else {"Mean"}
    model <- if (LMRSE < m_sig){
      model
    } else {"Mean"}
    
    #Record values
    fithe[nrow(fithe)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMp, 
                               LM_sig, meanhe, m_sig, model)
  }
  
  return(fithe)
}

#' Builds models for top-ht height allometry dynamics of surveyed Species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @return dataframe
#' @export

htDyn <- function(dat, thres = 5, pnts = 10, p = 0.05) {
  
  dat <- datClean(veg = dat,  base = "base", top = "top", he = "he", ht = "ht")
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"),
           bRat = ht/top)
  
  priorList <- unique(dat$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitht <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                      'linear' = character(0), 'Mean' = character(0), 'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    studySpecies <- dat %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$bRat)
    
    if (length(unique(x, incomparables = FALSE))>2 & length(unique(y, incomparables = FALSE))>1) {
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x)) && length(unique(x, incomparables = FALSE))>1) {
      LM<-lm(y ~ x)
      LMSum <- base::summary(LM)
      LMa <- LMSum$coefficients[2]
      LMb <- LMSum$coefficients[1]
      LMRSE <- LMSum$sigma
      LMp <- if (LMRSE == 0) {
        0
      } else {
        round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
      }
      LM_sig <- if (berryFunctions::is.error(LMp)) {
        ""
      }
      else if (LMp < 0.001) {
        "***"    
      } else if (LMp < 0.01) {
        "**"    
      } else if (LMp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(LM)
    } else {
      LMa <- NA
      LMb <- NA
      LRSE <- NA
      LMp <- 1
      LM_sig <- ""
    }
    }
    #Summary stats
    meanht <- round(mean(y, na.rm = TRUE),1)
    m_sig <- round(mRSE(dat = y),3)
    model <- if (LMp < p) {
      "Linear"
    } else {"Mean"}
    model <- if (LMRSE < m_sig){
      model
    } else {"Mean"}
    
    #Record values
    fitht[nrow(fitht)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMp, 
                               LM_sig, meanht, m_sig, model)
  }
  
  return(fitht)
}


#' Builds models for top-w height allometry dynamics of surveyed Species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @param bTest Multiples of mean + mRSE for which Burr & quadratic models can predict 
#' beyond the observed mean + standard deviation
#' @param maxiter The maximum number of iterations for model fitting
#' @return dataframe
#' @export

wDyn <- function(dat, width = "width", top = "top", 
                 thres = 5, pnts = 10, p = 0.05, bTest = 10, maxiter = 1000) {

  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"),
           Rat = as.numeric(width)/top)  
  
  # Find missing data
  entries <- which(is.na(dat[width]))
  if (length(entries)>0) {
    cat(" datClean removed these rows as they were missing crown widths", "\n", entries, "\n", "\n")
    dat <- dat[-entries,] 
  }
  
  priorList <- unique(dat$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitw <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_Rsq' = numeric(0), 'lin_p' = numeric(0),
                     'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_Rsq' = numeric(0), 'NE_p' = numeric(0),
                     'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_Rsq' = numeric(0), 'B_p' = numeric(0),
                     'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_Rsq' = numeric(0), 'Bin_p' = numeric(0),
                     'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'q_Rsq' = numeric(0), 'Q_p' = numeric(0), 
                     'Mean' = character(0), 'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=maxiter, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- dat %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$Rat)
    
    if (length(unique(x, incomparables = FALSE))>2 & length(unique(y, incomparables = FALSE))>1) {
      
      #Linear
      if (!berryFunctions::is.error(lm(y ~ x))) {
        LM<-lm(y ~ x)
        LMSum <- base::summary(LM)
        LMa <- LMSum$coefficients[2]
        LMb <- LMSum$coefficients[1]
        LMRSE <- LMSum$sigma
        LMRsq <- LMSum$r.squared
        LMp <- if (LMRSE == 0) {
          0
        } else {
          round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
        }
        rm(LM)
      } else {
        LMa <- NA
        LMb <- NA
        LMRSE <- 100
        LMRsq <- 0
        LMp <- 1
      }
      
      #Negative exponential
      init1<-c(k=50,r=0.5)
      if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T))) {
        NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
        NESum <- base::summary(NE)
        k <- NESum$coefficients[1]
        r <- NESum$coefficients[2]
        NERSE <- NESum$sigma
        NERsq <- cor(predict(NE, newdata=x), y)**2
        NEp <- round(max(NESum$coefficients[7],NESum$coefficients[8]),5)
        rm(NE)
      } else {
        k <- NA
        r <- NA
        NERSE <- 100
        NERsq <- 0
        NEp <- 1
      }
      
      #Burr
      init2<-c(a=3,b=2)
      if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
        Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
        BSum <- base::summary(Burr)
        Ba <- BSum$coefficients[1]
        Bb <- BSum$coefficients[2]
        
        #Added control for Burr
        f <- function(x){Ba*Bb*((0.1*x^(Ba-1))/((1+(0.1*x)^Ba)^Bb+1))}
        if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
          BRSE <- 100
          BRsq <- 0
          Bp <- 1
        } else if (bTest * (mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE)) > (optimize(f = f, interval=c(0, 150), maximum=FALSE))$objective) {
          BRSE <- 100
          BRsq <- 0
          Bp <- 1
        } else {
          BRSE <- BSum$sigma
          BRsq <- cor(predict(Burr, newdata=x), y)**2
          Bp <- round(max(BSum$coefficients[7],BSum$coefficients[8]),5)
        }
        rm(Burr)
      } else {
        Ba <- NA
        Bb <- NA
        BRSE <- 100
        BRsq <- 0
        Bp <- 1
      }
      
      #Binomial
      init3<-c(s=1,sd=3, m=20)
      if (!berryFunctions::is.error(nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control))) {
        Bin<-nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control)
        BinSum <- base::summary(Bin)
        Bs <- BinSum$coefficients[1]
        Bsd <- BinSum$coefficients[2]
        Bm <- BinSum$coefficients[3]
        BinRSE <- BinSum$sigma
        BinRsq <- cor(predict(Bin, newdata=x), y)**2
        Binp <- round(max(BinSum$coefficients[10],BinSum$coefficients[11],BinSum$coefficients[12]),5)
        rm(Bin)
      } else {
        Bs <- NA
        Bsd <- NA
        Bm <- NA
        BinRSE <- 100
        BinRsq <- 0
        Binp <- 1
      }
      
      #Quadratic
      init4 <- c(a = 1, b = 2, c = 0)
      if (!berryFunctions::is.error(nls(y ~ a*x^2 + b*x + c, data = studySpecies, 
                                        start = init4, trace = T, control = control))) {
        q <- nls(y ~ a*x^2 + b*x + c, data = studySpecies, start = init4, 
                 trace = T, control = control)
        qSum <- base::summary(q)
        qa <- qSum$coefficients[1]
        qb <- qSum$coefficients[2]
        qc <- qSum$coefficients[3]
        
        #Added control for Quadratic
        f <- function(x){qa*x^2 + qb*x + qc}
        if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
          qRSE <- 100
          qRsq <- 0
          qp <- 1
        } else {
          qRSE <- qSum$sigma
          qRsq <- cor(predict(q, newdata=x), y)**2
          qp <- round(max(qSum$coefficients[10], qSum$coefficients[11], 
                          qSum$coefficients[12]), 5)
        }
        rm(q)
      } else {
        qa <- NA
        qb <- NA
        qc <- NA
        qRSE <- 100
        qRsq <- 0
        qp <- 1
      }
      
      #Summary stats
      meanw <- round(mean(y, na.rm = TRUE),1)
      m_sig <- round(mRSE(dat = y),3)
      
      #Choose the best model
      listRsq <- c(LMRsq, NERsq, BRsq, BinRsq, qRsq)
      listRSE <- c(LMRSE, NERSE, BRSE, BinRSE, qRSE)
      listMod <- c("Linear", "NegExp", "Burr", "Binomial", "Quadratic")
      
      if (listRSE[which(listRsq == max(listRsq))]<=m_sig) {
        model <- listMod[which(listRsq == max(listRsq))]
      } else {
        model <- "Mean"
      }
      
    } else {
      LMa <- NA
      LMb <- NA
      LMRSE <- 100
      LMRsq <- 0
      LMp <- 1
      k <- NA
      r <- NA
      NERSE <- 100
      NERsq <- 0
      NEp <- 1
      Ba <- NA
      Bb <- NA
      BRSE <- 100
      BRsq <- 0
      Bp <- 1
      Bs <- NA
      Bsd <- NA
      Bm <- NA
      BinRSE <- 100
      BinRsq <- 0
      Binp <- 1
      qa <- NA
      qb <- NA
      qc <- NA
      qRSE <- 100
      qRsq <- 0
      qp <- 1
      model <- "Mean"      
      #Summary stats
      meanw <- round(mean(y, na.rm = TRUE),1)
      m_sig <- round(mRSE(dat = y),3)
    }
    
    #Record values
    fitw[nrow(fitw)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMRsq, LMp, k, r, NERSE, NERsq, NEp,
                             Ba, Bb, BRSE, BRsq, Bp, Bs, Bsd, Bm, BinRSE, BinRsq, Binp, qa, qb, qc, qRSE, qRsq, qp, 
                             meanw, m_sig, model)
  }
  
  return(fitw)
}


#' Builds the plant growth table used in dynamic modelling
#' 
#' Models plant growth from field data in a series of models,
#' selecting the best model and providing error statistics.
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @param bTest 
#' @param cTest 
#' @param Sr Rate of increase for surface litter in a negative exponential curve
#' @param Sk Asymptote for surface litter in a negative exponential curve
#' @param Sa 
#' @param Sb 
#' @param Sc 
#' @param NSr Rate of increase for NS fuels in a negative exponential curve
#' @param NSk Asymptote for NS fuels in a negative exponential curve
#' @param NSa 
#' @param NSb 
#' @param NSc 
#' @param maxiter The maximum number of iterations for model fitting
#'
#' @return dataframe
#' @export

floraDynamics <- function(dat, thres = 5, pnts = 10, p = 0.01, bTest  = 2, cTest  = 10, maxiter = 1000,
                          Sr = 0, Sk = 0, Sa = 0, Sb = 0, Sc = 0, 
                          NSr = 0, NSk = 0, NSa = 0, NSb = 0, NSc = 0){
  
  coverChange <- coverDyn(dat, thres = thres, pnts = pnts, p = p, bTest  = cTest, maxiter = maxiter)
  topChange <- topDyn(dat, thres = thres, pnts = pnts, p = p, bTest  = bTest, maxiter = maxiter)
  baseChange <- baseDyn(dat, thres = thres, pnts = pnts, p = p, bTest  = bTest, maxiter = maxiter)
  he_Change <- heDyn(dat, thres = thres, pnts = pnts, p = p)
  ht_Change <- htDyn(dat, thres = thres, pnts = pnts, p = p)
  w_Change <- wDyn(dat, thres = thres, pnts = pnts, p = p, bTest  = bTest, maxiter = maxiter)
  
  # Collect models
  # Cover
  a <- rep(NA, length(coverChange$Species))
  for (sp in 1:length(coverChange$Species)) {
    if (coverChange$Model[sp] == 'Linear') {
      a[sp] <- coverChange$lin_a[sp]
    } else
      if (coverChange$Model[sp] == 'NegExp') {
        a[sp] <- coverChange$k[sp]
      } else
        if (coverChange$Model[sp] == 'Burr') {
          a[sp] <- coverChange$Ba[sp]
        } else
          if (coverChange$Model[sp] == 'Binomial') {
            a[sp] <- coverChange$scale[sp]
          } else
            if (coverChange$Model[sp] == 'Quadratic') {
              a[sp] <- coverChange$Qa[sp]
            } else {a[sp] <- 0}
  }
  b <- rep(NA, length(coverChange$Species))
  for (sp in 1:length(coverChange$Species)) {
    if (coverChange$Model[sp] == 'Linear') {
      b[sp] <- coverChange$lin_b[sp]
    } else
      if (coverChange$Model[sp] == 'NegExp') {
        b[sp] <- coverChange$r[sp]
      } else
        if (coverChange$Model[sp] == 'Burr') {
          b[sp] <- coverChange$Bb[sp]
        } else
          if (coverChange$Model[sp] == 'Binomial') {
            b[sp] <- coverChange$sd[sp]
          } else
            if (coverChange$Model[sp] == 'Quadratic') {
              b[sp] <- coverChange$Qb[sp]
            } else {b[sp] <- coverChange$Mean[sp]}
  }
  c <- rep(NA, length(coverChange$Species))
  for (sp in 1:length(coverChange$Species)) {
    if (coverChange$Model[sp] == 'Linear') {
      c[sp] <- NA
    } else
      if (coverChange$Model[sp] == 'NegExp') {
        c[sp] <- NA
      } else
        if (coverChange$Model[sp] == 'Burr') {
          c[sp] <- NA
        } else
          if (coverChange$Model[sp] == 'Binomial') {
            c[sp] <- coverChange$Binm[sp]
          } else
            if (coverChange$Model[sp] == 'Quadratic') {
              c[sp] <- coverChange$Qc[sp]
            } else {c[sp] <- NA}
  }
  RSE <- rep(NA, length(coverChange$Species))
  for (sp in 1:length(coverChange$Species)) {
    if (coverChange$Model[sp] == 'Linear') {
      RSE[sp] <- coverChange$lin_Sigma[sp]
    } else
      if (coverChange$Model[sp] == 'NegExp') {
        RSE[sp] <- coverChange$NE_sigma[sp]
      } else
        if (coverChange$Model[sp] == 'Burr') {
          RSE[sp] <- coverChange$B_sigma[sp]
        } else
          if (coverChange$Model[sp] == 'Binomial') {
            RSE[sp] <- coverChange$Bin_sigma[sp]
          } else
            if (coverChange$Model[sp] == 'Quadratic') {
              RSE[sp] <- coverChange$Q_sigma[sp]
            } else {RSE[sp] <- coverChange$Mean_sigma[sp]}
  }
  Rsq <- rep(NA, length(coverChange$Species))
  for (sp in 1:length(coverChange$Species)) {
    if (coverChange$Model[sp] == 'Linear') {
      Rsq[sp] <- coverChange$lin_Rsq[sp]
    } else
      if (coverChange$Model[sp] == 'NegExp') {
        Rsq[sp] <- coverChange$NE_Rsq[sp]
      } else
        if (coverChange$Model[sp] == 'Burr') {
          Rsq[sp] <- coverChange$B_Rsq[sp]
        } else
          if (coverChange$Model[sp] == 'Binomial') {
            Rsq[sp] <- coverChange$Bin_Rsq[sp]
          } else
            if (coverChange$Model[sp] == 'Quadratic') {
              Rsq[sp] <- coverChange$q_Rsq[sp]
            } else {Rsq[sp] <- NA}
  }
  CovMean <- rep(NA, length(coverChange$Species))
  for (sp in 1:length(coverChange$Species)) {
    CovMean[sp] <- coverChange$Mean[sp]
  }
  
  # Height
  aa <- rep(NA, length(topChange$Species))
  for (sp in 1:length(topChange$Species)) {
    if (topChange$Model[sp] == 'Linear') {
      aa[sp] <- topChange$lin_a[sp]
    } else
      if (topChange$Model[sp] == 'NegExp') {
        aa[sp] <- topChange$k[sp]
      } else
        if (topChange$Model[sp] == 'Burr') {
          aa[sp] <- topChange$Ba[sp]
        } else
          if (topChange$Model[sp] == 'Binomial') {
            aa[sp] <- topChange$scale[sp]
          } else
            if (topChange$Model[sp] == 'Quadratic') {
              aa[sp] <- topChange$Qa[sp]
            } else {aa[sp] <- 0}
  }
  bb <- rep(NA, length(topChange$Species))
  for (sp in 1:length(topChange$Species)) {
    if (topChange$Model[sp] == 'Linear') {
      bb[sp] <- topChange$lin_b[sp]
    } else
      if (topChange$Model[sp] == 'NegExp') {
        bb[sp] <- topChange$r[sp]
      } else
        if (topChange$Model[sp] == 'Burr') {
          bb[sp] <- topChange$Bb[sp]
        } else
          if (topChange$Model[sp] == 'Binomial') {
            bb[sp] <- topChange$sd[sp]
          } else
            if (topChange$Model[sp] == 'Quadratic') {
              bb[sp] <- topChange$Qb[sp]
            } else {bb[sp] <- topChange$Mean[sp]}
  }
  cc <- rep(NA, length(topChange$Species))
  for (sp in 1:length(topChange$Species)) {
    if (topChange$Model[sp] == 'Linear') {
      cc[sp] <- NA
    } else
      if (topChange$Model[sp] == 'NegExp') {
        cc[sp] <- NA
      } else
        if (topChange$Model[sp] == 'Burr') {
          cc[sp] <- NA
        } else
          if (topChange$Model[sp] == 'Binomial') {
            cc[sp] <- topChange$Binm[sp]
          } else
            if (topChange$Model[sp] == 'Quadratic') {
              cc[sp] <- topChange$Qc[sp]
            } else {cc[sp] <- NA}
  }
  RSE1 <- rep(NA, length(topChange$Species))
  for (sp in 1:length(topChange$Species)) {
    if (topChange$Model[sp] == 'Linear') {
      RSE1[sp] <- topChange$lin_Sigma[sp]
    } else
      if (topChange$Model[sp] == 'NegExp') {
        RSE1[sp] <- topChange$NE_sigma[sp]
      } else
        if (topChange$Model[sp] == 'Burr') {
          RSE1[sp] <- topChange$B_sigma[sp]
        } else
          if (topChange$Model[sp] == 'Binomial') {
            RSE1[sp] <- topChange$Bin_sigma[sp]
          } else
            if (topChange$Model[sp] == 'Quadratic') {
              RSE1[sp] <- topChange$Q_sigma[sp]
            } else {RSE1[sp] <- topChange$Mean_sigma[sp]}
  }
  Rsq1 <- rep(NA, length(topChange$Species))
  for (sp in 1:length(topChange$Species)) {
    if (topChange$Model[sp] == 'Linear') {
      Rsq1[sp] <- topChange$lin_Rsq[sp]
    } else
      if (topChange$Model[sp] == 'NegExp') {
        Rsq1[sp] <- topChange$NE_Rsq[sp]
      } else
        if (topChange$Model[sp] == 'Burr') {
          Rsq1[sp] <- topChange$B_Rsq[sp]
        } else
          if (topChange$Model[sp] == 'Binomial') {
            Rsq1[sp] <- topChange$Bin_Rsq[sp]
          } else
            if (topChange$Model[sp] == 'Quadratic') {
              Rsq1[sp] <- topChange$q_Rsq[sp]
            } else {Rsq1[sp] <- NA}
  }
  hMean <- rep(NA, length(topChange$Species))
  for (sp in 1:length(topChange$Species)) {
    hMean[sp] <- topChange$Mean[sp]
  }
  
  # Base
  aaa <- rep(NA, length(baseChange$Species))
  for (sp in 1:length(baseChange$Species)) {
    if (baseChange$Model[sp] == 'Linear') {
      aaa[sp] <- baseChange$lin_a[sp]
    } else
      if (baseChange$Model[sp] == 'NegExp') {
        aaa[sp] <- baseChange$k[sp]
      } else
        if (baseChange$Model[sp] == 'Burr') {
          aaa[sp] <- baseChange$Ba[sp]
        } else
          if (baseChange$Model[sp] == 'Binomial') {
            aaa[sp] <- baseChange$scale[sp]
          } else
            if (baseChange$Model[sp] == 'Quadratic') {
              aaa[sp] <- baseChange$Qa[sp]
            } else {aaa[sp] <- 0}
  }
  bbb <- rep(NA, length(baseChange$Species))
  for (sp in 1:length(baseChange$Species)) {
    if (baseChange$Model[sp] == 'Linear') {
      bbb[sp] <- baseChange$lin_b[sp]
    } else
      if (baseChange$Model[sp] == 'NegExp') {
        bbb[sp] <- baseChange$r[sp]
      } else
        if (baseChange$Model[sp] == 'Burr') {
          bbb[sp] <- baseChange$Bb[sp]
        } else
          if (baseChange$Model[sp] == 'Binomial') {
            bbb[sp] <- baseChange$sd[sp]
          } else
            if (baseChange$Model[sp] == 'Quadratic') {
              bbb[sp] <- baseChange$Qb[sp]
            } else {bbb[sp] <- baseChange$Mean[sp]}
  }
  ccc <- rep(NA, length(baseChange$Species))
  for (sp in 1:length(baseChange$Species)) {
    if (baseChange$Model[sp] == 'Linear') {
      ccc[sp] <- NA
    } else
      if (baseChange$Model[sp] == 'NegExp') {
        ccc[sp] <- NA
      } else
        if (baseChange$Model[sp] == 'Burr') {
          ccc[sp] <- NA
        } else
          if (baseChange$Model[sp] == 'Binomial') {
            ccc[sp] <- baseChange$Binm[sp]
          } else
            if (baseChange$Model[sp] == 'Quadratic') {
              ccc[sp] <- baseChange$Qc[sp]
            } else {ccc[sp] <- NA}
  }
  RSE2 <- rep(NA, length(baseChange$Species))
  for (sp in 1:length(baseChange$Species)) {
    if (baseChange$Model[sp] == 'Linear') {
      RSE2[sp] <- baseChange$lin_Sigma[sp]
    } else
      if (baseChange$Model[sp] == 'NegExp') {
        RSE2[sp] <- baseChange$NE_sigma[sp]
      } else
        if (baseChange$Model[sp] == 'Burr') {
          RSE2[sp] <- baseChange$B_sigma[sp]
        } else
          if (baseChange$Model[sp] == 'Binomial') {
            RSE2[sp] <- baseChange$Bin_sigma[sp]
          } else
            if (baseChange$Model[sp] == 'Quadratic') {
              RSE2[sp] <- baseChange$Q_sigma[sp]
            } else {RSE2[sp] <- baseChange$Mean_sigma[sp]}
  }
  Rsq2 <- rep(NA, length(baseChange$Species))
  for (sp in 1:length(baseChange$Species)) {
    if (baseChange$Model[sp] == 'Linear') {
      Rsq2[sp] <- baseChange$lin_Rsq[sp]
    } else
      if (baseChange$Model[sp] == 'NegExp') {
        Rsq2[sp] <- baseChange$NE_Rsq[sp]
      } else
        if (baseChange$Model[sp] == 'Burr') {
          Rsq2[sp] <- baseChange$B_Rsq[sp]
        } else
          if (baseChange$Model[sp] == 'Binomial') {
            Rsq2[sp] <- baseChange$Bin_Rsq[sp]
          } else
            if (baseChange$Model[sp] == 'Quadratic') {
              Rsq2[sp] <- baseChange$q_Rsq[sp]
            } else {Rsq2[sp] <- NA}
  }
  bMean <- rep(NA, length(baseChange$Species))
  for (sp in 1:length(baseChange$Species)) {
    bMean[sp] <- baseChange$Mean[sp]
  }
  
  # He
  aA <- rep(NA, length(he_Change$Species))
  for (sp in 1:length(he_Change$Species)) {
    if (he_Change$Model[sp] == 'Linear') {
      aA[sp] <- he_Change$lin_a[sp]
    } else
      if (he_Change$Model[sp] == 'NegExp') {
        aA[sp] <- he_Change$k[sp]
      } else
        if (he_Change$Model[sp] == 'Burr') {
          aA[sp] <- he_Change$Ba[sp]
        } else
          if (he_Change$Model[sp] == 'Binomial') {
            aA[sp] <- he_Change$scale[sp]
          } else
            if (he_Change$Model[sp] == 'Quadratic') {
              aA[sp] <- he_Change$Qa[sp]
            } else {aA[sp] <- 0}
  }
  bB <- rep(NA, length(he_Change$Species))
  for (sp in 1:length(he_Change$Species)) {
    if (he_Change$Model[sp] == 'Linear') {
      bB[sp] <- he_Change$lin_b[sp]
    } else
      if (he_Change$Model[sp] == 'NegExp') {
        bB[sp] <- he_Change$r[sp]
      } else
        if (he_Change$Model[sp] == 'Burr') {
          bB[sp] <- he_Change$Bb[sp]
        } else
          if (he_Change$Model[sp] == 'Binomial') {
            bB[sp] <- he_Change$sd[sp]
          } else
            if (he_Change$Model[sp] == 'Quadratic') {
              bB[sp] <- he_Change$Qb[sp]
            } else {bB[sp] <- he_Change$Mean[sp]}
  }
  cC <- rep(NA, length(he_Change$Species))
  for (sp in 1:length(he_Change$Species)) {
    if (he_Change$Model[sp] == 'Linear') {
      cC[sp] <- NA
    } else
      if (he_Change$Model[sp] == 'NegExp') {
        cC[sp] <- NA
      } else
        if (he_Change$Model[sp] == 'Burr') {
          cC[sp] <- NA
        } else
          if (he_Change$Model[sp] == 'Binomial') {
            cC[sp] <- he_Change$Binm[sp]
          } else
            if (he_Change$Model[sp] == 'Quadratic') {
              cC[sp] <- he_Change$Qc[sp]
            } else {cC[sp] <- NA}
  }
  RSE3 <- rep(NA, length(he_Change$Species))
  for (sp in 1:length(he_Change$Species)) {
    if (he_Change$Model[sp] == 'Linear') {
      RSE3[sp] <- he_Change$lin_Sigma[sp]
    } else
      if (he_Change$Model[sp] == 'NegExp') {
        RSE3[sp] <- he_Change$NE_sigma[sp]
      } else
        if (he_Change$Model[sp] == 'Burr') {
          RSE3[sp] <- he_Change$B_sigma[sp]
        } else
          if (he_Change$Model[sp] == 'Binomial') {
            RSE3[sp] <- he_Change$Bin_sigma[sp]
          } else
            if (he_Change$Model[sp] == 'Quadratic') {
              RSE3[sp] <- he_Change$Q_sigma[sp]
            } else {RSE3[sp] <- he_Change$Mean_sigma[sp]}
  }
  heMean <- rep(NA, length(he_Change$Species))
  for (sp in 1:length(he_Change$Species)) {
    heMean[sp] <- he_Change$Mean[sp]
  }
  
  # Ht
  aT <- rep(NA, length(ht_Change$Species))
  for (sp in 1:length(ht_Change$Species)) {
    if (ht_Change$Model[sp] == 'Linear') {
      aT[sp] <- ht_Change$lin_a[sp]
    } else
      if (ht_Change$Model[sp] == 'NegExp') {
        aT[sp] <- ht_Change$k[sp]
      } else
        if (ht_Change$Model[sp] == 'Burr') {
          aT[sp] <- ht_Change$Ba[sp]
        } else
          if (ht_Change$Model[sp] == 'Binomial') {
            aT[sp] <- ht_Change$scale[sp]
          } else
            if (ht_Change$Model[sp] == 'Quadratic') {
              aT[sp] <- ht_Change$Qa[sp]
            } else {aT[sp] <- 0}
  }
  bT <- rep(NA, length(ht_Change$Species))
  for (sp in 1:length(ht_Change$Species)) {
    if (ht_Change$Model[sp] == 'Linear') {
      bT[sp] <- ht_Change$lin_b[sp]
    } else
      if (ht_Change$Model[sp] == 'NegExp') {
        bT[sp] <- ht_Change$r[sp]
      } else
        if (ht_Change$Model[sp] == 'Burr') {
          bT[sp] <- ht_Change$Bb[sp]
        } else
          if (ht_Change$Model[sp] == 'Binomial') {
            bT[sp] <- ht_Change$sd[sp]
          } else
            if (ht_Change$Model[sp] == 'Quadratic') {
              bT[sp] <- ht_Change$Qb[sp]
            } else {bT[sp] <- ht_Change$Mean[sp]}
  }
  cT <- rep(NA, length(ht_Change$Species))
  for (sp in 1:length(ht_Change$Species)) {
    if (ht_Change$Model[sp] == 'Linear') {
      cT[sp] <- NA
    } else
      if (ht_Change$Model[sp] == 'NegExp') {
        cT[sp] <- NA
      } else
        if (ht_Change$Model[sp] == 'Burr') {
          cT[sp] <- NA
        } else
          if (ht_Change$Model[sp] == 'Binomial') {
            cT[sp] <- ht_Change$Binm[sp]
          } else
            if (ht_Change$Model[sp] == 'Quadratic') {
              cT[sp] <- ht_Change$Qc[sp]
            } else {cT[sp] <- NA}
  }
  RSEt <- rep(NA, length(ht_Change$Species))
  for (sp in 1:length(ht_Change$Species)) {
    if (ht_Change$Model[sp] == 'Linear') {
      RSEt[sp] <- ht_Change$lin_Sigma[sp]
    } else
      if (ht_Change$Model[sp] == 'NegExp') {
        RSEt[sp] <- ht_Change$NE_sigma[sp]
      } else
        if (ht_Change$Model[sp] == 'Burr') {
          RSEt[sp] <- ht_Change$B_sigma[sp]
        } else
          if (ht_Change$Model[sp] == 'Binomial') {
            RSEt[sp] <- ht_Change$Bin_sigma[sp]
          } else
            if (ht_Change$Model[sp] == 'Quadratic') {
              RSEt[sp] <- ht_Change$Q_sigma[sp]
            } else {RSEt[sp] <- ht_Change$Mean_sigma[sp]}
  }
  htMean <- rep(NA, length(ht_Change$Species))
  for (sp in 1:length(ht_Change$Species)) {
    htMean[sp] <- ht_Change$Mean[sp]
  }
  
  # Width
  aw <- rep(NA, length(w_Change$Species))
  for (sp in 1:length(w_Change$Species)) {
    if (w_Change$Model[sp] == 'Linear') {
      aw[sp] <- w_Change$lin_a[sp]
    } else
      if (w_Change$Model[sp] == 'NegExp') {
        aw[sp] <- w_Change$k[sp]
      } else
        if (w_Change$Model[sp] == 'Burr') {
          aw[sp] <- w_Change$Ba[sp]
        } else
          if (w_Change$Model[sp] == 'Binomial') {
            aw[sp] <- w_Change$scale[sp]
          } else
            if (w_Change$Model[sp] == 'Quadratic') {
              aw[sp] <- w_Change$Qa[sp]
            } else {aw[sp] <- 0}
  }
  bw <- rep(NA, length(w_Change$Species))
  for (sp in 1:length(w_Change$Species)) {
    if (w_Change$Model[sp] == 'Linear') {
      bw[sp] <- w_Change$lin_b[sp]
    } else
      if (w_Change$Model[sp] == 'NegExp') {
        bw[sp] <- w_Change$r[sp]
      } else
        if (w_Change$Model[sp] == 'Burr') {
          bw[sp] <- w_Change$Bb[sp]
        } else
          if (w_Change$Model[sp] == 'Binomial') {
            bw[sp] <- w_Change$sd[sp]
          } else
            if (w_Change$Model[sp] == 'Quadratic') {
              bw[sp] <- w_Change$Qb[sp]
            } else {bw[sp] <- w_Change$Mean[sp]}
  }
  cw <- rep(NA, length(w_Change$Species))
  for (sp in 1:length(w_Change$Species)) {
    if (w_Change$Model[sp] == 'Linear') {
      cw[sp] <- NA
    } else
      if (w_Change$Model[sp] == 'NegExp') {
        cw[sp] <- NA
      } else
        if (w_Change$Model[sp] == 'Burr') {
          cw[sp] <- NA
        } else
          if (w_Change$Model[sp] == 'Binomial') {
            cw[sp] <- w_Change$Binm[sp]
          } else
            if (w_Change$Model[sp] == 'Quadratic') {
              cw[sp] <- w_Change$Qc[sp]
            } else {cw[sp] <- NA}
  }
  RSEw <- rep(NA, length(w_Change$Species))
  for (sp in 1:length(w_Change$Species)) {
    if (w_Change$Model[sp] == 'Linear') {
      RSEw[sp] <- w_Change$lin_Sigma[sp]
    } else
      if (w_Change$Model[sp] == 'NegExp') {
        RSEw[sp] <- w_Change$NE_sigma[sp]
      } else
        if (w_Change$Model[sp] == 'Burr') {
          RSEw[sp] <- w_Change$B_sigma[sp]
        } else
          if (w_Change$Model[sp] == 'Binomial') {
            RSEw[sp] <- w_Change$Bin_sigma[sp]
          } else
            if (w_Change$Model[sp] == 'Quadratic') {
              RSEw[sp] <- w_Change$Q_sigma[sp]
            } else {RSEw[sp] <- w_Change$Mean_sigma[sp]}
  }
  Rsqw <- rep(NA, length(w_Change$Species))
  for (sp in 1:length(w_Change$Species)) {
    if (w_Change$Model[sp] == 'Linear') {
      Rsqw[sp] <- w_Change$lin_Rsq[sp]
    } else
      if (w_Change$Model[sp] == 'NegExp') {
        Rsqw[sp] <- w_Change$NE_Rsq[sp]
      } else
        if (w_Change$Model[sp] == 'Burr') {
          Rsqw[sp] <- w_Change$B_Rsq[sp]
        } else
          if (w_Change$Model[sp] == 'Binomial') {
            Rsqw[sp] <- w_Change$Bin_Rsq[sp]
          } else
            if (w_Change$Model[sp] == 'Quadratic') {
              Rsqw[sp] <- w_Change$q_Rsq[sp]
            } else {Rsqw[sp] <- NA}
  }
  wMean <- rep(NA, length(w_Change$Species))
  for (sp in 1:length(w_Change$Species)) {
    wMean[sp] <- w_Change$Mean[sp]
  }
  
  models <- data.frame('Species'=coverChange$Species, 
                       'Cover' = coverChange$Model, 'C_a' = a, 'C_b' = b, 'C_c' = c, 'C_RSE' = RSE, 'C_Rsq' = Rsq, 'meanCover' = CovMean,
                       'Top' = topChange$Model, 'T_a' = aa, 'T_b' = bb, 'T_c' = cc, 'T_RSE' = RSE1, 'T_Rsq' = Rsq1, 'meanTop' = hMean,
                       'Base' = baseChange$Model, 'B_a' = aaa, 'B_b' = bbb, 'B_c' = ccc, 'B_RSE' = RSE2, 'B_Rsq' = Rsq2, 'meanBase' = bMean,
                       'He' = he_Change$Model, 'He_a' = aA, 'He_b' = bB, 'He_c' = cC, 'He_RSE' = RSE3, 'meanHe' = heMean,
                       'Ht' = ht_Change$Model, 'Ht_a' = aT, 'Ht_b' = bT, 'Ht_c' = cT, 'Ht_RSE' = RSEt, 'meanHt' = htMean,
                       'Width' = w_Change$Model, 'w_a' = aw, 'w_b' = bw, 'w_c' = cw, 'w_RSE' = RSEw, 'w_Rsq' = Rsqw, 'meanW' = wMean)
  # Dead biomass
  litter <- case_when(
    Sr != 0  ~ "NegExp",
    Sa != 0 ~ "Quadratic",
    TRUE ~ as.character("Set")
  ) 
  suspNS <- case_when(
    NSr != 0  ~ "NegExp",
    NSa != 0 ~ "Quadratic",
    TRUE ~ as.character("Set")
  )
  models[length(models$Species)+1,1] <- "Litter"
  models$Dead[length(models$Species)] <- litter
  models$d_Max[length(models$Species)] <- Sk
  models$d_Rate[length(models$Species)] <- Sr
  models$d_a[length(models$Species)] <- Sa
  models$d_b[length(models$Species)] <- Sb
  models$d_c[length(models$Species)] <- Sc
  models[length(models$Species)+1,1] <- "suspNS"
  models$Dead[length(models$Species)] <- suspNS
  models$d_Max[length(models$Species)] <- NSk
  models$d_Rate[length(models$Species)] <- NSr
  models$d_a[length(models$Species)] <- NSa
  models$d_b[length(models$Species)] <- NSb
  models$d_c[length(models$Species)] <- NSc
  return(models)
}


#' Predicts proportion cover at a given age
#' 
#' Selects model from tabled models per species
#' @param mods A table of models fit to species, format as per modCollector
#' @param sp The name of a species in the table mods
#' @param Age The age of the plant (years since fire)
#' @return value
#' @export
#' 
pCover <- function(mods, sp, Age = 10){
  n <- as.numeric(which(mods$Species == sp))
  if (length(n)>1) {
    cat("There is more than one entry for a species", "\n")
  }
  if (mods$Cover[n] == "NegExp") {
    c <- as.numeric(mods$C_a[n]) * (1 - exp(-as.numeric(mods$C_b[n]) * Age))
  } else if (mods$Cover[n] == "Burr") {
    c <- as.numeric(mods$C_a[n]) * as.numeric(mods$C_b[n]) * ((0.1 * Age ^ (as.numeric(mods$C_a[n]) - 1)) / ((1 + (0.1 * Age) ^ as.numeric(mods$C_a[n])) ^ as.numeric(mods$C_b[n]) + 1))
  } else if (mods$Cover[n] == "Linear") {
    c <- as.numeric(mods$C_a[n]) * Age + as.numeric(mods$C_b[n])
  } else if (mods$Cover[n] == "Binomial") {
    c <- as.numeric(mods$C_a[n])*(1/(as.numeric(mods$C_b[n])*sqrt(2*pi)))*exp(-((Age-as.numeric(mods$C_c[n]))^2)/(2*as.numeric(mods$C_b[n])^2))
  } else if (mods$Cover[n] == "Quadratic") {
    c <- as.numeric(mods$C_a[n])*Age^2 + as.numeric(mods$C_b[n])*Age + as.numeric(mods$C_c[n])
  } else if (mods$Cover[n] == "Mean") {
    c <- as.numeric(mods$meanCover[n])
  }
  c <- min(max(0,c),100)
  return(c)
}

#' Predicts plant top height at a given age
#' 
#' Selects model from tabled models per species
#' @param mods A table of models fit to species, format as per modCollector
#' @param sp The name of a species in the table mods
#' @param Age The age of the plant (years since fire)
#' @return value
#' @export
#' 
pTop <- function(mods, sp, Age = 10){
  n <- as.numeric(which(mods$Species == sp))
  if (length(n)>1) {
    cat("There is more than one entry for a species", "\n")
  }
  if (mods$Top[n] == "NegExp") {
    c <- as.numeric(mods$T_a[n]) * (1 - exp(-as.numeric(mods$T_b[n]) * Age))
  } else if (mods$Top[n] == "Burr") {
    c <- as.numeric(mods$T_a[n]) * as.numeric(mods$T_b[n]) * ((0.1 * Age ^ (as.numeric(mods$T_a[n]) - 1)) / ((1 + (0.1 * Age) ^ as.numeric(mods$T_a[n])) ^ as.numeric(mods$T_b[n]) + 1))
  } else if (mods$Top[n] == "Linear") {
    c <- as.numeric(mods$T_a[n]) * Age + as.numeric(mods$T_b[n])
  } else if (mods$Top[n] == "Binomial") {
    c <- as.numeric(mods$T_a[n])*(1/(as.numeric(mods$T_b[n])*sqrt(2*pi)))*exp(-((Age-as.numeric(mods$T_c[n]))^2)/(2*as.numeric(mods$T_b[n])^2))
  } else if (mods$Top[n] == "Quadratic") {
    c <- as.numeric(mods$T_a[n])*Age^2 + as.numeric(mods$T_b[n])*Age + as.numeric(mods$T_c[n])
  } else if (mods$Top[n] == "Mean") {
    c <- as.numeric(mods$meanTop[n])
  }
  c <- max(0,c)
  return(c)
}

#' Predicts plant base height at a given age
#' 
#' Selects model from tabled models per species
#' @param mods A table of models fit to species, format as per modCollector
#' @param sp The name of a species in the table mods
#' @param Age The age of the plant (years since fire)
#' @return value
#' @export
#' 
pBase <- function(mods, sp, Age = 10){
  n <- as.numeric(which(mods$Species == sp))
  if (length(n)>1) {
    cat("There is more than one entry for a species", "\n")
  }
  if (mods$Base[n] == "NegExp") {
    c <- as.numeric(mods$B_a[n]) * (1 - exp(-as.numeric(mods$B_b[n]) * Age))
  } else if (mods$Base[n] == "Burr") {
    c <- as.numeric(mods$B_a[n]) * as.numeric(mods$B_b[n]) * ((0.1 * Age ^ (as.numeric(mods$B_a[n]) - 1)) / ((1 + (0.1 * Age) ^ as.numeric(mods$B_a[n])) ^ as.numeric(mods$B_b[n]) + 1))
  } else if (mods$Base[n] == "Linear") {
    c <- as.numeric(mods$B_a[n]) * Age + as.numeric(mods$B_b[n])
  } else if (mods$Base[n] == "Binomial") {
    c <- as.numeric(mods$B_a[n])*(1/(as.numeric(mods$B_b[n])*sqrt(2*pi)))*exp(-((Age-as.numeric(mods$B_c[n]))^2)/(2*as.numeric(mods$B_b[n])^2))
  } else if (mods$Base[n] == "Quadratic") {
    c <- as.numeric(mods$B_a[n])*Age^2 + as.numeric(mods$B_b[n])*Age + as.numeric(mods$B_c[n])
  } else if (mods$Base[n] == "Mean") {
    c <- as.numeric(mods$meanBase[n])
  }
  c <- min(max(0,c),1)
  return(c)
}

#' Predicts plant He at a given age
#' 
#' Selects model from tabled models per species
#' @param mods A table of models fit to species, format as per modCollector
#' @param sp The name of a species in the table mods
#' @param Age The age of the plant (years since fire)
#' @return value
#' @export
#' 
pHe <- function(mods, sp, Age = 10){
  n <- as.numeric(which(mods$Species == sp))
  if (length(n)>1) {
    cat("There is more than one entry for a species", "\n")
  }
  if (mods$He[n] == "NegExp") {
    c <- as.numeric(mods$He_a[n]) * (1 - exp(-as.numeric(mods$He_b[n]) * Age))
  } else if (mods$He[n] == "Burr") {
    c <- as.numeric(mods$He_a[n]) * as.numeric(mods$He_b[n]) * ((0.1 * Age ^ (as.numeric(mods$He_a[n]) - 1)) / ((1 + (0.1 * Age) ^ as.numeric(mods$He_a[n])) ^ as.numeric(mods$He_b[n]) + 1))
  } else if (mods$He[n] == "Linear") {
    c <- as.numeric(mods$He_a[n]) * Age + as.numeric(mods$He_b[n])
  } else if (mods$He[n] == "Binomial") {
    c <- as.numeric(mods$He_a[n])*(1/(as.numeric(mods$He_b[n])*sqrt(2*pi)))*exp(-((Age-as.numeric(mods$He_c[n]))^2)/(2*as.numeric(mods$He_b[n])^2))
  } else if (mods$He[n] == "Quadratic") {
    c <- as.numeric(mods$He_a[n])*Age^2 + as.numeric(mods$He_b[n])*Age + as.numeric(mods$He_c[n])
  } else if (mods$He[n] == "Mean") {
    c <- as.numeric(mods$meanHe[n])
  }
  c <- max(0,c)
  return(c)
}

#' Predicts plant Ht at a given age
#' 
#' Selects model from tabled models per species
#' @param mods A table of models fit to species, format as per modCollector
#' @param sp The name of a species in the table mods
#' @param Age The age of the plant (years since fire)
#' @return value
#' @export
#' 
pHt <- function(mods, sp, Age = 10){
  n <- as.numeric(which(mods$Species == sp))
  if (length(n)>1) {
    cat("There is more than one entry for a species", "\n")
  }
  if (mods$Ht[n] == "NegExp") {
    c <- as.numeric(mods$Ht_a[n]) * (1 - exp(-as.numeric(mods$Ht_b[n]) * Age))
  } else if (mods$Ht[n] == "Burr") {
    c <- as.numeric(mods$Ht_a[n]) * as.numeric(mods$Ht_b[n]) * ((0.1 * Age ^ (as.numeric(mods$Ht_a[n]) - 1)) / ((1 + (0.1 * Age) ^ as.numeric(mods$Ht_a[n])) ^ as.numeric(mods$Ht_b[n]) + 1))
  } else if (mods$Ht[n] == "Linear") {
    c <- as.numeric(mods$Ht_a[n]) * Age + as.numeric(mods$Ht_b[n])
  } else if (mods$Ht[n] == "Binomial") {
    c <- as.numeric(mods$Ht_a[n])*(1/(as.numeric(mods$Ht_b[n])*sqrt(2*pi)))*exp(-((Age-as.numeric(mods$Ht_c[n]))^2)/(2*as.numeric(mods$Ht_b[n])^2))
  } else if (mods$Ht[n] == "Quadratic") {
    c <- as.numeric(mods$Ht_a[n])*Age^2 + as.numeric(mods$Ht_b[n])*Age + as.numeric(mods$Ht_c[n])
  } else if (mods$Ht[n] == "Mean") {
    c <- as.numeric(mods$meanHt[n])
  }
  c <- max(0,c)
  return(c)
}

#' Predicts plant crown width at a given age
#' 
#' Selects model from tabled models per species
#' @param mods A table of models fit to species, format as per modCollector
#' @param sp The name of a species in the table mods
#' @param Age The age of the plant (years since fire)
#' @return value
#' @export
#' 
pWidth <- function(mods, sp, Age = 10){
  n <- as.numeric(which(mods$Species == sp))
  if (length(n)>1) {
    cat("There is more than one entry for a species", "\n")
  }
  if (mods$Width[n] == "NegExp") {
    c <- as.numeric(mods$w_a[n]) * (1 - exp(-as.numeric(mods$w_b[n]) * Age))
  } else if (mods$Width[n] == "Burr") {
    c <- as.numeric(mods$w_a[n]) * as.numeric(mods$w_b[n]) * ((0.1 * Age ^ (as.numeric(mods$w_a[n]) - 1)) / ((1 + (0.1 * Age) ^ as.numeric(mods$w_a[n])) ^ as.numeric(mods$w_b[n]) + 1))
  } else if (mods$Width[n] == "Linear") {
    c <- as.numeric(mods$w_a[n]) * Age + as.numeric(mods$w_b[n])
  } else if (mods$Width[n] == "Binomial") {
    c <- as.numeric(mods$w_a[n])*(1/(as.numeric(mods$w_b[n])*sqrt(2*pi)))*exp(-((Age-as.numeric(mods$w_c[n]))^2)/(2*as.numeric(mods$w_b[n])^2))
  } else if (mods$Width[n] == "Quadratic") {
    c <- as.numeric(mods$w_a[n])*Age^2 + as.numeric(mods$w_b[n])*Age + as.numeric(mods$w_c[n])
  } else if (mods$Width[n] == "Mean") {
    c <- as.numeric(mods$meanW[n])
  }
  c <- max(0,c)
  return(c)
}

#' Stratum test
#' FAULTY, DON'T USE
#'
#' @param clust 

stratTest <- function(clust) {
  
  clust <- clust %>%
    mutate(mid = (base+top+he+ht)/4)
  sTab <- clust %>%
    group_by(cluster)%>%
    summarise_if(is.numeric, mean)
  o<- sTab[wrapr::orderv(sTab[,11]),]
  
  o$test <- 0
  for (n in 2:nrow(o)) {
    o$test[n] <- as.numeric(o$mid[n]<sum(o$top[1]:o$top[n-1])) # Using a sequence of this form adds numbers in increments of 1
  }
  
  out <- sum(o$test)
  return(out)
}

#' Arranges survey data into strata using k-means clustering
#' 
#' Data are stratified into 2-4 strata, then the largest number of strata 
#' are chosen where p<0.01. If none qualify, then the most significant
#' division is chosen.
#' 
#' Strata are sorted by top height and appended to the original data
#' 
#' 
#' @param veg A dataframe listing plant species with columns describing crown dimensions
#' @param pN The number of the point in the transect
#' @param spName Name of the field with the species name
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @param he Name of the field with dimension he
#' @param mStrat Maximum number of strata
#' @param sepSig p value to define significant stratum separation
#' @param ht Name of the field with dimension ht
#'
#' @return Dataframe
#' @export

frameStratify <- function(veg, pN ="Point", spName ="Species", base = "base", top = "top", he = "he", ht = "ht", mStrat = 4, sepSig = 0.001)
{
  
  veg <- datClean(veg = veg, base, top, he, ht)
  veg_subset <- veg %>% dplyr::select(all_of(c(pN, spName, base, top, he, ht)))
  veg_subset <- veg_subset[complete.cases(veg_subset), ] # Omit NAs in relevant columns
  veg_subset <- veg_subset %>% #log-scale dimensions for stratification
    mutate(base = case_when(veg_subset[,3] == 0 ~ 0.001, TRUE ~ veg_subset[,3]),
           lBase = log(veg_subset[,3]),
           lBase = case_when(is.infinite(lBase) ~ -6.9, TRUE ~ veg_subset[,3]),
           lTop = log(veg_subset[,4]),
           he = case_when(veg_subset[,5] == 0 ~ 0.001, TRUE ~ veg_subset[,5]),
           lhe = log(veg_subset[,5]),
           lhe = case_when(is.infinite(lhe) ~ -6.9, TRUE ~ veg_subset[,5]),
           lht = log(veg_subset[,6]))
  df <- scale(veg_subset[, c(7,8,9,10)])
  
  # Find the best division of strata
  sig <- vector()
  set.seed(123)
  if (!is.error(kmeans(df, centers = 2, nstart = 25))) {
    for (nstrat in 2:mStrat) {
      set.seed(123)
      if (!is.error(kmeans(df, centers = nstrat, nstart = 25))){
        km.res <- kmeans(df, centers = nstrat, nstart = 25)
        clust <- cbind(veg_subset, cluster = km.res$cluster)
        #testa <- frame:::stratTest(clust) Test is faulty
        testa <- 0
        test <- aov(cluster ~ base * top * he * ht, data = clust)
        sig[nstrat] <- base::summary(test)[[1]][["Pr(>F)"]][[3]]+testa
      }
    }
    if (length(which(sig < sepSig)) > 0) {
      nstrat <- as.numeric(max(which(sig < sepSig)))
    } else {
      if (length(sig[!is.na(sig)])>0) {
        nstrat <- as.numeric(min(which(sig == min(sig, na.rm = TRUE))))
      } else {
        nstrat <- 1
      }
    }
    rm(list=".Random.seed", envir=globalenv())
    set.seed(123)
    km.res <- kmeans(df, centers = nstrat, nstart = 25)
    clust <- cbind(veg_subset, cluster = km.res$cluster)
    
    # Summarise strata and order by mean height
    h <- clust %>% 
      mutate(mid = (base+top+he+ht)/4)%>%
      group_by(cluster) %>% 
      summarise_if(is.numeric, mean)
    h <- h[wrapr::orderv(h[,11]),] %>% 
      mutate(Stratum = 1:nstrat) %>% 
      select(cluster, Stratum)
    
    strat <- left_join(clust, h, by = "cluster") %>% 
      dplyr::select(all_of(c(pN, spName, top, "Stratum")))
    veg <- left_join(veg, strat, by = c(pN, spName, top))
  } else {
    veg$Stratum <- 1
  }
  rm(list=".Random.seed", envir=globalenv())
  return(veg)
}


#' Finds the distribution of species richness at a point
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @return dataframe
#' @export

rich <- function(dat, thres = 5, pnts = 10) {
  
  # Group minor species
  spCov <- frame::specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- suppressMessages(left_join(dat, spCov))%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"))
  
  y <- suppressMessages(dat %>%
                          group_by(Site, Point) %>%
                          summarise(n_distinct(Species)))
  
  #DATA ANALYSIS
  fitr <- data.frame('Mean' = character(0), 'SD' = character(0), 'Min' = character(0), 'Max' = character(0), stringsAsFactors=F)
  
  #Summary stats
  meanw <- round(mean(y$`n_distinct(Species)`, na.rm = TRUE),1)
  sdw <- round(sd(y$`n_distinct(Species)`, na.rm = TRUE), 2)
  minw <- as.numeric(min(y$`n_distinct(Species)`, na.rm = TRUE))
  maxw <- as.numeric(max(y$`n_distinct(Species)`, na.rm = TRUE))
  
  #Record values
  fitr[nrow(fitr)+1,] <- c(meanw, sdw, minw, maxw)
  
  return(fitr)
}



#' Finds the distribution of species richness per stratum
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' base - base height of each species
#' top - top height of each species
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param pN The number of the point in the transect
#' @param spName Name of the field with the species name
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @param he Name of the field with dimension he
#' @param ht Name of the field with dimension ht
#' @return dataframe

richS <- function(dat, thres = 0, pnts = 10, 
                  pN ="Point",  spName ="Species",  base = "base", 
                  top = "top", he = "he", ht = "ht", sepSig = 0.001) {
  
  spCov <- frame::specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- suppressMessages(left_join(dat, spCov))%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"))
  
  datS <- frame::frameStratify(veg = dat, pN = pN, spName = spName, base = base,
                          top = top, he = he, ht = ht, sepSig = sepSig)
  
  out <- suppressMessages(datS %>%
                            group_by(Stratum) %>%
                            summarise(n_distinct(Species)))
  out$Richness <- as.numeric(out$`n_distinct(Species)`)
  out <- out %>% select(Stratum, Richness)
  
  return(out)
}


#' Divides site data into consecutively numbered strata with
#' base and top heights
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' A field with the base height of the plants
#' A field with the top height of the plants
#' he
#' ht
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param pN The number of the point in the transect
#' @param spName Name of the field with the species name
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @param he Name of the field with dimension he
#' @param ht Name of the field with dimension ht
#' @return dataframe
#' @export
#' 
stratSite <- function(dat, thres = 0, pN ="Point",  spName ="Species",  base = "base", 
                      top = "top", he = "he", ht = "ht", sepSig = 0.001)  {
  pnts <- nrow(dat)
  strataDet <- data.frame(Stratum = numeric(0), Cover = numeric(0), 
                          Base = numeric(0), Top = numeric(0), stringsAsFactors = F)
  strat <- frame::frameStratify(veg = dat, pN = pN, spName = spName, base = base,
                                top = top, he = he, ht = ht, sepSig = sepSig)
  for (st in 1:as.numeric(max(strat$Stratum))) {
    stratSub <- strat %>% filter(Stratum == st)
    spnts <- unique(stratSub$Point, incomparables = FALSE)
    co <- length(spnts)/pnts
    b <- mean(stratSub$base)
    t <- mean(stratSub$top)
    strataDet[nrow(strataDet) + 1, ] <- c(as.numeric(st), as.numeric(co),
                                          as.numeric(b), as.numeric(t))
  }
  return(strataDet)
}


#' Calculates the range, mean and sd of 
#' species richness in each plant stratum
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' A field with the base height of the plants
#' A field with the top height of the plants
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param pN The number of the point in the transect
#' @param spName Name of the field with the species name
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @param he Name of the field with dimension he
#' @param ht Name of the field with dimension ht
#' @return dataframe
#' @export
#' 
stratRich <- function(dat, thres = 5, pnts = 10, 
                      pN ="Point",  spName ="Species",  base = "base", 
                      top = "top", he = "he", ht = "ht") {
  
  richList <- data.frame('S1' = numeric(0), 'S2' = numeric(0),
                         'S3' = numeric(0), 'S4' = numeric(0))
  slist <- unique(dat$Site, incomparables = FALSE)
  
  for (s in slist) {
    datSite <- filter(dat, Site == s)
    sRich <- richS(dat = datSite, thres = thres, pnts = pnts, 
                   pN = pN,  spName = spName,  base = base, 
                   top = top, he = he, ht = ht)
    nstrat <- as.numeric(max(sRich$Stratum))
    # Record values
    richList[which(slist == s), 1] <- sRich$Richness[1]
    if (nstrat > 1) {
      richList[which(slist == s),2] <- sRich$Richness[2]
    }
    if (nstrat > 2) {
      richList[which(slist == s),3] <- sRich$Richness[3]
    }
    if (nstrat > 3) {
      richList[which(slist == s),4] <- sRich$Richness[4]
    }
  }
  
  #DATA ANALYSIS
  fitr <- data.frame('Stratum' = numeric(0), 'Mean' = numeric(0), 'SD' = numeric(0), 
                     'Min' = numeric(0), 'Max' = numeric(0), stringsAsFactors=F)
  fitr[1,1] <- 1
  fitr[1,2] <- mean(richList$S1, na.rm = TRUE)
  fitr[1,3] <- sd(richList$S1, na.rm = TRUE)
  fitr[1,4] <- min(richList$S1, na.rm = TRUE)
  fitr[1,5] <- max(richList$S1, na.rm = TRUE)
  fitr[2,1] <- 2
  fitr[2,2] <- mean(richList$S2, na.rm = TRUE)
  fitr[2,3] <- sd(richList$S2, na.rm = TRUE)
  fitr[2,4] <- min(richList$S2, na.rm = TRUE)
  fitr[2,5] <- max(richList$S2, na.rm = TRUE)
  fitr[3,1] <- 3
  fitr[3,2] <- mean(richList$S3, na.rm = TRUE)
  fitr[3,3] <- sd(richList$S3, na.rm = TRUE)
  fitr[3,4] <- min(richList$S3, na.rm = TRUE)
  fitr[3,5] <- max(richList$S3, na.rm = TRUE)
  fitr[4,1] <- 4
  fitr[4,2] <- mean(richList$S4, na.rm = TRUE)
  fitr[4,3] <- sd(richList$S4, na.rm = TRUE)
  fitr[4,4] <- min(richList$S4, na.rm = TRUE)
  fitr[4,5] <- max(richList$S4, na.rm = TRUE)
  
  
  return(fitr)
}

#' Constructs the table F_flora from formatted survey data
#'
#' @param veg The dataframe containing the input data
#' @param pN The number of the point in the transect
#' @param spName Name of the field with the species name
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @param he Name of the field with dimension he
#' @param ht Name of the field with dimension ht
#' @param wid Name of the field with dimension width
#' @param rec Name of the field with the record number
#' @param sN Optional field with a site name
#' @param surf Weight of surface litter in t/ha
#' @return dataframe
#' @export
#'
#'

buildFlora <- function(veg, pN ="Point",  spName ="Species", base = "base", top = "top", he = "he", ht = "ht",
                       wid = "width", rec = "Site", sN = "SiteName", surf = 20, sepSig = 0.001) {

  vegA <- frame::frameStratify(veg = veg, pN = pN, spName = spName,
                          base = base, top = top, he = he, ht = ht, sepSig = sepSig)
  
  # Summarise species
  spCount <- vegA %>%
    dplyr::count(Stratum, Species, name = "comp")
  suppressMessages(spM <- vegA %>%
                     group_by(Stratum, Species) %>%
                     summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))))
  # Remove faulty data
  entries <- which(is.na(spM[wid]))
  if (length(entries)>0) {
    cat(length(entries), "Species were removed due to faulty crown width data", "\n")
    spM <- spM[-entries,]
  }
  
  suppressMessages(spSD <- vegA %>%
                     group_by(Stratum, Species) %>%
                     summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE))))
  if (length(entries)>0) {
    spSD <- spSD[-entries,]
  }
  
  # Correct empty entries
  singles <- which(is.na(spSD[top]))
  if (length(singles)>0) {
    spSD[top][is.na(spSD[top])]<-0.01
  }
  
  suppressMessages(spMax <- vegA %>%
                     group_by(Stratum, Species) %>%
                     summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))))
  if (length(entries)>0) {
    spMax <- spMax[-entries,]
  }
  suppressMessages(spMin <- vegA %>%
                     group_by(Stratum, Species) %>%
                     summarise(across(where(is.numeric), ~ min(.x, na.rm = TRUE))))
  if (length(entries)>0) {
    spMin <- spMin[-entries,]
  }
  tab <- (left_join(spMin, spCount, by = c("Stratum", "Species")))
  
  # Collate into table
  ns <- vegA %>% dplyr::select(all_of(c(rec, sN)))
  record <- matrix(nrow = length(spMin$Species))
  florTab <- data.frame(record)
  if (hasArg(rec)) {
    florTab$record <- ns$Site[1]
  } else {
    print("record field has not been named")
  }
  
  if (hasArg(sN)) {
    florTab$site <- ns$SiteName[1]
  } else {
    florTab$site <- NA
    print("site field has not been named")
  }
  
  florTab$species <- tab$Species
  florTab$stratum <- tab$Stratum
  florTab$comp <- tab$comp
  florTab$base <- round(spM$base,2)
  florTab$he <- round(spM$he,2)
  florTab$ht <- round(spM$ht,2)
  florTab$top <- round(spM$top,2)
  florTab$w <- round(spM$width,2)
  florTab$Hs <- round(spSD$top,2)
  florTab$Hr <- round(spMax$top - spMin$top,2)
  florTab$weight <- NA
  florTab$diameter <- NA
  s <- c(florTab$record[1], florTab$site[1], "Litter", NA, NA, NA, NA, NA, NA, NA, NA, NA, surf, 0.005)
  florTab <- rbind(florTab, s)
  
  return(florTab)
}

#' Constructs the table F_structure from formatted survey data
#'
#' @param veg The dataframe containing the input data
#' @param pN The number of the point in the transect
#' @param spName Name of the field with the species name
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @param he Name of the field with dimension he
#' @param ht Name of the field with dimension ht
#' @param rec Name of the field with the record number
#' @param sN Optional field with a site name
#' @return dataframe
#' @export
#'

buildStructure <- function(veg, pN ="Point", spName ="Species", base = "base", top = "top", 
                           he = "he", ht = "ht", rec = "Site", sN = "SiteName", sepSig = 0.001) {
  
  # 1. Horizontal relationships  
  vegA <- frame::frameStratify(veg = veg, pN = pN, spName = spName, base = base,
                          top = top, he = he, ht = ht, sepSig = sepSig)
  suppressMessages(StratC <- vegA %>%
                     select(Point, Stratum)%>%
                     group_by(Stratum, Point) %>%
                     summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))))
  suppressMessages(StratW <- vegA %>%
                     select(Stratum, width)%>%
                     group_by(Stratum) %>%
                     summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))))
  
  sepTab <- as.data.frame(table(StratC$Stratum))
  pnts <- n_distinct(StratC$Point)
  sep <- vector()
  for (s in 1:max(StratC$Stratum)) {
    sep[s] <- max(sqrt((as.numeric(StratW[s,2])^2)/(sepTab[s,2]/pnts)), as.numeric(StratW[s,2]))
  }
  
  # 2. Vertical relationships
  StratO <- as.data.frame(table(StratC))
  ns_e <- vector()
  ns_m <- vector()
  e_m <- vector()
  e_c <- vector()
  m_c <- vector()
  for (pt in unique(StratO$Point, incomparables = FALSE)) {
    point <- filter(StratO, Point == pt)
    if (max(StratC$Stratum) > 2) {
      ns_e[pt] <- point$Freq[1]*point$Freq[2]
      if (max(StratC$Stratum) == 4) {
        ns_m[pt] <- point$Freq[1]*point$Freq[3]
        e_m[pt] <- point$Freq[2]*point$Freq[3]
        e_c[pt] <- point$Freq[2]*point$Freq[4]
        m_c[pt] <- point$Freq[3]*point$Freq[4]
      } else {
        e_c[pt] <- point$Freq[2]*point$Freq[3]
      }
    }
  }
  
  # 3. Species richness
  suppressMessages(StratR <- vegA %>%
                     select(Stratum, Species)%>%
                     group_by(Stratum) %>%
                     summarise(across(everything(), n_distinct)))
  
  # Collate into table
  ns <- vegA %>% dplyr::select(all_of(c(rec, sN)))
  record <- matrix(nrow = 1)
  strucTab <- data.frame(record)
  if (hasArg(rec)) {
    strucTab$record <- ns$Site[1]
  } else {
    print("record field has not been named")
  }
  
  if (hasArg(sN)) {
    strucTab$site <- ns$SiteName[1]
  } else {
    strucTab$site <- NA
    print("site field has not been named")
  }
  ## Separation
  strucTab$NS <- NA
  strucTab$El <- NA
  strucTab$Mid <- NA
  strucTab$Can <- NA
  strucTab$NS <- round(sep[1],2)
  if (max(StratC$Stratum) > 2) {
    strucTab$El <- round(sep[2],2)
    if (max(StratC$Stratum) == 4) {
      strucTab$Mid <- round(sep[3],2)
      strucTab$Can <- round(sep[4],2)
    } else {
      strucTab$Can <- round(sep[3],2)
    }
  } else {
    strucTab$Can <- round(sep[2],2)
  }
  ## Overlap
  strucTab$ns_e <- NA
  strucTab$ns_m <- NA
  strucTab$e_m <- NA
  strucTab$e_c <- NA
  strucTab$m_c <- NA
  if (max(StratC$Stratum) > 2) {
    strucTab$ns_e <- sum(ns_e)>=as.numeric(pnts/2)
    if (max(StratC$Stratum) == 4) {
      strucTab$ns_m <- sum(ns_m)>=as.numeric(pnts/2)
      strucTab$e_m <- sum(e_m)>=as.numeric(pnts/2)
      strucTab$e_c <- sum(e_c)>=as.numeric(pnts/2)
      strucTab$m_c <- sum(m_c)>=as.numeric(pnts/2)
    } else {
      strucTab$e_c <- sum(e_c)>=as.numeric(pnts/2)
    }
  } 
  ## Richness
  strucTab$nsR <- NA
  strucTab$eR <- NA
  strucTab$mR <- NA
  strucTab$cR <- NA
  strucTab$nsR <- StratR$Species[1]
  if (max(StratC$Stratum) > 2) {
    strucTab$eR <- StratR$Species[2]
    if (max(StratC$Stratum) == 4) {
      strucTab$mR <- StratR$Species[3]
      strucTab$cR <- StratR$Species[4]
    } else {
      strucTab$cR <- StratR$Species[3]
    }
  } else {
    strucTab$cR <- StratR$Species[2]
  }
  
  return(strucTab)
}


#' Finds the height of near-surface litter
#' 
#' @param default.species.params Plant traits database
#' @param density Wood density (kg.m-3)
#' @param cover Percent cover of suspended layer
#' @param wNS Width of NS patches (m)
#' @param age Years since last fire
#' @param aQ Parameter for a quadratic trend; leave as NA if trend is negative exponential
#' @param bQ Parameter for a quadratic trend; leave as NA if trend is negative exponential
#' @param cQ Parameter for a quadratic trend; leave as NA if trend is negative exponential
#' @param dec Logical - TRUE allows for quadratic model to decline as vegetation thins,
#' @param maxNS Asymptote for negative exponential increase in NS
#' @param rateNS Rate for negative exponential increase in NS
#'
#' @return List
#' @export

susp <- function(default.species.params, density = 300, cover = 0.8,
                 age = 10, aQ = NA, bQ = NA, cQ = NA, maxNS = NA, rateNS = NA, dec = TRUE)
{
  suspDat <- filter(default.species.params, name == "suspNS")
  
  # Model packing
  if (dec == TRUE) {
    suspNS <- if (!is.na(maxNS)) {
      maxNS*(1-exp(-rateNS*age))
    } else if (!is.na(aQ)) {
      pmax(0,(aQ * age^2 + bQ * age + cQ))
    } else {
      0.1
    }
  } else {
    preLit <- vector()
    for (x in 1:age) {
      preLit[x] <- if (!is.na(maxNS)) {
        maxNS*(1-exp(-rateNS*x))
      } else if (!is.na(aQ)) {
        pmax(0,(aQ * x^2 + bQ * x + cQ))
      } else {
        0.1
      }
    }
    suspNS <- max(preLit)
  }
  
  nSticks <- (0.1*suspNS/cover)/(pi*(suspDat$leafThickness/2)^2*density) #Number of sticks
  top <- max(0.1,floor(nSticks*suspDat$leafSeparation) - 1) * suspDat$leafSeparation # Sets the lowest layer on the ground, finds height as separation * number of layers
  if (length(top) == 0) {
    top <- 0
  }
  
  return(list(top, suspNS))
}


#' Models the weight of surface litter from time since fire 
#' using either an Olson negative exponential function or a Burr curve
#'
#' @param negEx Value determining the model used. 
#' 1 = olson, 2 = Burr
#' @param max Maximum weight (t/ha)
#' @param rate Growth rate for a negative exponential function
#' @param a Parameter in the Burr equation
#' @param b Parameter in the Burr equation
#' @param age The number of years since last fire
#' @return value
#' @export

litter <- function(negEx = 1, max = 54.22, rate = 0.026, a = 3.35, b = 0.832, age = 10)
{
  if (negEx == 1) {
    surf <- max(4, max*(1-exp(-rate*age)))
  } else {
    surf <- max(4, a*b*((0.1*age^(a-1))/((1+(0.1*age)^a)^b+1)))
  }
  
  return(surf)
}


#' Combines multiple transects into one
#'
#' @param alldata Raw survey data with one or more transects
#'
#' @return dataframe
#' @export
#'

transectLong <- function(alldata){
  Sites <- unique(alldata$Site)
  maxPoint <- max(dplyr::filter(alldata, Site == Sites[1])$Point)
  out <- dplyr::filter(alldata, Site == Sites[1])
  
  # Rename consecutive sites
  for (s in Sites) {
    if (s > Sites[1]) {
      outA <- dplyr::filter(alldata, Site == Sites[s]) %>%
        mutate(Site = Sites[1],
               Point = Point + maxPoint)
      maxPoint <- max(outA$Point)
      out <- rbind(out, outA)
    }
  }
  return(out)
}

#' Processes field survey data into tables formatted for fire modelling
#'
#' @param dat The dataframe containing the formatted field survey data
#' @param default.species.params Plant traits database
#' @param pN The number of the point in the transect
#' @param spName Name of the field with the species name
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @param he Name of the field with dimension he
#' @param ht Name of the field with dimension ht
#' @param wid Name of the field with dimension width
#' @param rec Name of the field with the record number
#' @param sN Optional field with a site name
#' @param max Maximum weight of surface litter (t/ha)
#' @param rate Growth rate for surface litter in a negative exponential function
#' @param age Optional field with site age 
#' @param surf Weight of surface litter in t/ha
#' @param density Wood density (kg.m-3)
#' @param cover Percent cover of suspended layer
#' @param aQ Parameter for a quadratic trend; leave as NA if trend is negative exponential
#' @param bQ Parameter for a quadratic trend; leave as NA if trend is negative exponential
#' @param cQ Parameter for a quadratic trend; leave as NA if trend is negative exponential
#' @param maxNS Parameter for a negative exponential trend; leave as NA if trend is quadratic
#' @param rateNS Parameter for a negative exponential trend; leave as NA if trend is quadratic
#' @param thin Logical - TRUE uses plant self-thinning models in Dynamics,
#' otherwise plant cover remains at the highest point it has reached by that age
#' @param sLit Logical - TRUE allows surface litter to decline if the model does so, otherwise
#' the maximum value to that age is maintained
#' @param dec Logical - TRUE allows near surface surface litter to decline if the model does so, otherwise
#' @param negEx Value determining the model used. 1 = olson, 2 = Burr 
#' @param a 
#' @param b 
#' @param wNS Width of near surface patches (m)
#' @param sepSig 
#' the maximum value to that age is maintained
#' @return list
#' @export
#' 

frameSurvey <- function(dat, default.species.params, pN ="Point", spName ="Species", base = "base", top = "top", he = "he", ht = "ht",
                        wid = "width", rec = "Site", sN = "SiteName", negEx = 1, max = 54.22, rate = 0.026, a = 3.35, b = 0.832, age = NA, 
                        surf = 10, density = 300, cover = 0.8, aQ = NA, bQ = NA, cQ = NA, maxNS = NA, rateNS = NA, wNS = 1,
                        thin = TRUE, sLit  = TRUE, dec = TRUE, sepSig = 0.001) {
  
  # Find missing data
  entries <- which(is.na(dat[top]))
  if (length(entries)>0) {
    cat(" These rows were removed as they were missing top heights", "\n", entries, "\n", "\n")
    dat <- dat[-entries,] 
  }
  
  # Fill empty dimensions
  entries <- which(is.na(dat[ht])|is.na(dat[he]))
  if (length(entries)>0) {
    cat(" Empty values of ht & he were filled with top and base heights for these rows", "\n", entries, "\n", "\n")
    dat[ht][is.na(dat[ht])]<-dat[top][is.na(dat[ht])]
    dat[he][is.na(dat[he])]<-dat[base][is.na(dat[he])]
  }
  
  # Remove faulty data
  entries <- which(dat[ht]<dat[he]|dat[top]<dat[base])
  if (length(entries)>0) {
    cat(" These rows were removed as upper and lower heights conflicted", "\n", entries)
    dat <- dat[-entries,]
  }
  
  # Loop through records
  silist <- unique(dat$Site, incomparables = FALSE)
  Structure <- data.frame()
  Flora <- data.frame()
  
  for (rec in silist) {
    veg <- filter(dat, Site == rec)
    print(rec)
    # Find surface litter
    if (!is.na(age)) {
      AGE <- veg[1,age]
      
      # Control self-thinning
      (if (thin == TRUE && sLit == TRUE) {
        surf <- round(litter(negEx, max, rate, a, b, AGE),0)
      } else {
        preLit <- vector()
        for (x in 1:AGE) {
          preLit[x] <- litter(negEx, max, rate, a, b, x)
        }
        surf <- max(preLit)
      })
    }
    
    # Add suspended litter
    if (cover != 0) {
      if (thin == TRUE && dec == TRUE) {
        decline <- TRUE
      } else {
        decline <- FALSE
      }
      suspNS <- susp(default.species.params, density = density, cover = cover,
           age = AGE, aQ = aQ, bQ = bQ, cQ = cQ, maxNS = maxNS, rate = rateNS, dec = decline)
      top <- suspNS[[1]]
      #Update tables
      if (top > 0) {
        row <- nrow(veg)
        rows <- round(cover*length(unique(veg$Point)),0)
        for (r in 1:rows) {
          veg[row+r,1] <- AGE
          veg[row+r,2] <- "suspNS"
          veg[row+r,3] <- 0
          veg[row+r,4] <- top
          veg[row+r,5] <- 0
          veg[row+r,6] <- top
          veg[row+r,7] <- wNS
          veg[row+r,8] <- AGE
          veg[row+r,9] <- rec
        }
      }
    }
    
    Struct <- buildStructure(veg, pN ="Point", spName ="Species", base = "base", top = "top", 
                             rec = "Site", sN = "SiteName", sepSig = sepSig)
    Flor <- buildFlora(veg, pN ="Point",  spName ="Species", base = "base", top = "top", he = "he", ht = "ht",
                       wid = "width", rec = "Site", sN = "SiteName", surf = surf, sepSig = sepSig)
    Structure <- rbind(Structure, Struct)
    Flora <- rbind(Flora, Flor)
  }
  
  return(list(Flora, Structure))
}

#' Grows a table of species to a given age, listing cover and variability
#' 
#' Provides options to control for growth, self-thinning and self-pruning
#'
#' @param Dynamics Table of growth models as output by floraDynamics
#' @param Age Age of the site in years
#' @param growth Logical - TRUE uses plant growth models in Dynamics,
#' otherwise mean values are used and plants are not grown
#' @param thin Logical - TRUE uses plant self-thinning models in Dynamics,
#' otherwise plant cover remains at the highest point it has reached by that age
#' @param prune Logical - TRUE uses plant self-pruning models in Dynamics,
#' otherwise plants are not self-pruned and lower canopy heights remain at their lowest points
#' @return dataframe
#' @export

growPlants <- function(Dynamics, Age = 10, growth = TRUE, thin = TRUE, prune = TRUE) {
  Contenders <- data.frame('Species' = character(0), 'Cover' = numeric(0), 'cRSE' = numeric(0),
                           'base' = numeric(0), 'he' = numeric(0), 'ht' = numeric(0), 
                           'top' = numeric(0), 'tRSE' = numeric(0), 'w' = numeric(0))
  
  for (sp in 1:(nrow(Dynamics)-2)) {
    Contenders[sp,1] <- max(Dynamics$Species[sp],0)
    Contenders[sp,3] <- as.numeric(Dynamics$C_RSE[sp])/100
    
    # Control growth
    if (growth == TRUE) {
      Contenders[sp,7] <- max(pTop(mods = Dynamics, sp=Dynamics$Species[sp], Age = Age),0)
      Contenders[sp,8] <- as.numeric(Dynamics$T_RSE[sp])/100 
      Contenders[sp,9] <- max(pWidth(mods = Dynamics, sp=Dynamics$Species[sp], Age = Age),0)
    } else {
      Contenders[sp,7] <- max(as.numeric(Dynamics$meanTop[sp]),0)
      Contenders[sp,8] <- 0.01 
      Contenders[sp,9] <- max(as.numeric(Dynamics$meanW[sp]),0) 
    }
    
    # Plant Ht follows top height in all circumstances
    Contenders[sp,6] <- max(pHt(mods = Dynamics, sp=Dynamics$Species[sp], Age = Age),0)
    
    # Control self-pruning by keeping base heights at minimum values if prune == FALSE
    if (prune == TRUE) {
      Contenders[sp,5] <- max(min(pHe(mods = Dynamics, sp=Dynamics$Species[sp], Age = Age), Contenders[sp,6]),0)
      Contenders[sp,4] <- max(min(pBase(mods = Dynamics, sp=Dynamics$Species[sp], Age = Age), Contenders[sp,7]),0)
    } else {
      preHe <- vector()
      preBase <- vector()
      for (x in 1:Age) {
        preHe[x] <- max(min(pHe(mods = Dynamics, sp=Dynamics$Species[sp], Age = x), Contenders[sp,6]),0)
        preBase[x] <- max(min(pBase(mods = Dynamics, sp=Dynamics$Species[sp], Age = x), Contenders[sp,7]),0)
      }
      Contenders[sp,5] <- min(preHe)
      Contenders[sp,4] <- min(preBase)
    }
    
    # Assign no cover to plants with 0 height or foliage
    covTest <- if ((Contenders[sp,7] < 0.05) 
                   || (Contenders[sp,4] == Contenders[sp,7] & Contenders[sp,6] == Contenders[sp,5]) 
                   || (Contenders[sp,9] <= 0.05)) {
      0
    } else {
      1
    }
    
    # Control self-thinning
    (if (thin == TRUE) {
      Contenders[sp,2] <- covTest * pCover(mods = Dynamics, sp=Dynamics$Species[sp], Age = Age)
    } else {
      preCov <- vector()
      for (x in 1:Age) {
        preCov[x] <- pCover(mods = Dynamics, sp=Dynamics$Species[sp], Age = x)
      }
      Contenders[sp,2] <- covTest * max(preCov)
    })
  }
  # Remove species with no cover
  entries <- which(Contenders$Cover == 0)
  if (length(entries)>0) {
    Contenders <- Contenders[-entries,] 
  }
  
  return(Contenders)
}

#' Creates a dataset of plants with expected size and cover
#' randomly varied within natural ranges
#'
#' @param Dynamics Table of growth models as output by floraDynamics
#' @param pointRich table of species richness likelihood for a point
#' as output by function rich
#' @param default.species.params Plant traits database
#' @param pnts number of points in the transect
#' @param Age Age of the site in years
#' @param growth Logical - TRUE uses plant growth models in Dynamics,
#' otherwise mean values are used and plants are not grown
#' @param thin Logical - TRUE uses plant self-thinning models in Dynamics,
#' otherwise mean values are used and plant cover remains constant
#' @param prune Logical - TRUE uses plant self-pruning models in Dynamics,
#' otherwise plants are not self-pruned and lower canopy heights remain at their lowest points
#' @return dataframe
#' @export

pseudoTransect <- function(Dynamics, pointRich, default.species.params,
                           pnts = 10, Age = 10, growth = TRUE, thin = TRUE, prune = TRUE) 
{
  # Build species choice function
  chooseSp <- function(L, spList, drops) {
    temp <- L
    xRow <- data.frame()
    Lentry <- as.numeric(nrow(xRow))
    tot <- sum(temp)
    while (Lentry == 0 & tot >= nSpecies) {
      entry <-
        spList[which(temp == Rfast::nth(temp, 1, descending = TRUE)),]
      if (Rfast::nth(temp, 1, descending = TRUE) > 100 * runif(1)*(drops>0)) {
        xRow <- rbind(xRow, entry)
      } else {
        temp[which(temp == Rfast::nth(temp, 1, descending = TRUE))] <-
          0
      }
      Lentry <- as.numeric(nrow(xRow))
      tot <- sum(temp > 0)
    }
    return(xRow)
  }
  # Potential species
  spList <- growPlants(Dynamics, Age = Age, growth = growth, thin = thin, prune = prune) %>%
    mutate(Cover = 100-(pmin(1,5/top)*(100-Cover))) # Weight cover by height above 5m to reflect angle to canopy
  
  out <- data.frame('Point' = numeric(0), 'Species' = character(0), 'base' = numeric(0),
                    'top' = numeric(0), 'he' = numeric(0), 'ht' = numeric(0), 'width' = numeric(0),
                    'Age' = numeric(0), 'Site' = numeric(0), 'SiteName' = character(0), stringsAsFactors = FALSE)
  
  # Build pseudo-transect 
  for (r in 1:pnts) {
    # Find the number of species at this point
    spN  <-  round(rtnorm(n = 1, mean = as.numeric(pointRich$Mean[1]),
                          sd = as.numeric(pointRich$SD[1]),
                          a = as.numeric(pointRich$Min[1]),
                          b = as.numeric(pointRich$Max[1])),0)
    
    # Calculate likelihood for each species
    L <- vector()
    for (sp in 1:length(spList$Species)) {
      L[sp] <- rtnorm(n = 1, mean = as.numeric(spList$Cover[sp]),
                      sd = as.numeric(spList$cRSE[sp]), a = 0, b = 100)
    }
    
    # Choose species and record dimensions
    drops <- length(spList$Species) - spN
    for (nSpecies in 1:spN) {
      count <- nrow(out)+1
      entry <- chooseSp(L, spList, drops)
      if (nrow(entry) > 0) {
        # Remove chosen species from next option
        L[which(spList$Species == entry$Species)] <- 0
        out[count, 1] <- r
        out[count, 2] <- entry$Species[1]
        out[count, 4] <- round(rtnorm(n = 1, mean = as.numeric(entry$top[1]),
                                      sd = as.numeric(entry$tRSE[1]), a = 0, b = 150),2)
        out[count, 3] <- round(out[count, 4] * as.numeric(entry$base[1]),2)
        out[count, 5] <- round(out[count, 4] * as.numeric(entry$he[1]),2)
        out[count, 6] <- round(out[count, 4] * as.numeric(entry$ht[1]),2)
        out[count, 7] <- round(out[count, 4] * as.numeric(entry$w[1]),2)
        out[count, 8] <- Age
        out[count, 9] <- Age
        out[count, 10] <- Age
      } else {
        drops <- drops - 1
      }
    }
  }
  return(out)
}

#' Performs basic data checks on survey data, 
#' fixes or removes data and lists changes
#'
#' @param dat The dataframe to be cleaned
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @param he Name of the field with dimension he
#' @param ht Name of the field with dimension ht
#' @return dataframe
#' @export

datClean <- function(veg,  base = "base", top = "top", he = "he", ht = "ht") {
  
  # Find missing data
  entries <- which(is.na(veg[top]))
  if (length(entries)>0) {
    cat(" Removed these rows as they were missing top heights", "\n", entries, "\n", "\n")
    veg <- veg[-entries,] 
  }
  
  # Fill empty dimensions
  entries <- which(is.na(veg[ht])|is.na(veg[he]))
  if (length(entries)>0) {
    cat(" Filled empty values of ht & he with top and base heights for these rows", "\n", entries, "\n", "\n")
    veg[ht][is.na(veg[ht])]<-veg[top][is.na(veg[ht])]
    veg[he][is.na(veg[he])]<-veg[base][is.na(veg[he])]
  }
  
  # Remove plants <1cm
  entries <- which(veg[top]<0.01)
  if (length(entries)>0) {
    cat(" Removed these rows as species were too small to model", "\n", entries, "\n", "\n")
    veg <- veg[-entries,]
  }
  
  # Set low bases to ground level
  entries <- which(veg[he]<0.01 | veg[base]<0.01)
  if (length(entries)>0) {
    cat(" Set one or both base values for these rows to zero", "\n", entries, "\n", "\n")
    veg <- veg %>%
      mutate(he = case_when(he < 0.01 ~ 0, TRUE ~ he)) %>%
      mutate(base = case_when(base < 0.01 ~ 0, TRUE ~ base))
  }
  
  # Remove faulty data
  entries <- which(veg[ht]<veg[he]|veg[top]<veg[base])
  if (length(entries)>0) {
    cat(" Removed these rows as upper and lower heights conflicted", "\n", entries)
    veg <- veg[-entries,]
  }
  return(veg)
}


#' Removes leaf trait diversity
#' 
#' Initially pulls out species names "suspNS" and "Log"
#'
#' @param default.species.params Plant traits database
#' @return dataframe
#' @export

ctrlDiversity <- function(default.species.params){
  
  no.succession.params <- filter(default.species.params, name != "suspNS" & name != "Log")
  no.succession.params$leafForm <- "Flat"
  means <- no.succession.params %>%
    summarise_if(is.numeric, mean)
  for (name in colnames(means)) {
    no.succession.params[,name] <- means[1,name]
  }
  
  sus <- filter(default.species.params, name == "suspNS" | name == "Log")
  no.succession.params <- rbind(no.succession.params, sus)
  
  return(no.succession.params)
}


#' Finds percent cover of surveyed Species and groups minor Species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' Species - name of the surveyed Species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be kept single
#' @param pnts The number of points measured in a transect
#' @return dataframe
#' @export

specCover <- function(dat, thres = 5, pnts = 10) {
  
  #List Species and ages
  spList <- unique(dat$Species, incomparables = FALSE)
  ages <- unique(dat$Age, incomparables = FALSE)
  
  #Create empty summary dataframe
  spCover <- data.frame('Species' = character(0), 'Age' = numeric(0), 'Cover' = numeric(0), stringsAsFactors=F)
  
  #DATA COLLECTION
  for (sp in 1:length(spList)) {
    for (age in ages) {
      spName <- dat %>% filter(Species == spList[sp])
      spAge <- spName %>% filter(Age == age)
      
      #Percent cover
      sppnts <- unique(spAge$Point, incomparables = FALSE)
      covSp <- as.numeric(length(sppnts))*(100/pnts)
      
      #Record values
      spCover[nrow(spCover)+1,] <- c(as.character(spList[sp]), as.numeric(age), as.numeric(covSp))
    }
  }
  
  #Group minor Species
  spCover$Cover <- as.numeric(as.character(spCover$Cover))
  spShort <- spCover %>%
    group_by(Species) %>%
    summarise_if(is.numeric, mean)
  #List minor Species, then rename in dataset
  minor <- spShort %>% filter(Cover < thres)
  minList <- unique(minor$Species, incomparables = FALSE)
  if (length(minList)>0) {
    for (snew in 1:length(minList)) {
      spCover[spCover == minor$Species[snew]] <- "Minor Species"
    }
  }
  return(spCover)
}