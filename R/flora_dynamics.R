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
  
spCov <- specCover(dat = dat, thres = thres, pnts = pnts)
  priorList <- unique(spCov$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitCov <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                       'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_p' = numeric(0),
                       'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_p' = numeric(0),
                       'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_p' = numeric(0),
                       'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'Q_p' = numeric(0),
                       'linear' = character(0), 'NegExp' = character(0), 'Burr' = character(0), 'Binomial' = character(0), 'Quadratic' = character(0), 'Mean' = character(0),
                       'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
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
      LMp <- if (LMRSE == 0) {
        0
      } else {
        round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
      }
      LM_sig <- if (LMp < 0.001) {
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
    
    #Negative exponential
    init1<-c(k=50,r=0.5)
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init2,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
      NESum <- base::summary(NE)
      k <- NESum$coefficients[1]
      r <- NESum$coefficients[2]
      NERSE <- NESum$sigma
      NEp <- round(max(NESum$coefficients[7],NESum$coefficients[8]),5)
      NE_sig <- if (NEp < 0.001) {
        "***"    
      } else if (NEp < 0.01) {
        "**"    
      } else if (NEp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(NE)
    } else {
      k <- NA
      r <- NA
      NERSE <- NA
      NEp <- 1
      NE_sig <- ""
    }
    
    #Burr
    init2<-c(a=3,b=2)
    if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
      Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
      BSum <- base::summary(Burr)
      Ba <- BSum$coefficients[1]
      Bb <- BSum$coefficients[2]
      BRSE <- BSum$sigma
      #Added control for Burr
      f <- function(x){Ba*Bb*((0.1*x^(Ba-1))/((1+(0.1*x)^Ba)^Bb+1))}
      Bp <- if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        1
      } else if (bTest * (mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE)) > (optimize(f = f, interval=c(0, 150), maximum=FALSE))$objective) {
        1
      } else {round(max(BSum$coefficients[7],BSum$coefficients[8]),5)}
      B_sig <- if (Bp < 0.001) {
        "***"    
      } else if (Bp < 0.01) {
        "**"    
      } else if (Bp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(Burr)
    } else {
      Ba <- NA
      Bb <- NA
      BRSE <- NA
      Bp <- 1
      B_sig <- ""
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
      Binp <- round(max(BinSum$coefficients[10],BinSum$coefficients[11],BinSum$coefficients[12]),5)
      Bin_sig <- if (Binp < 0.001) {
        "***"    
      } else if (Binp < 0.01) {
        "**"    
      } else if (Binp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(Bin)
    } else {
      Bs <- NA
      Bsd <- NA
      Bm <- NA
      BinRSE <- NA
      Binp <- 1
      Bin_sig <- ""
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
      qRSE <- qSum$sigma
      #Added control for Quadratic
      f <- function(x){qa*x^2 + qb*x + qc}
      qp <- if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        1
      } else {round(max(qSum$coefficients[10], qSum$coefficients[11], 
                        qSum$coefficients[12]), 5)}
      q_sig <- if (qp < 0.001) {
        "***"
      }
      else if (qp < 0.01) {
        "**"
      }
      else if (qp < 0.05) {
        "*"
      }
      else {
        ""
      }
      rm(q)
    }
    else {
      qa <- NA
      qb <- NA
      qc <- NA
      qRSE <- NA
      qp <- 1
      q_sig <- ""
    }
    
    #Summary stats
    meanCov <- round(mean(y, na.rm = TRUE),1)
    m_sig <- round(mRSE(dat = y),3)
    model <- if (min(Binp,Bp,NEp,LMp, qp, na.rm = TRUE)<p) {
      if (Bp<=min(LMp,NEp,Binp, qp, na.rm = TRUE)) {
        "Burr"
      } else {
        if (LMp<=min(Bp,NEp,Binp, qp, na.rm = TRUE)) {
          "Linear"
        } else {
          if (NEp<=min(LMp,Bp,Binp, qp, na.rm = TRUE)) {
            "NegExp"
          } else {
            if (Binp<=min(LMp,Bp,NEp, qp, na.rm = TRUE)) {
              "Binomial"
            } else {
              if (qp<=min(LMp,Bp,NEp, Binp, na.rm = TRUE)) {
                "Quadratic"
              }
            }
          }
        }
      }
    } else {"Mean"}
    
    #Record values
    fitCov[nrow(fitCov)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMp, k, r, NERSE, NEp,
                                 Ba, Bb, BRSE, Bp, Bs, Bsd, Bm, BinRSE, Binp, qa, qb, qc, qRSE, qp, 
                                 LM_sig, NE_sig, B_sig, Bin_sig, q_sig, meanCov, m_sig, model)
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

topDyn <- function(dat, thres = 5, pnts = 10, p = 0.05, bTest = 10, maxiter = 1000) {
  
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"))
  
  priorList <- unique(dat$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitTop <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                       'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_p' = numeric(0),
                       'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_p' = numeric(0),
                       'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_p' = numeric(0),
                       'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'Q_p' = numeric(0),
                       'linear' = character(0), 'NegExp' = character(0), 'Burr' = character(0), 'Binomial' = character(0), 'Quadratic' = character(0), 'Mean' = character(0),
                       'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=maxiter, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- dat %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$top)
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
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
      LM_sig <- if (LMp < 0.001) {
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
    
    #Negative exponential
    init1<-c(k=50,r=0.5)
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init2,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
      NESum <- base::summary(NE)
      k <- NESum$coefficients[1]
      r <- NESum$coefficients[2]
      NERSE <- NESum$sigma
      NEp <- round(max(NESum$coefficients[7],NESum$coefficients[8]),5)
      NE_sig <- if (NEp < 0.001) {
        "***"    
      } else if (NEp < 0.01) {
        "**"    
      } else if (NEp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(NE)
    } else {
      k <- NA
      r <- NA
      NERSE <- NA
      NEp <- 1
      NE_sig <- ""
    }
    
    #Burr
    init2<-c(a=3,b=2)
    if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
      Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
      BSum <- base::summary(Burr)
      Ba <- BSum$coefficients[1]
      Bb <- BSum$coefficients[2]
      BRSE <- BSum$sigma
      #Added control for Burr
      f <- function(x){Ba*Bb*((0.1*x^(Ba-1))/((1+(0.1*x)^Ba)^Bb+1))}
      Bp <- if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        1
      } else if (bTest * (mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE)) > (optimize(f = f, interval=c(0, 150), maximum=FALSE))$objective) {
        1
      } else {round(max(BSum$coefficients[7],BSum$coefficients[8]),5)}
      B_sig <- if (Bp < 0.001) {
        "***"    
      } else if (Bp < 0.01) {
        "**"    
      } else if (Bp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(Burr)
    } else {
      Ba <- NA
      Bb <- NA
      BRSE <- NA
      Bp <- 1
      B_sig <- ""
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
      Binp <- round(max(BinSum$coefficients[10],BinSum$coefficients[11],BinSum$coefficients[12]),5)
      Bin_sig <- if (Binp < 0.001) {
        "***"    
      } else if (Binp < 0.01) {
        "**"    
      } else if (Binp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(Bin)
    } else {
      Bs <- NA
      Bsd <- NA
      Bm <- NA
      BinRSE <- NA
      Binp <- 1
      Bin_sig <- ""
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
      qRSE <- qSum$sigma
      #Added control for Quadratic
      f <- function(x){qa*x^2 + qb*x + qc}
      qp <- if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        1
      } else {round(max(qSum$coefficients[10], qSum$coefficients[11], 
                        qSum$coefficients[12]), 5)}
      q_sig <- if (qp < 0.001) {
        "***"
      }
      else if (qp < 0.01) {
        "**"
      }
      else if (qp < 0.05) {
        "*"
      }
      else {
        ""
      }
      rm(q)
    }
    else {
      qa <- NA
      qb <- NA
      qc <- NA
      qRSE <- NA
      qp <- 1
      q_sig <- ""
    }
    
    #Summary stats
    meanTop <- round(mean(y, na.rm = TRUE),1)
    m_sig <- round(mRSE(dat = y),3)
    model <- if (min(Binp,Bp,NEp,LMp, qp, na.rm = TRUE)<p) {
      if (Bp<=min(LMp,NEp,Binp, qp, na.rm = TRUE)) {
        "Burr"
      } else {
        if (LMp<=min(Bp,NEp,Binp, qp, na.rm = TRUE)) {
          "Linear"
        } else {
          if (NEp<=min(LMp,Bp,Binp, qp, na.rm = TRUE)) {
            "NegExp"
          } else {
            if (Binp<=min(LMp,Bp,NEp, qp, na.rm = TRUE)) {
              "Binomial"
            } else {
              if (qp<=min(LMp,Bp,NEp, Binp, na.rm = TRUE)) {
                "Quadratic"
              }
            }
          }
        }
      }
    } else {"Mean"}
    
    #Record values
    fitTop[nrow(fitTop)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMp, k, r, NERSE, NEp,
                                 Ba, Bb, BRSE, Bp, Bs, Bsd, Bm, BinRSE, Binp, qa, qb, qc, qRSE, qp, 
                                 LM_sig, NE_sig, B_sig, Bin_sig, q_sig, meanTop, m_sig, model)
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

baseDyn <- function(dat, thres = 5, pnts = 10, p = 0.05, bTest = 10, maxiter = 1000) {
  
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"),
           bRat = base/top)
  
  priorList <- unique(dat$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitBase <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                        'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_p' = numeric(0),
                        'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_p' = numeric(0),
                        'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_p' = numeric(0),
                        'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'Q_p' = numeric(0),
                        'linear' = character(0), 'NegExp' = character(0), 'Burr' = character(0), 'Binomial' = character(0), 'Quadratic' = character(0), 'Mean' = character(0),
                        'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=maxiter, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- dat %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$bRat)
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
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
    
    #Negative exponential
    init1<-c(k=50,r=0.5)
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init2,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
      NESum <- base::summary(NE)
      k <- NESum$coefficients[1]
      r <- NESum$coefficients[2]
      NERSE <- NESum$sigma
      NEp <- round(max(NESum$coefficients[7],NESum$coefficients[8]),5)
      NE_sig <- if (NEp < 0.001) {
        "***"    
      } else if (NEp < 0.01) {
        "**"    
      } else if (NEp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(NE)
    } else {
      k <- NA
      r <- NA
      NERSE <- NA
      NEp <- 1
      NE_sig <- ""
    }
    
    #Burr
    init2<-c(a=3,b=2)
    if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
      Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
      BSum <- base::summary(Burr)
      Ba <- BSum$coefficients[1]
      Bb <- BSum$coefficients[2]
      BRSE <- BSum$sigma
      #Added control for Burr
      f <- function(x){Ba*Bb*((0.1*x^(Ba-1))/((1+(0.1*x)^Ba)^Bb+1))}
      Bp <- if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        1
      } else if (bTest * (mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE)) > (optimize(f = f, interval=c(0, 150), maximum=FALSE))$objective) {
        1
      } else {round(max(BSum$coefficients[7],BSum$coefficients[8]),5)}
      B_sig <- if (Bp < 0.001) {
        "***"    
      } else if (Bp < 0.01) {
        "**"    
      } else if (Bp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(Burr)
    } else {
      Ba <- NA
      Bb <- NA
      BRSE <- NA
      Bp <- 1
      B_sig <- ""
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
      Binp <- round(max(BinSum$coefficients[10],BinSum$coefficients[11],BinSum$coefficients[12]),5)
      Bin_sig <- if (Binp < 0.001) {
        "***"    
      } else if (Binp < 0.01) {
        "**"    
      } else if (Binp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(Bin)
    } else {
      Bs <- NA
      Bsd <- NA
      Bm <- NA
      BinRSE <- NA
      Binp <- 1
      Bin_sig <- ""
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
      qRSE <- qSum$sigma
      #Added control for Quadratic
      f <- function(x){qa*x^2 + qb*x + qc}
      qp <- if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        1
      } else {round(max(qSum$coefficients[10], qSum$coefficients[11], 
                        qSum$coefficients[12]), 5)}
      q_sig <- if (qp < 0.001) {
        "***"
      }
      else if (qp < 0.01) {
        "**"
      }
      else if (qp < 0.05) {
        "*"
      }
      else {
        ""
      }
      rm(q)
    }
    else {
      qa <- NA
      qb <- NA
      qc <- NA
      qRSE <- NA
      qp <- 1
      q_sig <- ""
    }
    
    #Summary stats
    meanBase <- round(mean(y, na.rm = TRUE),1)
    m_sig <- round(mRSE(dat = y),3)
    model <- if (min(Binp,Bp,NEp,LMp, qp, na.rm = TRUE)<p) {
      if (Bp<=min(LMp,NEp,Binp, qp, na.rm = TRUE)) {
        "Burr"
      } else {
        if (LMp<=min(Bp,NEp,Binp, qp, na.rm = TRUE)) {
          "Linear"
        } else {
          if (NEp<=min(LMp,Bp,Binp, qp, na.rm = TRUE)) {
            "NegExp"
          } else {
            if (Binp<=min(LMp,Bp,NEp, qp, na.rm = TRUE)) {
              "Binomial"
            } else {
              if (qp<=min(LMp,Bp,NEp, Binp, na.rm = TRUE)) {
                "Quadratic"
              }
            }
          }
        }
      }
    } else {"Mean"}
    model <- if (min(LMRSE, NERSE, BRSE, BinRSE, qRSE, na.rm = TRUE) < m_sig){
      model
    } else {"Mean"}
    
    #Record values
    fitBase[nrow(fitBase)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMp, k, r, NERSE, NEp,
                                   Ba, Bb, BRSE, Bp, Bs, Bsd, Bm, BinRSE, Binp, qa, qb, qc, qRSE, qp, 
                                   LM_sig, NE_sig, B_sig, Bin_sig, q_sig, meanBase, m_sig, model)
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

heDyn <- function(dat, thres = 5, pnts = 10, p = 0.05) {
  
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
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
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
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
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

wDyn <- function(dat, thres = 5, pnts = 10, p = 0.05, bTest = 10, maxiter = 1000) {
  
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"),
           Rat = as.numeric(width)/top)
  
  priorList <- unique(dat$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitw <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                     'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_p' = numeric(0),
                     'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_p' = numeric(0),
                     'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_p' = numeric(0),
                     'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'Q_p' = numeric(0),
                     'linear' = character(0), 'NegExp' = character(0), 'Burr' = character(0), 'Binomial' = character(0), 'Quadratic' = character(0), 
                     'Mean' = character(0), 'Mean_sigma' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=maxiter, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- dat %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$Rat)
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
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
    
    #Negative exponential
    init1<-c(k=50,r=0.5)
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init2,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
      NESum <- base::summary(NE)
      k <- NESum$coefficients[1]
      r <- NESum$coefficients[2]
      NERSE <- NESum$sigma
      NEp <- round(max(NESum$coefficients[7],NESum$coefficients[8]),5)
      NE_sig <- if (NEp < 0.001) {
        "***"    
      } else if (NEp < 0.01) {
        "**"    
      } else if (NEp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(NE)
    } else {
      k <- NA
      r <- NA
      NERSE <- NA
      NEp <- 1
      NE_sig <- ""
    }
    
    #Burr
    init2<-c(a=3,b=2)
    if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
      Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
      BSum <- base::summary(Burr)
      Ba <- BSum$coefficients[1]
      Bb <- BSum$coefficients[2]
      BRSE <- BSum$sigma
      #Added control for Burr
      f <- function(x){Ba*Bb*((0.1*x^(Ba-1))/((1+(0.1*x)^Ba)^Bb+1))}
      Bp <- if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        1
      } else if (bTest * (mean(y, na.rm = TRUE) - sd(y, na.rm = TRUE)) > (optimize(f = f, interval=c(0, 150), maximum=FALSE))$objective) {
        1
      } else {round(max(BSum$coefficients[7],BSum$coefficients[8]),5)}
      B_sig <- if (Bp < 0.001) {
        "***"    
      } else if (Bp < 0.01) {
        "**"    
      } else if (Bp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(Burr)
    } else {
      Ba <- NA
      Bb <- NA
      BRSE <- NA
      Bp <- 1
      B_sig <- ""
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
      Binp <- round(max(BinSum$coefficients[10],BinSum$coefficients[11],BinSum$coefficients[12]),5)
      Bin_sig <- if (Binp < 0.001) {
        "***"    
      } else if (Binp < 0.01) {
        "**"    
      } else if (Binp < 0.05) {
        "*" 
      } else {
        ""
      }
      rm(Bin)
    } else {
      Bs <- NA
      Bsd <- NA
      Bm <- NA
      BinRSE <- NA
      Binp <- 1
      Bin_sig <- ""
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
      qRSE <- qSum$sigma
      #Added control for Quadratic
      f <- function(x){qa*x^2 + qb*x + qc}
      qp <- if (bTest * (sd(y, na.rm = TRUE) + mean(y, na.rm = TRUE)) < (optimize(f = f, interval=c(0, 150), maximum=TRUE))$objective) {
        1
      } else {round(max(qSum$coefficients[10], qSum$coefficients[11], 
                        qSum$coefficients[12]), 5)}
      q_sig <- if (qp < 0.001) {
        "***"
      }
      else if (qp < 0.01) {
        "**"
      }
      else if (qp < 0.05) {
        "*"
      }
      else {
        ""
      }
      rm(q)
    } else {
      qa <- NA
      qb <- NA
      qc <- NA
      qRSE <- NA
      qp <- 1
      q_sig <- ""
    }
    
    #Summary stats
    meanw <- round(mean(y, na.rm = TRUE),1)
    m_sig <- round(mRSE(dat = y),3)
    model <- if (min(Binp,Bp,NEp,LMp, qp, na.rm = TRUE)<p) {
      if (Bp<=min(LMp,NEp,Binp, qp, na.rm = TRUE)) {
        "Burr"
      } else {
        if (LMp<=min(Bp,NEp,Binp, qp, na.rm = TRUE)) {
          "Linear"
        } else {
          if (NEp<=min(LMp,Bp,Binp, qp, na.rm = TRUE)) {
            "NegExp"
          } else {
            if (Binp<=min(LMp,Bp,NEp, qp, na.rm = TRUE)) {
              "Binomial"
            } else {
              if (qp<=min(LMp,Bp,NEp, Binp, na.rm = TRUE)) {
                "Quadratic"
              }
            }
          }
        }
      }
    } else {"Mean"}
    model <- if (min(LMRSE, NERSE, BRSE, BinRSE, qRSE, na.rm = TRUE) < m_sig){
      model
    } else {"Mean"}
    
    #Record values
    fitw[nrow(fitw)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMp, k, r, NERSE, NEp,
                             Ba, Bb, BRSE, Bp, Bs, Bsd, Bm, BinRSE, Binp, qa, qb, qc, qRSE, qp, 
                             LM_sig, NE_sig, B_sig, Bin_sig, q_sig,
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
#' @param maxiter The maximum number of iterations for model fitting
#' @return dataframe
#' @export

floraDynamics <- function(dat, thres = 5, pnts = 10, p = 0.01, bTest  = 10, maxiter = 1000,
                          Sr = 0, Sk = 0, Sa = 0, Sb = 0, Sc = 0, 
                          NSr = 0, NSk = 0, NSa = 0, NSb = 0, NSc = 0){
  
  coverChange <- coverDyn(dat, thres = thres, pnts = pnts, p = p, bTest  = bTest, maxiter = maxiter)
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
  
  models <- data.frame('Species'=coverChange$Species, 
                       'Cover' = coverChange$Model, 'C_a' = a, 'C_b' = b, 'C_c' = c, 'C_RSE' = RSE,
                       'Top' = topChange$Model, 'T_a' = aa, 'T_b' = bb, 'T_c' = cc, 'T_RSE' = RSE1,
                       'Base' = baseChange$Model, 'B_a' = aaa, 'B_b' = bbb, 'B_c' = ccc, 'B_RSE' = RSE2, 
                       'He' = he_Change$Model, 'He_a' = aA, 'He_b' = bB, 'He_c' = cC, 'He_RSE' = RSE3, 
                       'Ht' = ht_Change$Model, 'Ht_a' = aT, 'Ht_b' = bT, 'Ht_c' = cT, 'Ht_RSE' = RSEt,  
                       'Width' = w_Change$Model, 'w_a' = aw, 'w_b' = bw, 'w_c' = cw, 'w_RSE' = RSEw)
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
    c <- as.numeric(mods$C_b[n])
  }
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
    c <- as.numeric(mods$T_b[n])
  }
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
    c <- as.numeric(mods$B_b[n])
  }
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
    c <- as.numeric(mods$He_b[n])
  }
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
    c <- as.numeric(mods$Ht_b[n])
  }
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
    c <- as.numeric(mods$w_b[n])
  }
  return(c)
}

#' Arranges survey data into strata using k-means clustering
#' 
#' @param veg A dataframe listing plant species with columns describing crown dimensions
#' @param cols A list of the columns that will be used for classification, and the species name
#' @param nstrat Number of strata. Defaults to 4
#' @return Dataframe
#' @export

stratify <- function(veg, cols, nstrat = 4)
{
  veg_subset <- veg[ , cols]
  veg_subset <- veg_subset[complete.cases(veg_subset), ] # Omit NAs in relevant columns
  veg_subset <- veg_subset %>%
    mutate(base = case_when(veg_subset[,3] == 0 ~ 0.001, TRUE ~ veg_subset[,3]),
           top = log(veg_subset[2]),
           base = log(base))
  df <- scale(veg_subset[, c(4,5)])
  set.seed(123)
  km.res <- kmeans(df, centers = nstrat, nstart = 25)
  clust <- cbind(veg_subset, cluster = km.res$cluster)
  h <- clust[order(clust[,2]),] %>% 
    group_by(cluster) %>% 
    summarise_if(is.numeric, mean)
  
  h <- h[order(h[,2]),] %>% 
    mutate(Stratum = 1:nstrat) %>% 
    select(cluster, Stratum)
  strat <- left_join(clust, h, by = "cluster") %>% 
    select(Species, Stratum) 
  veg <- left_join(veg, strat, by = "Species")
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
#' @param p The maximum allowable p value for a model
#' @param bTest Multiples of mean + mRSE for which Burr & quadratic models can predict 
#' beyond the observed mean + standard deviation
#' @param maxiter The maximum number of iterations for model fitting
#' @return dataframe
#' @export

rich <- function(dat, thres = 5, pnts = 10, p = 0.05) {
  
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
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param cols A list of the columns to be used for stratification, including base & top
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @param nstrat The maximum number of strata
#' @return dataframe

richS <- function(dat, thres = 5, cols, pnts = 10, nstrat = 4) {
  
  spCov <- frame::specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- suppressMessages(left_join(dat, spCov))%>%
    mutate(Species = replace(Species, which(Cover < thres), "Minor Species"))
  
  datS <- stratify(veg = dat, cols = cols, nstrat = nstrat)
  
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
#' 
#' Species that are less common than the set threshold are combined as "Minor Species"
#'
#' @param dat The dataframe containing the input data
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param base Name of the field with the base height
#' @param top Name of the field with the top height
#' @return dataframe
#' @export
#' 
stratSite <- function(dat, thres = 0, pnts = 10, base = "base", top = "top")  {
  strataDet <- data.frame(Stratum = numeric(0), Cover = numeric(0), 
                          Base = numeric(0), Top = numeric(0), stringsAsFactors = F)
  r <- rich(dat, thres = thres, pnts = pnts)
  nstrat <- round(min(4, as.numeric(r$Mean)),0)
  strat <- stratify(dat, cols = c(3:6), nstrat = nstrat)
  for (st in 1:nstrat) {
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
#' @param cols A list of column numbers used for stratifying
#' @param thres The minimum percent cover (0-100) of a Species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param nstrat Maximum number of strata
#' @return dataframe
#' @export
#' 
stratRich <- function(dat, cols, thres = 5, pnts = 10, nstrat = 4) {
  
  richList <- data.frame('S1' = numeric(0), 'S2' = numeric(0),
                         'S3' = numeric(0), 'S4' = numeric(0))
  slist <- unique(dat$Site, incomparables = FALSE)
  
  for (s in slist) {
    datSite <- filter(dat, Site == s)
    sRich <- richS(dat = datSite, cols = cols, thres = thres, 
                   pnts = pnts, nstrat = 4)
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

#' Finds % cover of surveyed Species and groups minor Species
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