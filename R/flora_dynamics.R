#' Finds % cover of surveyed species and groups minor species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' species - name of the surveyed species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a species that will be kept single
#' @param pnts The number of points measured in a transect
#' @return dataframe
#' @export

specCover <- function(dat, thres = 5, pnts = 10) {
  
  #List species and ages
  spList <- unique(dat$species, incomparables = FALSE)
  ages <- unique(dat$Age, incomparables = FALSE)
  
  #Create empty summary dataframe
  spCover <- data.frame('Species' = character(0), 'Age' = numeric(0), 'Cover' = numeric(0), stringsAsFactors=F)
  
  #DATA COLLECTION
  for (sp in 1:length(spList)) {
    for (age in ages) {
      spName <- dat %>% filter(species == spList[sp])
      spAge <- spName %>% filter(Age == age)
      
      #Percent cover
      sppnts <- unique(spAge$Point, incomparables = FALSE)
      covSp <- as.numeric(length(sppnts))*(100/pnts)
      
      #Record values
      spCover[nrow(spCover)+1,] <- c(as.character(spList[sp]), as.numeric(age), as.numeric(covSp))
    }
  }
  
  #Group minor species
  spCover$Cover <- as.numeric(as.character(spCover$Cover))
  spShort <- spCover %>%
    group_by(Species) %>%
    summarise_if(is.numeric, mean)
  #List minor species, then rename in dataset
  minor <- spShort %>% filter(Cover < thres)
  minList <- unique(minor$Species, incomparables = FALSE)
  
  for (snew in 1:length(minList)) {
    spCover[spCover==minor$Species[snew]]<-"Minor species"
  }
  return(spCover)
}


#' Builds models for cover dynamics of surveyed species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' species - name of the surveyed species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @return dataframe
#' @export

coverDyn <- function(dat, thres = 5, pnts = 10, p = 0.05) {
  
spCov <- specCover(dat = dat, thres = thres, pnts = pnts)
  priorList <- unique(spCov$Species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitCov <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                       'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_p' = numeric(0),
                       'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_p' = numeric(0),
                       'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_p' = numeric(0),
                       'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'Q_p' = numeric(0),
                       'linear' = character(0), 'NegExp' = character(0), 'Burr' = character(0), 'Binomial' = character(0), 'Quadratic' = character(0), 'Mean' = character(0),
                       'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=100000, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- spCov %>% filter(Species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$Cover)
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
      LM<-lm(y ~ x)
      LMSum <- base::summary(LM)
      LMa <- round(LMSum$coefficients[2],2)
      LMb <- round(LMSum$coefficients[1],2)
      LMRSE <- round(LMSum$sigma,2)
      LMp <- round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
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
      LMp <- NA
      LM_sig <- ""
    }
    
    #Negative exponential
    init1<-c(k=50,r=0.5)
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init2,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
      NESum <- base::summary(NE)
      k <- round(NESum$coefficients[1],2)
      r <- round(NESum$coefficients[2],2)
      NERSE <- round(NESum$sigma,2)
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
      NEp <- NA
      NE_sig <- ""
    }
    
    #Burr
    init2<-c(a=3,b=2)
    if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
      Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
      BSum <- base::summary(Burr)
      Ba <- round(BSum$coefficients[1],2)
      Bb <- round(BSum$coefficients[2],2)
      BRSE <- round(BSum$sigma,2)
      Bp <- round(max(BSum$coefficients[7],BSum$coefficients[8]),5)
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
      Bp <- NA
      B_sig <- ""
    }
    
    #Binomial
    init3<-c(s=1,sd=3, m=20)
    if (!berryFunctions::is.error(nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control))) {
      Bin<-nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control)
      BinSum <- base::summary(Bin)
      Bs <- round(BinSum$coefficients[1],2)
      Bsd <- round(BinSum$coefficients[2],2)
      Bm <- round(BinSum$coefficients[3],2)
      BinRSE <- round(BinSum$sigma,2)
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
      Binp <- NA
      Bin_sig <- ""
    }
    
    #Quadratic
    init4 <- c(a = 1, b = 2, c = 0)
    if (!berryFunctions::is.error(nls(y ~ a*x^2 + b*x + c, data = studySpecies, 
                                      start = init4, trace = T, control = control))) {
      q <- nls(y ~ a*x^2 + b*x + c, data = studySpecies, start = init4, 
               trace = T, control = control)
      qSum <- base::summary(q)
      qa <- round(qSum$coefficients[1], 2)
      qb <- round(qSum$coefficients[2], 2)
      qc <- round(qSum$coefficients[3], 2)
      qRSE <- round(qSum$sigma, 2)
      qp <- round(max(qSum$coefficients[10], qSum$coefficients[11], 
                      qSum$coefficients[12]), 5)
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
      qp <- NA
      q_sig <- ""
    }
    
    #Summary stats
    meanCov <- round(mean(y),1)
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
                                 LM_sig, NE_sig, B_sig, Bin_sig, q_sig,
                                 meanCov, model)
  }
  
  return(fitCov)
}



#' Builds models for top height dynamics of surveyed species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' species - name of the surveyed species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @return dataframe
#' @export

topDyn <- function(dat, thres = 5, pnts = 10, p = 0.05) {
  
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(species = replace(species, which(Cover < thres), "Minor species"))
  
  priorList <- unique(dat$species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitTop <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                       'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_p' = numeric(0),
                       'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_p' = numeric(0),
                       'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_p' = numeric(0),
                       'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'Q_p' = numeric(0),
                       'linear' = character(0), 'NegExp' = character(0), 'Burr' = character(0), 'Binomial' = character(0), 'Quadratic' = character(0), 'Mean' = character(0),
                       'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=100000, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- dat %>% filter(species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$top)
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
      LM<-lm(y ~ x)
      LMSum <- base::summary(LM)
      LMa <- round(LMSum$coefficients[2],2)
      LMb <- round(LMSum$coefficients[1],2)
      LMRSE <- round(LMSum$sigma,2)
      LMp <- round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
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
      LMp <- NA
      LM_sig <- ""
    }
    
    #Negative exponential
    init1<-c(k=50,r=0.5)
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init2,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
      NESum <- base::summary(NE)
      k <- round(NESum$coefficients[1],2)
      r <- round(NESum$coefficients[2],2)
      NERSE <- round(NESum$sigma,2)
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
      NEp <- NA
      NE_sig <- ""
    }
    
    #Burr
    init2<-c(a=3,b=2)
    if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
      Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
      BSum <- base::summary(Burr)
      Ba <- round(BSum$coefficients[1],2)
      Bb <- round(BSum$coefficients[2],2)
      BRSE <- round(BSum$sigma,2)
      Bp <- round(max(BSum$coefficients[7],BSum$coefficients[8]),5)
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
      Bp <- NA
      B_sig <- ""
    }
    
    #Binomial
    init3<-c(s=1,sd=3, m=20)
    if (!berryFunctions::is.error(nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control))) {
      Bin<-nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control)
      BinSum <- base::summary(Bin)
      Bs <- round(BinSum$coefficients[1],2)
      Bsd <- round(BinSum$coefficients[2],2)
      Bm <- round(BinSum$coefficients[3],2)
      BinRSE <- round(BinSum$sigma,2)
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
      Binp <- NA
      Bin_sig <- ""
    }
    
    #Quadratic
    init4 <- c(a = 1, b = 2, c = 0)
    if (!berryFunctions::is.error(nls(y ~ a*x^2 + b*x + c, data = studySpecies, 
                                      start = init4, trace = T, control = control))) {
      q <- nls(y ~ a*x^2 + b*x + c, data = studySpecies, start = init4, 
               trace = T, control = control)
      qSum <- base::summary(q)
      qa <- round(qSum$coefficients[1], 2)
      qb <- round(qSum$coefficients[2], 2)
      qc <- round(qSum$coefficients[3], 2)
      qRSE <- round(qSum$sigma, 2)
      qp <- round(max(qSum$coefficients[10], qSum$coefficients[11], 
                      qSum$coefficients[12]), 5)
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
      qp <- NA
      q_sig <- ""
    }
    
    #Summary stats
    meanTop <- round(mean(y),1)
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
                                 LM_sig, NE_sig, B_sig, Bin_sig, q_sig,
                                 meanTop, model)
  }
  
  return(fitTop)
}


#' Builds models for top-base height allometry dynamics of surveyed species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' species - name of the surveyed species
#' Age - age of the site since the triggering disturbance
#' 
#' Species that are less common than the set threshold are combined as "Minor species"
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) of a species that will be analysed
#' @param pnts The number of points measured in a transect
#' @param p The maximum allowable p value for a model
#' @return dataframe
#' @export

baseDyn <- function(dat, thres = 5, pnts = 10, p = 0.05) {
  
  spCov <- specCover(dat = dat, thres = 0, pnts = pnts)%>%
    group_by(Species)%>%
    summarise_if(is.numeric, mean)
  dat <- left_join(dat, spCov)%>%
    mutate(species = replace(species, which(Cover < thres), "Minor species"),
           bRat = base/top)
  
  priorList <- unique(dat$species, incomparables = FALSE)
  
  #DATA ANALYSIS
  fitBase <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                        'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_p' = numeric(0),
                        'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_p' = numeric(0),
                        'scale' = numeric(0), 'sd' = numeric(0), 'Binm' = numeric(0), 'Bin_sigma' = numeric(0), 'Bin_p' = numeric(0),
                        'Qa' = numeric(0), 'Qb' = numeric(0), 'Qc' = numeric(0), 'Q_sigma' = numeric(0), 'Q_p' = numeric(0),
                        'linear' = character(0), 'NegExp' = character(0), 'Burr' = character(0), 'Binomial' = character(0), 'Quadratic' = character(0), 'Mean' = character(0),
                        'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(priorList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=100000, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- dat %>% filter(species == priorList[SpeciesNumber])
    x <- as.numeric(studySpecies$Age)
    y <- as.numeric(studySpecies$bRat)
    
    #Linear
    if (!berryFunctions::is.error(lm(y ~ x))) {
      LM<-lm(y ~ x)
      LMSum <- base::summary(LM)
      LMa <- round(LMSum$coefficients[2],2)
      LMb <- round(LMSum$coefficients[1],2)
      LMRSE <- round(LMSum$sigma,2)
      LMp <- round(max(LMSum$coefficients[7],LMSum$coefficients[8]),5)
      LM_sig <- if (!berryFunctions::is.error(LMp)) {
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
      LMp <- NA
      LM_sig <- ""
    }
    
    #Negative exponential
    init1<-c(k=50,r=0.5)
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init2,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=studySpecies,start=init1,trace=T)
      NESum <- base::summary(NE)
      k <- round(NESum$coefficients[1],2)
      r <- round(NESum$coefficients[2],2)
      NERSE <- round(NESum$sigma,2)
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
      NEp <- NA
      NE_sig <- ""
    }
    
    #Burr
    init2<-c(a=3,b=2)
    if (!berryFunctions::is.error(nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control))) {
      Burr<-nls(y~a*b*((0.1*x^(a-1))/((1+(0.1*x)^a)^b+1)),data=studySpecies,start=init2,trace=T, control = control)
      BSum <- base::summary(Burr)
      Ba <- round(BSum$coefficients[1],2)
      Bb <- round(BSum$coefficients[2],2)
      BRSE <- round(BSum$sigma,2)
      Bp <- round(max(BSum$coefficients[7],BSum$coefficients[8]),5)
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
      Bp <- NA
      B_sig <- ""
    }
    
    #Binomial
    init3<-c(s=1,sd=3, m=20)
    if (!berryFunctions::is.error(nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control))) {
      Bin<-nls(y~s*(1/(sd*sqrt(2*pi)))*exp(-((x-m)^2)/(2*sd^2)),data=studySpecies,start=init3,trace=T, control = control)
      BinSum <- base::summary(Bin)
      Bs <- round(BinSum$coefficients[1],2)
      Bsd <- round(BinSum$coefficients[2],2)
      Bm <- round(BinSum$coefficients[3],2)
      BinRSE <- round(BinSum$sigma,2)
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
      Binp <- NA
      Bin_sig <- ""
    }
    
    #Quadratic
    init4 <- c(a = 1, b = 2, c = 0)
    if (!berryFunctions::is.error(nls(y ~ a*x^2 + b*x + c, data = studySpecies, 
                                      start = init4, trace = T, control = control))) {
      q <- nls(y ~ a*x^2 + b*x + c, data = studySpecies, start = init4, 
               trace = T, control = control)
      qSum <- base::summary(q)
      qa <- round(qSum$coefficients[1], 2)
      qb <- round(qSum$coefficients[2], 2)
      qc <- round(qSum$coefficients[3], 2)
      qRSE <- round(qSum$sigma, 2)
      qp <- round(max(qSum$coefficients[10], qSum$coefficients[11], 
                      qSum$coefficients[12]), 5)
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
      qp <- NA
      q_sig <- ""
    }
    
    #Summary stats
    meanBase <- round(mean(y),1)
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
    fitBase[nrow(fitBase)+1,] <- c(as.character(priorList[SpeciesNumber]), LMa, LMb, LMRSE, LMp, k, r, NERSE, NEp,
                                   Ba, Bb, BRSE, Bp, Bs, Bsd, Bm, BinRSE, Binp, qa, qb, qc, qRSE, qp, 
                                   LM_sig, NE_sig, B_sig, Bin_sig, q_sig,
                                   meanBase, model)
  }
  
  return(fitBase)
}