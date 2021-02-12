#' Builds models for cover dynamics of surveyed species
#'
#' Input table requires the following fields:
#' Point - numbered point in a transect
#' species - name of the surveyed species
#' Age - age of the site since the triggering disturbance
#'
#' @param dat The dataframe containing the input data,
#' @param thres The minimum percent cover (0-100) for following analyses
#' @param pnts The number of points measured in a transect
#' @return dataframe
#' @export

coverDyn <- function(dat, thres = 5, pnts = 10) {
  
  #List species and ages
  spList <- unique(dat$species, incomparables = FALSE)
  ages <- unique(dat$Age, incomparables = FALSE)
  
  #Create empty summary dataframe
  spCov <- data.frame('Species' = character(0), 'Age' = numeric(0), 'Cover' = numeric(0), stringsAsFactors=F)
  
  #DATA COLLECTION
  for (sp in 1:length(spList)) {
    for (age in ages) {
      spName <- dat %>% filter(species == spList[sp])
      spAge <- spName %>% filter(Age == age)
      
      #Percent cover
      sppnts <- unique(spAge$Point, incomparables = FALSE)
      covSp <- as.numeric(length(sppnts))*(100/pnts)
      
      #Record values
      spCov[nrow(spCov)+1,] <- c(as.character(spList[sp]), age, covSp)
    }
  }
  
  #DATA ANALYSIS
  fitCov <- data.frame('Species' = character(0), 'lin_a' = numeric(0), 'lin_b' = numeric(0),'lin_Sigma' = numeric(0), 'lin_p' = numeric(0),
                       'k' = numeric(0), 'r' = numeric(0), 'NE_sigma' = numeric(0), 'NE_p' = numeric(0),
                       'Ba' = numeric(0), 'Bb' = numeric(0), 'B_sigma' = numeric(0), 'B_p' = numeric(0),
                       'linear' = character(0), 'NegExp' = character(0), 'Burr' = character(0), 'Mean' = character(0),
                       'Status' = character(0), 'Model' = character(0), stringsAsFactors=F)
  
  for (sp in 1:length(spList)) {
    
    SpeciesNumber <- sp
    control=nls.control(maxiter=100000, tol=1e-7, minFactor = 1/999999999)
    studySpecies <- spCov %>% filter(Species == spList[SpeciesNumber])
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
    if (!berryFunctions::is.error(nls(y~k * (1-exp(-r*x)),data=Dat,start=init2,trace=T))) {
      NE<-nls(y~k * (1-exp(-r*x)),data=Dat,start=init1,trace=T)
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
    
    #Summary stats
    meanCov <- round(mean(y),1)
    status <- if (meanCov >= thres) {
      "common"
    } else {
      ""
    }
    model <- if (meanCov >= thres) {
      if (min(Bp,NEp,LMp, na.rm = TRUE)<0.05) {
        if (Bp<=min(LMp,NEp, na.rm = TRUE)) {
          "Burr"
        } else {
          if (LMp<=min(Bp,NEp, na.rm = TRUE)) {
            "Linear"
          } else {
            if (NEp<=min(LMp,Bp, na.rm = TRUE)) {
              "NegExp"
            }
          }
        }
      } else {"Mean"}
    } else {""}
    
    #Record values
    fitCov[nrow(fitCov)+1,] <- c(as.character(spList[SpeciesNumber]), LMa, LMb, LMRSE, LMp, k, r, NERSE, NEp,
                                 Ba, Bb, BRSE, Bp, LM_sig, NE_sig, B_sig, meanCov, status, model)
  }
  
  return(fitCov)
}