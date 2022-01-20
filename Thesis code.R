#One-way RM MANOVA and two-way ANOVA with scores2.csv
#Chris Lienaerts
#U880576

#install and load the packages needed
install.packages("rstatix")
install.packages("jtools")
library(rstatix)
library(gridExtra)
library(ggplot2)
library(jtools)

#set the working directory
setwd("/Users/chris/Documents/Studie/4e jaar/Thesis/data")

#read in the data
data <- read.csv(file = "scores2.csv", colClasses = c("factor", "factor", "factor", "numeric", "numeric"))
df <- data.frame(data)

#code needed for Figure 3 in the results section, retrieved from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

#summarize the data using summarySEwithin
#summary of pres.score
summ_data_pres <- summarySEwithin(df, measurevar="pres.score", withinvars="mod",
                             idvar="part.no", na.rm=FALSE, conf.interval=.95)

#summary of spat.score
summ_data_spat <- summarySEwithin(df, measurevar="spat.score", withinvars="mod",
                                  idvar="part.no", na.rm=FALSE, conf.interval=.95)

# Make separate graphs for pres.score and spat.score, with the 95% confidence interval
#pres.score
ggplot(summ_data_pres, aes(x = factor(mod, level = c('laptop', 'hmd')), y = pres.score, group=1)) +
  geom_line() +
  geom_errorbar(width=.1, aes(ymin=pres.score-ci, ymax=pres.score+ci)) +
  geom_point(shape=21, size=3, fill="white") +
  ylim(90,170) + 
  theme_apa(x.font.size = 25, y.font.size = 25, ) + 
  theme(axis.text = element_text(size=20)) +
  xlab("Modality") + ylab("Presence score")

#spat.score
ggplot(summ_data_spat, aes(x = factor(mod, level = c('laptop', 'hmd')), y = spat.score, group=1)) +
  geom_line() +
  geom_errorbar(width=.1, aes(ymin=spat.score-ci, ymax=spat.score+ci)) +
  geom_point(shape=21, size=3, fill="white") +
  ylim(25,100) + 
  theme_apa(x.font.size = 25, y.font.size = 25, ) + 
  theme(axis.text = element_text(size=20)) + 
  xlab("Modality") + ylab("Spatial score")



#checking the assumptions for the MANOVA
#visualize the data in two separate boxplots, one for presence score and one for spatial score
bxp_pres <- ggplot(data, aes(x=mod, y=pres.score)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
bxp_pres

bxp_spat <- ggplot(data, aes(x=mod, y=spat.score)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)
bxp_spat

#check for outliers
mahalanobis_distance(data = data[, c("pres.score", "spat.score")])$is.outlier

#check for normality
  #method 1
shapiro.test(data$pres.score)
shapiro.test(data$spat.score)

#plot frequency distributions for visually checking normality
hist(df$pres.score, main = "Frequency Distribution Presence score", xlab = "Presence score")
hist(df$spat.score, main = "Frequency Distribution Spatial score", xlab = "Spatial score")

#check for linearity
p1 <- data  %>% group_by(mod) %>% filter(mod == "laptop") %>% ggplot(aes(x = pres.score, y = spat.score)) + geom_point() + ggtitle("Laptop modality")
p2 <- data  %>% group_by(mod) %>% filter(mod == "hmd") %>% ggplot(aes(x = pres.score, y = spat.score)) + geom_point() + ggtitle("HMD modality") 
grid.arrange(p1, p2, ncol=2)

#check for multicollinearity
cor.test(df$pres.score, df$spat.score) #Pearson correlation test

#plot pres.score against spat.score to visually inspect correlation
ggplot(df, aes(x = pres.score, y = spat.score)) + xlab("Presence score") + ylab("Spatial score") + geom_point() + geom_smooth(method = lm, se = FALSE) + theme_apa()

#perform and summarize the MANOVA
man <- manova(cbind(spat.score, pres.score) ~ mod, data = data)
sum_man <- summary(man)
sum_man
lapply(sum_man, eta_squared)
#univariate test statistics
univ <- summary.aov(man)
univ
#effect sizes
lapply(univ, eta_squared)


#two-way ANOVA (additive)
aov2_add <- aov(spat.score ~ mod + pres.score , data = data)
summary(aov2_add)