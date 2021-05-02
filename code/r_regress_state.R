#Get packages.
#install.packages('tidyverse')
#install.packages('rms')
#install.packages('patchwork')
library(tidyverse) #Various useful functions.
library(rms) #Statistics.
library(lattice) #Statistical plotting for exploratory goal.
library(patchwork) #Plotting.

#Masked functions: src, summarize, format.pval, units, backsolve

#Set parameters manually for interactive running.
roi <- 'Whole' #Roi set name.
width <- '62' #dFC window width label.
k <- '5' #k based on elbow.

#Set code as working directory.

# Loading and Formatting Data ---------------------------------------------

#Load in the matrices for state flexibility score and cognitive score per subject
#and add more meaningful column names.
statepath <- paste('../outputs/r_stateflex/',roi,'/',width,'/',k,'/stateflex.csv',sep='')
cogpath <- '../outputs/cognition/cogscores.csv'
stateflex <- read_csv(statepath,col_names=c('Subject','Flex_Score'))
cogscores <- read_csv(cogpath,col_names=T)
cogscores <- cogscores %>%
  rename(g=gPCA,
         NIH_Fluid=CogFluidComp_Unadj,
         PMAT_Fluid=PMAT24_A_CR,
         NIH_Crystal=CogCrystalComp_Unadj,
         PicVoc_Crystal=PicVocab_Unadj)

#Merges datasets by subject.
flexcog <- stateflex %>%
  inner_join(cogscores)

#Converts subject id to string and gender to a labelled factor.
flexcog <- flexcog %>%
  mutate(Subject=as.character(Subject),
         Gender=factor(Gender,labels=c('Female','Male')))

#Converts age ranges to single numbers for use in analyses in lieu of the actual
#values from HCP. These are a random number between the ranges and in between 
#36-40 for 36+ to fit with the other 4-value separated ranges. For use as a
#x-var predictor for cognitive scores. Sets the seed for reproducibility.
set.seed(12345)
flexpred <- flexcog %>%
  mutate(Age=case_when(
           Age=='36+'~'36-40',
           TRUE ~ Age)
        ) %>% #Converts unbounded to bounded.
  separate(Age,into=c('Low_Bound','High_Bound'),sep='-') %>% #Separates into bounds.
  rowwise() %>% #Makes next verbs act row-wise.
  mutate(Age=sample(Low_Bound:High_Bound,1),
         Low_Bound=NULL,High_Bound=NULL) %>% #Finds the means and drops the bounds.
  ungroup() %>% #Ungroups row-wise groups to prevent unexpected behavior.
  mutate( 
    Range=case_when( 
      Age<=25 ~ 'Younger (22-25)', 
      Age>=26 & Age<=30 ~ 'Middle (26-30)', 
      Age>=31 ~ 'Older (31-36+)') #Merges 31-35 and 36-40, only one subject 36+.
        ) %>% #Reforms the age ranges in addition for visualization later.
  mutate(
    Range=factor(Range,levels=c('Younger (22-25)',
                                'Middle (26-30)',
                                'Older (31-36+)'))) %>% #Converts ranges to factors.
  select(Subject,Flex_Score,Gender,Age,Range,everything()) #Reorders for readability.

# Finding the Models -------------------------------------------------------

#Produce data summaries for each variable which are useful to rms.
flexpred.dd <- datadist(flexpred)
options(datadist='flexpred.dd')

#1. Conduct ANOVA F-test of overall regression for models including all 
#possible predictors (state flexibility, age, gender, all interactions) chosen 
#to determine whether there exists significant coefficients within (statistical 
#test of the fit of your model compared to intercept only model, Ha = at least 
#one coefficient != 0). 
anova(ols(g~Flex_Score*Age*Gender,data=flexpred)) #g-factor, p ~0.2097
anova(ols(NIH_Fluid~Flex_Score*Age*Gender,data=flexpred)) #NIH Fluid, p ~0.0138
anova(ols(PMAT_Fluid~Flex_Score*Age*Gender,data=flexpred)) #PMAT Fluid, p ~0.0283
anova(ols(NIH_Crystal~Flex_Score*Age*Gender,data=flexpred)) #NIH Crystal, p ~0.3419
anova(ols(PicVoc_Crystal~Flex_Score*Age*Gender,data=flexpred)) #PicVoc Crystal, p ~0.2901

#Indicates that only fluid intelligence models are significant.

#2. For the significant full models, use the ANOVA overall regression test to 
#determine significance of variations of reduced models and use adjusted R2 
#(adj R2) to determine the model which best fits the data.

#NIH fluid intelligence.

#Create a list of reducing variations of ordinary least squares fitted models.
f.NIH.olslist <- list(
  ols(NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Flex_Score:Gender+Age:Gender+Flex_Score:Age:Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Flex_Score:Gender+Age:Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Flex_Score:Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Age:Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Gender+Age:Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Age,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Gender+Age:Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Flex_Score:Age,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Gender+Flex_Score:Gender,data=flexpred),
  ols(NIH_Fluid~Age+Gender+Age:Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age+Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Age,data=flexpred),
  ols(NIH_Fluid~Flex_Score+Gender,data=flexpred),
  ols(NIH_Fluid~Age+Gender,data=flexpred),
  ols(NIH_Fluid~Flex_Score,data=flexpred),
  ols(NIH_Fluid~Age,data=flexpred),
  ols(NIH_Fluid~Gender,data=flexpred)
)

#Apply ANOVA to the list of ols-fitted models and save the output to a list.
f.NIH.anovalist = list()
for (i in seq_along(f.NIH.olslist)) {
  f.NIH.anovalist[[i]] <- anova(f.NIH.olslist[[i]])
}

#Print out all the ANOVA overall regression tests from the list object to find 
#significance.
f.NIH.anovalist

#Print out all the adj R2 values from the list object to compare fit.
f.NIH.olslist

#Select the highest adj R2 model.

#p-value adj R2
#1 0.0138 0.107
#2 0.0107 0.107
#3 0.0420 0.067
#4 0.0051 0.116 <<< Flex_Score+Age+Gender+Flex_Score:Age+Age:Gender
#5 0.0217 0.083
#6 0.0208 0.076
#7 0.2007 0.021
#8 0.0158 0.082
#9 0.0453 0.051
#10 0.3472 0.003
#11 0.0066 0.091
#12 0.2035 0.017
#13 0.2774 0.006
#14 0.6063 -0.010
#15 0.1052 0.026
#16 0.7619 -0.009
#17 0.1106 0.016
#18 0.3648 -0.002

#PMAT24 fluid intelligence.

#Create a list of reducing variations of ordinary least squares fitted models.
f.PMAT.olslist <- list(
  ols(PMAT_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Flex_Score:Gender+Age:Gender+Flex_Score:Age:Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Flex_Score:Gender+Age:Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Flex_Score:Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Age:Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Gender+Flex_Score:Gender+Age:Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Gender+Flex_Score:Age,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Gender+Flex_Score:Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Gender+Age:Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Flex_Score:Age,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Gender+Flex_Score:Gender,data=flexpred),
  ols(PMAT_Fluid~Age+Gender+Age:Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age+Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Age,data=flexpred),
  ols(PMAT_Fluid~Flex_Score+Gender,data=flexpred),
  ols(PMAT_Fluid~Age+Gender,data=flexpred),
  ols(PMAT_Fluid~Flex_Score,data=flexpred),
  ols(PMAT_Fluid~Age,data=flexpred),
  ols(PMAT_Fluid~Gender,data=flexpred)
)

#Apply ANOVA to the list of ols-fitted models and save the output to a list.
f.PMAT.anovalist = list()
for (i in seq_along(f.PMAT.olslist)) {
  f.PMAT.anovalist[[i]] <- anova(f.PMAT.olslist[[i]])
}

#Print out all the ANOVA overall regression tests from the list object to find 
#significance.
f.PMAT.anovalist

#Print out all the adj R values from the list object to compare fit.
f.PMAT.olslist

#Select the highest adj R2 model.

#p-value adj R2
#1 0.0283 0.088
#2 0.0291 0.082
#3 0.0149 0.092
#4 0.0277 0.077
#5 0.0148 0.092
#6 0.0132 0.086
#7 0.0067 0.101 
#8 0.0153 0.083
#9 0.0242 0.065
#10 0.0035 0.104 <<< Flex_Score+Gender+Flex_Score:Gender
#11 0.0474 0.050
#12 0.0063 0.093
#13 0.0130 0.067
#14 0.0037 0.091
#15 0.0208 0.058
#16 0.0139 0.051
#17 0.075 0.022
#18 0.0115 0.054

#3. Analyze residual plots to check for assumptions for linear regression.
#Plots for exploratory use, will not save these plots for now.

#Save the ideal models into variables, for NIH and PMAT24 fluid intelligence.
f.NIH.ideal.ols <- ols(NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Age:Gender,
                       data=flexpred)
f.PMAT.ideal.ols <- ols(PMAT_Fluid~Flex_Score+Gender+Flex_Score:Gender,
                        data=flexpred)

#NIH fluid intelligence. 

#QQ plot: Data mostly fits line except for extreme values indicating normality.
qqnorm(resid(f.NIH.ideal.ols),main='Normal Q-Q Plot for Residuals')
qqline(resid(f.NIH.ideal.ols))

#Residuals vs fitted values: variability of residuals seem mostly consistent across
#fitted values aroudn the regression line, indicating constant variance.
xyplot(resid(f.NIH.ideal.ols)~fitted(f.NIH.ideal.ols),xlab='Fitted Values',
       ylab='Residuals',main="Residuals versus Fitted Values",type=c('p','r'))

#PMAT24 fluid intelligence.

#QQ plots: Data mostly fits line except for extreme values indicating normality
#for both models.
qqnorm(resid(f.PMAT.ideal.ols),main='Normal Q-Q Plot for Residuals')
qqline(resid(f.PMAT.ideal.ols))

#Residuals vs fitted values: variability of residuals seem mostly consistent across
#fitted values around the regression line for both models, indicating constant variance.
xyplot(resid(f.PMAT.ideal.ols)~fitted(f.PMAT.ideal.ols),xlab="Fitted Values",
       ylab="Residuals",main="Residuals versus Fitted Values",type=c('p','r'))

#Ideal models:
#NIH fluid intelligence: NIH_Fluid~Flex_Score+Age+Gender+Flex_Score:Age+Age:Gender
#PMAT24 fluid intelligence: PMAT_Fluid~Flex_Score+Gender+Flex_Score:Gender

# Plotting the Data -----------------------------------------------------------

#Plots to examine the data visually, sequentially add other information to
#separate and analyze the relationship between the cognitive score
#and state flexibility. Will format and save plots for later use. To save time
#and space, won't comment repeats of code which generally do the same thing.

#Set theme to light.
theme_set(theme_light())

#Plots to examine NIH Fluid Intelligence.
#1. Base plot for NIH fluid score-flexibility. 
#NIH fluid score does not seem to have a relationship with state flexibility by 
#itself based on the distance of points from the regression line and the lack of 
#a slope in the regression line nor data points. 
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Flex_Score)) + #Defines equation.
  geom_point(size=2) + #Points.
  stat_smooth(method='lm',formula=y~x,fullrange=T) + #Regression and error.
  scale_x_continuous(breaks=seq(0,250,by=50)) + #x-axis ticks for less crowding.
  coord_cartesian(ylim=c(75,150)) + #y-axis limits to encompass data.
  ylab('NIH Fluid Score') +
  xlab('State Flexibility Score') +
  ggtitle('NIH Fluid Score-State Flexibility') #Labels and title.
p1

#Save graph in specified location, in inch units, with width and heights of 8.
ggsave('../outputs/r_statefigures/NIH_Fluid~Flex.pdf',units='in',width=8,height=8)

#2. Base plot for NIH fluid score-age. 
#NIH fluid score does not seem to have a relationship with age by itself based 
#on the distance of points from the regression line and the lack of a slope 
#in the regression line or data points. 
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Age)) + 
  geom_point(size=2) + 
  stat_smooth(method='lm',formula=y~x,fullrange=T) + 
  scale_x_continuous(breaks=seq(21,41,by=6)) + 
  coord_cartesian(ylim=c(75,150)) + 
  ylab('NIH Fluid Score') +
  xlab('Age') +
  ggtitle('NIH Fluid Score-Age') 
p1
ggsave('../outputs/r_statefigures/NIH_Fluid~Age.pdf',units='in',width=8,height=8)

#3. Base plot for NIH fluid score-gender. 
#NIH fluid score does not seem to differ between the genders by itself based on 
#the lack of a difference between the values of the box plots. 
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Gender,fill=Gender)) + #Defines equation and color separation.
  geom_boxplot() + #Box plot.
  ylab('NIH Fluid Score') +
  xlab('Gender') +
  ggtitle('NIH Fluid Score-Gender') #Labels and title.
p1
ggsave('../outputs/r_statefigures/NIH_Fluid~Gender.pdf',units='in',width=8,height=8)

#4. NIH fluid score-flexibility separated by the three main age ranges in the data. 
#The relationship between NIH fluid score and state flexibility seems to change 
#with greater age - there is a negative relationship at younger ages, a 
#moderately positive one at the middle, and a positive one when older with
#the error for younger ages being farther from the other two indicating an actual
#difference. This supports the two-way interaction for age and state flexibility
#in the ideal model.
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Flex_Score,fill=Range,colour=Range)) + #Defines equation and color separation.
  geom_point(size=2) + #Points.
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + #Regression and error.
  xlab('State Flexibility Score') +
  ggtitle('NIH Fluid Score-State Flexibility',subtitle='Age Separated') + #Labels and title.
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='None') #Get rid of redundant text in axes and legend.
p1
p2 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Flex_Score,fill=Range,colour=Range)) + #Defines equation and color separation.
  geom_point(size=2) + #Points.
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + #Regression and error.
  facet_grid(cols=vars(Range)) + #Separates by variable into adjacent graphs.
  scale_x_continuous(breaks=seq(0,250,by=50)) + #x-axis ticks for less crowding.
  coord_cartesian(ylim=c(75,150)) + #y-axis limits to encompass data.
  theme(strip.background=element_rect(fill='grey90'), 
        strip.text=element_text(colour='black'), #Converts title strip aesthetics.
        legend.position='None', #Removes redundant label.
        axis.title.y=element_text(hjust=1.40)) + #Moves y-axis label to common area.
  xlab('State Flexibility Score') +
  ylab('NIH Fluid Score') #Label and title.
p2
p1/p2 #Combine graphs vertically.
ggsave('../outputs/r_statefigures/NIH_Fluid~Flex~AgeSep.pdf',units='in',width=8,height=8)

#5. NIH fluid score-flexibility separated by gender. 
#The relationship between NIH fluid score and state flexibility seems to only
#slightly differ between the genders - there is a positive relationship for 
#females, and a negative relationship for males, but the regression lines are 
#mostly within each others' errors indicating that the slopes are not different. 
#This might suggest that there is no two-way interaction for age and gender.
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Flex_Score,fill=Gender,color=Gender)) + 
  geom_point(size=2) + 
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + 
  scale_x_continuous(breaks=seq(0,250,by=50)) + 
  ggtitle('NIH Fluid Score-State Flexibility',subtitle='Gender Separated') + 
  theme(legend.position='None',
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 
p1
p2 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Flex_Score,fill=Gender,color=Gender)) + 
  geom_point(size=2) + 
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + 
  scale_x_continuous(breaks=seq(0,250,by=50)) +
  coord_cartesian(ylim=c(75,150)) +
  facet_wrap(~Gender) +
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black'),
        legend.position='None',
        axis.title.y=element_text(hjust=1.40)) +
  xlab('State Flexibility Score') +
  ylab('NIH Fluid Score') 
p2
p1/p2
ggsave('../outputs/r_statefigures/NIH_Fluid~Flex~GenSep.pdf',units='in',
       width=8,height=8)

#6. NIH fluid score-flexibility separated by the three age ranges in rows and 
#gender in columns. 
#The relationship between the NIH fluid score and state flexibility seems to 
#differ across age but not gender - from young to old there seems to be a 
#shift from a negative to a positive relationship in both genders, but there
#does not seem to be a difference between genders. This supports the existence
#of an flexibility-age interaction but not flexibility-gender, fitting with
#the ideal model. To make this more clear, we can collapse by gender and age to 
#compare the orientation of regression lines.
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Flex_Score,fill=Gender,colour=Gender)) + #Color separate by second variable.
  geom_point(size=2) + #Points.
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + #Regression and error.
  facet_grid(rows=vars(Range),cols=vars(Gender)) + #Separate by second and third variable.
  scale_x_continuous(breaks=seq(0,250,by=50)) + #x-axis ticks for less crowding.
  coord_cartesian(ylim=c(75,150)) + #y-axis limits to encompass data.
  xlab('State Flexibility Score') +
  ylab('NIH Fluid Intelligence Score') +
  ggtitle('NIH Fluid Score-State Flexibility',subtitle='Age-Gender Separated') + #Labels and title.
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black'),
        legend.position='None') #Converts title strip aesthetics and removes redundant legend.
p1
ggsave('../outputs/r_statefigures/NIH_Fluid~Flex~AgeGenSep.pdf',units='in',width=8,height=8)

#7. NIH fluid score-flexibility separated by the three age ranges in rows and gender 
#overlapped. 
#When overlapping gender regression lines, the lack of a difference in slope 
#between the genders becomes more evident, while maintaining the difference
#in slopes between ages - matching previous conclusions.
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Flex_Score,colour=Gender,fill=Gender)) + #Color separate by second variable.
  geom_point(size=2) + #Points.
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + #Regression and error.
  facet_grid(rows=vars(Range)) + #Separate by third variable.
  scale_x_continuous(breaks=seq(0,250,by=50)) + #x-axis ticks for less crowding.
  coord_cartesian(ylim=c(75,150)) + #y-axis limits to encompass data.
  xlab('State Flexibility Score') +
  ylab('NIH Fluid Intelligence Score') +
  ggtitle('NIH Fluid Score-State Flexibility',subtitle='Age-Gender Separated, Gender Overlapped') + #Labels and title.
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black')) #Converts title strip aesthetics.
p1
ggsave('../outputs/r_statefigures/NIH_Fluid~Flex~AgeSepGenOver.pdf',units='in',
       width=8,height=8)

#8. NIH fluid score-flexibility separated by gender and the three age ranges  
#overlapped. 
#When overlapping age regression lines, we again see that the same result as the
#previous plot in a different view. We can see that different age ranges have
#different slopes but different genders show the same slopes.
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Flex_Score,colour=Range,fill=Range)) + #Color separate by second variable.
  geom_point(size=2) + #Points.
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + #Regression and error.
  facet_grid(rows=vars(Gender)) + #Separate by third variable.
  scale_x_continuous(breaks=seq(0,250,by=50)) + #x-axis ticks for less crowding.
  coord_cartesian(ylim=c(75,150)) + #y-axis limits to encompass data.
  xlab('State Flexibility Score') +
  ylab('NIH Fluid Intelligence Score') +
  ggtitle('NIH Fluid Score-State Flexibility',subtitle='Age-Gender Separated, Age Overlapped') + #Labels and title.
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black')) #Converts title strip aesthetics.
p1
ggsave('../outputs/r_statefigures/NIH_Fluid~Flex~AgeOverGenSep.pdf',units='in',
       width=8,height=8)

#9. NIH fluid score-age separated by gender.
#In the model, there was an age-gender interaction term which might be of 
#peripheral interest. We can examine that with this graph. It seems that there 
#is a negative relationship between age and NIH fluid score for females, but 
#almost no relationship (slope ~0) between age and NIH fluid score for males. 
#This supports the existence of the age-gender interaction, fitting with the 
#ideal model.
p1 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Age,fill=Gender,color=Gender)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) +
  ggtitle('NIH Fluid Score-Age',subtitle='Gender Separated') +
  theme(legend.position='None',axis.title.x=element_blank(),
        axis.title.y=element_blank())
p1
p2 <- flexpred %>%
  ggplot(aes(y=NIH_Fluid,x=Age,fill=Gender,color=Gender)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + 
  coord_cartesian(ylim=c(75,150)) +
  facet_wrap(~Gender) +
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black'),
        legend.position='None',
        axis.title.y=element_text(hjust=1.40)) +
  xlab('Age') +
  ylab('NIH Fluid Score') 
p2
p1/p2
ggsave('../outputs/r_statefigures/NIH_Fluid~Age~GenSep.pdf',units='in',
       width=8,height=8)

#Plots to examine PMAT24 Fluid Intelligence.
#1. Base plot for PMAT24 fluid score-flexibility. 
#PMAT24 fluid score seems to have a positive relationship with state flexibility 
#by itself based on the distance of points from the regression line and the trend 
#of the regression line and data points. 
p1 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Flex_Score)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T) + 
  scale_x_continuous(breaks=seq(0,250,by=50)) +
  coord_cartesian(ylim=c(0,30)) +
  ylab('PMAT24 Fluid Score') +
  xlab('State Flexibility Score') +
  ggtitle('PMAT24 Fluid Score-State Flexibility')
p1
ggsave('../outputs/r_statefigures/PMAT24_Fluid~Flex.pdf',units='in',width=8,height=8)

#2. Base plot for PMAT fluid score-age. 
#PMAT24 fluid score does not seem to have a relationship with age by itself 
#based on the distance of points from the regression line and the trend of the 
#regression line and data points, though there seems to be a slight negative trend. 
p1 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Age)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T) + 
  scale_x_continuous(breaks=seq(21,41,by=6)) +
  coord_cartesian(ylim=c(0,30)) +
  ylab('PMAT24 Fluid Score') +
  xlab('Age') +
  ggtitle('PMAT24 Fluid Score-Age')
p1
ggsave('../outputs/r_statefigures/PMAT24_Fluid~Age.pdf',units='in',width=8,height=8)

#3. Base plot for NIH fluid score-gender. 
#NIH fluid score seems to differ between the genders, with the medians almost 
#being farther than the Q3-Q1 of each other (higher scores for males compared
#to females), and the distribution for male being left skewed compared to a 
#symmetric one for females. 
p1 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Gender,fill=Gender)) +
  geom_boxplot() +
  ylab('PMAT24 Fluid Score') +
  xlab('Gender') +
  ggtitle('PMAT24 Fluid Score-Gender')
p1
ggsave('../outputs/r_statefigures/PMAT24_Fluid~Gender.pdf',units='in',width=8,height=8)

#4. PMAT24 fluid score-flexibility separated by the three age ranges. 
#The relationship between PMAT24 fluid score and state flexibility seems to slightly
#change with greater age - the trend is similar to for NIH fluid score, but the 
#errors are much closer. This does not support a flexibility-age interaction as much,
#fitting the ideal model.
p1 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Flex_Score,fill=Range,colour=Range)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) +
  xlab('State Flexibility Score') +
  ggtitle('PMAT24 Fluid Score-State Flexibility',subtitle='Age Separated') +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='None')
p1
p2 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Flex_Score,fill=Range,colour=Range)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) +
  facet_grid(cols=vars(Range)) + 
  scale_x_continuous(breaks=seq(0,250,by=50)) +
  coord_cartesian(ylim=c(0,30)) +
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black'),
        legend.position='None',
        axis.title.y=element_text(hjust=1.35)) +
  xlab('State Flexibility Score') +
  ylab('PMAT24 Fluid Score') 
p2
p1/p2
ggsave('../outputs/r_statefigures/PMAT24_Fluid~Flex~AgeSep.pdf',units='in',width=8,height=8)

#5. PMAT24 fluid score-flexibility separated by gender. 
#The relationship between PMAT24 fluid score and state flexibility seems to differ 
#between the genders - there is greater positive relationship but a lower magnitude
#for females than males who seem to have a marginally positive relationship but a 
#higher magnitude, with errors being relatively far apart for lower flexibility
#scores. This is suggests a two-way flexibility-gender interaction, which supports
#the ideal model. 
p1 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Flex_Score,fill=Gender,color=Gender)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) +
  scale_x_continuous(breaks=seq(0,250,by=50)) +
  ggtitle('PMAT24 Fluid Score-State Flexibility',subtitle='Gender Separated') +
  theme(legend.position='None',axis.title.x=element_blank(),
        axis.title.y=element_blank())
p1
p2 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Flex_Score,fill=Gender,color=Gender)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) + 
  scale_x_continuous(breaks=seq(0,250,by=50)) +
  coord_cartesian(ylim=c(0,30)) +
  facet_wrap(~Gender) +
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black'),
        legend.position='None',
        axis.title.y=element_text(hjust=1.35)) +
  xlab('State Flexibility Score') +
  ylab('PMAT24 Fluid Score') 
p2
p1/p2
ggsave('../outputs/r_statefigures/PMAT24_Fluid~Flex~GenSep.pdf',units='in',width=8,height=8)

#6. PMAT24 fluid score-flexibility separated by the three age ranges in rows and 
#gender in columns. 
#When looking at the regression lines and error, there seems to be a difference 
#in slope between genders for the middle ages and not the other ages. This 
#suggests a three-way interaction, which doesn't fit the ideal model found.
#One explanation could be the lack of data for regression lines for graphs like
#this, when data is split into multiple graphs. Supporting this, there is the lack 
#of a gender difference for younger and older age ranges, rather than a different 
#set of gender slopes. To clarify, we could analyze the addition of this term 
#and collapse by gender and age to examine the orientation of regression lines.
p1 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Flex_Score,fill=Gender,colour=Gender)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) +
  facet_grid(rows=vars(Range),cols=vars(Gender)) + 
  scale_x_continuous(breaks=seq(0,250,by=50)) +
  coord_cartesian(ylim=c(0,30)) +
  xlab('State Flexibility Score') +
  ylab('PMAT24 Fluid Intelligence Score') +
  ggtitle('PMAT24 Fluid Score-State Flexibility',subtitle='Age-Gender Separated') +
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black'),
        legend.position='None')
p1
ggsave('../outputs/r_statefigures/PMAT24_Fluid~Flex~AgeGenSep.pdf',units='in',
       width=8,height=8)

#7. PMAT24 fluid score-flexibility separated by the three age ranges in rows and gender. 
#With the gender lines overlapped, it makes the same result clearer - it does look
#like there is an interaction where the middle age ranges show a greater 
#negative relationship for females compared to males with less overlapping errors 
#and no differences in other age ranges, despite what R2 adj says about 
#including the three-way interaction term.
p1 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Flex_Score,colour=Gender,fill=Gender)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) +
  facet_grid(rows=vars(Range)) + 
  scale_x_continuous(breaks=seq(0,250,by=50)) +
  coord_cartesian(ylim=c(0,30)) +
  xlab('State Flexibility Score') +
  ylab('PMAT24 Fluid Intelligence Score') +
  ggtitle('PMAT24 Fluid Score-State Flexibility',subtitle='Age-Gender Separated, Gender Overlapped') +
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black'))
p1
ggsave('../outputs/r_statefigures/PMAT24_Fluid~Flex~AgeSepGenOver.pdf',units='in',
       width=8,height=8)

#8. PMAT24 fluid score-flexibility separated by gender and the three age ranges  
#overlapped. 
#With the age lines overlapped, the gender difference for middle ages is also
#evident from this perspective - with the same conclusions. Getting more subjects'
#data might help show the true relationship.
p1 <- flexpred %>%
  ggplot(aes(y=PMAT_Fluid,x=Flex_Score,colour=Range,fill=Range)) +
  geom_point(size=2) +
  stat_smooth(method='lm',formula=y~x,fullrange=T,alpha=0.15) +
  facet_grid(rows=vars(Gender)) + 
  scale_x_continuous(breaks=seq(0,250,by=50)) +
  coord_cartesian(ylim=c(0,30)) +
  xlab('State Flexibility Score') +
  ylab('PMAT24 Fluid Intelligence Score') +
  ggtitle('PMAT24 Fluid Score-State Flexibility',subtitle='Age-Gender Separated, Age Overlapped') +
  theme(strip.background=element_rect(fill='grey90'),
        strip.text=element_text(colour='black'))
p1
ggsave('../outputs/r_statefigures/PMAT24_Fluid~Flex~AgeOverGenSep.pdf',units='in',width=8,height=8)

# Interpreting the Models --------------------------------------------------

#In the context of those graphs and the best models, we can put values on the
#relationship between fluid intelligence and state flexibility.

#ANOVA for state flexibility coefficients.

#For the model for NIH fluid intelligence, the p-value for the coefficient 
#for state flexibility is ~0.10 holding all other variables constant, which is 
#insignificant at a 0.05 level.
anova(f.NIH.ideal.ols)

#For the model for PMAT24 fluid intelligence, the p-value for the coefficient for
#state flexibility is ~0.027 holding all other variables constant, which is 
#significant at a 0.05 level.
anova(f.PMAT.ideal.ols)

#Effect for different levels.

#For NIH fluid intelligence. Flex_Score+Age+Gender+Flex_Score:Age+Age:Gender
#Specify x-axis range to scale effect value, specify middle of age range to
#represent that age range, need to specify due to interaction.

#For both males and females in the middle of the younger age range, for every 
#100 point increase in resting-state state flexibility score, there is a ~10.74 
#decrease in NIH fluid intelligence score. 
summary(f.NIH.ideal.ols,Flex_Score=c(0,100),Age=23.5) 

#For both males and females in the middle of the middle age range, for every 
#100 point increase in resting-state state flexibility score, there is a ~2.70 
#decrease in NIH fluid intelligence score.
summary(f.NIH.ideal.ols,Flex_Score=c(0,100),Age=28)

#For both males and females in the middle of the older age range, for every 
#100 point increase in resting-state state flexibility score, there is a ~7.12 
#increase in NIH fluid intelligence score.
summary(f.NIH.ideal.ols,Flex_Score=c(0,100),Age=33.5)

#For PMAT24 fluid intelligence. Flex_Score+Gender+Flex_Score:Gender
#Specify x-axis range to scale effect value, specify gender, need to specify due 
#to interaction.

#For males of all age ranges, for every 100 point increase in resting-state
#state flexibility score, there is a ~0.72 increase in PMAT24 fluid intelligence
#score.
summary(f.PMAT.ideal.ols,Flex_Score=c(0,100),Gender='Male')

#For females of all age ranges, for every 100 point increase in resting-state
#state flexibility score, there is a ~4.30 increase in PMAT24 fluid intelligence
#score.
summary(f.PMAT.ideal.ols,Flex_Score=c(0,100),Gender='Female')

