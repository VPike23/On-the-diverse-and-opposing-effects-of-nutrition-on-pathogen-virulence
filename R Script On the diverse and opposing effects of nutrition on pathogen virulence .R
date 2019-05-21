### Sample of R code used in our analysis 
### From: On the diverse and opposing effects of nutrition on pathogen virulence
### Authors: Victoria L. Pike, Katrina A. Lythgoe and Kayla C. King
### 19/2/19
### Refs: Assink, M., & Wibbelink, C. J. M. (2016). Fitting three-level meta-analytic models in R: A step-by-step tutorial. The Quantitative Methods for Psychology, 12(3), 154-174. DOI: 10.20982/tqmp.12.3.p154 


##### Setting up the data frame with variances and calculating cohens d/hedges g ----

#### Step 1 - loading required package ----

library(dplyr)
library(meta)
library(ggplot2)
library(dplyr)
library(metafor)
library(cowplot)


#### Step 2 - loading the data (Analysing Revision data) ----

QLmod <- read.csv("QualityRevisionDataCSV.csv", header = TRUE)


#### Step 2.5 adding column of effect size ID ----

QLmod$effect_size_ID <- 1:33
head(QLmod)


#### Step 3 - calculate Cohen's D ----

# calculate the numberator for the pooled standard deviation:

sd_pooled_n <- ((QLmod$N_for_M1-1)*(QLmod$M1_SD^2)) + ((QLmod$N_for_M2-1)*(QLmod$M2_SD^2))
str(sd_pooled_n)

# calculate the denominator for the pooled standard deviation:

sd_pooled_d <- QLmod$N_for_M1 + QLmod$N_for_M2 - 2
str(sd_pooled_d)

# calculate the pooled sample standard deviation:

sd_pooled <- sqrt(sd_pooled_n / sd_pooled_d)
str(sd_pooled)

# add the pooled standard deviation as a column (variable) in the QLmod dataframe:

QLmod$sd_pooled <- sd_pooled
str(QLmod)

# calculate Cohen's D:

QLmod$Cohen_D <- (QLmod$M2 - QLmod$M1) / QLmod$sd_pooled
str(QLmod)

#### Step 4 - calculate the variance and standard error of Cohen's D ----

QLmod$var_Cohen_D <- (QLmod$N_for_M1 + QLmod$N_for_M2)/(QLmod$N_for_M1 * QLmod$N_for_M2) + ((QLmod$Cohen_D^2)/(2*(QLmod$N_for_M1 + QLmod$N_for_M2)))
str(QLmod)  

QLmod$se_Cohen_D <- sqrt(QLmod$var_Cohen_D)  
str(QLmod)  

#### Step 5 - calculate Hedge's G ----

# calculating the J correction factor
# J = 1 - 3 / (4df - 1)
# df = n1 + n2 - 2

QLmod$df_pooled <- QLmod$N_for_M1 + QLmod$N_for_M2 - 2
QLmod$J <- 1 - 3/(4*QLmod$df_pooled - 1)
str(QLmod$J)

# all values of J should: 0 < J < 1

summary(QLmod$J)

# calculating Hedge's g:

QLmod$Hedges_G <- QLmod$J * QLmod$Cohen_D 
str(QLmod$Hedges_G)

#### Step 6 - calculate the variance and stanard error in Hedge's G ----

# variance of Hedge's G =  J^2 * Cohen_D:

QLmod$var_Hedges_G <- QLmod$J^2 * QLmod$var_Cohen_D
str(QLmod$var_Hedges_G)  

# standard error of Hedge's G:

QLmod$se_Hedges_G <- sqrt(QLmod$var_Hedges_G)
str(QLmod$se_Hedges_G)

QLmod

#### Calculating the Overall effect size ----

#### Estimate the overall effect by fitting an intercept-only model:

overall <- rma.mv(Hedges_G, var_Hedges_G, random = list(~ 1 | effect_size_ID, ~ 1 | Reference), tdist=
                    TRUE, data = QLmod)

#### Request a print of the results stored in the object ‘‘overall’’ in three digits:

summary(overall, digits=3)

#### Making a Graph of overall effect sizes ----

# First restructure data frame in terms of hedges g size: 

QLmod <- arrange(QLmod, Hedges_G) 


QLmod$Parasite <- factor(QLmod$Parasite, levels = rev(unique(QLmod$Parasite)), ordered=TRUE)


cols_subtype <- c("green","navyblue","black", "purple")

Fig1a <- ggplot(QLmod, aes(y=reorder(Parasite, Hedges_G ), x = Hedges_G))+
  geom_point(aes(colour = Type_of_parasite),size=2,show.legend = T) +
  geom_errorbarh(aes(xmin= Hedges_G - se_Hedges_G, 
                     xmax = Hedges_G + se_Hedges_G), 
                 height=.1) +
  scale_color_manual("Type_of_parasite",values=cols_subtype) +
  scale_x_continuous(name='Virulence effect size (Hedges g)') +
  ylab('Pathogen') +
  geom_vline(xintercept = c(0, 1.086), color= 'black', linetype=c('solid','dashed')) +
  annotate("rect", xmin = -0.263, xmax = 2.435, ymin = 0 , ymax = Inf,
           alpha = .1) +
  theme_classic()

Fig1a

### Moderator analysis ----

# Moderator 1:Treatment type is chosen as the reference category ---- 

TreatmentType <- rma.mv(Hedges_G, var_Hedges_G, mods = ~ Treatment_type, random = list(~ 1 | effect_size_ID, ~ 1 | Reference), tdist=TRUE, data=QLmod)
summary(TreatmentType, digits=3)

#### Making a plot of Treatment type (Code to get means and SE's): 

# create a data frame that contains the 'variable of interest' (in this case Treatment Type) and the Hedge's G values

TrT_df <- data.frame(QLmod$Treatment_type, QLmod$Hedges_G)

# check this data frame has been produced correctly

str(TrT_df)

# calculate the mean Hedge's G value for each of the different treatment types

TrT_df_mean <- aggregate(TrT_df[,2],list(TrT_df$QLmod.Treatment_type), mean)

# check this has been performed

TrT_df_mean

# give the data frame more meaningful and easy-to-read column names

TrT_colnames <- c("Treatment_Type", "Mean_Hedges_G")
colnames(TrT_df_mean) <- TrT_colnames
TrT_df_mean

# creating a function for calculating standard error
# make sure you change the 'z_TT' when you want to write a function for other variables (e.g. invert vs vert)

Standard_Error <- function(z_TrT){
  sd(z_TrT)/sqrt(length(z_TrT))
}

z_TrT <- TrT_df$QLmod.Hedges_G

# calculate the standard error of the mean (the mean Hedge's G) for each of the different treatment types

TrT_df_se <- aggregate(TrT_df[,2],list(TrT_df$QLmod.Treatment_type), Standard_Error)

# check this has been performed (remember standard errors are always positive - so a negative value indicates a mistake somewhere)

TrT_df_se

# give the data frame more eaningful and easy-to-read column names

TrT_se_colnames <- c("Treatment_Type", "SE")
colnames(TrT_df_se) <- TrT_se_colnames
TrT_df_se

# combine the treatment type, mean Hedge's G and standard error of the mean into a single data frame

TrT_graph_df <- data.frame(TrT_df_mean$Treatment_Type, TrT_df_mean$Mean_Hedges_G, TrT_df_se$SE)

# check the data frame has been created correctly 

TrT_graph_df

# plot a graph of effect size (y-axis) against treatment type (x-axis)
# remember - the x and y values (and the standard errors) are the column names from the new data frame (TT_graph_df)
# position_jitter to jitter the points 

TrT_plot <- ggplot(TrT_graph_df, aes(x = TrT_df_mean.Treatment_Type, y = TrT_df_mean.Mean_Hedges_G)) +
  geom_point() +
  geom_errorbar(aes(ymin = TrT_df_mean.Mean_Hedges_G - TrT_df_se.SE, 
                    ymax = TrT_df_mean.Mean_Hedges_G + TrT_df_se.SE), width=.1) +
  geom_point(data = QLmod, aes(x = Treatment_type, y = Hedges_G), alpha = 0.25, position=position_jitter(width =.05)) +
  geom_hline(yintercept = 0) +
  xlab('Treatment Type') +
  ylab('Hedges g') +
  theme_classic()

TrT_plot

# Host type is chosen as the reference category ---- 



# Removing rows of the data where the sample sizes are not large enough:

summary(QLmod$Type_of_host)

ModforHostType <- QLmod[ !(QLmod$Type_of_host %in% c('Coral ','Nematode','Zooplankton','Amphibian')), ]

ModforHostType$Type_of_host <- factor(ModforHostType$Type_of_host)

summary(ModforHostType$Type_of_host)

Hosttype <- rma.mv(Hedges_G, var_Hedges_G, mods = ~ Type_of_host, random = list(~ 1 | effect_size_ID, ~ 1 | Reference), tdist=TRUE, data=ModforHostType)

summary(Hosttype, digits=3)

#### Code to get means and SE's: 

# create a data frame that contains the 'variable of interest' (in this case Treatment Type) and the Hedge's G values

HoT_df <- data.frame(QLmod$Type_of_host, QLmod$Hedges_G)

# check this data frame has been produced correctly

str(HoT_df)

# calculate the mean Hedge's G value for each of the different treatment types

HoT_df_mean <- aggregate(HoT_df[,2],list(HoT_df$QLmod.Type_of_host), mean)

# check this has been performed

HoT_df_mean

# give the data frame more meaningful and easy-to-read column names

HoT_colnames <- c("Host_Type", "Mean_Hedges_G")
colnames(HoT_df_mean) <- HoT_colnames
HoT_df_mean

# creating a function for calculating standard error

Standard_Error <- function(z_HoT){
  sd(z_HoT)/sqrt(length(z_HoT))
}

z_HoT <- HoT_df$QLmod.Hedges_G

# calculate the standard error of the mean (the mean Hedge's G) for each of the different treatment types

HoT_df_se <- aggregate(HoT_df[,2],list(HoT_df$QLmod.Type_of_host), Standard_Error)

# check this has been performed 

HoT_df_se

# give the data frame more eaningful and easy-to-read column names

HoT_se_colnames <- c("Host_Type", "SE")
colnames(HoT_df_se) <- HoT_se_colnames
HoT_df_se

# combine the treatment type, mean Hedge's G and standard error of the mean into a single data frame

HoT_graph_df <- data.frame(HoT_df_mean$Host_Type, HoT_df_mean$Mean_Hedges_G, HoT_df_se$SE)

# check the data frame has been created correctly 

HoT_graph_df

# plot a graph of effect size (y-axis) against treatment type (x-axis)
# remember - the x and y values (and the standard errors) are the column names from the new data frame (TT_graph_df)

HoT_plot <- ggplot(HoT_graph_df, aes(x = HoT_df_mean.Host_Type, y = HoT_df_mean.Mean_Hedges_G)) +
  geom_point() +
  geom_errorbar(aes(ymin = HoT_df_mean.Mean_Hedges_G - HoT_df_se.SE, 
                    ymax = HoT_df_mean.Mean_Hedges_G + HoT_df_se.SE), width=.1) +
  geom_point(data = QLmod, aes(x = Type_of_host, y = Hedges_G), alpha = 0.25, position=position_jitter(width =.05)) +   
  xlab('Host') +
  ylab('Hedges g') +
  theme_classic()

HoT_plot


#### NOTE - the analysis above was carried out for the other moderator variables: Vertebrate vs invertebrate/ Host type / Pathogen type 




