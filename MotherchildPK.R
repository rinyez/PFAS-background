#******************************************************************************************************
# A Simple Pharmacokinetic Model of Prenatal and Postnatal Exposure to
# Perfluoroalkyl Substances (2016)
#
# Model coded in acslX by Marc-Andre Verner
# Translated into R by Yerin Jung
# Additional comments and code edits by Scott Bartell
#
# *Units*
# Time: hour; except "Age_" parameters & input time = [year]
# Volume: L
# Amounts: ug
#******************************************************************************************************

#install.package("deSolve")  
library(deSolve) 
library(dplyr)
library(tidyverse)

set.seed(271)
# Subject-specific characteristics
Age_conception = 30                     # Maternal Age at conception (years)
Age_delivery = Age_conception+0.769     # Maternal Age at delivery (years)
Weight_prepregnancy_M = 70              # Pre-pregnancy maternal body weight (kg)
Compound = 2                            # Compound (1:PFOS, 2:PFOA, 3:PFHxS, 4:PFNA)
Sex = 0                                 # Sex of child (Female=0; Male=1)
Total_bf = 6/12                         # Total duration of breastfeeding (years)
Total_pp = 12                           # Total postpartum duration for model predictions (years)
Total_pc = 0                            # Total pre-conception duration (years)
Dust = 0                                # Dust literature (0=Wu et al., 2015; 1= Winkens et al., 2018)
Diet = 0                                # Diet literature (0=Domingo et al., 2012; 1=Haug et al., 2010)

# Simulation time
Age_M =(Age_conception+seq(from=-Total_pc, to=0.769+Total_pp, by=1/365/24))    # Age of mother (years; by=1 hour)

# Median weight (kg) profiles (index: Age in years by=0.1 year)
Age_span=seq(from=0,to=100,by=0.1)
table_bw_females=list("years"=c(0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 13.0, 16.0, 20.0, 100.0),
                      "y"= c(3.4, 7.2, 9.5, 12.1, 18.0, 33.0, 46.0, 54.0, 58.0, 58.0))
table_bw_females = data.frame(approx(table_bw_females$years,table_bw_females$y, xout=round(Age_span,1), rule=2))
BW_F_func <- approxfun(table_bw_females$x, table_bw_females$y, rule = 2)
table_bw_males=list("years"=c(0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 13.0, 16.0, 20.0, 100.0),
                    "y"=c(3.5, 7.8, 10.3, 12.7, 18.0, 32.0, 46.0, 61.0, 71.0, 71.0))
table_bw_males = data.frame(approx(table_bw_males$years,table_bw_males$y, xout=round(Age_span,1), rule=2))
BW_M_func <- approxfun(table_bw_males$x, table_bw_males$y, rule = 2) 
table_fetal_weight=list("years"=c(0.0, 0.217, 0.353, 0.463, 0.551, 0.639, 0.768),
                        "y"=c(0.0, 0.0, 0.26, 0.69, 1.25, 2.02, 3.28))
table_fetal_weight = data.frame(approx(table_fetal_weight$years,table_fetal_weight$y, xout=round(seq(from=0,to=0.8,by=0.1),1), rule=2))
BW_C_func <- approxfun(table_fetal_weight$x, table_fetal_weight$y, rule = 2) 

# Children Dust ingestion rate by age [g/hour]
# indoor settled dust only
table_q_dust_child=list("years"=c(0.0, 0.5, 1.0, 6.0, 20.0),
                        "y"= c(0, 30, 60, 60, 30)/24/10^3)
table_q_dust_child = data.frame(approx(table_q_dust_child$years,table_q_dust_child$y, 
                                       xout=round(seq(from=0,to=20,by=0.1),1), rule=2))
q_dust_child_func <- approxfun(table_q_dust_child$x, table_q_dust_child$y, rule = 2)

# Children fraction of time spent indoors [unitless] (US EPA, 2011)
table_f_ind_child=list("years"=round(c(0.0, 1/12, 3/12, 6/12, 1, 2, 3, 6, 11, 16),1),
                        "y"= c(1440, 1432, 1414, 1301, 1353, 1316, 1278, 1244, 1260, 1248)/60/24)
table_f_ind_child = data.frame(approx(table_f_ind_child$years,table_f_ind_child$y, 
                                       xout=round(seq(from=0,to=20,by=0.1),1), rule=2))
f_ind_child_func <- approxfun(table_f_ind_child$x, table_f_ind_child$y, rule = 2)

# Adult fraction of time spent indoors [unitless] (US EPA, 2011)
q_dust_m = 30/24/10^3           # Dust ingestion rate [g/hour]; indoor settled dust only
f_ind_m = 1159/60/24            # Mean daily time spent indoors at the age of 18-64 y.o.
F_uptake = 0.8                  # The gastrointestinal uptake factor (Gebbink et al. 2015; Y. Wu et al., 2020 )

#******************************************************************************************************
# Compound
PFOS = 1*(Compound == 1)
PFOA = 1*(Compound == 2)
PFHXS = 1*(Compound == 3)
PFNA = 1*(Compound == 4)

#Dust concentration (ng/g)
PFOS_dust = 29*(1-Dust)+1860*Dust
PFOA_dust = 41.4*(1-Dust)+746*Dust
PFHXS_dust = 3.47*(1-Dust)+1410*Dust
PFNA_dust = 13.3*(1-Dust)+89.5*Dust
# Initial PFAS level in indoor dust (ug/g)
dust_conc = (PFOS_dust*PFOS + PFOA_dust*PFOA + PFNA_dust*PFNA + PFHXS_dust*PFHXS)/1000

# Daily intake through diet (ng/kg/d)
PFOS_diet = 2.26*(1-Diet) + 1.2*Diet
PFOA_diet = 6.37*(1-Diet) + 0.55*Diet
PFNA_diet = 1.47
PFHXS_diet = 0.06
diet_edi = PFOS_diet*PFOS + PFOA_diet*PFOA + PFNA_diet*PFNA + PFHXS_diet*PFHXS                           


# Volume of distribution (l/kg)
PFOS_VD_BW = 0.230    # PFOS
PFOA_VD_BW = 0.170    # PFOA
PFHXS_VD_BW = 0.213   # PFHXS
PFNA_VD_BW = 0.170    # PFNA (Zhang et al., 2013)
VD_BW = PFOS*PFOS_VD_BW+PFOA*PFOA_VD_BW+PFHXS*PFHXS_VD_BW+PFNA*PFNA_VD_BW

# Plasma:milk level ratio
PFOS_PMilk = 0.0138   # PFOS
PFOA_PMilk = 0.0577   # PFOA
PFHXS_PMilk = 0.0140  # PFHxS
PFNA_PMilk = 0.0056   # PFNA (Kim et al., 2011)
PMilk = PFOS*PFOS_PMilk+PFOA*PFOA_PMilk+PFHXS*PFHXS_PMilk+PFNA*PFNA_PMilk

# Cord:maternal plasma level ratio
PFOS_P_CM = 0.454     # PFOS
PFOA_P_CM = 0.783     # PFOA
PFHXS_P_CM = 0.556    # PFHxS
PFNA_P_CM = 0.468    # PFNA (Kim et al., 2011)
P_CM = PFOS*PFOS_P_CM+PFOA*PFOA_P_CM+PFHXS*PFHXS_P_CM+PFNA*PFNA_P_CM

# Half-lives (years)
PFOS_HL = 5.4         # PFOS
PFOA_HL = 2.3         # PFOA (Bartell et al., 2010)
PFHXS_HL = 8.5        # PFHxS
PFNA_HL = 3.9         # PFNA (Zhang et al., 2013)
# Half-lives (hours)
HALF_LIFE = (PFOS*PFOS_HL+PFOA*PFOA_HL+PFHXS*PFHXS_HL+PFNA*PFNA_HL)*365*24 

# placental transfer
ktrans1 = P_CM # Mother->fetus placental transfer  (l/h)
ktrans2 = 1.0  # Fetus->mother placental transfer  (l/h)
#******************************************************************************************************

## PBPK model
# Unit: age [y] input time [h] mass [ug] volume [l] unless noted
pbpkmodel <- function(Time, State, Parmeters) {
  with(as.list(c(State, Paras)), {
    #Switches for Time Events    #*********************************************************************
    if (Age_conception-Time/24/365>0) {            #Conception
      IO_CONCEPTION = 0
    } else {
      IO_CONCEPTION = 1
    }
    if (Age_delivery-Time/24/365>0) {              #Delivery
      IO_DELIVERY = 0
    } else {
      IO_DELIVERY = 1
    }
    if (Total_bf-(Time/24/365-Age_delivery)>0 & (Time/24/365-Age_delivery)>0 ) {  #Breasfeeding
      IO_BF = 1
    } else {
      IO_BF = 0
    }
    
    # Age of child (year)
    Age_C = IO_DELIVERY*(Time/24/365-Age_delivery)  #[y]
    if (Total_bf-Age_C>0){                   # Breastfeeding
      IO_END_TOTAL_BF = 0
    } else {
      IO_END_TOTAL_BF = 1
    }
    if (0.500-Age_C>0){                      # Child age > 6 m.o.
      IO_CHILD_6m = 0 
    } else {
      IO_CHILD_6m = 1
    }
    if (1.000-Age_C>0){                      # Child age > 12 m.o.
      IO_CHILD_12m = 0
    } else {
      IO_CHILD_12m = 1
    }
    if (2.500-Age_C>0){                      # Child age > 30 m.o.
      IO_CHILD_30m = 0
    } else {
      IO_CHILD_30m = 1
    }
    #******************************************************************************************************
    # weight (kg)
    BW_M = (Weight_prepregnancy_M/BW_F_func(Age_conception))*BW_F_func(Time/365/24)
    BW_C = ((IO_CONCEPTION-IO_DELIVERY)*BW_C_func(Time/365/24-Age_conception) + 
                     (IO_DELIVERY)*((1-Sex)*BW_F_func(Age_C)+Sex*BW_M_func(Age_C)))
    #******************************************************************************************************
    # Volume of distribution (l)
    VD_M = VD_BW*BW_M # Maternal volume of distribution (l)
    VD_C = VD_BW*BW_C+1E-03 # Child volume of distribution (l)
    
    C_MOTHER = A_MOTHER/VD_M # Concentration in maternal compartment (ug/l)
    C_CHILD = A_CHILD/VD_C # Concentration in child compartment (ug/l)

    # Placental diffusion
    DIFF_M_C = (IO_CONCEPTION-IO_DELIVERY)*ktrans1*C_MOTHER # Mother->fetus placental transfer rate (ug/h)
    DIFF_C_M = (IO_CONCEPTION-IO_DELIVERY)*ktrans2*C_CHILD # Fetus->mother placental transfer rate (ug/h)
    #******************************************************************************************************
    # Breast milk transfer
    VOLUME_MILK_1ST_YEAR = BW_C*(-0.312*(Age_C*365)+157.7)/24  # (g/h)
    VOLUME_MILK_2ND_YEAR = (-0.763*Age_C*365 + 720.3)/24       # (g/h)
    VOLUME_MILK = IO_BF*(IO_DELIVERY-IO_END_TOTAL_BF)* 
      ((1-IO_CHILD_12m)*VOLUME_MILK_1ST_YEAR 
       +(IO_CHILD_12m-IO_CHILD_30m)*VOLUME_MILK_2ND_YEAR)/1000 # Breast milk consumption rate (l/h)

    # Breast milk consumption rate (l/h)
    C_MILK = PMilk*C_MOTHER # Concentration in breast milk (ug/l)
    TRANSFER_BF = VOLUME_MILK*C_MILK # Mother-child lactational transfer (ug/h)
    #******************************************************************************************************
    # Elimination
    ELIMINATION_M = C_MOTHER*VD_M*log(2)/HALF_LIFE # Elimination from the maternal compartment (ug/h)
    ELIMINATION_C = C_CHILD*VD_C*log(2)/HALF_LIFE  # Elimination from the child compartment (ug/h)
    #******************************************************************************************************
    # Dosing
    INTAKE_M =diet_edi/1000/24*BW_M  +                                    # Maternal PFAS dietary intake (ug/h)
      dust_conc*q_dust_m*f_ind_m*F_uptake                                 # Maternal PFAS dust intake (ug/h)
    INTAKE_C = IO_CHILD_6m*((diet_edi/1000/24)/BW_M*BW_C +                # Child PFAS dietary intake (ug/h)
      dust_conc*q_dust_child_func(Age_C)*f_ind_child_func(Age_C)*F_uptake) # Child PFAS dust intake (ug/h)
    
    #******************************************************************************************************
    # Mass-balance differential equations
    # Maternal compartment
    RA_MOTHER = INTAKE_M-ELIMINATION_M-TRANSFER_BF-DIFF_M_C+DIFF_C_M # Rate of amount in maternal compartment (ug/h)
    # Child compartment
    RA_CHILD = INTAKE_C+TRANSFER_BF-ELIMINATION_C+DIFF_M_C-DIFF_C_M # Rate of amount in child compartment (ug/h)

    ## Mass balance
    Tmass = A_MOTHER+A_CHILD
    list(c(RA_MOTHER,RA_CHILD), 
         c(dose=INTAKE_M,Tmass=Tmass, serum_M=C_MOTHER, serum_C=C_CHILD,Age_C=Age_C) )
      })
  }

# For initial conditions, assume mom is at steady state at time of conception
C_SS <- (diet_edi/1000/24+(dust_conc*q_dust_m*f_ind_m*F_uptake)/Weight_prepregnancy_M) / 
  (log(2)/HALF_LIFE) / VD_BW   # maternal prepregnancy body weight cancelled out 
# (ug/h/kg)/(/h)/(l/kg) = ug/L
State <- c(A_MOTHER=C_SS*VD_BW*Weight_prepregnancy_M, A_CHILD=0) #[ug]
Paras <- c(ktrans1, ktrans2, HALF_LIFE)

# PBPK output
# input unit: initial state (ng), time (hour)
# output unit: serum [ug/L; ng/mL], age [year]
out <- ode(y = State,
           times = Age_M*365*24, 
           func = pbpkmodel, 
           parms = Paras) 

# extract results
preds = as.data.frame(out)
preds$time <- preds$time/24 #time unit: hour -> day
#write.csv(preds, file="pfna_Wu_Domingo.csv")    # as a csv.file

#Figure S1
plot(preds$time/365 - Age_conception, preds$serum_M, ylim=c(1,50), 
     log="y", xlab="Year after conception", ylab="Maternal PFOA serum level (ng/mL)", typ="l",xlim=c(0,3))
abline(v=0.769,lty=2)
text(0.769/2,10,"Gestation",col="red",cex=0.9,pos=1)
abline(v=0.769+Total_bf,lty=2)
text(0.769+Total_bf/2,16,"Breastfeeding",col="blue",cex=0.8,pos=1)