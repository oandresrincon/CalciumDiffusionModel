#-----------------------------------------------------------------------------
#cts.py: constants used in model of Ca+2 dynamics
#Description: set constants and calculate concentration of binding sites at
#rest in fibers type I, IIA, IIX and IIB. 
#-----------------------------------------------------------------------------
#import requiered libraries
import numpy as np
import pandas as pd
#-----------------------------------------------------------------------------
#resting [Ca2+] in the cytoplasm (uM)
Ca0=0.106
Ca0a=['Ca0 (\u03BCM)',Ca0,Ca0,Ca0,Ca0]
cts=Ca0a
#-----------------------------------------------------------------------------
#resting [Mg2+] in the cytoplasm (uM)
Mg0=780
Mg0a=['Mg0 (\u03BCM)',Mg0,Mg0,Mg0,Mg0]
cts = np.vstack((cts,Mg0a))
#-----------------------------------------------------------------------------
#resting [Ca2+] in the SR (uM)
Ca0sr1=1.14e+3 
Ca0sr2a=1.07e+3 
Ca0sr2x=1.07e+3
Ca0sr2b=1.07e+3
Ca0sra=['Ca0sr (\u03BCM)',Ca0sr1,Ca0sr2a,Ca0sr2x,Ca0sr2b]
cts=np.vstack((cts,Ca0sra))
#-----------------------------------------------------------------------------
#resting [Ca2+] in the MITO (uM)
Ca0mito=0.15
Ca0mitoa=['Ca0mito (\u03BCM)',Ca0mito,Ca0mito,Ca0mito,Ca0mito]
cts=np.vstack((cts,Ca0mitoa))
#-----------------------------------------------------------------------------
#total ATP concentration (uM)
ATPT1=5000
ATPT2a=8000
ATPT2x=8000
ATPT2b=8000
ATPTa=['ATPT (\u03BCM)',ATPT1,ATPT2a,ATPT2x,ATPT2b]
cts=np.vstack((cts,ATPTa))
#-----------------------------------------------------------------------------
#Ca dissociation constant for sites of ATP (uM)
KdCaATP=2200
KdCaATPa=['KdCaATP (\u03BCM)',KdCaATP,KdCaATP,KdCaATP,KdCaATP]
cts=np.vstack((cts,KdCaATPa))
#-----------------------------------------------------------------------------
#reverse rate constant for Ca binding to S (per ms)
kCafATP=0.1364e-1
kCafATPa=['kCafATP (ms⁻¹)',kCafATP,kCafATP,kCafATP,kCafATP]
cts=np.vstack((cts,kCafATPa))
#-----------------------------------------------------------------------------
#reverse rate constant for Ca binding to S (per ms)
kCarATP=kCafATP*KdCaATP
kCarATPa=['kCarATP (ms⁻¹)',kCarATP,kCarATP,kCarATP,kCarATP]
cts=np.vstack((cts,kCarATPa))
#-----------------------------------------------------------------------------
#fraction of ATP bound with Ca at rest
CaATP01=ATPT1*Ca0/(Ca0+KdCaATP)
CaATP02a=ATPT2a*Ca0/(Ca0+KdCaATP)
CaATP02x=ATPT2x*Ca0/(Ca0+KdCaATP)
CaATP02b=ATPT2b*Ca0/(Ca0+KdCaATP)
CaATP0a=['CaATP0 (ms⁻¹)',CaATP01,CaATP02a,CaATP02x,CaATP02b]
cts=np.vstack((cts,CaATP0a))
#-----------------------------------------------------------------------------
#total Dye concentration (uM)
DyeT=229.1
DyeTa=['DyeT (\u03BCM)',DyeT,DyeT,DyeT,DyeT]
cts=np.vstack((cts,DyeTa))
#-----------------------------------------------------------------------------
cts=pd.DataFrame(cts)
cts.columns =['','I','IIA','IIX','IIB']
cts.to_excel('outputs/cts.xlsx',index=False)
#-----------------------------------------------------------------------------
#Ca dissociation constant for site S (uM)
KdCa_Dye=721.07427 
#-----------------------------------------------------------------------------
#forward rate constant for Ca binding to Dye (per uM per ms)
kCaf_Dye=1e+1
#-----------------------------------------------------------------------------
#reverse rate constant for Ca binding to Dye (per ms)
kCar_Dye=kCaf_Dye*KdCa_Dye 
#-----------------------------------------------------------------------------
#fraction of Dye bound with Ca at rest
fCaDye0=Ca0/(Ca0+KdCa_Dye)  
#-----------------------------------------------------------------------------
#total Troponin concentration (uM)
TropT1=120 
TropT2a=240 
TropT2x=240 
TropT2b=240 

Ca0mitoa=['Ca0mito (\u03BCM)',TropT1,TropT2a,TropT2x,TropT2b]
cts=np.vstack((cts,Ca0mitoa))
#-----------------------------------------------------------------------------
KdCa1_Trop=8.723#dissociation constant for binding of 1st Ca (uM)
kCaf1_Trop=1.77e-1 #forward rate constant for binding of 1 st Ca (uM)
kCar1_Trop=kCaf1_Trop*KdCa1_Trop #reverse rate constant for 1st bound (per ms)
KdCa2_Trop=0.194    #dissociation constant for binding of 2nd Ca (uM)
kCaf2_Trop=0.885e-1 #forward rate constant for binding of 2nd Ca (uM)
kCar2_Trop=kCaf2_Trop*KdCa2_Trop #reverse rate constant for 2nd bound Ca (per
# ms)
#-----------------------------------------------------------------------------
#reverse rate constant for Ca binding to Dye (per ms)
#kCar_Trop=kCaf_Trop*KdCa_Trop
#
##now calculate the fractional states of trop at rest
##-----------------------------------------------------------------------------
##fraction of Dye bound with Ca at rest
#fCaTrop0=Ca0/(Ca0+KdCa_Trop)  
##-----------------------------------------------------------------------------
##total Parvalbumin concentration (uM)
#Parv=1500   
##-----------------------------------------------------------------------------
##Ca dissociation constant for site S (uM)
#KdCa_Parv=0.012 
##-----------------------------------------------------------------------------
##forward rate constant for Ca binding to S (per uM per ms)
#kCaf_Parv = 0.417e-1 
##-----------------------------------------------------------------------------
##reverse rate constant for Ca binding to S (per ms)
#kCar_Parv = kCaf_Parv*KdCa_Parv 
##-----------------------------------------------------------------------------
##Mg dissociation constant for site S (uM)
#KdMg_Parv=90.9 
##-----------------------------------------------------------------------------
##forward rate constant for Mg binding to S (per uM per ms)
#kMgf_Parv = 3.3e-5 
##-----------------------------------------------------------------------------
##reverse rate constant for Ca binding to S (per ms)
#kMgr_Parv = kMgf_Parv*KdMg_Parv 
##-----------------------------------------------------------------------------
##next calculate the fractional states of parvalbumin at rest
##Ca-binding to Parv (resting)
#fCaParv0=Ca0/(Ca0+KdCa_Parv*(1+Mg0/KdMg_Parv)) 
##-----------------------------------------------------------------------------
##Mg-binding to Parv (resting)
#fMgParv0=Mg0/(Mg0+KdMg_Parv*(1+Ca0/KdCa_Parv)) 
##-----------------------------------------------------------------------------
##maximum flux of NCX (uM per ms)
#VNCX=5.46 
##-----------------------------------------------------------------------------
##membrane potential in resting state (mV)
#dVme=180 
##-----------------------------------------------------------------------------
##dissociation constant of Ca+2 to NCX molecules (uM)
#KdCa_NCX=140 
##-----------------------------------------------------------------------------
##dissociation constant of Na+ to NCX molecules (uM)
#KdNa_NCX=14e+3 
##-----------------------------------------------------------------------------
##extracelular [Na+] at rest (uM)
#Na_ex=140e+3 
##-----------------------------------------------------------------------------
##cytosolic [Na+] at rest (uM)
#Na_cy=10e+3 
##-----------------------------------------------------------------------------
##extracelular [Ca+2] at rest (uM)
#Ca_ex=1.0e+3 
##-----------------------------------------------------------------------------
##faraday constant (J per mol per mV)
#Fc=96.484 
##-----------------------------------------------------------------------------
##gas constant (J per mol per K)
#Rc=8.314 
##-----------------------------------------------------------------------------
##tempeture (K)
#T=296.15 
##-----------------------------------------------------------------------------
## next calculate the NCX flux at rest
#CaNCX0=((VNCX)\
#*(np.exp(+(((0.5)*dVme*Fc)/(Rc*T)))*((((Na_ex**(3))*Ca0))/((KdNa_NCX**(3))\
#*KdCa_NCX)))-(np.exp(-(((0.5)*dVme*Fc)/(Rc*T)))*((((Na_cy**(3))*Ca_ex))\
#/((KdNa_NCX**(3))*KdCa_NCX))))/(1+((Na_ex**(3))/(KdNa_NCX**(3)))+(Ca0\
#/KdCa_NCX)+(((Na_ex**(3))*Ca0)/((KdNa_NCX**(3))*KdCa_NCX))+((Na_cy**(3))\
#/(KdNa_NCX**(3)))+(Ca_ex/KdCa_NCX)+(((Na_cy**(3))*Ca_ex)/((KdNa_NCX**(3))\
#*KdCa_NCX)))
##-----------------------------------------------------------------------------
##maximum flux of SOCE (uM per ms)
#VSOCE=35e-3                                                           
##-----------------------------------------------------------------------------
##SOCE hill coefficient (A.U.)
#hSOCE=4.7  
##-----------------------------------------------------------------------------
##dissociation constant of Ca+2 to SOCE molecules (uM)
#KdSOCE=0.3e+3      
##-----------------------------------------------------------------------------
## next calculate the SOCE flux at rest (uM)
#CaSOCE0=(VSOCE*KdSOCE**hSOCE)/((Ca0_SR**hSOCE)+(KdSOCE**hSOCE))
##-----------------------------------------------------------------------------
##total mitochindrial buffer concentration (uM)
#B=2 
##-----------------------------------------------------------------------------
##forward rate constant for Ca2+ binding to B (per uM per ms)
#kCaf_B=0.8e-3 
##-----------------------------------------------------------------------------
##reverse rate constant for Ca2+ binding to B (per ms)
#kCar_B=0.192e-3 
##-----------------------------------------------------------------------------
##Ca2+ dissociation constant for site B (uM)
#kdCa_B=kCar_B/kCaf_B 
##-----------------------------------------------------------------------------
##maximum flux of NCE (uM per ms).
#VNCE=2.8e-3 
##-----------------------------------------------------------------------------
##mitochondrial membrane potential in resting state (mV)
#dVmito=190 
##-----------------------------------------------------------------------------
##dissociation constant of Ca+2 to NCX molecules (uM)
#KdCa_NCE=1.1 
##-----------------------------------------------------------------------------
##dissociation constant of Na+ to NCX molecules (uM)
#KdNa_NCE=8.2e+3 
##-----------------------------------------------------------------------------
##cytosolic [Na+] at rest (uM)
#Na_MITO=5e+3 
##-----------------------------------------------------------------------------
##maximum flux of MCU (uM per ms)
#VMCU=24.5e-3
##-----------------------------------------------------------------------------
##MCU hill coefficient (A.U.)
#hMCU=2 
##-----------------------------------------------------------------------------
##dissociation constant of Ca+2 to MCU molecules (uM)
#KdMCU=1.2 
##-----------------------------------------------------------------------------
##[B] bound with Ca2+ at rest (uM)
#CaB0=Ca0_MITO/(Ca0_MITO+kdCa_B) 
##-----------------------------------------------------------------------------
## next calculate the MCU flux at rest
#CaMCU0=(VMCU*Ca0**hMCU)/((Ca0**hMCU)+(KdMCU**hMCU)) 
##-----------------------------------------------------------------------------
## next calculate the NCE flux at rest
#CaNCE0=((VNCE)\
#*(np.exp(+(((0.5)*dVmito*Fc)/(Rc*T)))*((((Na_cy**(3))*Ca0_MITO))/((KdNa_NCE\
#**(3))*KdCa_NCE)))-(np.exp(-(((0.5)*dVmito*Fc)/(Rc*T)))*((((Na_MITO**(3))*Ca0\
#))/((KdNa_NCE**(3))*KdCa_NCE))))/(1+((Na_cy**(3))/(KdNa_NCE**(3)))+(Ca0_MITO\
#/KdCa_NCE)+(((Na_cy**(3))*Ca0_MITO)/((KdNa_NCE**(3))*KdCa_NCE))+((Na_MITO**\
#(3))/(KdNa_NCE**(3)))+(Ca0/KdCa_NCE)+(((Na_MITO**(3))*Ca0)/((KdNa_NCE**(3))*\
#KdCa_NCE)))                                                          
##-----------------------------------------------------------------------------
