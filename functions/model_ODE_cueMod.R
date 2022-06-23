library(rootSolve) #stode
library(hydroGOF) #NSE
library(ggplot2) #histogram plots
library(reshape2) #melt
library(scales) #trans_breaks, pretty colors for figures
library(gridExtra) #multi-panel ggplot figure
library(grid) #multi-panel ggplot figure
library(cowplot) #multi-panel ggplot figure
library(deSolve)


model_ODE_cueMod <- function(data, times){
  FXEQ <- function(t, y, pars) {
    with (as.list(c(y, pars)),{
      
      #Carbon fluxes
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flows to and from MIC_1
      LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)#LIT_1)   #MIC_1 decomp of MET lit
      LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)#LIT_2)   #MIC_1 decomp of STRUC lit
      SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)#SOM_3)   #decomp of SOMa by MIC_1
      MICtrn[1] = MIC_1^densDep * turnover[1]  * fPHYS[1] #MIC_1 turnover to PHYSICAL SOM
      MICtrn[2] = MIC_1^densDep * turnover[1]  * fCHEM[1] #MIC_1 turnover to CHEMICAL SOM
      MICtrn[3] = MIC_1^densDep * turnover[1]  * fAVAI[1] #MIC_1 turnover to AVAILABLE SOM
      
      #Flows to and from MIC_2
      LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)#LIT_1)   #decomp of MET litter
      LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)#LIT_2)   #decomp of SRUCTURAL litter
      SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)#SOM_3)   #decomp of SOMa by MIC_2
      MICtrn[4] = MIC_2^densDep * turnover[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM
      MICtrn[5] = MIC_2^densDep * turnover[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM
      MICtrn[6] = MIC_2^densDep * turnover[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM
      
      DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
      OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +#SOM_2)) +
                     (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))#SOM_2)))  #oxidation of C to A
      
      #Nitrogen fluxes
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flows to and from MIC_1
      LITminN[1] =  LITmin[1]*LIT_1_N/(LIT_1+1e-100)#LITmin[1]/CN_m
      LITminN[2] =  LITmin[2]*LIT_2_N/(LIT_2++1e-100)#LITmin[2]/CN_s
      SOMminN[1] =  SOMmin[1]*(SOM_3_N/(SOM_3++1e-100))#SOMmin[1]*(SOM_3_N/SOM_3)#*relRateN
      MICtrnN[1] =  MICtrn[1]*MIC_1_N/(MIC_1+1e-100)#MICtrn[1]/CN_r
      MICtrnN[2] =  MICtrn[2]*MIC_1_N/(MIC_1+1e-100)#MICtrn[2]/CN_r
      MICtrnN[3] =  MICtrn[3]*MIC_1_N/(MIC_1+1e-100)#MICtrn[3]/CN_r
      
      #Flows to and from MIC_2
      LITminN[3] =  LITmin[3]*LIT_1_N/(LIT_1+1e-100)#LITmin[3]/CN_m
      LITminN[4] =  LITmin[4]*LIT_2_N/(LIT_2+1e-100)#LITmin[4]/CN_s
      SOMminN[2] =  SOMmin[2]*(SOM_3_N/(SOM_3+1e-100))#SOMmin[2]*(SOM_3_N/SOM_3)#*relRateN
      MICtrnN[4] =  MICtrn[4]*MIC_2_N/(MIC_2+1e-100)#MICtrn[4]/CN_K
      MICtrnN[5] =  MICtrn[5]*MIC_2_N/(MIC_2+1e-100)#MICtrn[5]/CN_K
      MICtrnN[6] =  MICtrn[6]*MIC_2_N/(MIC_2+1e-100)#MICtrn[6]/CN_K
      
      DEsorbN    =  DEsorb[1]*(SOM_1_N/(SOM_1+1e-100))#*relRateN
      OXIDATN    =  OXIDAT[1]*(SOM_2_N/(SOM_2+1e-100))#*relRateN
      DINup[1]   = (1-Nleak)*DIN*MIC_1/(MIC_1+MIC_2+1e-100) #Partitions DIN between microbial pools based on relative biomass
      DINup[2]   = (1-Nleak)*DIN*MIC_2/(MIC_1+MIC_2+1e-100)
      
      #####
      upMIC_1    = CUE[1]*(LITmin[1] + SOMmin[1]) + CUE[2]*(LITmin[2])
      upMIC_1_N  = NUE*(LITminN[1] + SOMminN[1]) + NUE*(LITminN[2]) + DINup[1]
      CNup[1]    = (upMIC_1)/(upMIC_1_N+1e-100) #Trying to optimize run speed here by avoiding /0
      Overflow[1] = (upMIC_1) - (upMIC_1_N)*min(CN_r, CNup[1])
      Nspill[1]   = (upMIC_1_N) - (upMIC_1)/max(CN_r, CNup[1])
      
      upMIC_2    = CUE[3]*(LITmin[3] + SOMmin[2]) + CUE[4]*(LITmin[4])
      upMIC_2_N  = NUE*(LITminN[3] + SOMminN[2]) + NUE*(LITminN[4]) + DINup[2]
      CNup[2]    = (upMIC_2)/(upMIC_2_N+1e-100)
      Overflow[2] = (upMIC_2) - (upMIC_2_N)*min(CN_K, CNup[2])
      Nspill[2]   = (upMIC_2_N) - (upMIC_2)/max(CN_K, CNup[2])
      ######
      
      dLIT_1 = Inputs[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
      dLIT_2 = Inputs[2]*(1-FI[2]) - LITmin[2] - LITmin[4]
      dMIC_1 = CUE[1]*(LITmin[1] + SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3]) - Overflow[1]
      dMIC_2 = CUE[3]*(LITmin[3] + SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6]) - Overflow[2] 
      dSOM_1 = Inputs[1]*FI[1] + MICtrn[1] + MICtrn[4] - DEsorb 
      dSOM_2 = Inputs[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
      dSOM_3 = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
      
      dLIT_1_N = Inputs[1]*(1-FI[1])/CN_m - LITminN[1] - LITminN[3]
      dLIT_2_N = Inputs[2]*(1-FI[2])/CN_s - LITminN[2] - LITminN[4]
      dMIC_1_N = NUE*(LITminN[1] + SOMminN[1]) + NUE*(LITminN[2]) - sum(MICtrnN[1:3]) + DINup[1] - Nspill[1]
      dMIC_2_N = NUE*(LITminN[3] + SOMminN[2]) + NUE*(LITminN[4]) - sum(MICtrnN[4:6]) + DINup[2] - Nspill[2]
      dSOM_1_N = Inputs[1]*FI[1]/CN_m + MICtrnN[1] + MICtrnN[4] - DEsorbN
      dSOM_2_N = Inputs[2]*FI[2]/CN_s + MICtrnN[2] + MICtrnN[5] - OXIDATN
      dSOM_3_N = MICtrnN[3] + MICtrnN[6] + DEsorbN + OXIDATN - SOMminN[1] - SOMminN[2]
      
      dDIN = (1-NUE)*(LITminN[1] + LITminN[2] + SOMminN[1]) +  #Inputs from r decomp
        (1-NUE)*(LITminN[3] + LITminN[4] + SOMminN[2]) +  #Inputs from K decomp
        Nspill[1] + Nspill[2] - DINup[1] - DINup[2]    #Uptake to microbial pools and spillage
      LeachingLoss = Nleak*DIN
      dDIN = dDIN-LeachingLoss #N leaching losses    
      
      list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3,dLIT_1_N, dLIT_2_N, dMIC_1_N, dMIC_2_N, dSOM_1_N, dSOM_2_N, dSOM_3_N,dDIN))
    })
  }
  
  
  TFXEQ <- function(t, y, pars) {
    with (as.list(c(y, pars)),{
      
      #Carbon fluxes
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flows to and from MIC_1
      LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)#LIT_1)   #MIC_1 decomp of MET lit
      LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)#LIT_2)   #MIC_1 decomp of STRUC lit
      SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)#SOM_3)   #decomp of SOMa by MIC_1
      MICtrn[1] = MIC_1^densDep * turnover[1]  * fPHYS[1] #MIC_1 turnover to PHYSICAL SOM
      MICtrn[2] = MIC_1^densDep * turnover[1]  * fCHEM[1] #MIC_1 turnover to CHEMICAL SOM
      MICtrn[3] = MIC_1^densDep * turnover[1]  * fAVAI[1] #MIC_1 turnover to AVAILABLE SOM
      
      #Flows to and from MIC_2
      LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)#LIT_1)   #decomp of MET litter
      LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)#LIT_2)   #decomp of SRUCTURAL litter
      SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)#SOM_3)   #decomp of SOMa by MIC_2
      MICtrn[4] = MIC_2^densDep * turnover[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM
      MICtrn[5] = MIC_2^densDep * turnover[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM
      MICtrn[6] = MIC_2^densDep * turnover[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM
      
      DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
      OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +#SOM_2)) +
                     (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))#SOM_2)))  #oxidation of C to A
      
      #Nitrogen fluxes
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Flows to and from MIC_1
      LITminN[1] =  LITmin[1]*LIT_1_N/(LIT_1+1e-100)#LITmin[1]/CN_m
      LITminN[2] =  LITmin[2]*LIT_2_N/(LIT_2++1e-100)#LITmin[2]/CN_s
      SOMminN[1] =  SOMmin[1]*(SOM_3_N/(SOM_3++1e-100))#SOMmin[1]*(SOM_3_N/SOM_3)#*relRateN
      MICtrnN[1] =  MICtrn[1]*MIC_1_N/(MIC_1+1e-100)#MICtrn[1]/CN_r
      MICtrnN[2] =  MICtrn[2]*MIC_1_N/(MIC_1+1e-100)#MICtrn[2]/CN_r
      MICtrnN[3] =  MICtrn[3]*MIC_1_N/(MIC_1+1e-100)#MICtrn[3]/CN_r
      
      #Flows to and from MIC_2
      LITminN[3] =  LITmin[3]*LIT_1_N/(LIT_1+1e-100)#LITmin[3]/CN_m
      LITminN[4] =  LITmin[4]*LIT_2_N/(LIT_2+1e-100)#LITmin[4]/CN_s
      SOMminN[2] =  SOMmin[2]*(SOM_3_N/(SOM_3+1e-100))#SOMmin[2]*(SOM_3_N/SOM_3)#*relRateN
      MICtrnN[4] =  MICtrn[4]*MIC_2_N/(MIC_2+1e-100)#MICtrn[4]/CN_K
      MICtrnN[5] =  MICtrn[5]*MIC_2_N/(MIC_2+1e-100)#MICtrn[5]/CN_K
      MICtrnN[6] =  MICtrn[6]*MIC_2_N/(MIC_2+1e-100)#MICtrn[6]/CN_K
      
      DEsorbN    =  DEsorb[1]*(SOM_1_N/(SOM_1+1e-100))#*relRateN
      OXIDATN    =  OXIDAT[1]*(SOM_2_N/(SOM_2+1e-100))#*relRateN
      DINup[1]   = (1-Nleak)*DIN*MIC_1/(MIC_1+MIC_2+1e-100) #Partitions available DIN between microbial pools based on relative biomass
      DINup[2]   = (1-Nleak)*DIN*MIC_2/(MIC_1+MIC_2+1e-100)
      
      #####
      upMIC_1    = CUE[1]*(LITmin[1] + SOMmin[1]) + CUE[2]*(LITmin[2])
      upMIC_1_N  = NUE*(LITminN[1] + SOMminN[1]) + NUE*(LITminN[2]) + DINup[1]
      CNup[1]    = (upMIC_1)/(upMIC_1_N+1e-100) #Trying to optimize run speed here by avoiding /0
      Overflow[1] = (upMIC_1) - (upMIC_1_N)*min(CN_r, CNup[1])
      Nspill[1]   = (upMIC_1_N) - (upMIC_1)/max(CN_r, CNup[1])
      
      upMIC_2    = CUE[3]*(LITmin[3] + SOMmin[2]) + CUE[4]*(LITmin[4])
      upMIC_2_N  = NUE*(LITminN[3] + SOMminN[2]) + NUE*(LITminN[4]) + DINup[2]
      CNup[2]    = (upMIC_2)/(upMIC_2_N+1e-100)
      Overflow[2] = (upMIC_2) - (upMIC_2_N)*min(CN_K, CNup[2])
      Nspill[2]   = (upMIC_2_N) - (upMIC_2)/max(CN_K, CNup[2])
      ######
      
      dLIT_1 = Inputs[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
      dLIT_2 = Inputs[2]*(1-FI[2]) - LITmin[2] - LITmin[4]
      dMIC_1 = CUE[1]*(LITmin[1] + SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3]) - Overflow[1]
      dMIC_2 = CUE[3]*(LITmin[3] + SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6]) - Overflow[2] 
      dSOM_1 = Inputs[1]*FI[1] + MICtrn[1] + MICtrn[4] - DEsorb 
      dSOM_2 = Inputs[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
      dSOM_3 = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
      
      dLIT_1_N = Inputs[1]*(1-FI[1])/CN_m - LITminN[1] - LITminN[3]
      dLIT_2_N = Inputs[2]*(1-FI[2])/CN_s - LITminN[2] - LITminN[4]
      dMIC_1_N = NUE*(LITminN[1] + SOMminN[1]) + NUE*(LITminN[2]) - sum(MICtrnN[1:3]) + DINup[1] - Nspill[1]
      dMIC_2_N = NUE*(LITminN[3] + SOMminN[2]) + NUE*(LITminN[4]) - sum(MICtrnN[4:6]) + DINup[2] - Nspill[2]
      dSOM_1_N = Inputs[1]*FI[1]/CN_m + MICtrnN[1] + MICtrnN[4] - DEsorbN
      dSOM_2_N = Inputs[2]*FI[2]/CN_s + MICtrnN[2] + MICtrnN[5] - OXIDATN
      dSOM_3_N = MICtrnN[3] + MICtrnN[6] + DEsorbN + OXIDATN - SOMminN[1] - SOMminN[2]
      
      dDIN = (1-NUE)*(LITminN[1] + LITminN[2] + SOMminN[1]) +  #Inputs from r decomp
        (1-NUE)*(LITminN[3] + LITminN[4] + SOMminN[2]) +  #Inputs from K decomp
        Nspill[1] + Nspill[2] - DINup[1] - DINup[2]    #Uptake to microbial pools and spillage
      LeachingLoss = Nleak*DIN
      dDIN = dDIN-LeachingLoss #N leaching losses  
      
      #Tracer tier calculations
      # III. Carbon, tracer
      TLITmin[1] = TLIT_1/LIT_1*LITmin[1]   #MIC_1 decomp of MET lit
      TLITmin[2] = TLIT_2/LIT_2*LITmin[2]   #MIC_1 decomp of STRUC lit
      TSOMmin[1] = TSOM_3/SOM_3*SOMmin[1]   #Decomp of SOMa by MIC_1
      TMICtrn[1] = TMIC_1/MIC_1*MICtrn[1]   #MIC_1 turnover to SOMp
      TMICtrn[2] = TMIC_1/MIC_1*MICtrn[2]   #MIC_1 turnover to SOMc
      TMICtrn[3] = TMIC_1/MIC_1*MICtrn[3]   #MIC_1 turnover to SOMa
      
      #Flows to and from MIC_2
      TLITmin[3] = TLIT_1/LIT_1*LITmin[3]   #decomp of MET litter
      TLITmin[4] = TLIT_2/LIT_2*LITmin[4]   #decomp of SRUCTURAL litter
      TSOMmin[2] = TSOM_3/SOM_3*SOMmin[2]   #decomp of PHYSICAL SOM by MIC_1
      TMICtrn[4] = TMIC_2/MIC_2*MICtrn[4]   #MIC_2 turnover to SOMp
      TMICtrn[5] = TMIC_2/MIC_2*MICtrn[5]   #MIC_2 turnover to SOMc
      TMICtrn[6] = TMIC_2/MIC_2*MICtrn[6]   #MIC_2 turnover to SOMa
      
      TDEsorb[1]    = TSOM_1/SOM_1*DEsorb[1]#SOM_1[2] * desorb  #* (MIC_1 + MIC_2)  #desorbtion of PHYS to AVAIL (function of fCLAY)
      TOXIDAT[1]    = TSOM_2/SOM_2*OXIDAT[1]
      
      #Flows to and from MIC_1
      TLITminN[1] =  TLITmin[1]*TLIT_1_N/(TLIT_1+1e-100)
      TLITminN[2] =  TLITmin[2]*TLIT_2_N/(TLIT_2+1e-100)
      TMICtrnN[1] =  TMICtrn[1]*TMIC_1_N/(TMIC_1+1e-100)
      TMICtrnN[2] =  TMICtrn[2]*TMIC_1_N/(TMIC_1+1e-100)
      TMICtrnN[3] =  TMICtrn[3]*TMIC_1_N/(TMIC_1+1e-100)
      TSOMminN[1] =  TSOMmin[1]*TSOM_3_N/(TSOM_3+1e-100)  
      
      #Flows to and from MIC_2
      TLITminN[3] =  TLITmin[3]*TLIT_1_N/(TLIT_1+1e-100)
      TLITminN[4] =  TLITmin[4]*TLIT_2_N/(TLIT_2+1e-100)
      TMICtrnN[4] =  TMICtrn[4]*TMIC_2_N/(TMIC_2+1e-100)
      TMICtrnN[5] =  TMICtrn[5]*TMIC_2_N/(TMIC_2+1e-100)
      TMICtrnN[6] =  TMICtrn[6]*TMIC_2_N/(TMIC_2+1e-100)
      TSOMminN[2] =  TSOMmin[2]*TSOM_3_N/(TSOM_3+1e-100)
      
      TDEsorbN =  TDEsorb*TSOM_1_N/(TSOM_1+1e-100)
      TOXIDATN =  TOXIDAT*TSOM_2_N/(TSOM_2+1e-100)
      TDINup[1] = TDIN/(DIN+1e-100)*DINup[1]
      TDINup[2] = TDIN/(DIN+1e-100)*DINup[2]
      
      #####
      TNimport[1] = fracNImportr*LeachingLoss
      TNimport[2] = fracNImportK*LeachingLoss
      
      upTMIC_1    = CUE[1]*(TLITmin[1] + TSOMmin[1]) + CUE[2]*TLITmin[2]
      upTMIC_1_N  = NUE*(TLITminN[1] + TSOMminN[1]) + NUE*TLITminN[2] + TDINup[1] + TNimport[1]
      TCNup[1]    = upTMIC_1/(upTMIC_1_N+1e-100)
      TOverflow[1] = upTMIC_1 - upTMIC_1_N*min(CN_r, TCNup[1])
      TNspill[1]   = upTMIC_1_N - upTMIC_1/max(CN_r, TCNup[1])
      
      upTMIC_2    = CUE[3]*(TLITmin[3] + TSOMmin[2]) + CUE[4]*TLITmin[4]
      upTMIC_2_N  = NUE*(TLITminN[3] + TSOMminN[2]) + NUE*TLITminN[4] + TDINup[2] + TNimport[2]
      TCNup[2]    = upTMIC_2/(upTMIC_2_N+1e-100)
      TOverflow[2] = upTMIC_2 - upTMIC_2_N*min(CN_K, TCNup[2])
      TNspill[2]   = upTMIC_2_N - upTMIC_2/max(CN_K, TCNup[2])
      ######
      
      dTLIT_1 = - TLITmin[1] - TLITmin[3]
      dTLIT_2 = - TLITmin[2] - TLITmin[4]
      dTMIC_1 = CUE[1]*(TLITmin[1] + TSOMmin[1]) + CUE[2]*TLITmin[2] - sum(TMICtrn[1:3]) - TOverflow[1]
      dTMIC_2 = CUE[3]*(TLITmin[3] + TSOMmin[2]) + CUE[4]*TLITmin[4] - sum(TMICtrn[4:6]) - TOverflow[2]
      dTSOM_1 = TMICtrn[1] + TMICtrn[4] - TDEsorb[1]
      dTSOM_2 = TMICtrn[2] + TMICtrn[5] - TOXIDAT[1]
      dTSOM_3 = TMICtrn[3] + TMICtrn[6] + TDEsorb[1] + TOXIDAT[1] - TSOMmin[1] - TSOMmin[2] 
      
      dTLIT_1_N = - TLITminN[1] - TLITminN[3]
      dTLIT_2_N = - TLITminN[2] - TLITminN[4]
      dTMIC_1_N = NUE*TLITminN[1] + NUE*TLITminN[2] - sum(TMICtrnN[1:3]) + TDINup[1] + NUE*TSOMminN[1] + TNimport[1] - TNspill[1]
      dTMIC_2_N = NUE*TLITminN[3] + NUE*TLITminN[4] - sum(TMICtrnN[4:6]) + TDINup[2] + NUE*TSOMminN[2] + TNimport[2] - TNspill[2]
      dTSOM_1_N = TMICtrnN[1] + TMICtrnN[4] - TDEsorbN[1]
      dTSOM_2_N = TMICtrnN[2] + TMICtrnN[5] - TOXIDATN[1]
      dTSOM_3_N = TMICtrn[3]/CN_r + TMICtrn[6]/CN_K + TDEsorbN[1] + TOXIDATN[1] - TSOMminN[1] - TSOMminN[2]
      
      dTDIN = (1-NUE)*(TLITminN[1] + TLITminN[2] + TSOMminN[1]) +  #Inputs from r decomp
        (1-NUE)*(TLITminN[3] + TLITminN[4] + TSOMminN[2]) +  #Inputs from K decomp
        TNspill[1] + TNspill[2] - TDINup[1] - TDINup[2]    #Uptake to microbial pools and spillage
      dTDIN = dTDIN-(Nleak)*TDIN #N leaching losses
      
      list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3, dLIT_1_N, 
             dLIT_2_N, dMIC_1_N, dMIC_2_N, dSOM_1_N, dSOM_2_N, dSOM_3_N, dDIN,
             dTLIT_1, dTLIT_2, dTMIC_1, dTMIC_2, dTSOM_1, dTSOM_2, dTSOM_3, dTLIT_1_N, 
             dTLIT_2_N, dTMIC_1_N, dTMIC_2_N, dTSOM_1_N, dTSOM_2_N, dTSOM_3_N, dTDIN))
    })
  }
  
  
  SSoutput <<- data.frame(matrix(NA,nrow=14, ncol=19))
  colnames(SSoutput) = c("Site", "MICr","MICK","LITm","LITs","SOMp","SOMc","SOMa",
                         "MICrN","MICKN","LITmN","LITsN","SOMpN","SOMcN","SOMaN","InorgN",
                         "NminTot","NminNet","Resp")
  
  
  LITtype  <<- c('TRAEf', 'PIREf','THPLf','ACSAf','QUPRf','DRGLf')
  bagMET   <<- c(10.6, 36.2, 37.4, 56.8, 37.1, 49.3) #from Gordon's LitterCharacteristics.txt
  bagLIG   <<- c(16.2, 19.2, 26.7, 15.9, 23.5, 10.9) # % from Gordon's LitterCharacteristics.txt
  bagN     <<- c(0.38, 0.59, 0.62, 0.81, 1.03, 1.97) # %N 
  bagCN    <<- c(133.3,92.7, 83.1, 61.8, 50.5, 24.2)
  calcN    <<- (1 / bagCN) / 2.5 * 100
  calcMET  <<- 0.85 - 0.013 * bagLIG/calcN #as calculated in DAYCENT
  bagMET   <<- calcMET
  
  LTERdata <- data
  #LTERdata = read.csv("LTER_SITE_test.csv") #site level forcing variables
  
  ANPP_C  <<- LTERdata$ANPP      # convert to gC/m2/y from g/m2/y
  strSite <<- as.character(LTERdata$Site)  #convert site names to string
  nsites  <<- length(strSite)
  npts   <<- 6*10*14 #6 litter * 10 years * 14 sites
  clay  <<- LTERdata$CLAY2/100 
  d1 <<- LTERdata$depth
  tsoi  <<- LTERdata$MAT
  nsites <<- length(LTERdata$Site)
  lig   <<- LTERdata$LIG/100
  Nnew  <<- 1/LTERdata$CN/2.5             #N in litter additions
  fMET1 <<- 0.85 - 0.013 * lig / Nnew    #as partitioned in Daycent
  
  
  siteSpecificParameters = function() {
    
    #Parameters related to inputs
    EST_LIT_in  <<- ANPP_C[s] / (365*24)   #gC/m2/h (from g/m2/y, Knapp et al. Science 2001)
    BAG_LIT_in  <<- 100      #gC/m2/h
    soilDepth       <<- d1[s]
    h2y         <<- 24*365
    MICROtoECO  <<- soilDepth * 1e4 * 1e6 / 1e6   #mgC/cm3 to kgC/km2
    EST_LIT     <<- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h 
    BAG_LIT     <<- BAG_LIT_in  * 1e3 / 1e4    #mgC/cm2/h
    fMET        <<- fMET1[s]
    Inputs        <<- array(NA, dim=2)              #Litter inputs to MET/STR
    Inputs[1]     <<- (EST_LIT / soilDepth) * fMET      #partitioned to layers
    Inputs[2]     <<- (EST_LIT / soilDepth) * (1-fMET)
    FI       <<- c(0.05,0.3)#c(0.05, 0.05)#
    
    BAG      <<- array(NA, dim=c(6,2))              #litter BAG inputs to MET/STR
    for (i in 1:6) {
      BAG[i,1]   <<- (BAG_LIT / soilDepth) * bagMET[i]      #partitioned to layers
      BAG[i,2]   <<- (BAG_LIT / soilDepth) * (1-bagMET[i])
    }
    
    #Parameters related to stabilization mechanisms
    fCLAY       <<- clay[s]
    fPHYS    <<- 0.1 * c(.15 * exp(1.3*fCLAY), 0.1 * exp(0.8*fCLAY)) #Sulman et al. 2018
    fCHEM    <<- 3*c(0.1 * exp(-3*fMET) * 1, 0.3 * exp(-3*fMET) * 1) 	#Sulman et al. 2018 #fraction to SOMc
    fAVAI    <<- 1-(fPHYS + fCHEM)
    desorb   <<- 2e-5      * exp(-4.5*(fCLAY)) #Sulman et al. 2018
    desorb   <<- 0.05*desorb
    Nleak   <<- 0.2#.1   #This is the proportion N lost from DIN pool at each hourly time step.
    
    #Parameters related to microbial physiology and pool stoichiometry
    CUE        <<- c(LTERdata$CUE1[s], LTERdata$CUE2[s], 
                     LTERdata$CUE3[s], LTERdata$CUE4[s])  #for LITm and LITs entering MICr and MICK, respectively
    NUE        <<- .85         #Nitrogen stoichiometry of fixed pools
    CN_m        <<- 15
    CN_s        <<- (LTERdata$CN[s]-CN_m*fMET)/(1-fMET)
    CN_s_BAG    <<-  (bagCN-CN_m*bagMET)/(1-bagMET)
    CN_r        <<-6#5
    CN_K        <<-10#8
    
    turnover      <<- c(5.2e-4*exp(0.3*(fMET)), 2.4e-4*exp(0.1*(fMET)))	#WORKS BETTER FOR N_RESPONSE RATIO
    turnover_MOD1 <<- sqrt(ANPP_C[s]/100)  #basicaily standardize against NWT
    turnover_MOD1[turnover_MOD1 < 0.6] <<- 0.6 # correction not used in LIDET resutls 
    turnover_MOD1[turnover_MOD1 > 1.3] <<- 1.3      #Sulman et al. 2018
    turnover      <<- turnover * turnover_MOD1
    turnover <<- turnover/2.2
    turnover <<- turnover^2*0.55/(.45*Inputs)
    densDep <<- 2#1 #This exponent controls the density dependence of microbial turnover. Currently anything other than 1 breaks things.
    
    fracNImportr  <<-  0 #No N import for r strategists
    fracNImportK  <<-  0.2 #Only K strategists can import N
    
    #Parameters related to temperature-sensitive enzyme kinetics
    TSOI        <<- tsoi[s]   
    #Calculate Vmax & (using parameters from German 2012)
    #from "gamma" simulations "best", uses max Vslope, min Kslope
    Vslope   <<- array(NA,dim=6)
    Vslope[1]<<- 0.043 #META LIT to MIC_1
    Vslope[2]<<- 0.043 #STRU LIT to MIC_1 
    Vslope[3]<<- 0.063 #AVAI SOM to MIC_1 
    Vslope[4]<<- 0.043 #META LIT to MIC_2 
    Vslope[5]<<- 0.063 #STRU LIT to MIC_2 
    Vslope[6]<<- 0.063 #AVAI SOM to MIC_2 
    Vint     <<- 5.47
    aV       <<- 8e-6
    aV       <<- aV*.06 #Forward
    Vmax     <<- exp(TSOI * Vslope + Vint) * aV
    
    Kslope   <<- array(NA,dim=6)
    Kslope[1]<<- 0.017 #META LIT to MIC_1
    Kslope[2]<<- 0.027 #STRU LIT to MIC_1 
    Kslope[3]<<- 0.017 #AVAI SOM to MIC_1 
    Kslope[4]<<- 0.017 #META LIT to MIC_2
    Kslope[5]<<- 0.027 #STRU LIT to MIC_2
    Kslope[6]<<- 0.017 #AVAI SOM to MIC_2
    Kint     <<- 3.19
    aK       <<- 10
    aK       <<- aK/20 #Forward
    Km       <<- exp(Kslope * TSOI + Kint) * aK
    
    #Enzyme kinetic modifiers:
    k        <<- 2.0    #2.0			#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
    a        <<- 2.0    #2.2			#increased from 4.0 to 4.5
    cMAX     <<- 1.4                    #ORIG 1.4 Maximum CHEM SOM scalar w/   0% Clay 
    cMIN     <<- 1.2                    #ORIG 1.4 Minimum CHEM SOM scalar w/ 100% Clay 
    cSLOPE   <<- cMIN - cMAX            #Slope of linear function of cSCALAR for CHEM SOM  
    pSCALAR  <<- a * exp(-k*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
    
    #------------!!MODIFIERS AS IN MIMICS2_b!!---------------
    MOD1     <<- c(10, 2*.75, 10, 3, 3*.75, 2) 
    MOD2     <<- c( 8, 2 ,4 * pSCALAR, 2, 4, 6 * pSCALAR) 	
    
    VMAX     <<- Vmax * MOD1
    KM       <<- Km / MOD2
    KO       <<- c(6,6)     #Values from Sulman et al. 2018
  }
  
  
  initializePools = function() {
    LIT_1  <<- 1
    LIT_2  <<- 1
    MIC_1  <<- 1
    MIC_2  <<- 1
    SOM_1  <<- 1
    SOM_2  <<- 1
    SOM_3  <<- 1
    
    LIT_1_N  <<- .1
    LIT_2_N  <<- .1
    MIC_1_N  <<- .1
    MIC_2_N  <<- .1
    SOM_1_N  <<- .1
    SOM_2_N  <<- .1
    SOM_3_N  <<- .1
    DIN      <<- .1
    
    LITmin  <<- array(NA, dim=4)
    MICtrn  <<- array(NA, dim=6)
    SOMmin  <<- array(NA, dim=2)
    DEsorb  <<- array(NA, dim=1)
    OXIDAT  <<- array(NA, dim=1)
    LITminN   <<- array(NA, dim=4)
    MICtrnN   <<- array(NA, dim=6)
    SOMminN   <<- array(NA, dim=2)
    DEsorbN   <<- array(NA, dim=1)
    OXIDATN   <<- array(NA, dim=1)
    DINup     <<- array(NA, dim=2)
    Overflow  <<- array(NA, dim=2)
    Nspill    <<- array(NA, dim=2)
    CNup      <<- array(NA, dim=2)
    upMIC_1   <<-  array(NA, dim=1)
    upMIC_1_N <<-  array(NA, dim=1)
    upMIC_2   <<-  array(NA, dim=1)
    upMIC_2_N <<-  array(NA, dim=1)
  }
  
  findSteadyState = function() {
    Tpars <<- c( Inputs = Inputs, VMAX = VMAX, KM = KM, CUE = CUE, 
                 fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                 turnover = turnover, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                 desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT,
                 LITminN = LITminN, SOMminN = SOMminN, MICtrnN = MICtrnN,
                 DEsorbN = DEsorbN, OXIDATN = OXIDATN,
                 KO = KO,
                 CNup=CNup, DINup=DINup,Nspill=Nspill, Overflow=Overflow, 
                 upMIC_1=upMIC_1, upMIC_1_N=upMIC_1_N,
                 upMIC_2=upMIC_2, upMIC_2_N=upMIC_2_N,
                 NUE=NUE, CN_m=CN_m, CN_s=CN_s, CN_r=CN_r, CN_K=CN_K,Nleak=Nleak,densDep=densDep)
    Ty    <<- c(LIT_1 = LIT_1, LIT_2 = LIT_2, 
                MIC_1 = MIC_1, MIC_2 = MIC_2, 
                SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3,
                LIT_1_N = LIT_1_N, LIT_2_N = LIT_2_N, 
                MIC_1_N = MIC_1_N, MIC_2_N = MIC_2_N, 
                SOM_1_N = SOM_1_N, SOM_2_N = SOM_2_N, SOM_3_N = SOM_3_N,
                DIN = DIN)
    test  <<- stode(y = Ty, time = 1E8, fun = FXEQ, parms = Tpars, positive = TRUE)
    test[[1]]
  }
  
  
  
  ode_steadyState <- function(){
    Tpars <<- c( Inputs = Inputs, VMAX = VMAX, KM = KM, CUE = CUE, 
                 fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                 turnover = turnover, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                 desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT,
                 LITminN = LITminN, SOMminN = SOMminN, MICtrnN = MICtrnN,
                 DEsorbN = DEsorbN, OXIDATN = OXIDATN,
                 KO = KO,
                 CNup=CNup, DINup=DINup,Nspill=Nspill, Overflow=Overflow, 
                 upMIC_1=upMIC_1, upMIC_1_N=upMIC_1_N,
                 upMIC_2=upMIC_2, upMIC_2_N=upMIC_2_N,
                 NUE=NUE, CN_m=CN_m, CN_s=CN_s, CN_r=CN_r, CN_K=CN_K,Nleak=Nleak,densDep=densDep)
    Ty    <<- c(LIT_1 = LIT_1, LIT_2 = LIT_2, 
                MIC_1 = MIC_1, MIC_2 = MIC_2, 
                SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3,
                LIT_1_N = LIT_1_N, LIT_2_N = LIT_2_N, 
                MIC_1_N = MIC_1_N, MIC_2_N = MIC_2_N, 
                SOM_1_N = SOM_1_N, SOM_2_N = SOM_2_N, SOM_3_N = SOM_3_N,
                DIN = DIN)
    
    t_1 <- 1:times
    
    ode_solve <- ode(y = Ty, times = t_1, func = FXEQ, parms = Tpars)
    
    return(ode_solve)
  }
  
  
  outputSteadyState <- function(){
    LIT_1    <<- test[[1]][[1]]
    LIT_2    <<- test[[1]][[2]]
    MIC_1    <<- test[[1]][[3]]
    MIC_2    <<- test[[1]][[4]]
    SOM_1    <<- test[[1]][[5]]
    SOM_2    <<- test[[1]][[6]]
    SOM_3    <<- test[[1]][[7]]
    
    LIT_1_N    <<- test[[1]][[8]]
    LIT_2_N    <<- test[[1]][[9]]
    MIC_1_N    <<- test[[1]][[10]]
    MIC_2_N    <<- test[[1]][[11]]
    SOM_1_N    <<- test[[1]][[12]]
    SOM_2_N    <<- test[[1]][[13]]
    SOM_3_N    <<- test[[1]][[14]]
    DIN        <<- test[[1]][[15]]
    Site <<- strSite[s]
    
    #Carbon fluxes
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #decomp of SOMa by MIC_1
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of SOMa by MIC_2
    
    #Nitrogen fluxes
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Flows to and from MIC_1
    LITminN[1] =  LITmin[1]*LIT_1_N/(LIT_1+1e-100)#LITmin[1]/CN_m
    LITminN[2] =  LITmin[2]*LIT_2_N/(LIT_2++1e-100)#LITmin[2]/CN_s
    SOMminN[1] =  SOMmin[1]*(SOM_3_N/(SOM_3++1e-100))#SOMmin[1]*(SOM_3_N/SOM_3)#*relRateN
    LITminN[3] =  LITmin[3]*LIT_1_N/(LIT_1+1e-100)#LITmin[3]/CN_m
    LITminN[4] =  LITmin[4]*LIT_2_N/(LIT_2+1e-100)#LITmin[4]/CN_s
    SOMminN[2] =  SOMmin[2]*(SOM_3_N/(SOM_3+1e-100))#SOMmin[2]*(SOM_3_N/SOM_3)#*relRateN
    DINup[1]   = (1-Nleak)*DIN*MIC_1/(MIC_1+MIC_2+1e-100) #Partitions available DIN between microbial pools based on relative biomass
    DINup[2]   = (1-Nleak)*DIN*MIC_2/(MIC_1+MIC_2+1e-100)
    
    #####
    upMIC_1    = CUE[1]*(LITmin[1] + SOMmin[1]) + CUE[2]*(LITmin[2])
    upMIC_1_N  = NUE*(LITminN[1] + SOMminN[1]) + NUE*(LITminN[2]) + DINup[1]
    CNup[1]    = (upMIC_1)/(upMIC_1_N+1e-100) #Trying to optimize run speed here by avoiding /0
    Overflow[1] = (upMIC_1) - (upMIC_1_N)*min(CN_r, CNup[1])
    Nspill[1]   = (upMIC_1_N) - (upMIC_1)/max(CN_r, CNup[1])
    
    upMIC_2    = CUE[3]*(LITmin[3] + SOMmin[2]) + CUE[4]*(LITmin[4])
    upMIC_2_N  = NUE*(LITminN[3] + SOMminN[2]) + NUE*(LITminN[4]) + DINup[2]
    CNup[2]    = (upMIC_2)/(upMIC_2_N+1e-100)
    Overflow[2] = (upMIC_2) - (upMIC_2_N)*min(CN_K, CNup[2])
    Nspill[2]   = (upMIC_2_N) - (upMIC_2)/max(CN_K, CNup[2])
    ######
    
    dLIT_1 = Inputs[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dLIT_2 = Inputs[2]*(1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_1 = CUE[1]*(LITmin[1] + SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3]) - Overflow[1]
    dMIC_2 = CUE[3]*(LITmin[3] + SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6]) - Overflow[2] 
    dSOM_1 = Inputs[1]*FI[1] + MICtrn[1] + MICtrn[4] - DEsorb 
    dSOM_2 = Inputs[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
    dSOM_3 = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    dLIT_1_N = Inputs[1]*(1-FI[1])/CN_m - LITminN[1] - LITminN[3]
    dLIT_2_N = Inputs[2]*(1-FI[2])/CN_s - LITminN[2] - LITminN[4]
    dMIC_1_N = NUE*(LITminN[1] + SOMminN[1]) + NUE*(LITminN[2]) - sum(MICtrnN[1:3]) + DINup[1] - Nspill[1]
    dMIC_2_N = NUE*(LITminN[3] + SOMminN[2]) + NUE*(LITminN[4]) - sum(MICtrnN[4:6]) + DINup[2] - Nspill[2]
    dSOM_1_N = Inputs[1]*FI[1]/CN_m + MICtrnN[1] + MICtrnN[4] - DEsorbN
    dSOM_2_N = Inputs[2]*FI[2]/CN_s + MICtrnN[2] + MICtrnN[5] - OXIDATN
    dSOM_3_N = MICtrnN[3] + MICtrnN[6] + DEsorbN + OXIDATN - SOMminN[1] - SOMminN[2]
    
    dDIN = (1-NUE)*(LITminN[1] + LITminN[2] + SOMminN[1]) +  #Inputs from r decomp
      (1-NUE)*(LITminN[3] + LITminN[4] + SOMminN[2]) +  #Inputs from K decomp
      Nspill[1] + Nspill[2] - DINup[1] - DINup[2]    #Uptake to microbial pools and spillage
    LeachingLoss = Nleak*DIN
    dDIN = dDIN-LeachingLoss #N leaching losses    
    
    Resp=(1-CUE[1])*(LITmin[1]+ SOMmin[1]) + (1-CUE[2])*(LITmin[2]) +
      (1-CUE[3])*(LITmin[3]+ SOMmin[2]) + (1-CUE[4])*(LITmin[4]) + Overflow[1] + Overflow[2]
    NminTot=(1-NUE)*(LITminN[1] + LITminN[2] + SOMminN[1]) +  #Inputs from r decomp
      (1-NUE)*(LITminN[3] + LITminN[4] + SOMminN[2]) +  #Inputs from K decomp
      Nspill[1] + Nspill[2]
    NminNet=NminTot-DINup[1]-DINup[2]
    
    newrow <<- cbind.data.frame(Site, MIC_1,MIC_2,LIT_1,LIT_2,SOM_1,SOM_2,SOM_3,
                                MIC_1_N,MIC_2_N,LIT_1_N,LIT_2_N,SOM_1_N,SOM_2_N,SOM_3_N,DIN, 
                                NminTot,NminNet,Resp,
                                stringsAsFactors = FALSE)
    SSoutput[s,] <<- newrow
  }
  
  x1 <- range(1,length(LTERdata$Site))
  t1 <- as.list(x1[1]:x1[2])
  valiSites <- t1
  out <- data.frame(matrix(nrow = length(valiSites), ncol = 4))
  colnames(out) <- c('SOM_1', 'SOM_2', 'SOM_3')
  
  for(s in valiSites){
    #Find steady state solution, plot to examine and output values to table
    siteSpecificParameters()
    initializePools()
    
    new_out_ODE <- ode_steadyState()
    assign(paste0('meanODE', s), data.frame(new_out_ODE))
    temp <- data.frame(new_out_ODE)
    len <- length(temp$SOM_1)
    out$SOM_1[s] <- temp$SOM_1[len]
    out$SOM_2[s] <- temp$SOM_2[len]
    out$SOM_3[s] <- temp$SOM_3[len]
    
    
  }
  
  return(out)
}


