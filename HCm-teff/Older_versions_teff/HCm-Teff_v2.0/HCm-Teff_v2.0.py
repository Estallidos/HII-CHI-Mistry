
# Filename: HCm_Teff_v01.py


print ' ---------------------------------------------------------------------'
print 'This is HII-CHI-mistry_Teff v. 2.1'
print ' See Perez-Montero, E. et al. (in prep.) for details'
print  ' Insert the name of your input text file with following columns:'
print '12+log(O/H), 3727 [OII], 5007 [OIII], 6725 [SII], 9069 [SIII]'
print 'elative to Beta and their corresponding error.'
print '---------------------------------------------------------------------'
print ''

import string
import numpy as np

# Reading of models grids. These can be changed

grid = np.loadtxt('C13_WMb_Teff_30-55.dat')

# Input file reading

input0 = np.loadtxt(raw_input('Insert input file name:'))
output = []



if input0.shape == (10,):
   input1 = [0,0,0,0,0,0,0,0,0,0,input0[0],input0[1],input0[2],input0[3],input0[4],input0[5],input0[6],input0[7],input0[8],input0[9]]
   input = np.reshape(input1,(2,10))
else:
   input = input0


print 'Reading grids ....'
print ''
print ''
print '---------------------------------------'
print '(%)   12+log(O/H)  T_eff(K)    log(U)'
print '---------------------------------------'


# Beginning of loop of calculation

count = 0
for tab in input:
   count = count + 1


   n = 100
   OH_mc = []
   Teff_mc = []
   logU_mc = []
   eOH_mc = []
   eTeff_mc = []
   elogU_mc = []

   for monte in range(0,n,1):


      OH_p = 0
      logU_p = 0
      Teff_p = 0
      den_OH = 0
      den_Teff = 0
      OH_e = 0
      Teff_e = 0
      logU_e = 0
      den_OH_e = 0
      den_Teff_e = 0
      tol_max = 1e3


      if tab[2] == 0:
         OII_3727_obs = 0 
      else:
         OII_3727_obs = np.random.normal(tab[2],tab[3]+1e-4)
         if OII_3727_obs <= 0: OII_3727_obs = 0
      if tab[4] == 0:
         OIII_5007_obs = 0
      else:
         OIII_5007_obs = np.random.normal(tab[4],tab[5]+1e-4)
         if OIII_5007_obs <= 0: OIII_5007_obs = 0
      if tab[6] == 0:
         SII_6725_obs = 0
      else:
         SII_6725_obs = np.random.normal(tab[6],tab[7]+1e-4)
         if SII_6725_obs <= 0: SII_6725_obs = 0
      if tab[8] == 0:
         SIII_9069_obs = 0
      else:   
         SIII_9069_obs = np.random.normal(tab[8],tab[9]+1e-4)
         if SIII_9069_obs <= 0: SIII_9069_obs = 0
      if OII_3727_obs == 0 or OIII_5007_obs == 0:
         O2O3_obs = -10
         R23_obs = -10
      else:
         O2O3_obs = np.log10(OII_3727_obs / OIII_5007_obs) 
         R23_obs = np.log10(OII_3727_obs + OIII_5007_obs )
      if SII_6725_obs == 0 or SIII_9069_obs == 0:
         S2S3_obs = -10
         S23_obs = -10
      else:
         S2S3_obs = np.log10(SII_6725_obs / SIII_9069_obs )
         S23_obs = (SII_6725_obs + SIII_9069_obs )
      if O2O3_obs == -10 or S2S3_obs == -10:
         eta_obs = -10
         ieta_obs = -10
      else:
         eta_obs = O2O3_obs - S2S3_obs
         ieta_obs = S2S3_obs - np.log10(OIII_5007_obs/OII_3727_obs)

      

# Interpolation of grid at specific O/H


      if tab[0] > 0:
         OH = np.random.normal(tab[0],tab[1]+1e-3)
         OH_mc.append(OH)

         grid_T0 = []
         if OH <= 7.1:
            OH = 7.1
            i0 = 0
            i1 = 40
         elif OH >= 7.1 and OH < 7.4:
            i0 = 0
            i1 = 40
         elif OH >= 7.4 and OH < 7.7:
            i0 = 40
            i1 = 80
         elif OH >= 7.7 and OH < 8.0:
            i0 = 80
            i1 = 120
         elif OH >= 8.0 and OH < 8.3:
            i0 = 120
            i1 = 160
         elif OH >= 8.3 and OH < 8.6:
            i0 = 160
            i1 = 200
         elif OH >= 8.6 and OH < 8.9:
            i0 = 200
            i1 = 240
         elif OH >= 8.9:
            OH = 8.9
            i0 = 200
            i1 = 240
       
         for x in range(0,40):
            for y in range(0,7):
               grid_T0.append(grid[i0+x,y]*np.abs(0.3-OH+grid[i0,0])/0.3+grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)
            
      #         grid_T0.append(grid[i0+x,y]*np.abs(0.3-grid[i0,0]+OH)/0.3 + grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)


         grid_T = np.reshape(grid_T0,(40,7))

      else:
         OH = 0
         OH_mc.append(OH)
         grid_T = grid
   



# Calculation of T and log U



      if S2S3_obs == -10 and O2O3_obs == -10:
         Teff = 0
         logU = 0
      else:
         CHI_O2O3 = 0
         CHI_S2S3 = 0
         CHI_S23 = 0
         CHI_ETA = 0
         CHI_iETA = 0

         for index in grid_T:
            if index[5] == 0 or index[6] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            elif S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            else:
               CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/np.log10(index[5]/index[6])
               CHI_S23 = (index[5]+index[6]-S23_obs)**2/((index[5]+index[6]))
            if index[3] == 0 or index[4] == 0:
               CHI_O2O3 = tol_max
            elif O2O3_obs == -10:
               CHI_O2O3 = 0
            else:
               CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/np.log10(index[3]/index[4])
            if index[5] == 0 or index[6] == 0 or index[3] == 0 or index[4] == 0:
               CHI_ETA = tol_max
               CHI_iETA = tol_max    
            elif eta_obs == -10:
               CHI_ETA = 0
               CHI_iETA = 0
            else:
               CHI_ETA = ((np.log10(index[3]/index[4])-np.log10(index[5]/index[6])) - eta_obs)**2 /\
                  np.abs(np.log10(index[3]/index[4])-np.log10(index[5]/index[6]))
               CHI_iETA = ((np.log10(index[5]/index[6])-np.log10(index[4]/index[3])) - ieta_obs)**2 /\
                  np.abs(np.log10(index[5]/index[6])-np.log10(index[4]/index[3]))



            if eta_obs == -10:
               CHI_Teff = (CHI_O2O3**2 + CHI_S2S3**2 )**0.5
            else:
               CHI_Teff = (CHI_ETA**2 + CHI_iETA**2 )**0.5

         
         
            Teff_p = index[1] *(1/CHI_Teff)**2 + Teff_p
            logU_p = index[2] *(1/CHI_Teff)**2 + logU_p         
            den_Teff = (1/CHI_Teff)**2 + den_Teff
         Teff = Teff_p / den_Teff 
         logU = logU_p / den_Teff



# Calculation of T and log U errors


      if S2S3_obs == -10 and O2O3_obs == -10:
         eTeff = 0
         elogU = 0
      else:
         CHI_O2O3 = 0
         CHI_S2S3 = 0
         CHI_S23 = 0
         CHI_ETA = 0
         CHI_iETA = 0


         for index in grid_T:
            if index[5] == 0 or index[6] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            elif S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            else:
               CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/np.log10(index[5]/index[6])
               CHI_S23 = (index[5]+index[6]-S23_obs)**2/((index[5]+index[6]))
            if index[3] == 0 or index[4] == 0:
               CHI_O2O3 = tol_max
            elif O2O3_obs == -10:
               CHI_O2O3 = 0
            else:
               CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/np.log10(index[3]/index[4])
            if index[5] == 0 or index[6] == 0 or index[3] == 0 or index[4] == 0:
               CHI_ETA = tol_max
               CHI_iETA = tol_max    
            elif eta_obs == -10:
               CHI_ETA = 0
               CHI_iETA = 0
            else:
               CHI_ETA = ((np.log10(index[3]/index[4])-np.log10(index[5]/index[6])) - eta_obs)**2 /\
                  np.abs(np.log10(index[3]/index[4])-np.log10(index[5]/index[6]))
               CHI_iETA = ((np.log10(index[5]/index[6])-np.log10(index[4]/index[3])) - ieta_obs)**2 /\
                  np.abs(np.log10(index[5]/index[6])-np.log10(index[4]/index[3]))



            if eta_obs == -10:
               CHI_Teff = (CHI_O2O3**2 + CHI_S2S3**2 )**0.5
            else:
               CHI_Teff = (CHI_ETA**2 + CHI_iETA**2 )**0.5
         
         
            Teff_e = np.abs(index[1] - Teff) * (1/CHI_Teff)**2 + Teff_e
            logU_e = np.abs(index[2] - logU) * (1/CHI_Teff)**2 + logU_e         
            den_Teff_e = 1 * (1/CHI_Teff**2) + den_Teff_e

         eTeff = Teff_e / den_Teff_e 
         elogU = logU_e / den_Teff_e

      Teff_mc.append(Teff)
      logU_mc.append(logU)
      eTeff_mc.append(eTeff)
      elogU_mc.append(elogU)


   

   if tab[0] > 0:
      OHf = tab[0]
      eOHf = tab[1]
   else: 
      OHf = 0
      eOHf = 0
   Tefff = np.mean(Teff_mc)
   eTefff = (np.std(Teff_mc)**2 + np.mean(eTeff_mc)**2)**0.5
   logUf = np.mean(logU_mc)
   elogUf = (np.std(logU_mc)**2+np.mean(elogU_mc)**2)**0.5

   
   output.append(tab[2])
   output.append(tab[3])
   output.append(tab[4])
   output.append(tab[5])
   output.append(tab[6])
   output.append(tab[7])
   output.append(tab[8])
   output.append(tab[9])

   output.append(OHf)
   output.append(eOHf)
   output.append(Tefff)
   output.append(eTefff)
   output.append(logUf)
   output.append(elogUf)  

   if input0.shape >= (10,) and count == 1: continue

   print round(100*(count)/float(len(input)),1),'%','', round(OHf,2), round(eOHf,2), 100*int(Tefff/100),100*int(eTefff/100),round(logUf,2),round(elogUf,2) 


out = np.reshape(output,(len(input),14))

np.savetxt('output.dat',out,fmt='%.4f')
print '________________________________'
print 'Results are stored in output.dat'



















