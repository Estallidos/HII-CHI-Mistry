# Filename: HCm_Teff_v3.0.py


print ' ---------------------------------------------------------------------'
print 'This is HII-CHI-mistry_Teff v. 3.0'
print ' See Perez-Montero et al (2019) for details'
print  ' Insert the name of your input text file with the following columns:'
print '12+log(O/H), 3727 [OII], 5007 [OIII], 6725 [SII], 9069 [SIII]'
print 'with their corresponding errors in adjacent columns'
print 'relative to Hbeta or 0 for missing information.'
print '---------------------------------------------------------------------'
print ''

import string
import numpy as np


# Input file reading

input00 = raw_input('Insert input file name:')
input0 = np.loadtxt(input00)

output = []


# Iterations for Montecarlo error derivation
n = 25



# Reading of models grids. These can be changed


print ''
question = True
while question:
   print'-------------------------------------------------'
   print '(1) Plane-parallel geometry'
   print '(2) Spherical geometry'
   print'-------------------------------------------------'
   geo = int(raw_input('Choose geometry of the models:'))

   if geo==1:
      geo_type = 'Plane-parallel geometry'
      grid = np.loadtxt('C17_WMb_Teff_30-55_pp.dat')
      print ''
      question = False
   elif geo==2:
      geo_type = 'Spherical geometry'
      grid = np.loadtxt('C17_WMb_Teff_30-55_sph.dat')
      print ''
      question = False


   print ''




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
      tol_max = 1e2


      OII_3727_obs = np.random.normal(tab[2],tab[3]+1e-3)
      if OII_3727_obs <= 0 : OII_3727_obs = 0
      OIII_5007_obs = np.random.normal(tab[4],tab[5]+1e-3)
      if OIII_5007_obs <= 0: OIII_5007_obs = 0
      SII_6725_obs = np.random.normal(tab[6],tab[7]+1e-3)
      if SII_6725_obs <= 0: SII_6725_obs = 0
      SIII_9069_obs = np.random.normal(tab[8],tab[9]+1e-3)
      if SIII_9069_obs <= 0: SIII_9069_obs = 0
      if tab[2] == 0 or tab[4] == 0 or OII_3727_obs == 0 or OIII_5007_obs == 0:
         O2O3_obs = -10
         R23_obs = -10
      else:
         O2O3_obs = np.log10(OII_3727_obs / OIII_5007_obs) 
         R23_obs = np.log10(OII_3727_obs + OIII_5007_obs )
      if tab[6] == 0 or tab[8] == 0 or SII_6725_obs == 0 or SIII_9069_obs == 0:
         S2S3_obs = -10
         S23_obs = -10
      else:
         S2S3_obs = np.log10(SII_6725_obs / SIII_9069_obs )
         S23_obs = (SII_6725_obs + SIII_9069_obs )

      

# Interpolation of grid at specific O/H


      if tab[0] > 0:
         OH = np.random.normal(tab[0],tab[1]+1e-3)
         OH_mc.append(OH)


         grid_T0 = []
         if OH <= 7.1:
            OH = 7.1
            i0 = 0
            i1 = 72
         elif OH >= 7.1 and OH < 7.4:
            i0 = 0
            i1 = 72
         elif OH >= 7.4 and OH < 7.7:
            i0 = 72
            i1 = 144
         elif OH >= 7.7 and OH < 8.0:
            i0 = 144
            i1 = 216
         elif OH >= 8.0 and OH < 8.3:
            i0 = 216
            i1 = 288
         elif OH >= 8.3 and OH < 8.6:
            i0 = 288
            i1 = 360
         elif OH >= 8.6 and OH < 8.9:
            i0 = 360
            i1 = 432
         elif OH >= 8.9:
            OH = 8.9
            i0 = 360
            i1 = 432
       
         for x in range(0,72):
            for y in range(0,7):
               grid_T0.append(grid[i0+x,y]*np.abs(0.3-OH+grid[i0,0])/0.3+grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)
            
      #         grid_T0.append(grid[i0+x,y]*np.abs(0.3-grid[i0,0]+OH)/0.3 + grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)


         grid_T = np.reshape(grid_T0,(72,7))

      else:
         OH = 0
         OH_mc.append(OH)
         grid_T = grid
   


#      np.savetxt('int_models.dat',grid_T)

# Calculation of T and log U


      if S2S3_obs == -10 and O2O3_obs == -10:
         Teff = 0
         logU = 0
      else:
         CHI_O2O3 = 0
         CHI_S2S3 = 0
         CHI_S23 = 0

         for index in grid_T:
            if index[5] == 0 or index[6] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            elif S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            else:
               CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/S2S3_obs
               CHI_S23 = (index[5]+index[6]-S23_obs)**2/S23_obs
            if index[3] == 0 or index[4] == 0:
               CHI_O2O3 = tol_max
            elif O2O3_obs == -10:
               CHI_O2O3 = 0
            else:
               CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/O2O3_obs


            CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 )**0.5
         
            Teff_p = index[1]*(1/CHI_Teff)**2 + Teff_p
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

         for index in grid_T:
            if index[5] == 0 or index[6] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            elif S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            else:
               CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/S2S3_obs
               CHI_S23 = (index[5]+index[6]-S23_obs)**2/S23_obs
            if index[3] == 0 or index[4] == 0:
               CHI_O2O3 = tol_max
            elif O2O3_obs == -10:
               CHI_O2O3 = 0
            else:
               CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/O2O3_obs


            CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 )**0.5
         
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

lineas_header = [' HII-CHI-mistry v.3.0 output file',' Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+geo_type,'','O2Hb eO2Hb O3Hb  eO3Hb  eO3Hb S2Hb  eS2Hb S3Hb  eS3HbO/H   eO/H  Teff  eTeff logU elogU']





header = '\n'.join(lineas_header)



np.savetxt('output.dat',out,fmt=' '.join(['%.3f']*8+['%.2f']*2+['%.i']*2+['%.2f']*2),header=header)
print '________________________________'
print 'Results are stored in output.dat'


