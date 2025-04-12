# Filename: HCm_UV_v3.0.py

print ' ---------------------------------------------------------------------'
print ' This is HII-CHI-mistry for UV version 3.0'
print ' See Perez-Montero, & Amorin (2017) for details'
print  ' Insert the name of your input text file with following columns:'
print ' Lya 1216, CIV 1549, OIII 1665, CIII 1909, Hb 4861, OIII 5007'
print 'in arbitrary units and reddening corrected. Each column must be'
print 'followed by its corresponding flux error.'
print '---------------------------------------------------------------------'



import string
import numpy as np

# Iterations for Montecarlo error derivation
n = 50


# Reading of models grids. These can be changed

grid1 = np.loadtxt('C13_cha_1Myr_v2.0_uv.dat')
grid2 = np.loadtxt('C13_cha_1Myr_logU_adapted_emp_v2.0_uv.dat')
grid3 = np.loadtxt('C13_cha_1Myr_logU-CO_adapted_emp_v2.0_uv.dat')

print ''
print 'The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O'
print ''
res_CO = 0.125


# Input file reading

input0 = np.loadtxt(raw_input('Insert input file name:'))
output = []


if input0.shape == (12,):
   input1 = [input0[0],input0[1],input0[2],input0[3],input0[4],input0[5],input0[6],input0[7],input0[8],input0[9],input0[10],input0[11],0,0,0,0,0,0,0,0,0,0,0,0]
   input = np.reshape(input1,(2,12))
else:
   input = input0


print 'Reading grids ....'
print ''
print ''
print '----------------------------------------------------------------'
print '(%)   Grid  12+log(O/H)  log(C/O)    log(U)'
print '-----------------------------------------------------------------'

# Beginning of loop of calculation

count = 0
for tab in input:
   count = count + 1


   OH_mc = []
   CO_mc = []
   logU_mc = []
   OHe_mc = []
   COe_mc = []
   logUe_mc = []  


   output.append(tab[0])
   output.append(tab[1])
   output.append(tab[2])
   output.append(tab[3])
   output.append(tab[4])
   output.append(tab[5])
   output.append(tab[6])
   output.append(tab[7])
   output.append(tab[8])
   output.append(tab[9])
   output.append(tab[10])
   output.append(tab[11])

   for monte in range(0,n,1):


      OH_p = 0
      logU_p = 0
      CO_p = 0
      den_OH = 0
      den_CO = 0
      OH_e = 0
      CO_e = 0
      logU_e = 0
      den_OH_e = 0
      den_CO_e = 0
      tol_max = 1e2

      if tab[0] == 0:
         Lya_1216_obs = 0
      else:
         Lya_1216_obs = np.random.normal(tab[0],tab[1]+1e-3)
         if Lya_1216_obs <= 0: Lya_1216_obs = 0
      if tab[2] == 0:
         CIV_1549_obs = 0
      else:
         CIV_1549_obs = np.random.normal(tab[2],tab[3]+1e-3)
         if CIV_1549_obs <= 0: CIV_15_obs = 0
      if tab[4] == 0:
         OIII_1665_obs = 0
      else:
         OIII_1665_obs = np.random.normal(tab[4],tab[5]+1e-3)
         if OIII_1665_obs <= 0: OIII_1665_obs = 0
      if tab[6] == 0:
         CIII_1909_obs = 0
      else:
         CIII_1909_obs = np.random.normal(tab[6],tab[7]+1e-3)
         if CIII_1909_obs <= 0: CIII_1909_obs = 0
      if tab[8] == 0:
         Hb_4861_obs = 0
      else:
         Hb_4861_obs = np.random.normal(tab[8],tab[9]+1e-3)
         if Hb_4861_obs <= 0: Hb_4861_obs = 0
      if tab[10] == 0:
         OIII_5007_obs = 0
      else:
         OIII_5007_obs = np.random.normal(tab[10],tab[1]+1e-3)
         if OIII_5007_obs <= 0: OIII_5007_obs = 0
      if OIII_1665_obs == 0 or OIII_5007_obs == 0:
         ROIII_obs = 0
      else:
         ROIII_obs = OIII_5007_obs/OIII_1665_obs
      if Lya_1216_obs == 0 or CIII_1909_obs == 0:
         C34_obs = 0
      else:
         C34_obs = (CIII_1909_obs + CIV_1549_obs) / (Lya_1216_obs)
      if CIII_1909_obs == 0 or OIII_1665_obs == 0:
         C3O3_obs = 0
      else:   
         C3O3_obs = (CIII_1909_obs + CIV_1549_obs) / (OIII_1665_obs)
      if CIII_1909_obs == 0 or CIV_1549_obs == 0:
         C3C4_obs = 0
      else:
         C3C4_obs = CIII_1909_obs/CIV_1549_obs
      if OIII_1665_obs == 0 or Lya_1216_obs == 0:
         O3_obs = 0
      else:
         O3_obs = OIII_1665_obs/Lya_1216_obs
      if CIII_1909_obs == 0 or Hb_4861_obs == 0:
         C34Hb_obs = 0
      else:
         C34Hb_obs = (CIII_1909_obs + CIV_1549_obs) / Hb_4861_obs
      if (OIII_1665_obs == 0 and OIII_5007_obs == 0) or Hb_4861_obs == 0:
         O3Hb_obs = 0
      else:
         O3Hb_obs = OIII_5007_obs / Hb_4861_obs
              
         



# Selection of grid
   
      if ROIII_obs > 0:
         grid = grid1
         if monte == n-1: output.append(1)
         grid_type = 1
      elif C3O3_obs > 0:
         grid = grid2
         if monte == n-1: output.append(2)
         grid_type = 2
      else:
         grid = grid3
         if monte == n-1: output.append(3)
         grid_type = 3



  # Calculation of C/O 

  
      if C3O3_obs == 0:
         CO = -10
      else:
         CHI_ROIII = 0
         CHI_C3O3 = 0   

         for index in grid:
            if ROIII_obs == 0:
               CHI_ROIII = 0
            elif index[4] == 0:
               CHI_ROIII = tol_max
            else:
               CHI_ROIII = (index[5]/index[4] - ROIII_obs)**2/(index[5]/index[4])
            if C3O3_obs == 0:
               CHI_C3O3 = 0          
            elif index[4] == 0 or index[6] == 0:
               CHI_C3O3 = tol_max
            else:
               CHI_C3O3 =(((index[6]+index[7])/index[4]) - C3O3_obs)**2/((index[6]+index[7])/index[4])

            CHI_CO = (CHI_ROIII**2 + CHI_C3O3**2 )**0.5
            CO_p = index[1] / CHI_CO + CO_p
            den_CO = 1 / CHI_CO + den_CO
         CO = CO_p / den_CO 

  # Calculation of C/O error

  
      if C3O3_obs == 0:
         eCO = 0
      else:
         CHI_ROIII = 0
         CHI_C3O3 = 0   

         for index in grid:
            if ROIII_obs == 0:
               CHI_ROIII = 0
            elif index[4] == 0:
               CHI_ROIII = tol_max
            else:
               CHI_ROIII = (index[5]/index[4] - ROIII_obs)**2/(index[5]/index[4])
            if C3O3_obs == 0:
               CHI_C3O3 = 0          
            elif index[4] == 0 or index[6] == 0:
               CHI_C3O3 = tol_max
            else:
               CHI_C3O3 =(((index[6]+index[7])/index[4]) - C3O3_obs)**2/((index[6]+index[7])/index[4])

            CHI_CO = (CHI_ROIII**2 + CHI_C3O3**2 )**0.5
            CO_e = (index[1] - CO)**2 / CHI_CO + CO_e
            den_CO_e = 1 / CHI_CO + den_CO_e  
         eCO = CO_e / den_CO_e 

      
# Calculation of O/H and log U

      if C34_obs == 0 and O3_obs == 0 and ROIII_obs == 0 and C34Hb_obs == 0 and O3Hb_obs == 0:
         OH = 0
         logU = 0
      else:
         CHI_ROIII = 0 
         CHI_C3C4 = 0
         CHI_C34 = 0
         CHI_C34Hb = 0
         CHI_O3 = 0
         CHI_O3Hb = 0

         CHI_OH = 0
         for index in grid:
            if CO > -10 and np.abs(index[1] - CO) > np.abs(eCO+res_CO):
               continue
            else:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index[4] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index[5]/index[4] - ROIII_obs)**2/(index[5]/index[4])
               if C34_obs == 0:
                  CHI_C34 = 0
               elif index[4] == 0:
                  CHI_C34 = tol_max
               else:
                  CHI_C34 = ((index[6]+index[7])/index[3] - C34_obs)**2/((index[6]+index[7])/index[3])
               if C34Hb_obs == 0:
                  CHI_C34Hb = 0
               elif index[4] == 0:
                  CHI_C34Hb = tol_max
               else:
                  CHI_C34Hb = (index[6]+index[7] - C34_obs)**2/(index[6]+index[7])
               if C3C4_obs == 0:
                  CHI_C3C4 = 0
               elif index[7] == 0 or index[6] == 0:
                  CHI_C3C4 = tol_max
               else:
                  CHI_C3C4 = (index[6]/index[7] - C3C4_obs)**2/(index[6]/index[7])
               if O3_obs == 0:
                  CHI_O3 = 0
               elif index[4] == 0 or index[3] == 0:
                  CHI_O3 = tol_max
               else:
                  CHI_O3 = (index[4]/index[3] - O3_obs)**2/(index[4]/index[3])
               if O3Hb_obs == 0:
                  CHI_O3Hb = 0
               elif index[4] == 0:
                  CHI_O3Hb = tol_max
               else:
                  CHI_O3Hb = (index[4] - O3Hb_obs)**2/index[4]


               if C34Hb_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34Hb**2 + CHI_O3Hb*2 + CHI_C3C4**2)**0.5
               else:	
                  CHI_OH = (CHI_ROIII**2 + CHI_C34**2 + CHI_O3*2 + CHI_C3C4**2 )**0.5


            OH_p = index[0] / CHI_OH + OH_p
            logU_p = index[2] / CHI_OH + logU_p
            den_OH = 1 / CHI_OH + den_OH
         
         OH = OH_p / den_OH
         logU = logU_p / den_OH

# Calculation of error of O/H and logU

      if C34_obs == 0 and O3_obs == 0 and ROIII_obs == 0 and C34Hb_obs == 0 and O3Hb_obs == 0:
         eOH = 0
         elogU = 0
      else:
         CHI_ROIII = 0 
         CHI_C3C4 = 0
         CHI_C34 = 0
         CHI_C34Hb = 0
         CHI_O3 = 0
         CHI_O3Hb = 0
         CHI_OH = 0

         for index in grid:
            if CO > -10 and np.abs(index[1] - CO) > np.abs(eCO+res_CO):
               continue
            else:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index[4] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index[5]/index[4] - ROIII_obs)**2/(index[5]/index[4])
               if C34_obs == 0:
                  CHI_C34 = 0
               elif index[4] == 0:
                  CHI_C34 = tol_max
               else:
                  CHI_C34 = ((index[6]+index[7])/index[3] - C34_obs)**2/((index[6]+index[7])/index[3])
               if C34Hb_obs == 0:
                  CHI_C34Hb = 0
               elif index[4] == 0:
                  CHI_C34Hb = tol_max
               else:
                  CHI_C34Hb = (index[6]+index[7] - C34_obs)**2/(index[6]+index[7])
               if C3C4_obs == 0:
                  CHI_C3C4 = 0
               elif index[7] == 0 or index[6] == 0:
                  CHI_C3C4 = tol_max
               else:
                  CHI_C3C4 = (index[6]/index[7] - C3C4_obs)**2/(index[6]/index[7])
               if O3_obs == 0:
                  CHI_O3 = 0
               elif index[4] == 0 or index[3] == 0:
                  CHI_O3 = tol_max
               else:
                  CHI_O3 = (index[4]/index[3] - O3_obs)**2/(index[4]/index[3])
               if O3Hb_obs == 0:
                  CHI_O3Hb = 0
               elif index[4] == 0:
                  CHI_O3Hb = tol_max
               else:
                  CHI_O3Hb = (index[4] - O3Hb_obs)**2/index[4]


               if C34Hb_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34Hb**2 + CHI_O3Hb*2 + CHI_C3C4**2)**0.5
               else:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34**2 + CHI_O3*2 + CHI_C3C4**2 )**0.5

            if CHI_OH == 0:
               OH_e = OH_e
               logU_e = logU_e
               den_OH_e = den_OH_e
            else:
               OH_e = (index[0] - OH)**2 / CHI_OH + OH_e
               logU_e = (index[2] - logU)**2 / CHI_OH + logU_e
               den_OH_e = 1 / CHI_OH + den_OH_e 

         eOH = OH_e / den_OH_e
         elogU = logU_e / den_OH_e 


      OH_mc.append(OH)
      CO_mc.append(CO)
      logU_mc.append(logU)
      OHe_mc.append(eOH)
      COe_mc.append(eCO)
      logUe_mc.append(elogU)


   OHff = np.mean(OH_mc)
   eOHff = (np.std(OH_mc)**2+np.mean(OHe_mc)**2)**0.5
   COff = np.mean(CO_mc)
   eCOff = (np.std(CO_mc)**2+np.mean(COe_mc)**2)**0.5
   logUff = np.mean(logU_mc)
   elogUff = (np.std(logU_mc)**2+np.std(logUe_mc)**2)**0.5

         
   output.append(OHff)
   output.append(eOHff)
   output.append(COff)
   output.append(eCOff)
   output.append(logUff)
   output.append(elogUff)



   if input0.shape >= (12,) and count == 0: continue   
   print round(100*(count)/float(len(input)),1),'%',grid_type,'', round(OHff,3), round(eOHff,3),'',round(COff,3), round(eCOff,3), '',round(logUff,3), round(elogUff,3)


out = np.reshape(output,(len(input),19))

np.savetxt('output.dat',out,fmt='%.4f')
print '________________________________'
print 'Results are stored in output.dat'
            

























