# Filename: HII-CHI-mistry_v02.1.py

print ' ---------------------------------------------------------------------'
print ' This is HII-CHI-mistry version 2.1'
print ' See Perez-Montero, E. (2014) for details'
print  ' Insert the name of your input text file with following columns:'
print ' 3727 [OII], 4363 [OIII], 5007 [OIII], 6584 [NII], 6725 [SII]'
print '---------------------------------------------------------------------'


import string
import numpy as np

# Reading of models grids. These can be changed

print ''
question = True
while question:
   inter = int(raw_input('Choose models [0] No interpolated [1] Interpolated: '))
   if inter == 0:
      grid1 = np.loadtxt('C13_cha_1Myr_v2.0.dat')
      grid2 = np.loadtxt('C13_cha_1Myr_logU_adapted_emp_v2.0.dat')
      grid3 = np.loadtxt('C13_cha_1Myr_logU-NO_adapted_emp_v2.0.dat')
      print 'No interpolation for the models are going to be used.'
      print 'The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O'
      print ''
      question = False
      res_NO = 0.125
   elif inter == 1:
      grid1 = np.loadtxt('C13_cha_1Myr_v2.0int.dat')
      grid2 = np.loadtxt('C13_cha_1Myr_logU_adapted_emp_v2.0int.dat')
      grid3 = np.loadtxt('C13_cha_1Myr_logU-NO_adapted_emp_v2.0int.dat')
      print ''
      print 'The grid has a resolution 0f 0.02dex for O/H and 0.025dex for N/O'
      print ''
      question = False
      res_NO = 0.125


# Input file reading

input0 = np.loadtxt(raw_input('Insert input file name:'))
output = []

if input0.shape == (5,):
   input1 = [input0[0],input0[1],input0[2],input0[3],input0[4],0,0,0,0,0]
   input = np.reshape(input1,(2,5))
else:
   input = input0


print 'Reading grids ....'
print ''
print ''
print '----------------------------------------------------------------'
print '(%)   Grid  12+log(O/H)  log(N/O)    log(U)'
print '-----------------------------------------------------------------'


# Beginning of loop of calculation

for tab in range(len(input)):


   OH_p = 0
   logU_p = 0
   NO_p = 0
   den_OH = 0
   den_NO = 0
   OH_e = 0
   NO_e = 0
   logU_e = 0
   den_OH_e = 0
   den_NO_e = 0
   tol_max = 1e2

   OII_3727_obs = input[tab,0]
   output.append(OII_3727_obs)
   OIII_4363_obs = input[tab,1]
   output.append(OIII_4363_obs)
   OIII_5007_obs = input[tab,2]
   output.append(OIII_5007_obs)
   if OIII_4363_obs == 0  or OIII_5007_obs == 0:
      ROIII_obs = 0
   else:
      ROIII_obs = OIII_5007_obs / OIII_4363_obs
   NII_6584_obs = input[tab,3]
   output.append(NII_6584_obs)
   SII_6725_obs = input[tab,4]
   output.append(SII_6725_obs)
   if NII_6584_obs == 0 or OII_3727_obs == 0:
      N2O2_obs = -10
   else:
      N2O2_obs = np.log10(NII_6584_obs / OII_3727_obs)
   if NII_6584_obs == 0 or SII_6725_obs == 0:
      N2S2_obs = -10
   else:   
      N2S2_obs = np.log10(NII_6584_obs / SII_6725_obs)
   if OII_3727_obs == 0 or OIII_5007_obs == 0:
      O2O3_obs = 0
      R23_obs = -10
   else:
      O2O3_obs = (OII_3727_obs / OIII_5007_obs )
      R23_obs = np.log10(OII_3727_obs + OIII_5007_obs )
   if OIII_5007_obs == 0 or NII_6584_obs == 0:
      O3N2_obs = -10
   else:
      O3N2_obs = np.log10( OIII_5007_obs / NII_6584_obs )


# Selection of grid

   if ROIII_obs > 0:
      grid = grid1
      output.append(1)
      grid_type = 1
   elif N2O2_obs > -10 or N2S2_obs > -10:
      grid = grid2
      output.append(2)
      grid_type = 2
   else:
      grid = grid3
      output.append(3)
      grid_type = 3
   

# Calculation of array of chis

   chi = []
   for index in range(len(grid)):
      chi.append(grid[index,0])
      chi.append(grid[index,1])
      chi.append(grid[index,2])
      if ROIII_obs == 0: 
         CHI_ROIII = 0
      elif grid[index,4] == 0:
         CHI_ROIII = tol_max
      else:
         CHI_ROIII = (grid[index,5]/grid[index,4]- ROIII_obs)**2/(grid[index,5]/grid[index,4])
      chi.append(CHI_ROIII)
      if N2O2_obs == -10:
         CHI_N2O2 = 0
      elif grid[index,3] == 0 or grid[index,6] == 0:
         CHI_N2O2 = tol_max
      else:
         CHI_N2O2 =(np.log10(grid[index,6]/grid[index,3]) - N2O2_obs)**2/(abs(np.log10(grid[index,6]/grid[index,3])+1e-5))
      chi.append(CHI_N2O2)
      if N2S2_obs == -10: 
         CHI_N2S2 = 0
      elif grid[index,6] == 0 or grid[index,7] == 0:
         CHI_N2S2 = tol_max
      else:
         CHI_N2S2 =(np.log10(grid[index,6]/grid[index,7]) - N2S2_obs)**2/(abs(np.log10(grid[index,6]/grid[index,7])+1e-5))
      
      chi.append(CHI_N2S2)
      if NII_6584_obs == 0:
         CHI_NII = 0
      elif grid[index,6] == 0:
         CHI_NII = tol_max
      else:
         CHI_NII = (grid[index,6] - NII_6584_obs)**2/grid[index,6]
      chi.append(CHI_NII)
      if OII_3727_obs == 0:
         CHI_OII = 0
      elif grid[index,3] == 0:
         CHI_OII = tol_max
      else:
         CHI_OII = (grid[index,3] - OII_3727_obs)**2/grid[index,3]
      chi.append(CHI_OII)
      if OIII_5007_obs == 0:
         CHI_OIII = 0
      elif grid[index,5] == 0:
         CHI_OIII = tol_max
      else:
         CHI_OIII = (grid[index,5] - OIII_5007_obs)**2/grid[index,5]
      chi.append(CHI_OIII)
      if OII_3727_obs == 0 or OIII_5007_obs == 0:
         CHI_O2O3 = 0
         CHI_R23 = 0
      elif grid[index,3] == 0 or grid[index,5] == 0:
         CHI_O2O3 = tol_max
         CHI_R23 = tol_max
      else:
         CHI_O2O3 = (grid[index,3]/grid[index,5] - O2O3_obs)**2/(grid[index,3]/grid[index,5])
         CHI_R23 = (np.log10(grid[index,3]+grid[index,5])-R23_obs)**2/(np.abs(np.log10(grid[index,3]+grid[index,5]+1e-5)))
      chi.append(CHI_O2O3)
      chi.append(CHI_R23)
      if OIII_5007_obs == 0 or NII_6584_obs == 0:
         CHI_O3N2 = 0
      elif grid[index,5] == 0 or grid[index,6] == 0:
         CHI_O3N2 = tol_max
      else:
         CHI_O3N2 = (np.log10(grid[index,5]/grid[index,6]) - O3N2_obs)**2/(np.abs(np.log10(grid[index,5]/grid[index,6]+1e-5)))
      chi.append(CHI_O3N2)

   chi_mat = np.reshape(chi,(len(grid),12))
   
# Calculation of N/O and errors
   
   if N2O2_obs == -10 and N2S2_obs == -10:
      NO = -10
      eNO = 0
   else:
      for index in range(len(chi_mat)):
         if chi_mat[index,4] == 0 and chi_mat[index,5] == 0:
            continue
         else:
            CHI_NO = (chi_mat[index,3]**2 + chi_mat[index,4]**2 + chi_mat[index,5]**2)**0.5
            NO_p = chi_mat[index,1] / CHI_NO + NO_p
            den_NO = 1 / CHI_NO + den_NO
         NO = NO_p / den_NO 
      for index in range(len(chi_mat)):
         if chi_mat[index,4] == 0 and chi_mat[index,5] == 0:
            continue
         else:
            CHI_NOe = (chi_mat[index,3]**2 + chi_mat[index,4]**2 + chi_mat[index,5]**2)**0.5
            NO_e = (chi_mat[index,1] - NO)**2 / CHI_NOe + NO_e
            den_NO_e = 1 / CHI_NOe + den_NO_e  
         eNO = NO_e / den_NO_e 
               

# Calculation of O/H and error

   if R23_obs == -10 and NII_6584_obs == 0:
      OH = 0
      eOH = 0
      logU = 0
      elogU = 0
   else:
      for index in range(len(chi_mat)):
         if NO > -10 and np.abs(chi_mat[index,1] - NO) > np.abs(eNO+res_NO):
            continue
         else:
            if ROIII_obs > 0:
               CHI_OH = (chi_mat[index,3]**2 + chi_mat[index,6]**2 + chi_mat[index,7]**2 + chi_mat[index,8]**2)**0.5
            elif NII_6584_obs > 0 and OII_3727_obs > 0:
               CHI_OH = (chi_mat[index,6]**2 + chi_mat[index,9]**2 + chi_mat[index,10]**2)**0.5 
            elif NII_6584_obs > 0 and OII_3727_obs == 0:
               CHI_OH = (chi_mat[index,6]**2 + chi_mat[index,9]**2 + chi_mat[index,11]**2)**0.5
            elif NII_6584_obs == 0:
               CHI_OH = (chi_mat[index,9]**2 + chi_mat[index,10]**2)**0.5

         OH_p = chi_mat[index,0] / CHI_OH + OH_p
         logU_p = chi_mat[index,2] / CHI_OH + logU_p
         den_OH = 1 / CHI_OH + den_OH
         
      OH = OH_p / den_OH
      logU = logU_p / den_OH

      for index in range(len(chi_mat)):
         if NO > -10 and np.abs(chi_mat[index,1] - NO) > np.abs(eNO+res_NO):
            continue
         else:
            if ROIII_obs > 0:
               CHI_OHe = (chi_mat[index,3]**2 + chi_mat[index,6]**2 + chi_mat[index,7]**2 + chi_mat[index,8]**2)**0.5
            elif NII_6584_obs > 0 and OII_3727_obs > 0:
               CHI_OHe = (chi_mat[index,6]**2 + chi_mat[index,9]**2 + chi_mat[index,10]**2)**0.5 
            elif NII_6584_obs > 0 and OII_3727_obs == 0:
               CHI_OHe = (chi_mat[index,6]**2 + chi_mat[index,9]**2 + chi_mat[index,11]**2)**0.5
            elif NII_6584_obs == 0:
               CHI_OHe = (chi_mat[index,9]**2 + chi_mat[index,10]**2)**0.5
         OH_e = (chi_mat[index,0] - OH)**2 / CHI_OHe + OH_e
         logU_e = (chi_mat[index,2] - logU)**2 / CHI_OHe + logU_e
         den_OH_e = 1 / CHI_OHe + den_OH_e 
      eOH = OH_e / den_OH_e
      elogU = logU_e / den_OH_e 

   output.append(OH)
   output.append(eOH)
   output.append(NO)
   output.append(eNO)
   output.append(logU)
   output.append(elogU)
         
   if input0.shape >= (5,) and tab == 1: continue
   print round(100*(tab+1)/float(len(input)),1),'%',grid_type,'', round(OH,3), round(eOH,3),'',round(NO,3), round(eNO,3), '',round(logU,3), round(elogU,3)

out = np.reshape(output,(len(input),12))

np.savetxt('output.dat',out,fmt='%.4f')
print '________________________________'
print 'Results are stored in output.dat'
                           
      

      
