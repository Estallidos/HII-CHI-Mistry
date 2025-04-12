# Filename: HII-CHI-mistry_v02.0_UV.py

print ' ---------------------------------------------------------------------'
print ' This is HII-CHI-mistry for UV version 2.0'
print ' See Perez-Montero & Amorin (2016) for details'
print  ' Insert the name of your input text file with following columns:'
print ' Lya 1216, CIV 1549, OIII 1665, CIII 1909, Hb 4861, OIII 5007'
print '---------------------------------------------------------------------'
print ''


import string
import numpy as np

# Reading of models grids. These can be changed

grid1 = np.loadtxt('C13_cha_1Myr_v2.0_uv.dat')
grid2 = np.loadtxt('C13_cha_1Myr_logU_adapted_emp_v2.0_uv.dat')
grid3 = np.loadtxt('C13_cha_1Myr_logU-CO_adapted_emp_v2.0_uv.dat')
res_CO = 0.125


# Input file reading

input0 = np.loadtxt(raw_input('Insert input file name:'))
output = []


if input0.shape == (6,):
   input1 = [input0[0],input0[1],input0[2],input0[3],input0[4],input0[5],0,0,0,0,0,0]
   input = np.reshape(input1,(2,6))
else:
   input = input0


print 'Reading grids ....'
print ''
print ''
print '----------------------------------------------------------------'
print '(%)   Grid  12+log(O/H)  log(C/O)    log(U)'
print '-----------------------------------------------------------------'

# Beginning of loop of calculation


for tab in range(len(input)):


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

   Lya_1216_obs = input[tab,0]
   output.append(Lya_1216_obs)
   CIV_1549_obs = input[tab,1]
   output.append(CIV_1549_obs)
   OIII_1665_obs = input[tab,2]
   output.append(OIII_1665_obs)
   CIII_1909_obs = input[tab,3]
   output.append(CIII_1909_obs)
   Hb_4861_obs = input[tab,4]
   output.append(Hb_4861_obs)
   OIII_5007_obs = input[tab,5]
   output.append(OIII_5007_obs)
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
      output.append(1)
      grid_type = 1
   elif C3O3_obs > 0:
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
      if C3O3_obs == 0:
         CHI_C3O3 = 0
      elif grid[index,4] == 0 or grid[index,6] == 0:
         CHI_C3O3 = tol_max
      else:
         CHI_C3O3 = ((grid[index,6]+grid[index,7])/grid[index,4] - C3O3_obs)**2/((grid[index,6]+grid[index,7])/grid[index,4])
      chi.append(CHI_C3O3)
      if C3C4_obs == 0:
         CHI_C3C4 = 0
      elif grid[index,6] == 0 or grid[index,7] == 0:
         CHI_C3C4 = tol_max
      else:
         CHI_C3C4 = (grid[index,6]/grid[index,7] - C3C4_obs)**2/(grid[index,6]/grid[index,7])
      chi.append(CHI_C3C4)
      if C34_obs == 0:
         CHI_C34 = 0
      elif grid[index,6] == 0 or grid[index,3] == 0:
         CHI_C34 = tol_max
      else:
         CHI_C34 = ((grid[index,6]+grid[index,7])/grid[index,3] - C34_obs)**2/((grid[index,6]+grid[index,7])/grid[index,3])
      chi.append(CHI_C34)
      if O3_obs == 0:
         CHI_O3 = 0
      elif grid[index,4] == 0 or grid[index,3] == 0:
         CHI_O3 = tol_max
      else:
         CHI_O3 = (grid[index,4]/grid[index,3] - O3_obs)**2/(grid[index,4]/grid[index,3])
      chi.append(CHI_O3)
      if C34Hb_obs == 0:
         CHI_C34Hb = 0
      elif grid[index,6] == 0:
         CHI_C34Hb = tol_max
      else:
         CHI_C34Hb = (grid[index,6]+grid[index,7] - C34Hb_obs)**2 / (grid[index,6] + grid[index,7])
      chi.append(CHI_C34Hb)
      if O3Hb_obs == 0:
         CHI_O3Hb = 0
      elif grid[index,4] == 0:
         CHI_O3Hb = tol_max
      else:
         CHI_O3Hb = (grid[index,4] - O3Hb_obs)**2/grid[index,4]
      chi.append(CHI_O3Hb)


   chi_mat = np.reshape(chi,(len(grid),10))
   
# Calculation of C/O and errors

   if C3O3_obs == 0:
      CO = -10
      eCO = 0
   else:
      for index in range(len(chi_mat)):
         if chi_mat[index,4] == 0:
            continue
         else:
            CHI_CO = (chi_mat[index,3]**2 + chi_mat[index,4]**2)**0.5
            CO_p = chi_mat[index,1] / CHI_CO + CO_p
            den_CO = 1 / CHI_CO + den_CO
         CO = CO_p / den_CO 
      for index in range(len(chi_mat)):
         if chi_mat[index,4] == 0:
            continue
         else:
            CHI_COe = (chi_mat[index,3]**2 + chi_mat[index,4]**2)**0.5
            CO_e = (chi_mat[index,1] - CO)**2 / CHI_COe + CO_e
            den_CO_e = 1 / CHI_COe + den_CO_e  
         eCO = CO_e / den_CO_e 

# Calculation of O/H and error

   if C34_obs == 0 and O3_obs == 0 and ROIII_obs == 0 and C34Hb_obs == 0 and O3Hb_obs == 0:
      OH = 0
      eOH = 0
      logU = 0
      elogU = 0
   else:
      for index in range(len(chi_mat)):
         if CO > -10 and np.abs(chi_mat[index,1] - CO) > np.abs(eCO+res_CO):
            continue
         else:
            CHI_OH = (chi_mat[index,3]**2 + chi_mat[index,5]**2 + chi_mat[index,6]**2 + chi_mat[index,7]**2
               + chi_mat[index,8]**2 + chi_mat[index,9])**0.5

         OH_p = chi_mat[index,0] / CHI_OH + OH_p
         logU_p = chi_mat[index,2] / CHI_OH + logU_p
         den_OH = 1 / CHI_OH + den_OH
         
      OH = OH_p / den_OH
      logU = logU_p / den_OH

      for index in range(len(chi_mat)):
         if CO > -10 and np.abs(chi_mat[index,1] - CO) > np.abs(eCO+res_CO):
            continue
         else:
            CHI_OHe = (chi_mat[index,3]**2 + chi_mat[index,5]**2 + chi_mat[index,6]**2 + chi_mat[index,7]**2 + 
               chi_mat[index,8]**2 + chi_mat[index,9]**2)**0.5
         OH_e = (chi_mat[index,0] - OH)**2 / CHI_OHe + OH_e
         logU_e = (chi_mat[index,2] - logU)**2 / CHI_OHe + logU_e
         den_OH_e = 1 / CHI_OHe + den_OH_e 
      eOH = OH_e / den_OH_e
      elogU = logU_e / den_OH_e 

   output.append(OH)
   output.append(eOH)
   output.append(CO)
   output.append(eCO)
   output.append(logU)
   output.append(elogU)

   if input0.shape >= (5,) and tab == 1: continue
   print round(100*(tab+1)/float(len(input)),1),'%',grid_type,'', round(OH,3), round(eOH,3),'',round(CO,3), round(eCO,3), '',round(logU,3), round(elogU,3)

out = np.reshape(output,(len(input),13))

np.savetxt('output.dat',out,fmt='%.4f')
print '________________________________'
print 'Results are stored in output.dat'
            

























