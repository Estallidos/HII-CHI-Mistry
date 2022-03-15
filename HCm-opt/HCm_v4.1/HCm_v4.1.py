# Filename: HII-CHCm_v 4.1.py



import string
import numpy as np
import sys



print (' ---------------------------------------------------------------------')
print ('This is HII-CHI-mistry v. 4.1')
print (' See Perez-Montero, E. (2014) for details')
print  (' Insert the name of your input text file with the following columns:')
print (' 3727 [OII], 3868 [NeIII], 4363 [OIII], 5007 [OIII], 6584 [NII], 6725 [SII]')
print ('with their corresponding errors in adjacent columns')
print ('with 0 for missing information.')
print ('---------------------------------------------------------------------')





# Input file reading


if len(sys.argv) == 1:
   if int(sys.version[0]) < 3:
      input00 = raw_input('Insert input file name:')
   else:
      input00 = input('Insert input file name:')
else:
   input00 = str(sys.argv[1])
try:
   input0 = np.loadtxt(input00)
   if (input0.ndim == 1 and input0.shape[0] != 12) or (input0.ndim > 1 and input0.shape[1] != 12):
      print ('The input file does not have 12 columns. Please check')
      sys.exit()
   print ('The input file is:'+input00)
except:
   print ('Input file error: It does not exist or has wrong format')
   sys.exit()
print ('')


output = []


# Iterations for Montecarlo error derivation
if len(sys.argv) < 3:
   n = 25
else:
   n = int(sys.argv[2])
print ('The number of iterations for MonteCarlo simulation is: ',n)
print ('')



# Reading of models grids. These can be changed

print ('')
question = True
while question:
   print('-------------------------------------------------')
   print ('(1) POPSTAR with Chabrier IMF, age = 1 Myr')
   print ('(2) AGN, double component, a(OX) = -0.8, a(UV) = -1.0')
   print ('(3) AGN, double component, a(OX) = -1.2, a(UV) = -1.0')
   print('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      sed = raw_input('Choose SED of the models:')
   else:
      sed = input('Choose SED of the models:')
   if sed == '1' or sed == '2' or sed == '3': question = False 
print ('')

question = True
while question:
   if int(sys.version[0]) < 3:
      inter = raw_input('Choose models [0] No interpolated [1] Interpolated: ')
   else:
      inter = input('Choose models [0] No interpolated [1] Interpolated: ')
   if inter == '0' or inter == '1': question = False
print ('')

sed = int(sed)
inter = int(inter)

if inter == 0 and sed==1:
   sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation'
   grid1 = np.loadtxt('C17_cha_1Myr_v4.0.dat')
   grid2 = np.loadtxt('C17_cha_1Myr_logU_adapted_emp_v4.0.dat')
   grid3 = np.loadtxt('C17_cha_1Myr_logU-NO_adapted_emp_v4.0.dat')
   print ('No interpolation for the POPSTAR models is going to be used.')
   print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
   print ('')

   res_NO = 0.125
elif inter == 1 and sed==1:
   sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. interpolation'
   grid1 = np.loadtxt('C17_cha_1Myr_v4.0.dat')
   grid2 = np.loadtxt('C17_cha_1Myr_logU_adapted_emp_v4.0.dat')
   grid3 = np.loadtxt('C17_cha_1Myr_logU-NO_adapted_emp_v4.0.dat')
   grid1i = np.loadtxt('C17_cha_1Myr_v4.0int.dat')
   grid2i = np.loadtxt('C17_cha_1Myr_logU_adapted_emp_v4.0int.dat')
   grid3i = np.loadtxt('C17_cha_1Myr_logU-NO_adapted_emp_v4.0int.dat')
   print ('Interpolation for the POPSTAR models is going to be used.')
   print ('The grid has a resolution of 0.02dex for O/H and 0.025dex for N/O')
   print ('')

   res_NO = 0.125
elif inter == 0 and sed==2:
   sed_type = 'Double composite AGN, a(OX) = -0.8. No interpolation'
   grid1 = np.loadtxt('C17_agn_v4.0.dat')
   grid2 = np.loadtxt('C17_agn_v4.0.dat')
   grid3 = np.loadtxt('C17_agn_NO_adapted_emp_v4.0.dat')
   print ('No interpolation for the AGN a(ox) = -0.8 models is going to be used.')
   print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
   print ('')

   res_NO = 0.125
elif inter == 1 and sed==2:
   sed_type = 'Double composite AGN, a(OX) = -0.8. Interpolation'
   grid1 = np.loadtxt('C17_agn_v4.0.dat')
   grid2 = np.loadtxt('C17_agn_v4.0.dat')
   grid3 = np.loadtxt('C17_agn_NO_adapted_emp_v4.0.dat')
   grid1i = np.loadtxt('C17_agn_v4.0int.dat')
   grid2i = np.loadtxt('C17_agn_v4.0int.dat')
   grid3i = np.loadtxt('C17_agn_NO_adapted_emp_v4.0int.dat')
   print ('Interpolation for the AGN a(ox) = -0.8 models is going to be used.')
   print ('The grid has a resolution of 0.02dex for O/H and 0.025dex for N/O')
   print ('')

   res_NO = 0.125
elif inter == 0 and sed==3:
   sed_type = 'Double composite AGN, a(OX) = -1.2. No interpolation'
   grid1 = np.loadtxt('C17_agn_a12_v4.0.dat')
   grid2 = np.loadtxt('C17_agn_a12_v4.0.dat')
   grid3 = np.loadtxt('C17_agn_a12_NO_adapted_emp_v4.0.dat')
   print ('No interpolation for the AGN a(ox) = -1.2 models is going to be used.')
   print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
   print ('')

   res_NO = 0.125
elif inter == 1 and sed==3:
   sed_type = 'Double composite AGN, a(OX) = -1.2. Interpolation'
   grid1 = np.loadtxt('C17_agn_a12_v4.0.dat')
   grid2 = np.loadtxt('C17_agn_a12_v4.0.dat')
   grid3 = np.loadtxt('C17_agn_a12_NO_adapted_emp_v4.0.dat')
   grid1i = np.loadtxt('C17_agn_a12_v4.0int.dat')
   grid2i = np.loadtxt('C17_agn_a12_v4.0int.dat')
   grid3i = np.loadtxt('C17_agn_a12_NO_adapted_emp_v4.0int.dat')
   print ('Interpolation for the AGN a(ox) = -1.2 models is going to be used.')
   print ('The grid has a resolution of 0.02dex for O/H and 0.025dex for N/O')
   print ('')

   res_NO = 0.125


# Input file reading



if input0.shape == (12,):
   input1 = [0,0,0,0,0,0,0,0,0,0, 0, 0,input0[0],input0[1],input0[2],input0[3],input0[4],input0[5],input0[6],input0[7],input0[8],input0[9],input0[10],input0[11]]
   input = np.reshape(input1,(2,12))
else:
   input = input0


print ('Reading grids ....')
print ('')
print ('')
print ('----------------------------------------------------------------')
print ('(%)   Grid  12+log(O/H)  log(N/O)    log(U)')
print ('-----------------------------------------------------------------')


# Beginning of loop of calculation

count = 0
for tab in input:
   count = count + 1


   OH_mc = []
   NO_mc = []
   logU_mc = []
   OHe_mc = []
   NOe_mc = []
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
      NO_p = 0
      den_OH = 0
      den_NO = 0
      OH_e = 0
      NO_e = 0
      logU_e = 0
      den_OH_e = 0
      den_NO_e = 0
      tol_max = 1e2

      if tab[0] == 0:
         OII_3727_obs = 0
      else:
         OII_3727_obs = np.random.normal(tab[0],tab[1]+1e-5)
         if OII_3727_obs <= 0: OII_3727_obs = 0
      if tab[2] == 0:
         NeIII_3868_obs = 0
      else:
         NeIII_3868_obs = np.random.normal(tab[2],tab[3]+1e-5)
         if NeIII_3868_obs <= 0: NeIII_3868_obs = 0
      if tab[4] == 0:
         OIII_4363_obs = 0
      else:
         OIII_4363_obs = np.random.normal(tab[4],tab[5]+1e-5)
         if OIII_4363_obs <= 0: OIII_4363_obs = 0
      if tab[6] == 0:
         OIII_5007_obs = 0
      else:
         OIII_5007_obs = np.random.normal(tab[6],tab[7]+1e-5)
         if OIII_5007_obs <= 0: OIII_5007_obs = 0
      if OIII_4363_obs == 0  or OIII_5007_obs == 0:
         ROIII_obs = 0
      else:
         ROIII_obs = OIII_5007_obs / OIII_4363_obs
      if tab[8] == 0:
         NII_6584_obs = 0
      else:
         NII_6584_obs = np.random.normal(tab[8],tab[9]+1e-3)
         if NII_6584_obs <= 0: NII_6584_obs = 0
      if tab[10] == 0:
            SII_6725_obs = 0
      else:
         SII_6725_obs = np.random.normal(tab[10],tab[11]+1e-3)
         if SII_6725_obs <= 0: SII_6725_obs = 0
      if NII_6584_obs == 0 or OII_3727_obs == 0:
         N2O2_obs = -10
      else:
         N2O2_obs = np.log10(NII_6584_obs / OII_3727_obs)
      if NII_6584_obs == 0 or SII_6725_obs == 0:
         N2S2_obs = -10
      else:   
         N2S2_obs = np.log10(NII_6584_obs / SII_6725_obs)
      if OII_3727_obs == 0 or OIII_5007_obs== 0:
         O2O3_obs = 0
         R23_obs = -10
      else:
         R23_obs = np.log10(OII_3727_obs + OIII_5007_obs )
         O2O3_obs = (OII_3727_obs / OIII_5007_obs )
      if OII_3727_obs == 0 or NeIII_3868_obs== 0:
         O2Ne3_obs = 0
         R2Ne3_obs = -10
      else:
         O2Ne3_obs = (OII_3727_obs / NeIII_3868_obs )
         R2Ne3_obs = np.log10(OII_3727_obs + NeIII_3868_obs )
      if OIII_5007_obs == 0 or NII_6584_obs == 0:
         O3N2_obs = -10
      else:
         O3N2_obs = np.log10( OIII_5007_obs / NII_6584_obs )


# Selection of grid


   
      if ROIII_obs > 0:
         grid = grid1
         if monte == n-1: output.append(1)
         grid_type = 1
      elif N2O2_obs > -10 or N2S2_obs > -10:
         grid = grid2
         if monte == n-1: output.append(2)
         grid_type = 2
      else:
         grid = grid3
         if monte == n-1: output.append(3)
         grid_type = 3
      

# Calculation of N/O and errors

      if N2O2_obs == -10 and N2S2_obs == -10:
         NO = -10
      else:
         CHI_ROIII = 0
         CHI_N2O2 = 0
         CHI_N2S2 = 0

         for index in grid:
            if ROIII_obs == 0: 
               CHI_ROIII = 0
            elif index[5] == 0:
               CHI_ROIII = tol_max
            else:   
               CHI_ROIII = (index[6]/index[5]- ROIII_obs)**2/(index[6]/index[5])
            if N2O2_obs == -10:
               CHI_N2O2 = 0
            elif index[3] == 0 or index[7] == 0:
               CHI_N2O2 = tol_max
            else:
               CHI_N2O2 =(np.log10(index[7]/index[3]) - N2O2_obs)**2/(abs(np.log10(index[7]/index[3])+1e-5))
            if N2S2_obs == -10: 
               CHI_N2S2 = 0
            elif index[7] == 0 or index[8] == 0:
               CHI_N2S2 = tol_max
            else:
               CHI_N2S2 =(np.log10(index[7]/index[8]) - N2S2_obs)**2/(abs(np.log10(index[7]/index[8])+1e-5))

            CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2)**0.5
            NO_p = index[1] / CHI_NO**2 + NO_p
            den_NO = 1 / CHI_NO**2 + den_NO
         NO = NO_p / den_NO 


# Calculation of N/O error

   
      if N2O2_obs == -10 and N2S2_obs == -10:
         eNO = 0
      else:   
         CHI_ROIII = 0
         CHI_N2O2 = 0
         CHI_N2S2 = 0
         CHI_NO = 0
      
         for index in grid:
            if ROIII_obs == 0: 
               CHI_ROIII = 0
            elif index[5] == 0:
               CHI_ROIII = tol_max
            else:   
               CHI_ROIII = (index[6]/index[5]- ROIII_obs)**2/(index[6]/index[5])
            if N2O2_obs == -10:
               CHI_N2O2 = 0
            elif index[3] == 0 or index[7] == 0:
               CHI_N2O2 = tol_max
            else:
               CHI_N2O2 =(np.log10(index[7]/index[3]) - N2O2_obs)**2/(abs(np.log10(index[7]/index[3])+1e-5))
            if N2S2_obs == -10: 
               CHI_N2S2 = 0
            elif index[7] == 0 or index[8] == 0:
               CHI_N2S2 = tol_max
            else:
               CHI_N2S2 =(np.log10(index[7]/index[8]) - N2S2_obs)**2/(abs(np.log10(index[7]/index[8])+1e-5))

            CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2)**0.5
            NO_e = (index[1] - NO)**2 / CHI_NO**2 + NO_e
            den_NO_e = 1 / CHI_NO**2 + den_NO_e  
         eNO = NO_e / den_NO_e 
      


# Calculation of O/H and logU

      if R23_obs == -10 and NII_6584_obs == 0 and ROIII_obs == 0 and R2Ne3_obs == -10:
         OH = 0
         logU = 0
      else:
         CHI_ROIII = 0
         CHI_NII = 0
         CHI_OIII = 0
         CHI_OII = 0
         CHI_O2O3 = 0
         CHI_R23 = 0
         CHI_O2Ne3 = 0
         CHI_R2Ne3 = 0
         CHI_O3N2 = 0
         for index in grid:

            if NO > -10 and np.abs(index[1] - NO) > np.abs(eNO+res_NO):
               continue
            else:
               if ROIII_obs == 0: 
                  CHI_ROIII = 0
               elif index[5] == 0:
                  CHI_ROIII = tol_max
               else:   
                  CHI_ROIII = (index[6]/index[5]- ROIII_obs)**2/(index[6]/index[5])
               if OIII_5007_obs == 0:
                  CHI_OIII = 0
               elif index[6] == 0:
                  CHI_OIII = tol_max
               else:
                  CHI_OIII = (index[6] - OIII_5007_obs)**2/index[6]
               if OII_3727_obs == 0:
                  CHI_OII = 0
               elif index[3] == 0:
                  CHI_OII = tol_max
               else:
                  CHI_OII = (index[3] - OII_3727_obs)**2/index[3]
               if NII_6584_obs == 0:
                  CHI_NII = 0
               elif index[7] == 0:
                  CHI_NII = tol_max
               else:
                  CHI_NII = (index[7] - NII_6584_obs)**2/index[7]
               if OII_3727_obs == 0 or OIII_5007_obs == 0:
                  CHI_O2O3 = 0
                  CHI_R23 = 0
               elif index[3] == 0 or index[6] == 0:
                  CHI_O2O3 = tol_max
                  CHI_R23 = tol_max
               else:
                  CHI_O2O3 = (index[3]/index[6] - O2O3_obs)**2/(index[3]/index[6])
                  CHI_R23 = (np.log10(index[3]+index[6])-R23_obs)**2/   (np.abs(np.log10(index[3]+index[6]+1e-5)))
               if OII_3727_obs == 0 or NeIII_3868_obs == 0:
                  CHI_O2Ne3 = 0
                  CHI_R2Ne3 = 0
               elif index[3] == 0 or index[4] == 0:
                  CHI_O2Ne3 = tol_max
                  CHI_R2Ne3 = tol_max
               else:
                  CHI_O2Ne3 = (index[3]/index[4] - O2Ne3_obs)**2/(index[3]/index[4])
                  CHI_R2Ne3 = (np.log10(index[3]+index[4])-R2Ne3_obs)**2/   (np.abs(np.log10(index[3]+index[4]+1e-5)))
               if OIII_5007_obs == 0 or NII_6584_obs == 0:
                  CHI_O3N2 = 0
               elif index[6] == 0 or index[7] == 0:
                  CHI_O3N2 = tol_max
               else:
                  CHI_O3N2 = (np.log10(index[6]/index[7]) - O3N2_obs)**2/(np.abs(np.log10(index[6]/index[7]+1e-5)))


               if ROIII_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OIII**2 + CHI_OII**2 )**0.5
               elif NII_6584_obs > 0 and OII_3727_obs > 0:
                  CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O2Ne3**2 + CHI_R2Ne3**2)**0.5 
               elif NII_6584_obs > 0 and OII_3727_obs == 0:
                  CHI_OH = (CHI_NII**2 + CHI_O3N2**2)**0.5
               elif NII_6584_obs == 0:
                  CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O2Ne3**2 + CHI_R2Ne3**2)**0.5

            OH_p = index[0] / CHI_OH**2 + OH_p
            logU_p = index[2] / CHI_OH**2 + logU_p
            den_OH = 1 / CHI_OH**2 + den_OH
         
         OH = OH_p / den_OH
         logU = logU_p / den_OH


# Calculation of error of O/H and logU

      if R23_obs == -10 and NII_6584_obs == 0 and ROIII_obs == 0 and R2Ne3_obs == -10:
         eOH = 0
         elogU = 0
      else:
         CHI_ROIII = 0
         CHI_NII = 0
         CHI_OIII = 0
         CHI_OII = 0
         CHI_O2O3 = 0
         CHI_R23 = 0
         CHI_O2Ne3 = 0
         CHI_R2Ne3 = 0
         CHI_O3N2 = 0
         for index in grid:

            if NO > -10 and np.abs(index[1] - NO) > np.abs(eNO+res_NO):
               continue
            else:
               if ROIII_obs == 0: 
                  CHI_ROIII = 0
               elif index[5] == 0:
                  CHI_ROIII = tol_max
               else:   
                  CHI_ROIII = (index[6]/index[5]- ROIII_obs)**2/(index[6]/index[5])
               if OIII_5007_obs == 0:
                  CHI_OIII = 0
               elif index[6] == 0:
                  CHI_OIII = tol_max
               else:
                  CHI_OIII = (index[6] - OIII_5007_obs)**2/index[6]
               if OII_3727_obs == 0:
                  CHI_OII = 0
               elif index[3] == 0:
                  CHI_OII = tol_max
               else:
                  CHI_OII = (index[3] - OII_3727_obs)**2/index[3]
               if NII_6584_obs == 0:
                  CHI_NII = 0
               elif index[7] == 0:
                  CHI_NII = tol_max
               else:
                  CHI_NII = (index[7] - NII_6584_obs)**2/index[7]
               if OII_3727_obs == 0 or OIII_5007_obs == 0:
                  CHI_O2O3 = 0
                  CHI_R23 = 0
               elif index[3] == 0 or index[6] == 0:
                  CHI_O2O3 = tol_max
                  CHI_R23 = tol_max
               else:
                  CHI_O2O3 = (index[3]/index[6] - O2O3_obs)**2/(index[3]/index[6])
                  CHI_R23 = (np.log10(index[3]+index[6])-R23_obs)**2/   (np.abs(np.log10(index[3]+index[6]+1e-5)))
               if OII_3727_obs == 0 or NeIII_3868_obs == 0:
                  CHI_O2Ne3 = 0
                  CHI_R2Ne3 = 0
               elif index[3] == 0 or index[4] == 0:
                  CHI_O2Ne3 = tol_max
                  CHI_R2Ne3 = tol_max
               else:
                  CHI_O2Ne3 = (index[3]/index[4] - O2Ne3_obs)**2/(index[3]/index[4])
                  CHI_R2Ne3 = (np.log10(index[3]+index[4])-R2Ne3_obs)**2/   (np.abs(np.log10(index[3]+index[4]+1e-5)))
               if OIII_5007_obs == 0 or NII_6584_obs == 0:
                  CHI_O3N2 = 0
               elif index[6] == 0 or index[7] == 0:
                  CHI_O3N2 = tol_max
               else:
                  CHI_O3N2 = (np.log10(index[6]/index[7]) - O3N2_obs)**2/(np.abs(np.log10(index[6]/index[7]+1e-5)))


               if ROIII_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OIII*2 + CHI_OII**2)**0.5
               elif NII_6584_obs > 0 and OII_3727_obs > 0:
                  CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O2Ne3**2 + CHI_R2Ne3**2)**0.5 
               elif NII_6584_obs > 0 and OII_3727_obs == 0:
                  CHI_OH = (CHI_NII**2 + CHI_O3N2**2)**0.5
               elif NII_6584_obs == 0:
                  CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O2Ne3**2 + CHI_R2Ne3**2)**0.5


            if CHI_OH == 0:
               OH_e = OH_e
               logU_e = logU_e
               den_OH_e = den_OH_e
            else:
               OH_e = (index[0] - OH)**2 / CHI_OH**2 + OH_e
               logU_e = (index[2] - logU)**2 / CHI_OH**2 + logU_e
               den_OH_e = 1 / CHI_OH**2 + den_OH_e 

         eOH = OH_e / den_OH_e
         elogU = logU_e / den_OH_e 


# Iterations for interpolated models

      if inter == 0:
         NOf = NO
         OHf = OH
         logUf = logU

      
      elif inter == 1:
         if grid_type == 1:
            igrid = grid1i
         elif grid_type == 2:
            igrid = grid2i
         elif grid_type == 3:
            igrid = grid3i
         iN = igrid[:,1]
         cond_N = (iN > NO-0.125) & (iN < NO + 0.125)
         if NO > -10:
            igrid_N = igrid[cond_N]
         else:
            igrid_N = igrid    
         iO = igrid_N[:,0]
         cond_O = (iO < OH + 0.1) & (iO > OH - 0.1)
         igrid_O = igrid_N[cond_O]
      
         CHI_ROIII = 0
         CHI_N2O2 = 0
         CHI_N2S2 = 0
         CHI_OIII = 0
         CHI_OII = 0
         CHI_NII = 0
         CHI_O2O3 = 0
         CHI_R23 = 0
         CHI_O3N2 = 0
         CHI_O2Ne3 = 0
         CHI_R2Ne3 = 0
         CHI_NO = 0
         CHI_OH = 0
         NO_p = 0
         den_NO = 0
         OH_p = 0
         logU_p = 0
         den_OH = 0
      

         for index in igrid_O:
            if ROIII_obs == 0: 
               CHI_ROIII = 0
            elif index[5] == 0:
               CHI_ROIII = tol_max
            else:   
               CHI_ROIII = (index[6]/index[5]- ROIII_obs)**2/(index[6]/index[5])
               if OIII_5007_obs == 0:
                  CHI_OIII = 0
               elif index[6] == 0:
                  CHI_OIII = tol_max
               else:
                  CHI_OIII = (index[6] - OIII_5007_obs)**2/index[6]
               if OII_3727_obs == 0:
                  CHI_OII = 0
               elif index[3] == 0:
                  CHI_OII = tol_max
               else:
                  CHI_OII = (index[3] - OII_3727_obs)**2/index[3]
            if N2O2_obs == -10:
               CHI_N2O2 = 0
            elif index[3] == 0 or index[7] == 0:
               CHI_N2O2 = tol_max
            else:
               CHI_N2O2 =(np.log10(index[7]/index[3]) - N2O2_obs)**2/(abs(np.log10(index[7]/index[3])+1e-5))
            if N2S2_obs == -10: 
               CHI_N2S2 = 0
            elif index[7] == 0 or index[8] == 0:
               CHI_N2S2 = tol_max
            else:
               CHI_N2S2 =(np.log10(index[7]/index[8]) - N2S2_obs)**2/(abs(np.log10(index[7]/index[8])+1e-5))
            if NII_6584_obs == 0:
               CHI_NII = 0
            elif index[7] == 0:
               CHI_NII = tol_max
            else:
               CHI_NII = (index[7] - NII_6584_obs)**2/index[7]
            if OII_3727_obs == 0 or OIII_5007_obs == 0:
               CHI_O2O3 = 0
               CHI_R23 = 0
            elif index[3] == 0 or index[6] == 0:
               CHI_O2O3 = tol_max
               CHI_R23 = tol_max
            else:
               CHI_O2O3 = (index[3]/index[6] - O2O3_obs)**2/(index[3]/index[6])
               CHI_R23 = (np.log10(index[3]+index[6])-R23_obs)**2/(np.abs(np.log10(index[3]+index[6]+1e-5)))
            if OII_3727_obs == 0 or NeIII_3868_obs == 0:
               CHI_O2Ne3 = 0
               CHI_R2Ne3 = 0
            elif index[3] == 0 or index[4] == 0:
               CHI_O2Ne3 = tol_max
               CHI_R2Ne3 = tol_max
            else:
               CHI_O2Ne3 = (index[3]/index[4] - O2Ne3_obs)**2/(index[3]/index[4])
               CHI_R2Ne3 = (np.log10(index[3]+index[4])-R2Ne3_obs)**2/(np.abs(np.log10(index[3]+index[4]+1e-5)))

            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif index[6] == 0 or index[7] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index[6]/index[7]) - O3N2_obs)**2/(np.abs(np.log10(index[6]/index[7]+1e-5)))



            CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2)**0.5
            if CHI_NO == 0:
               NO_p = NO_p
               den_NO = den_NO
            else:
               NO_p = index[1] / CHI_NO**2 + NO_p
               den_NO = 1 / CHI_NO**2 + den_NO

            if ROIII_obs > 0:
               CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2)**0.5
            elif NII_6584_obs > 0 and OII_3727_obs > 0:
               CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O2Ne3**2 + CHI_R2Ne3**2)**0.5 
            elif NII_6584_obs > 0 and OII_3727_obs == 0:
               CHI_OH = (CHI_NII**2 + CHI_O3N2**2)**0.5
            elif NII_6584_obs == 0:
               CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O2Ne3**2 + CHI_R2Ne3**2 )**0.5

            OH_p = index[0] / CHI_OH**2 + OH_p
            logU_p = index[2] / CHI_OH**2 + logU_p
            den_OH = 1 / CHI_OH**2 + den_OH

         if NO == -10:
            NOf = NO
         else:
            NOf = NO_p / den_NO

         if OH == 0:
            OHf = OH
            logUf = logU
         else:         
            OHf = OH_p / den_OH
            logUf = logU_p / den_OH


      OH_mc.append(OHf)
      NO_mc.append(NOf)
      logU_mc.append(logUf)
      OHe_mc.append(eOH)
      NOe_mc.append(eNO)
      logUe_mc.append(elogU)
   
   OHff = np.mean(OH_mc)
   eOHff = (np.std(OH_mc)**2+np.mean(OHe_mc)**2)**0.5
   NOff = np.mean(NO_mc)
   eNOff = (np.std(NO_mc)**2+np.mean(NOe_mc)**2)**0.5
   logUff = np.mean(logU_mc)
   elogUff = (np.std(logU_mc)**2+np.std(logUe_mc)**2)**0.5

   




   
   output.append(OHff)
   output.append(eOHff)
   output.append(NOff)
   output.append(eNOff)
   output.append(logUff)
   output.append(elogUff)
         

   if input0.shape >= (12,) and count == 0: continue   
   print (round(100*(count)/float(len(input)),1),'%',grid_type,'', round(OHff,3), round(eOHff,3),'',round(NOff,3), round(eNOff,3), '',round(logUff,3), round(elogUff,3))


out = np.reshape(output,(len(input),19))



lineas_header = [' HII-CHI-mistry v.4.1 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','O2Hb eO2Hb  Ne3Hb  eNeHb O3aHb  eO3aHb O3nHb  eO3nHb N2Hb   eN2Hb  S2Hb   eS2Hb  i O/H     eO/H  N/O    eN/O  logU   elogU']

header = '\n'.join(lineas_header)

np.savetxt(input00+'_hcm-output.dat',out,fmt=' '.join(['%.4f']*12+['%i']+['%.3f']*6),header=header)
print ('________________________________')
print ('Results are stored in ' + input00 + '_hcm-output.dat')
                           
      

      
