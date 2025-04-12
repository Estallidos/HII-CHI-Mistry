# Filename: HII-CHCm-IR_v1.1.py



import string
import numpy as np
import sys
#sys.stderr = open('errorlog.txt', 'w')



#Function for interpolation of grids

def interpolate(grid,z,zmin,zmax,n):
   ncol = 15
   vec = []
   for col in range(ncol):
      inter = 0
      no_inter = 0
      for row in range(0,len(grid)):
         if grid[row,z] < zmin or grid[row,z] > zmax: continue
         if z == 2: x = 0; y = 1
         if z == 1: x = 0; y = 2
         if z == 0: x = 1; y = 2
         if row == (len(grid)-1):
            vec.append(grid[row,col])
            no_inter = no_inter + 1
         elif grid[row,x] < grid[row+1,x] or grid[row,y] < grid[row+1,y] :
            vec.append(grid[row,col])
            no_inter = no_inter + 1
         else:
            inter = inter + 1
            for index in range(0,n):
               i = grid[row,col]+(index)*(grid[row+1,col]-grid[row,col])/n
               vec.append(i)
   out = np.transpose(np.reshape(vec,(-1,n*inter+no_inter)))
   return out



print (' ---------------------------------------------------------------------')
print ('This is HII-CHI-mistry IR v. 1.1')
print (' See Fernandez-Ontiveros et al (2020) for details')
print  (' Insert the name of your input text file with the following columns:')
print (' HI 4.05m, HI 7.46m, [SIV] 10.5m, HI 12.4m, [NeII] 12.8m, [NeIII] 15.5m, [SIII] 18.7m, [SIII] 33.7m, [OIII] 52m, [NIII] 57m, [OIII] 88m and [NII] 122m')
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
   if (input0.ndim == 1 and input0.shape[0] != 24) or (input0.ndim > 1 and input0.shape[1] != 24):
      print ('The input file does not have 24 columns. Please check')
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
while question:
   if int(sys.version[0]) < 3:
      inter = raw_input('Choose models [0] No interpolated [1] Interpolated: ')
   else:
      inter = input('Choose models [0] No interpolated [1] Interpolated: ')
   if inter == '0' or inter == '1': question = False
print ('')

sed = 1
inter = int(inter)

if inter == 0 and sed==1:
   sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation'
   grid1 = np.loadtxt('C17_popstar_v1.1_ir.dat')
   grid2 = np.loadtxt('C17_popstar_logU_adapted_emp_v1.1_ir.dat')
   grid3 = np.loadtxt('C17_popstar_logU-NO_adapted_emp_v1.1_ir.dat')
   print ('No interpolation for the POPSTAR models is going to be used.')
   print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
   print ('')

   res_NO = 0.125
elif inter == 1 and sed==1:
   sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. interpolation'
   grid1 = np.loadtxt('C17_popstar_v1.1_ir.dat')
   grid2 = np.loadtxt('C17_popstar_logU_adapted_emp_v1.1_ir.dat')
   grid3 = np.loadtxt('C17_popstar_logU-NO_adapted_emp_v1.1_ir.dat')
   print ('Interpolation for the POPSTAR models is going to be used.')
   print ('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
   print ('')

   res_NO = 0.125
   res_NO = 0.125


# Input file reading



if input0.shape == (24,):
   input1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,input0[0],input0[1],input0[2],input0[3],input0[4],input0[5],input0[6],input0[7],input0[8],input0[9],input0[10],input0[11],input0[12],input0[13],input0[14],input0[15],input0[16],input0[17],input0[18],input0[19],input0[20],input0[21],input0[22],input0[23]]
   input = np.reshape(input1,(2,24))
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
   output.append(tab[12])
   output.append(tab[13])
   output.append(tab[14])
   output.append(tab[15])
   output.append(tab[16])
   output.append(tab[17])
   output.append(tab[18])
   output.append(tab[19])
   output.append(tab[20])
   output.append(tab[21])
   output.append(tab[22])
   output.append(tab[23])


# Selection of grid
   
   if tab[18] > 0:
      grid = grid2
      grid_type = 2
      output.append(2)
   else:
      grid = grid3
      grid_type = 3
      output.append(3)      

# Calculation of N/O

   if tab[18] == 0 and (tab[16] == 0 or tab[20] == 0):
      NOff = -10
      eNOff = 0
   else:
      for monte in range(0,n,1):
         NO_p = 0
         den_NO = 0
         NO_e = 0
         den_NO_e = 0
         tol_max = 1e2

         if tab[16] == 0:
            OIII_52_obs = 0
         else:
            OIII_52_obs = np.random.normal(tab[16],tab[17]+1e-5)
            if OIII_52_obs <= 0: OIII_52_obs = 0
         if tab[18] == 0:
            NIII_57_obs = 0
         else:
            NIII_57_obs = np.random.normal(tab[18],tab[19]+1e-5)
            if NIII_57_obs <= 0: NIII_57_obs = 0
         if tab[20] == 0:
            OIII_88_obs = 0
         else:
            OIII_88_obs = np.random.normal(tab[20],tab[21]+1e-5)
            if OIII_88_obs <= 0: OIII_88_obs = 0
         if OIII_52_obs == 0:
            N3O3a_obs = -10
         else:
            N3O3a_obs = np.log10(NIII_57_obs / OIII_52_obs)
         if OIII_88_obs == 0:
            N3O3b_obs = -10
         else:
            N3O3b_obs = np.log10(NIII_57_obs / OIII_88_obs)


         CHI_N3O3a = 0
         CHI_N3O3b = 0
         CHI_NO = 0

         for index in grid:
            if N3O3a_obs == -10: 
               CHI_N3O3a = 0
            elif index[11] == 0 or index[12] == 0:
               CHI_N3O3a = tol_max
            else:   
               CHI_N3O3a = (np.log10(index[12]/index[11])- N3O3a_obs)**2/np.log10(index[12]/index[11])
            if N3O3b_obs == -10: 
               CHI_N3O3b = 0
            elif index[13] == 0 or index[12] == 0:
               CHI_N3O3b = tol_max
            else:   
               CHI_N3O3b = (np.log10(index[12]/index[13])- N3O3b_obs)**2/np.log10(index[12]/index[13])

            CHI_NO = (CHI_N3O3a**2 + CHI_N3O3b**2 )**0.5

            NO_p = index[1] / (CHI_NO) + NO_p
            den_NO = 1 / (CHI_NO) + den_NO
         NO = NO_p / den_NO 


# Calculation of N/O error

   
         CHI_N3O3a = 0
         CHI_N3O3b = 0
         CHI_NO = 0

         for index in grid:
            if N3O3a_obs == -10: 
               CHI_N3O3a = 0
            elif index[11] == 0 or index[12] == 0:
               CHI_N3O3a = tol_max
            else:   
               CHI_N3O3a = (np.log10(index[12]/index[11])- N3O3a_obs)**2/np.log10(index[12]/index[11])
            if N3O3b_obs == -10: 
               CHI_N3O3b = 0
            elif index[13] == 0 or index[12] == 0:
               CHI_N3O3b = tol_max
            else:   
               CHI_N3O3b = (np.log10(index[12]/index[13])- N3O3b_obs)**2/np.log10(index[12]/index[13])

            CHI_NO = (CHI_N3O3a**2 + CHI_N3O3b**2 )**0.5

            NO_e = (index[1] - NO)**2 / (CHI_NO) + NO_e
            den_NO_e = 1 / (CHI_NO) + den_NO_e  
         eNO = NO_e / den_NO_e 



#Iterations for the interpolation mode

         if inter == 0 or NO == -10:
            NOf = NO
         elif inter == 1:
            igrid = grid[np.lexsort((grid[:,0],grid[:,2]))]
            igrid = interpolate(igrid,1,NO-eNO-0.125,NO+eNO,10)


            CHI_N3O3a = 0
            CHI_N3O3b = 0
            CHI_NO = 0
            NO_p = 0
            den_NO = 0

            for index in igrid:
               if N3O3a_obs == -10: 
                  CHI_N3O3a = 0
               elif index[11] == 0 or index[12] == 0:
                  CHI_N3O3a = tol_max
               else:   
                  CHI_N3O3a = (np.log10(index[12]/index[11])- N3O3a_obs)**2/np.log10(index[12]/index[11])
               if N3O3b_obs == -10: 
                  CHI_N3O3b = 0
               elif index[10] == 0 or index[9] == 0:
                  CHI_N3O3b = tol_max
               else:   
                  CHI_N3O3b = (np.log10(index[9]/index[10])- N3O3b_obs)**2/np.log10(index[9]/index[10])


               CHI_NO = (CHI_N3O3a**2 + CHI_N3O3b**2 )**0.5
               if CHI_NO == 0:
                  NO_p = NO_p
                  den_NO = den_NO
               else:
                  NO_p = index[1] / CHI_NO + NO_p
                  den_NO = 1 / CHI_NO + den_NO

            NOf = NO_p / den_NO



         NO_mc.append(NOf)
         NOe_mc.append(eNO)

      NOff = np.mean(NO_mc)
      eNOff = (np.std(NO_mc)**2+np.mean(NOe_mc)**2)**0.5



# Creation of a constrained grid on N/O

   if NOff == -10:
      grid_c = grid
   else:
      grid_mac = []
      for index in grid:
         if np.abs(index[1] - NOff) > np.abs(eNOff+res_NO):
            continue
         else:
            grid_mac.append(index[0])
            grid_mac.append(index[1])
            grid_mac.append(index[2])
            grid_mac.append(index[3])
            grid_mac.append(index[4])
            grid_mac.append(index[5])
            grid_mac.append(index[6])
            grid_mac.append(index[7])
            grid_mac.append(index[8])
            grid_mac.append(index[9])
            grid_mac.append(index[10])
            grid_mac.append(index[11])
            grid_mac.append(index[12])
            grid_mac.append(index[13])
            grid_mac.append(index[14])

         grid_c = np.reshape(grid_mac,(int(len(grid_mac)/15),15))



# Calculation of O/H and logU


   for monte in range(0,n,1):

      OH_p = 0
      logU_p = 0
      den_OH = 0
      OH_e = 0
      logU_e = 0
      den_OH_e = 0
      tol_max = 1e2

      if tab[0] == 0:
         HI_4m_obs = 0
      else:
         HI_4m_obs = np.random.normal(tab[0],tab[1]+1e-5)
         if HI_4m_obs <= 0: HI_4m_obs = 0
      if tab[2] == 0:
         HI_7m_obs = 0
      else:
         HI_7m_obs = np.random.normal(tab[2],tab[3]+1e-5)
         if HI_7m_obs <= 0: HI_7m_obs = 0
      if tab[4] == 0:
         SIV_10m_obs = 0
      else:
         SIV_10m_obs = np.random.normal(tab[4],tab[5]+1e-5)
         if SIV_10m_obs <= 0: SIV_10m_obs = 0
      if tab[6] == 0:
         HI_12m_obs = 0
      else:
         HI_12m_obs = np.random.normal(tab[6],tab[7]+1e-5)
         if HI_12m_obs <= 0: HI_12m_obs = 0
      if tab[8] == 0:
         NeII_12m_obs = 0
      else:
         NeII_12m_obs = np.random.normal(tab[8],tab[9]+1e-3)
         if NeII_12m_obs <= 0: NeII_12m_obs = 0
      if tab[10] == 0:
            NeIII_15m_obs = 0
      else:
         NeIII_15m_obs = np.random.normal(tab[10],tab[11]+1e-3)
         if NeIII_15m_obs <= 0: NeIII_15m_obs = 0
      if tab[12] == 0:
            SIII_18m_obs = 0
      else:
         SIII_18m_obs = np.random.normal(tab[12],tab[13]+1e-3)
         if SIII_18m_obs <= 0: SIII_18m_obs = 0
      if tab[14] == 0:
            SIII_33m_obs = 0
      else:
         SIII_33m_obs = np.random.normal(tab[14],tab[15]+1e-3)
         if SIII_33m_obs <= 0: SIII_33m_obs = 0
      if tab[16] == 0:
         OIII_52m_obs = 0
      else:
         OIII_52m_obs = np.random.normal(tab[16],tab[17]+1e-5)
         if OIII_52m_obs <= 0: OIII_52m_obs = 0
      if tab[18] == 0:
         NIII_57m_obs = 0
      else:
         NIII_57m_obs = np.random.normal(tab[18],tab[19]+1e-5)
         if NIII_57m_obs <= 0: NIII_57m_obs = 0
      if tab[20] == 0:
         OIII_88m_obs = 0
      else:
         OIII_88m_obs = np.random.normal(tab[20],tab[21]+1e-5)
         if OIII_88m_obs <= 0: OIII_88m_obs = 0
      if tab[22] == 0:
         NII_122m_obs = 0
      else:
         NII_122m_obs = np.random.normal(tab[22],tab[23]+1e-3)
         if NII_122m_obs <= 0: NII_122m_obs = 0
      if HI_4m_obs == 0 or NeII_12m_obs == 0 or NeIII_15m_obs == 0:
         Ne23a_obs = -10
      else:
         Ne23a_obs = np.log10((NeII_12m_obs + NeIII_15m_obs)/HI_4m_obs)
      if HI_7m_obs == 0 or NeII_12m_obs == 0 or NeIII_15m_obs == 0:
         Ne23b_obs = -10
      else:
         Ne23b_obs = np.log10((NeII_12m_obs + NeIII_15m_obs)/HI_7m_obs)
      if HI_12m_obs == 0 or NeII_12m_obs == 0 or NeIII_15m_obs == 0:
         Ne23c_obs = -10
      else:
         Ne23c_obs = np.log10((NeII_12m_obs + NeIII_15m_obs)/HI_12m_obs)
      if NeII_12m_obs == 0 or NeIII_15m_obs == 0:
         Ne2Ne3_obs = -10
      else:
         Ne2Ne3_obs = np.log10((NeII_12m_obs / NeIII_15m_obs))
      if HI_4m_obs == 0 or SIV_10m_obs == 0 or SIII_18m_obs == 0:
         S34a_obs = -10
      else:
         S34a_obs = np.log10((SIV_10m_obs + SIII_18m_obs)/HI_4m_obs)
      if HI_7m_obs == 0 or SIV_10m_obs == 0 or SIII_18m_obs == 0:
         S34b_obs = -10
      else:
         S34b_obs = np.log10((SIV_10m_obs + SIII_18m_obs)/HI_7m_obs)
      if HI_12m_obs == 0 or SIV_10m_obs == 0 or SIII_18m_obs == 0:
         S34c_obs = -10
      else:
         S34c_obs = np.log10((SIV_10m_obs + SIII_18m_obs)/HI_12m_obs)
      if HI_4m_obs == 0 or SIV_10m_obs == 0 or SIII_33m_obs == 0:
         S34d_obs = -10
      else:
         S34d_obs = np.log10((SIV_10m_obs + SIII_33m_obs)/HI_4m_obs)
      if HI_7m_obs == 0 or SIV_10m_obs == 0 or SIII_33m_obs == 0:
         S34e_obs = -10
      else:
         S34e_obs = np.log10((SIV_10m_obs + SIII_33m_obs)/HI_7m_obs)
      if HI_12m_obs == 0 or SIV_10m_obs == 0 or SIII_33m_obs == 0:
         S34f_obs = -10
      else:
         S34f_obs = np.log10((SIV_10m_obs + SIII_33m_obs)/HI_12m_obs)
      if SIV_10m_obs == 0  or SIII_18m_obs == 0:
         S3S4a_obs = -10
      else:
         S3S4a_obs = np.log10((SIII_18m_obs / SIV_10m_obs))
      if SIV_10m_obs == 0  or SIII_33m_obs == 0:
         S3S4b_obs = -10
      else:
         S3S4b_obs = np.log10(SIII_33m_obs / SIV_10m_obs)
      if HI_4m_obs == 0 or NIII_57m_obs == 0 or NII_122m_obs == 0:
         N23a_obs = -10
      else:
         N23a_obs = np.log10((NII_122m_obs + NIII_57m_obs)/HI_4m_obs)
      if HI_7m_obs == 0 or NIII_57m_obs == 0 or NII_122m_obs == 0:
         N23b_obs = -10
      else:
         N23b_obs = np.log10((NII_122m_obs + NIII_57m_obs)/HI_7m_obs)
      if HI_12m_obs == 0 or NIII_57m_obs == 0 or NII_122m_obs == 0:
         N23c_obs = -10
      else:
         N23c_obs = np.log10((NII_122m_obs + NIII_57m_obs)/HI_12m_obs)
      if NII_122m_obs == 0 or NIII_57m_obs == 0:
         N2N3_obs = -10
      else:
         N2N3_obs = np.log10((NII_122m_obs / NIII_57m_obs))
      if NII_122m_obs == 0 or OIII_52m_obs == 0:
         O3N2a_obs = -10
      else:
         O3N2a_obs = np.log10((OIII_52m_obs/NII_122m_obs))
      if NII_122m_obs == 0 or OIII_88m_obs == 0:
         O3N2b_obs = -10
      else:
         O3N2b_obs = np.log10((OIII_88m_obs/NII_122m_obs))

      if  Ne23a_obs == -10 and Ne23b_obs == -10 and Ne23c_obs == -10 and S34a_obs == -10 and S34b_obs == -10 and S34c_obs == -10 and S34d_obs == -10 and S34e_obs == -10 and S34f_obs == -10 and Ne2Ne3_obs == -10 and S3S4a_obs == -10 and S3S4b_obs == -10 and N23a_obs == -10 and N23b_obs == -10 and N23c_obs == -10 and N2N3_obs == -10 and O3N2a_obs == -10 and O3N2b_obs == -10: 
         OH = 0
         logU = 0
      else:
         CHI_Ne23a = 0
         CHI_Ne23b = 0
         CHI_Ne23c = 0
         CHI_Ne2Ne3 = 0
         CHI_S34a = 0
         CHI_S34b = 0
         CHI_S34c3 = 0
         CHI_S34d3 = 0
         CHI_S34e = 0
         CHI_S34f = 0
         CHI_S3S4a = 0
         CHI_S3S4b = 0
         CHI_N23a = 0
         CHI_N23b = 0
         CHI_N23c = 0
         CHI_N2N3 = 0
         CHI_O3N2a = 0
         CHI_O3N2b = 0
         CHI_OH = 0
         for index in grid_c:
            if Ne23a_obs == -10: 
               CHI_Ne23a = 0
            elif index[3] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23a = tol_max
            else:   
               CHI_Ne23a = (np.log10((index[7]+index[8])/index[3])- Ne23a_obs)**2/np.log10((index[7]+index[8])/index[3])
            if Ne23b_obs == -10: 
               CHI_Ne23b = 0
            elif index[4] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23b = tol_max
            else:   
               CHI_Ne23b = (np.log10((index[7]+index[8])/index[4])- Ne23b_obs)**2/np.log10((index[7]+index[8])/index[4])
            if Ne23c_obs == -10: 
               CHI_Ne23c = 0
            elif index[6] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23c = tol_max
            else:   
               CHI_Ne23c = (np.log10((index[7]+index[8])/index[6])- Ne23c_obs)**2/np.log10((index[7]+index[8])/index[6])
            if Ne2Ne3_obs == -10: 
               CHI_Ne2Ne3 = 0
            elif index[7] == 0 or index[8] == 0:
               CHI_Ne2Ne3 = tol_max
            else:   
               CHI_Ne2Ne3 = (np.log10((index[7]/index[8]))- Ne2Ne3_obs)**2/np.log10((index[7]/index[8]))
            if S34a_obs == -10: 
               CHI_S34a = 0
            elif index[3] == 0 or index[5] == 0 or index[9] == 0:
               CHI_S34a = tol_max
            else:   
               CHI_S34a = (np.log10((index[5]+index[9])/index[3])- S34a_obs)**2/np.log10((index[5]+index[9])/index[3])
            if S34b_obs == -10: 
               CHI_S34b = 0
            elif index[4] == 0 or index[5] == 0 or index[9] == 0:
               CHI_S34b = tol_max
            else:   
               CHI_S34b = (np.log10((index[5]+index[9])/index[4])- S34b_obs)**2/np.log10((index[5]+index[9])/index[4])
            if S34c_obs == -10: 
               CHI_S34c = 0
            elif index[6] == 0 or index[5] == 0 or index[9] == 0: 
               CHI_S34c = tol_max
            else:   
               CHI_S34c = (np.log10((index[5]+index[9])/index[6])- S34c_obs)**2/np.log10((index[5]+index[9])/index[6])
            if S34d_obs == -10: 
               CHI_S34d = 0
            elif index[3] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34d = tol_max
            else:   
               CHI_S34d = (np.log10((index[5]+index[10])/index[3])- S34d_obs)**2/np.log10((index[5]+index[10])/index[3])
            if S34e_obs == -10: 
               CHI_S34e = 0
            elif index[4] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34e = tol_max
            else:   
               CHI_S34d = (np.log10((index[5]+index[10])/index[4])- S34e_obs)**2/np.log10((index[5]+index[10])/index[4])
            if S34f_obs == -10: 
               CHI_S34f = 0
            elif index[6] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34f = tol_max
            else:   
               CHI_S34f = (np.log10((index[5]+index[10])/index[6])- S34f_obs)**2/np.log10((index[5]+index[10])/index[6])
            if S3S4a_obs == -10: 
               CHI_S3S4a = 0
            elif index[5] == 0 or index[9] == 0 : 
               CHI_S3S4a = tol_max
            else:   
               CHI_S3S4a = (np.log10(index[9]/index[5])- S3S4a_obs)**2/np.log10((index[9]/index[5]))
            if S3S4b_obs == -10: 
               CHI_S3S4b = 0
            elif index[5] == 0 or index[10] == 0 : 
               CHI_S3S4b = tol_max
            else:   
               CHI_S3S4b = (np.log10(index[10]/index[5])- S3S4b_obs)**2/np.log10((index[10]/index[5]))
            if N23a_obs == -10: 
               CHI_N23a = 0
            elif index[3] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23a = tol_max
            else:   
               CHI_N23a = (np.log10((index[12]+index[14])/index[3])- N23a_obs)**2/np.log10((index[12]+index[14])/index[3])
            if N23b_obs == -10: 
               CHI_N23b = 0
            elif index[4] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23b = tol_max
            else:   
               CHI_N23b = (np.log10((index[12]+index[14])/index[4])- N23b_obs)**2/np.log10((index[12]+index[14])/index[4])
            if N23c_obs == -10: 
               CHI_N23c = 0
            elif index[6] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23c = tol_max
            else:   
               CHI_N23c = (np.log10((index[12]+index[14])/index[6])- N23c_obs)**2/np.log10((index[12]+index[14])/index[6])
            if N2N3_obs == -10: 
               CHI_N2N3 = 0
            elif index[12] == 0 or index[14] == 0: 
               CHI_N2N3 = tol_max
            else:   
               CHI_N2N3 = (np.log10(index[14]/index[12])- N2N3_obs)**2/np.log10((index[14]/index[12]))
            if O3N2a_obs == -10: 
               CHI_O3N2a = 0
            elif index[11] == 0 or index[14] == 0: 
               CHI_O3N2a = tol_max
            else:   
               CHI_O3N2a = (np.log10(index[11]/index[14])- O3N2a_obs)**2/np.log10((index[11]/index[14]))
            if O3N2b_obs == -10: 
               CHI_O3N2b = 0
            elif index[13] == 0 or index[14] == 0: 
               CHI_O3N2b = tol_max
            else:   
               CHI_O3N2b = (np.log10(index[13]/index[14])- O3N2b_obs)**2/np.log10((index[13]/index[14]))

            CHI_OH = (CHI_Ne23a**2 + CHI_Ne23b**2 + CHI_Ne23c**2 + CHI_Ne2Ne3**2 + CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S3S4a**2 + CHI_S3S4b**2+CHI_N23a**2+CHI_N23b**2+CHI_N23c**2+CHI_N2N3**2+CHI_O3N2a**2+CHI_O3N2b**2)**0.5

            if CHI_OH == 0:
               OH_p = OH_p
               logU_p = logU_p
               den_OH = den_OH

            else:
               OH_p = index[0] / (CHI_OH) + OH_p
               logU_p = index[2] / (CHI_OH) + logU_p
               den_OH = 1 / (CHI_OH) + den_OH


         OH = OH_p / den_OH
         logU = logU_p / den_OH


#Calculation of error of O/H and logU

      if  Ne23a_obs == -10 and Ne23b_obs == -10 and Ne23c_obs == -10 and S34a_obs == -10 and S34b_obs == -10 and S34c_obs == -10 and S34d_obs == -10 and S34e_obs == -10 and S34f_obs == -10 and Ne2Ne3_obs == -10 and S3S4a_obs == -10 and S3S4b_obs == -10 and N23a_obs == -10 and N23b_obs == -10 and N23c_obs == -10 and N2N3_obs == -10 and O3N2a_obs == -10 and O3N2b_obs == -10: 
         eOH = 0
         elogU = 0
      else:
         CHI_Ne23a = 0
         CHI_Ne23b = 0
         CHI_Ne23c = 0
         CHI_Ne2Ne3 = 0
         CHI_S34a = 0
         CHI_S34b = 0
         CHI_S34c3 = 0
         CHI_S34d3 = 0
         CHI_S34e = 0
         CHI_S34f = 0
         CHI_S3S4a = 0
         CHI_S3S4b = 0
         CHI_N23a = 0
         CHI_N23b = 0
         CHI_N23c = 0
         CHI_N2N3 = 0
         CHI_O3N2a = 0
         CHI_O3N2b = 0
         CHI_OH = 0
         for index in grid_c:
            if Ne23a_obs == -10: 
               CHI_Ne23a = 0
            elif index[3] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23a = tol_max
            else:   
               CHI_Ne23a = (np.log10((index[7]+index[8])/index[3])- Ne23a_obs)**2/np.log10((index[7]+index[8])/index[3])
            if Ne23b_obs == -10: 
               CHI_Ne23b = 0
            elif index[4] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23b = tol_max
            else:   
               CHI_Ne23b = (np.log10((index[7]+index[8])/index[4])- Ne23b_obs)**2/np.log10((index[7]+index[8])/index[4])
            if Ne23c_obs == -10: 
               CHI_Ne23c = 0
            elif index[6] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23c = tol_max
            else:   
               CHI_Ne23c = (np.log10((index[7]+index[8])/index[6])- Ne23c_obs)**2/np.log10((index[7]+index[8])/index[6])
            if Ne2Ne3_obs == -10: 
               CHI_Ne2Ne3 = 0
            elif index[7] == 0 or index[8] == 0:
               CHI_Ne2Ne3 = tol_max
            else:   
               CHI_Ne2Ne3 = (np.log10(index[7]/index[8])- Ne2Ne3_obs)**2/np.log10((index[7]/index[8]))
            if S34a_obs == -10: 
               CHI_S34a = 0
            elif index[3] == 0 or index[5] == 0 or index[9] == 0:
               CHI_S34a = tol_max
            else:   
               CHI_S34a = (np.log10((index[5]+index[9])/index[3])- S34a_obs)**2/np.log10((index[5]+index[9])/index[3])
            if S34b_obs == -10: 
               CHI_S34b = 0
            elif index[4] == 0 or index[5] == 0 or index[9] == 0:
               CHI_S34b = tol_max
            else:   
               CHI_S34b = (np.log10((index[5]+index[9])/index[4])- S34b_obs)**2/np.log10((index[5]+index[9])/index[4])
            if S34c_obs == -10: 
               CHI_S34c = 0
            elif index[6] == 0 or index[5] == 0 or index[9] == 0: 
               CHI_S34c = tol_max
            else:   
               CHI_S34c = (np.log10((index[5]+index[9])/index[6])- S34c_obs)**2/np.log10((index[5]+index[9])/index[6])
            if S34d_obs == -10: 
               CHI_S34d = 0
            elif index[3] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34d = tol_max
            else:   
               CHI_S34d = (np.log10((index[5]+index[10])/index[3])- S34d_obs)**2/np.log10((index[5]+index[10])/index[3])
            if S34e_obs == -10: 
               CHI_S34e = 0
            elif index[4] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34e = tol_max
            else:   
               CHI_S34d = (np.log10((index[5]+index[10])/index[4])- S34e_obs)**2/np.log10((index[5]+index[10])/index[4])
            if S34f_obs == -10: 
               CHI_S34f = 0
            elif index[6] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34f = tol_max
            else:   
               CHI_S34f = (np.log10((index[5]+index[10])/index[6])- S34f_obs)**2/np.log10((index[5]+index[10])/index[6])
            if S3S4a_obs == -10: 
               CHI_S3S4a = 0
            elif index[5] == 0 or index[9] == 0 : 
               CHI_S3S4a = tol_max
            else:   
               CHI_S3S4a = (np.log10(index[9]/index[5])- S3S4a_obs)**2/np.log10((index[9]/index[5]))
            if S3S4b_obs == -10: 
               CHI_S3S4b = 0
            elif index[5] == 0 or index[10] == 0 : 
               CHI_S3S4b = tol_max
            else:   
               CHI_S3S4b = (np.log10(index[10]/index[5])- S3S4b_obs)**2/np.log10((index[10]/index[5]))
            if N23a_obs == -10: 
               CHI_N23a = 0
            elif index[3] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23a = tol_max
            else:   
               CHI_N23a = (np.log10((index[12]+index[14])/index[3])- N23a_obs)**2/np.log10((index[12]+index[14])/index[3])
            if N23b_obs == -10: 
               CHI_N23b = 0
            elif index[4] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23b = tol_max
            else:   
               CHI_N23b = (np.log10((index[12]+index[14])/index[4])- N23b_obs)**2/np.log10((index[12]+index[14])/index[4])
            if N23c_obs == -10: 
               CHI_N23c = 0
            elif index[6] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23c = tol_max
            else:   
               CHI_N23c = (np.log10((index[12]+index[14])/index[6])- N23c_obs)**2/np.log10((index[12]+index[14])/index[6])
            if N2N3_obs == -10: 
               CHI_N2N3 = 0
            elif index[12] == 0 or index[14] == 0: 
               CHI_N2N3 = tol_max
            else:   
               CHI_N2N3 = (np.log10(index[14]/index[12])- N2N3_obs)**2/np.log10((index[14]/index[12]))
            if O3N2a_obs == -10: 
               CHI_O3N2a = 0
            elif index[11] == 0 or index[14] == 0: 
               CHI_O3N2a = tol_max
            else:   
               CHI_O3N2a = (np.log10(index[11]/index[14])- O3N2a_obs)**2/np.log10((index[11]/index[14]))
            if O3N2b_obs == -10: 
               CHI_O3N2b = 0
            elif index[13] == 0 or index[14] == 0: 
               CHI_O3N2b = tol_max
            else:   
               CHI_O3N2b = (np.log10(index[13]/index[14])- O3N2b_obs)**2/np.log10((index[13]/index[14]))



            CHI_OH = (CHI_Ne23a**2 + CHI_Ne23b**2 + CHI_Ne23c**2 + CHI_Ne2Ne3**2 + CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S3S4a**2 + CHI_S3S4b**2+CHI_N23a**2+CHI_N23b**2+CHI_N23c**2+CHI_N2N3**2+CHI_O3N2a**2+CHI_O3N2b**2)**0.5

            if CHI_OH == 0:
               OH_e = OH_e
               logU_e = logU_e
               den_OH_e = den_OH_e
            else:
               OH_e = (index[0] - OH)**2 / (CHI_OH) + OH_e
               logU_e = (index[2] - logU)**2 / (CHI_OH) + logU_e
               den_OH_e = 1 / (CHI_OH) + den_OH_e 

         eOH = OH_e / den_OH_e
         elogU = logU_e / den_OH_e 


# Iterations for interpolated models

      if inter == 0 or OH == 0:
         OHf = OH
         logUf = logU
      elif inter == 1:
         igrid = interpolate(grid_c,2,logU-elogU-0.25,logU+elogU,10)
         igrid = igrid[np.lexsort((igrid[:,1],igrid[:,2]))]
         igrid = interpolate(igrid,0,OH-eOH-0.1,OH+eOH,10)
         igrid = igrid[np.lexsort((igrid[:,0],igrid[:,2]))]

         CHI_Ne23a = 0
         CHI_N23b = 0
         CHI_N23c = 0
         CHI_Ne2Ne3 = 0
         CHI_S34a = 0
         CHI_S34b = 0
         CHI_S34d = 0
         CHI_S34e = 0
         CHI_S34f = 0
         CHI_S3S4a = 0
         CHI_S3S4b = 0
         CHI_N23a = 0
         CHI_N23b = 0
         CHI_N23c = 0
         CHI_N2N3 = 0
         CHI_O3N2a = 0
         CHI_O3N2b = 0
         CHI_OH = 0
         OH_p = 0
         logU_p = 0
         den_OH = 0
      

         for index in igrid:
            if Ne23a_obs == -10: 
               CHI_Ne23a = 0
            elif index[3] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23a = tol_max
            else:   
               CHI_Ne23a = (np.log10((index[7]+index[8])/index[3])- Ne23a_obs)**2/np.log10((index[7]+index[8])/index[3])
            if Ne23b_obs == -10: 
               CHI_Ne23b = 0
            elif index[4] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23b = tol_max
            else:   
               CHI_Ne23b = (np.log10((index[7]+index[8])/index[4])- Ne23b_obs)**2/np.log10((index[7]+index[8])/index[4])
            if Ne23c_obs == -10: 
               CHI_Ne23c = 0
            elif index[6] == 0 or index[7] == 0 or index[8] == 0:
               CHI_Ne23c = tol_max
            else:   
               CHI_Ne23c = (np.log10((index[7]+index[8])/index[6])- Ne23c_obs)**2/np.log10((index[7]+index[8])/index[6])
            if Ne2Ne3_obs == -10: 
               CHI_Ne2Ne3 = 0
            elif index[7] == 0 or index[8] == 0:
               CHI_Ne2Ne3 = tol_max
            else:   
               CHI_Ne2Ne3 = (np.log10((index[7]/index[8]))- Ne2Ne3_obs)**2/np.log10((index[7]/index[8]))
            if S34a_obs == -10: 
               CHI_S34a = 0
            elif index[3] == 0 or index[5] == 0 or index[9] == 0:
               CHI_S34a = tol_max
            else:   
               CHI_S34a = (np.log10((index[5]+index[9])/index[3])- S34a_obs)**2/np.log10((index[5]+index[9])/index[3])
            if S34b_obs == -10: 
               CHI_S34b = 0
            elif index[4] == 0 or index[5] == 0 or index[9] == 0:
               CHI_S34b = tol_max
            else:   
               CHI_S34b = (np.log10((index[5]+index[9])/index[4])- S34b_obs)**2/np.log10((index[5]+index[9])/index[4])
            if S34c_obs == -10: 
               CHI_S34c = 0
            elif index[6] == 0 or index[5] == 0 or index[9] == 0: 
               CHI_S34c = tol_max
            else:   
               CHI_S34c = (np.log10((index[5]+index[9])/index[6])- S34c_obs)**2/np.log10((index[5]+index[9])/index[6])
            if S34d_obs == -10: 
               CHI_S34d = 0
            elif index[3] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34d = tol_max
            else:   
               CHI_S34d = (np.log10((index[5]+index[10])/index[3])- S34d_obs)**2/np.log10((index[5]+index[10])/index[3])
            if S34e_obs == -10: 
               CHI_S34e = 0
            elif index[4] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34e = tol_max
            else:   
               CHI_S34d = (np.log10((index[5]+index[10])/index[4])- S34e_obs)**2/np.log10((index[5]+index[10])/index[4])
            if S34f_obs == -10: 
               CHI_S34f = 0
            elif index[6] == 0 or index[5] == 0 or index[10] == 0: 
               CHI_S34f = tol_max
            else:   
               CHI_S34f = (np.log10((index[5]+index[10])/index[6])- S34f_obs)**2/np.log10((index[5]+index[10])/index[6])
            if S3S4a_obs == -10: 
               CHI_S3S4a = 0
            elif index[5] == 0 or index[9] == 0 : 
               CHI_S3S4a = tol_max
            else:   
               CHI_S3S4a = (np.log10(index[9]/index[5])- S3S4a_obs)**2/np.log10((index[9]/index[5]))
            if S3S4b_obs == -10: 
               CHI_S3S4b = 0
            elif index[5] == 0 or index[10] == 0 : 
               CHI_S3S4b = tol_max
            else:   
               CHI_S3S4b = (np.log10(index[10]/index[5])- S3S4b_obs)**2/np.log10((index[10]/index[5]))
            if N23a_obs == -10: 
               CHI_N23a = 0
            elif index[3] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23a = tol_max
            else:   
               CHI_N23a = (np.log10((index[12]+index[14])/index[3])- N23a_obs)**2/np.log10((index[12]+index[14])/index[3])
            if N23b_obs == -10: 
               CHI_N23b = 0
            elif index[4] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23b = tol_max
            else:   
               CHI_N23b = (np.log10((index[12]+index[14])/index[4])- N23b_obs)**2/np.log10((index[12]+index[14])/index[4])
            if N23c_obs == -10: 
               CHI_N23c = 0
            elif index[6] == 0 or index[12] == 0 or index[14] == 0: 
               CHI_N23c = tol_max
            else:   
               CHI_N23c = (np.log10((index[12]+index[14])/index[6])- N23c_obs)**2/np.log10((index[12]+index[14])/index[6])
            if N2N3_obs == -10: 
               CHI_N2N3 = 0
            elif index[12] == 0 or index[14] == 0: 
               CHI_N2N3 = tol_max
            else:   
               CHI_N2N3 = (np.log10(index[14]/index[12])- N2N3_obs)**2/np.log10((index[14]/index[12]))
            if O3N2a_obs == -10: 
               CHI_O3N2a = 0
            elif index[11] == 0 or index[14] == 0: 
               CHI_O3N2a = tol_max
            else:   
               CHI_O3N2a = (np.log10(index[11]/index[14])- O3N2a_obs)**2/np.log10((index[11]/index[14]))
            if O3N2b_obs == -10: 
               CHI_O3N2b = 0
            elif index[13] == 0 or index[14] == 0: 
               CHI_O3N2b = tol_max
            else:   
               CHI_O3N2b = (np.log10(index[13]/index[14])- O3N2b_obs)**2/np.log10((index[13]/index[14]))



            CHI_OH = (CHI_Ne23a**2 + CHI_Ne23b**2 + CHI_Ne23c**2 + CHI_Ne2Ne3**2 + CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S3S4a**2 + CHI_S3S4b**2+CHI_N23a**2+CHI_N23b**2+CHI_N23c**2+CHI_N2N3**2+CHI_O3N2a**2+CHI_O3N2b**2)**0.5


            OH_p = index[0] / CHI_OH**2 + OH_p
            logU_p = index[2] / CHI_OH**2 + logU_p
            den_OH = 1 / CHI_OH**2 + den_OH

         if OH == 0:
            OHf = OH
            logUf = logU
         else:         
            OHf = OH_p / den_OH
            logUf = logU_p / den_OH


      OH_mc.append(OHf)
      logU_mc.append(logUf)
      OHe_mc.append(eOH)
      logUe_mc.append(elogU)
   
   OHff = np.mean(OH_mc)
   eOHff = (np.std(OH_mc)**2+np.mean(OHe_mc)**2)**0.5
   logUff = np.mean(logU_mc)
   elogUff = (np.std(logU_mc)**2+np.mean(logUe_mc)**2)**0.5

   




   
   output.append(OHff)
   output.append(eOHff)
   output.append(NOff)
   output.append(eNOff)
   output.append(logUff)
   output.append(elogUff)
         

   if input0.shape >= (24,) and count == 1: continue   
   print (round(100*(count)/float(len(input)),1),'%',grid_type,'', round(OHff,3), round(eOHff,3),'',round(NOff,3), round(eNOff,3), '',round(logUff,3), round(elogUff,3))


out = np.reshape(output,(len(input),31))
if input0.shape == (24,): out = np.delete(out,obj=0,axis=0)



lineas_header = [' HII-CHI-mistry-IR v.1.1 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','Bra  eBra   Pfa   ePfa  S4     eS4   Hua   eHua Ne2-12  eNe2-12  Ne3-15  eNe3-15  S318  eS318  S3-33  eS3-33  O3-52 eO3-52  N3-57  eN3-57  O3-88  eO3-88   N2-122  eN2-122 i O/H     eO/H  N/O    eN/O  logU   elogU']

header = '\n'.join(lineas_header)

np.savetxt(input00+'_hcm-output.dat',out,fmt=' '.join(['%.4f']*24+['%i']+['%.3f']*6),header=header)
print ('________________________________')
print ('Results are stored in ' + input00 + '_hcm-output.dat')
                           
      

      
