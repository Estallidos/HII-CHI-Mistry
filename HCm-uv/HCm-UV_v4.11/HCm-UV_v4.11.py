# Filename: HCm_UV_v4.11.py

import string
import numpy as np
import sys
#sys.stderr = open('errorlog.txt', 'w')


#Function for interpolation of grids

def interpolate(grid,z,zmin,zmax,n):
   ncol = 9
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
print (' This is HII-CHI-mistry for UV version 4.11')
print (' See Perez-Montero, & Amorin (2017) for details')
print ( ' Insert the name of your input text file with some or all of the following columns:')
print (' Lya 1216, CIV 1549, HeII 1640, OIII 1665, CIII 1909, Hb 4861, OIII 5007')
print ('in arbitrary units and reddening corrected. Each column must be given')
print ('with labels and followed by its corresponding flux error.')
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
   input0 = np.genfromtxt(input00,dtype=None,names=True, encoding = 'ascii')
   print ('The input file is:'+input00)
except:
   print ('Input file error: It does not exist or has wrong format')

   sys.exit

print ('')

if input0.size == 1:
   input1 = np.stack((input0,input0))
else:
   input1 = input0


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
   print ('(2) BPASS v.2.1 a_IMF = 1.35, Mup = 300, age = 1Myr')
   print('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      sed = raw_input('Choose SED of the models:')
   else:
      sed = input('Choose SED of the models:')
   if sed == '1' or sed == '2' : question = False 
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



if sed==1 :
   grid1 = np.loadtxt('C17_popstar_uv_v4.0.dat')
   grid2 = np.loadtxt('C17_popstar_logU_adapted_emp_uv_v4.0.dat')
   grid3 = np.loadtxt('C17_popstar_logU-CO_adapted_emp_uv_v4.0.dat')
   if inter == 0:
      sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation'
      print ('No interpolation for the POPSTAR models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF  interpolated'
      print ('Interpolation for the POPSTAR models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
elif sed==2:
   grid1 = np.loadtxt('C17_bpass_uv_v4.1.dat')
   grid2 = np.loadtxt('C17_bpass_logU_adapted_emp_uv_v4.1.dat')
   grid3 = np.loadtxt('C17_bpass_logU-CO_adapted_emp_uv_v4.1.dat')
   if inter == 0:
      sed_type = 'BPASS a_IMF = 1.35, M_up = 300, age = 1Myr. No interpolation'
      print ('No interpolation for theBPASS models is going to be used.')
      print ('The grid has a resolution of 0.1 dex for O/H and 0.125 dex for N/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'BPASS a_IMF = 1.35, M_up = 300, age = 1Myr interpolated'
      print ('Interpolation for theBPASS models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
      res_CO = 0.125



grids = []
OHffs = []
eOHffs = []
COffs = []
eCOffs = []
logUffs = []
elogUffs = []

Label_ID = False
Label_Lya = False
Label_eLya = False
Label_CIV = False
Label_eCIV = False
Label_HeII = False
Label_eHeII = False
Label_OIII_1665 = False
Label_eOIII_1665 = False
Label_CIII = False
Label_eCIII = False
Label_OIII_5007 = False
Label_eOIII_5007 = False
Label_Hbeta = False
Label_eHbeta = False


for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == 'Lya_1216':
      Label_Lya = True
   if input1.dtype.names[col] == 'eLya_1216':
      Label_eLya = True
   if input1.dtype.names[col] == 'CIV_1549':
      Label_CIV = True
   if input1.dtype.names[col] == 'eCIV_1549':
      Label_eCIV = True
   if input1.dtype.names[col] == 'HeII_1640':
      Label_HeII = True
   if input1.dtype.names[col] == 'eHeII_1640':
      Label_eHeII = True
   if input1.dtype.names[col] == 'OIII_1665':
      Label_OIII_1665 = True
   if input1.dtype.names[col] == 'eOIII_1665':
      Label_eOIII_1665 = True
   if input1.dtype.names[col] == 'CIII_1909':
      Label_CIII = True
   if input1.dtype.names[col] == 'eCIII_1909':
      Label_eCIII = True
   if input1.dtype.names[col] == 'Hb_4861':
      Label_Hbeta = True
   if input1.dtype.names[col] == 'eHb_4861':
      Label_eHbeta = True
   if input1.dtype.names[col] == 'OIII_5007':
      Label_OIII_5007 = True
   if input1.dtype.names[col] == 'eOIII_5007':
      Label_eOIII_5007 = True


if Label_ID == False:
   Names = np.arange(1,input1.size+1,1)
else:
   Names = input1['ID']
if Label_Lya == False:
   Lya_1216 = np.zeros(input1.size)
else:
   Lya_1216 = input1['Lya_1216']
if Label_eLya == False:
   eLya_1216 = np.zeros(input1.size)
else:
   eLya_1216 = input1['eLya_1216']
if Label_CIV == False:
   CIV_1549 = np.zeros(input1.size)
else:
   CIV_1549 = input1['CIV_1549']
if Label_eCIV == False:
   eCIV_1549 = np.zeros(input1.size)
else:
   eCIV_1549 = input1['eCIV_1549']
if Label_HeII == False:
   HeII_1640 = np.zeros(input1.size)
else:
   HeII_1640 = input1['HeII_1640']
if Label_eHeII == False:
   eHeII_1640 = np.zeros(input1.size)
else:
   eHeII_1640 = input1['eHeII_1640']
if Label_OIII_1665 == False:
   OIII_1665 = np.zeros(input1.size)
else:
   OIII_1665 = input1['OIII_1665']
if Label_eOIII_1665 == False:
   eOIII_1665 = np.zeros(input1.size)
else:
   eOIII_1665 = input1['eOIII_1665']
if Label_CIII == False:
   CIII_1909 = np.zeros(input1.size)
else:
   CIII_1909 = input1['CIII_1909']
if Label_eCIII == False:
   eCIII_1909 = np.zeros(input1.size)
else:
   eCIII_1909 = input1['eCIII_1909']
if Label_Hbeta == False:
   Hb_4861 = np.zeros(len(input1))
else:
   Hb_4861 = input1['Hb_4861']
if Label_eHbeta == False:
   eHb_4861 = np.zeros(input1.size)
else:
   eHb_4861 = input1['eHb_4861']
if Label_OIII_5007 == False:
   OIII_5007 = np.zeros(input1.size)
else:
   OIII_5007 = input1['OIII_5007']
if Label_eOIII_5007 == False:
   eOIII_5007 = np.zeros(input1.size)
else:
   eOIII_5007 = input1['eOIII_5007']


output = np.zeros(input1.size, dtype=[('ID', 'U12'), ('Lya_1216', float),('eLya_1216', float),('CIV_1549', float),('eCIV_1549', float),('HeII_1640', float),('eHeII_1640', float),('OIII_1665', float),('eOIII_1665', float),('CIII_1909', float),('eCIII_1909', float),('Hb_4861', float),('eHb_4861', float),('OIII_5007', float),('eOIII_5007', float),('grid', int),('OH', float),('eOH', float),('CO', float),('eCO', float),('logU', float),('elogU', float)] )

output['ID'] = Names
output['Lya_1216'] = Lya_1216
output['eLya_1216'] = eLya_1216
output['CIV_1549'] = CIV_1549
output['eCIV_1549'] = eCIV_1549
output['HeII_1640'] = HeII_1640
output['eHeII_1640'] = eHeII_1640
output['OIII_1665'] = OIII_1665
output['eOIII_1665'] = eOIII_1665
output['CIII_1909'] = CIII_1909
output['eCIII_1909'] = eCIII_1909
output['Hb_4861'] = Hb_4861
output['eHb_4861'] = eHb_4861
output['OIII_5007'] = OIII_5007
output['eOIII_5007'] = eOIII_5007


print ('Reading grids ....')
print ('')
print ('')
print ('----------------------------------------------------------------')
print ('(%)   ID    Grid  12+log(O/H)  log(C/O)    log(U)')
print ('-----------------------------------------------------------------')

# Beginning of loop of calculation

count = 0
for tab in range(0,len(input1),1):

   count = count + 1


   OH_mc = []
   CO_mc = []
   logU_mc = []
   OHe_mc = []
   COe_mc = []
   logUe_mc = []  



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

      Lya_1216_obs = 0
      if Lya_1216[tab] == 0:
         Lya_1216_obs = 0
      else:
         while Lya_1216_obs <= 0:
            Lya_1216_obs = np.random.normal(Lya_1216[tab],eLya_1216[tab]+1e-5)
      CIV_1549_obs = 0
      if CIV_1549[tab] == 0:
         CIV_1549_obs = 0
      else:
         while CIV_1549_obs <= 0:
            CIV_1549_obs = np.random.normal(CIV_1549[tab],eCIV_1549[tab]+1e-5)
      HeII_1640_obs = 0
      if HeII_1640[tab] == 0:
         HeII_1640_obs = 0
      else:
         if HeII_1640_obs <= 0:
            HeII_1640_obs = np.random.normal(HeII_1640[tab],eHeII_1640[tab]+1e-5)
      OIII_1665_obs = 0
      if OIII_1665[tab] == 0:
         OIII_1665_obs = 0
      else:
         while OIII_1665_obs <= 0:
            OIII_1665_obs = np.random.normal(OIII_1665[tab],eOIII_1665[tab]+1e-5)
      CIII_1909_obs = 0
      if CIII_1909[tab] == 0:
         CIII_1909_obs = 0
      else:
         while CIII_1909_obs <= 0:
            CIII_1909_obs = np.random.normal(CIII_1909[tab],eCIII_1909[tab]+1e-5)
      Hb_4861_obs = 0
      if Hb_4861[tab] == 0:
         Hb_4861_obs = 0
      else:
         while Hb_4861_obs <= 0:
            Hb_4861_obs = np.random.normal(Hb_4861[tab],eHb_4861[tab]+1e-5)
      OIII_5007_obs = 0
      if OIII_5007[tab] == 0:
         OIII_5007_obs = 0
      else:
         while OIII_5007_obs <= 0:
            OIII_5007_obs = np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-5)
      if OIII_1665_obs == 0 or OIII_5007_obs == 0:
         ROIII_obs = 0
      else:
         ROIII_obs = OIII_5007_obs/OIII_1665_obs
      if Lya_1216_obs == 0 or CIII_1909_obs == 0:
         C34_obs = 0
      else:
         C34_obs = (CIII_1909_obs + CIV_1549_obs) / (Lya_1216_obs)
      if HeII_1640_obs == 0 or CIII_1909_obs == 0:
         C34He2_obs = 0
      else:
         C34He2_obs = (CIII_1909_obs + CIV_1549_obs) / (HeII_1640_obs)
      if CIII_1909_obs == 0 or OIII_1665_obs == 0:
         C3O3_obs = -10
      else:   
         C3O3_obs = np.log10((CIII_1909_obs) / (OIII_1665_obs))
      if CIII_1909_obs == 0 or CIV_1549_obs == 0:
         C3C4_obs = 0
      else:
         C3C4_obs = (CIII_1909_obs/CIV_1549_obs)
      if CIII_1909_obs == 0 or Hb_4861_obs == 0:
         C34Hb_obs = 0
      else:
         C34Hb_obs = (CIII_1909_obs + CIV_1549_obs) / Hb_4861_obs
              
         



# Selection of grid
   
      if OIII_1665[tab] > 0 and OIII_5007[tab] > 0:
         grid = grid1
         if monte == n-1: grids.append(1)
         grid_type = 1
      elif OIII_1665[tab] > 0 and CIII_1909[tab] > 0:
         grid = grid2
         if monte == n-1: grids.append(2)
         grid_type = 2
      else:
         grid = grid3
         if monte == n-1: grids.append(3)
         grid_type = 3



  # Calculation of C/O 

  
      if C3O3_obs == -10:
         CO = -10
      else:
         CHI_ROIII = 0
         CHI_C3O3 = 0   
         CHI_CO = 0

         for index in grid:
            if ROIII_obs == 0:
               CHI_ROIII = 0
            elif index[6] == 0 or index[8] == 0:
               CHI_ROIII = tol_max
            else:
               CHI_ROIII = (index[8]/index[6] - ROIII_obs)**2/(index[8]/index[6])
            if C3O3_obs == -10:
               CHI_C3O3 = 0          
            elif index[7] == 0 or index[6] == 0:
               CHI_C3O3 = tol_max
            else:
               CHI_C3O3 =(np.log10((index[7])/index[6]) - C3O3_obs)**2/np.log10((index[7])/(index[6]+1e-5))

            CHI_CO = (CHI_ROIII**2 + CHI_C3O3**2 )**0.5


            if CHI_CO == 0:
               CO_p = CO_p
               den_CO = den_CO
            else:
               CO_p = index[1] /np.exp(CHI_CO) + CO_p
               den_CO = 1 / np.exp(CHI_CO) + den_CO

         CO = CO_p / den_CO 

  # Calculation of C/O error

  
      if C3O3_obs == -10:
         eCO = 0
      else:
         CHI_ROIII = 0
         CHI_C3O3 = 0   
         CHI_CO = 0

         for index in grid:
            if ROIII_obs == 0:
               CHI_ROIII = 0
            elif index[6] == 0 or index[8] == 0:
               CHI_ROIII = tol_max
            else:
               CHI_ROIII = (index[8]/index[6] - ROIII_obs)**2/(index[8]/index[6])
            if C3O3_obs == -10:
               CHI_C3O3 = 0          
            elif index[7] == 0 or index[6] == 0:
               CHI_C3O3 = tol_max
            else:
               CHI_C3O3 =(np.log10((index[7])/index[6]) - C3O3_obs)**2/np.log10((index[7])/(index[6]+1e-5))

            CHI_CO = (CHI_ROIII**2 + CHI_C3O3**2 )**0.5


            if CHI_CO == 0:
               CO_e = CO_e
               den_CO_e = den_CO_e  
            else:
               CO_e = (index[1] - CO)**2 / np.exp(CHI_CO) + CO_e
               den_CO_e = 1 /np.exp(CHI_CO) + den_CO_e  


         eCO = CO_e / den_CO_e 

      
# Calculation of O/H and log U

      if C34_obs == 0 and ROIII_obs == 0 and C34Hb_obs == 0 and C34He2_obs == 0 :
         OH = 0
         logU = 0
      else:
         CHI_ROIII = 0 
         CHI_C3C4 = 0
         CHI_C34He2 = 0
         CHI_C34 = 0
         CHI_C34Hb = 0
         CHI_OH = 0
         for index in grid:
            if CO > -10 and np.abs(index[1] - CO) > np.abs(eCO+0.125):
               continue
            if CIV_1549_obs > 0 and index[4] == 0:
               continue
            if HeII_1640_obs > 0 and index[5] == 0:
               continue
            else:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index[6] == 0 or index[8] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index[8]/index[6] - ROIII_obs)**2/(index[8]/index[6])
               if C34_obs == 0:
                  CHI_C34 = 0
               elif index[3] == 0 or index[7] == 0:
                  CHI_C34 = tol_max
               else:
                  CHI_C34 = ((index[7]+index[4])/index[3] - C34_obs)**2/((index[7]+index[4])/index[3])
               if C34He2_obs == 0:
                  CHI_C34He2 = 0
               elif index[5] == 0 or index[7] == 0:
                  CHI_C34He2 = tol_max
               else:
                  CHI_C34He2 = ((index[7]+index[4])/index[5] - C34He2_obs)**2/((index[7]+index[4])/index[5])
               if C34Hb_obs == 0:
                  CHI_C34Hb = 0
               elif index[7] == 0:
                  CHI_C34Hb = tol_max
               else:
                  CHI_C34Hb = (index[7]+index[4] - C34Hb_obs)**2/(index[7]+index[4])
               if C3C4_obs == 0:
                  CHI_C3C4 = 0
               elif index[4] == 0 or index[7] == 0:
                  CHI_C3C4 = tol_max
               else:
                  CHI_C3C4 = (index[7]/index[4] - C3C4_obs)**2/(index[7]/index[4])

               if C34Hb_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34Hb**2  + CHI_C3C4**2)**0.5
               else:	
                  CHI_OH = (CHI_ROIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_C3C4**2 )**0.5



            if CHI_OH == 0:
               OH_p = OH_p
               logU_p = logU_p
               den_OH = den_OH
            else:
               OH_p = index[0] / np.exp(CHI_OH) + OH_p
               logU_p = index[2] / np.exp(CHI_OH) + logU_p
               den_OH = 1 /np.exp(CHI_OH) + den_OH


         OH = OH_p / den_OH
         logU = logU_p / den_OH

# Calculation of error of O/H and logU

      if C34_obs == 0 and ROIII_obs == 0 and C34Hb_obs == 0  and C34He2_obs == 0:
         eOH = 0
         elogU = 0
      else:
         CHI_ROIII = 0 
         CHI_C3C4 = 0
         CHI_C34 = 0
         CHI_C34He2 = 0
         CHI_C34Hb = 0
         CHI_OH = 0

         for index in grid:
            if CO > -10 and np.abs(index[1] - CO) > np.abs(eCO+res_CO):
               continue
            if CIV_1549_obs > 0 and index[4] == 0:
               continue
            if HeII_1640_obs > 0 and index[5] == 0:
               continue
            else:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index[6] == 0 or index[8] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index[8]/index[6] - ROIII_obs)**2/(index[8]/index[6])
               if C34_obs == 0:
                  CHI_C34 = 0
               elif index[3] == 0 or index[7] == 0:
                  CHI_C34 = tol_max
               else:
                  CHI_C34 = ((index[7]+index[4])/index[3] - C34_obs)**2/((index[7]+index[4])/index[3])
               if C34He2_obs == 0:
                  CHI_C34He2 = 0
               elif index[5] == 0 or index[7] == 0:
                  CHI_C34He2 = tol_max
               else:
                  CHI_C34He2 = ((index[7]+index[4])/index[5] - C34He2_obs)**2/((index[7]+index[4])/index[5])
               if C34Hb_obs == 0:
                  CHI_C34Hb = 0
               elif index[7] == 0:
                  CHI_C34Hb = tol_max
               else:
                  CHI_C34Hb = (index[7]+index[4] - C34Hb_obs)**2/(index[7]+index[4])
               if C3C4_obs == 0:
                  CHI_C3C4 = 0
               elif index[4] == 0 or index[7] == 0:
                  CHI_C3C4 = tol_max
               else:
                  CHI_C3C4 = (index[7]/index[4] - C3C4_obs)**2/(index[7]/index[4])


               if C34Hb_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34Hb**2 + CHI_C3C4**2)**0.5
               else:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_C3C4**2  )**0.5

            if CHI_OH == 0:
               OH_e = OH_e
               logU_e = logU_e
               den_OH_e = den_OH_e
            else:
               OH_e = (index[0] - OH)**2 /np.exp(CHI_OH) + OH_e
               logU_e = (index[2] - logU)**2 /np.exp(CHI_OH) + logU_e
               den_OH_e = 1 /np.exp(CHI_OH) + den_OH_e 

         eOH = OH_e / den_OH_e
         elogU = logU_e / den_OH_e 


# Iterations for interpolated models

      if inter == 0  or (OH == 0 and CO == -10):
         COf = CO
         OHf = OH
         logUf = logU
      elif inter == 1:
         if OH == 0:
            igrid = grid
         else:
            igrid = interpolate(grid,2,logU-elogU-0.25,logU+elogU+0.25,10)
            igrid = igrid[np.lexsort((igrid[:,1],igrid[:,2]))]
            igrid = interpolate(igrid,0,OH-eOH-0.1,OH+eOH+0.1,10)
         if CO == -10:
            igrid = igrid
         else:
            igrid = igrid[np.lexsort((igrid[:,0],igrid[:,2]))]
            igrid = interpolate(igrid,1,CO-eCO-0.125,CO+eCO+0.125,10)


         CHI_ROIII = 0
         CHI_C3O3 = 0   
         CHI_C3C4 = 0
         CHI_C34He2 = 0
         CHI_C34 = 0
         CHI_C34Hb = 0
         CHI_OH = 0
         CHI_CO = 0


         for index in igrid:
            if ROIII_obs == 0:
               CHI_ROIII = 0
            elif index[6] == 0 or index[8] == 0:
               CHI_ROIII = tol_max
            else:
               CHI_ROIII = (index[8]/index[6] - ROIII_obs)**2/(index[8]/index[6])
            if C3O3_obs == -10:
               CHI_C3O3 = 0          
            elif index[7] == 0 or index[6] == 0:
               CHI_C3O3 = tol_max
            else:
               CHI_C3O3 =(np.log10((index[7])/index[6]) - C3O3_obs)**2/np.log10((index[7])/(index[6]+1e-5))
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


               if C34Hb_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34Hb**2 + CHI_C3C4**2)**0.5
               else:	
                  CHI_OH = (CHI_ROIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_C3C4**2 )**0.5



            if CHI_OH == 0:
               OH_p = OH_p
               logU_p = logU_p
               den_OH = den_OH
            else:
               OH_p = index[0] /np.exp(CHI_OH) + OH_p
               logU_p = index[2] /np.exp(CHI_OH) + logU_p
               den_OH = 1 /np.exp(CHI_OH) + den_OH

            CHI_CO = (CHI_ROIII**2 + CHI_C3O3**2 )**0.5

            if CHI_CO == 0:
               CO_p = CO_p
               den_CO = den_CO
            else:
               CO_p = index[1] /np.exp(CHI_CO)**2 + CO_p
               den_CO = 1 /np.exp(CHI_CO)**2 + den_CO

         if CO == -10:
            COf = -10
         else:
            COf = CO_p / den_CO 

         if OH == 0:
            OHf = 0
            logUf = 0
         else:
            OHf = OH_p / den_OH
            logUf = logU_p / den_OH


      OH_mc.append(OHf)
      CO_mc.append(COf)
      logU_mc.append(logUf)
      OHe_mc.append(eOH)
      COe_mc.append(eCO)
      logUe_mc.append(elogU)


   OHff = np.mean(OH_mc)
   eOHff = (np.std(OH_mc)**2+np.mean(OHe_mc)**2)**0.5
   COff = np.mean(CO_mc)
   eCOff = (np.std(CO_mc)**2+np.mean(COe_mc)**2)**0.5
   logUff = np.mean(logU_mc)
   elogUff = (np.std(logU_mc)**2+np.mean(logUe_mc)**2)**0.5

         
   OHffs.append(OHff)
   eOHffs.append(eOHff)
   COffs.append(COff)
   eCOffs.append(eCOff)
   logUffs.append(logUff)
   elogUffs.append(elogUff)
         

   if input0.size == 1 and tab==0: continue

   print (round(100*(count)/float(input1.size),1),'%',Names[tab],grid_type,'', round(OHff,2), round(eOHff,2),'',round(COff,2), round(eCOff,2), '',round(logUff,2), round(elogUff,2))



output['grid'] = grids
output['OH'] = OHffs
output['eOH'] = eOHffs
output['CO'] = COffs
output['eCO'] = eCOffs
output['logU'] = logUffs
output['elogU'] = elogUffs

if input0.size == 1:  output = np.delete(output,obj=1,axis=0)


lineas_header = [' HII-CHI-mistry_UV v.4.11 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','ID.   Lya  eLya   1549   e1549 1640 e1640 1665  e1665  1909   e1909  Hbeta   eHbeta  5007    e5007  i O/H     eO/H  C/O    eC/O  logU   elogU']

header = '\n'.join(lineas_header)

np.savetxt(input00+'_hcm-uv-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*14+['%i']+['%.2f']*6),header=header)
print ('________________________________')
print ('Results are stored in '+input00+'_hcm-uv-output.dat')
            
