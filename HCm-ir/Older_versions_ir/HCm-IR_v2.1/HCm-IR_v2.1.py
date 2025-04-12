# Filename: HII-CHCm-IR_v2.1.py



import string
import numpy as np
import sys
#sys.stderr = open('errorlog.txt', 'w')



#Function for interpolation of grids

def interpolate(grid,z,zmin,zmax,n):
   ncol = 16
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
print ('This is HII-CHI-mistry IR v. 2.1')
print (' See Fernandez-Ontiveros et al (2021) for details')
print  (' Insert the name of your input text file with all or some of the following columns:')
print ('')
print (' HI 4.05m')
print ('HI 7.46m') 
print ('[SIV] 10.5m')
print ('HI 12.4m')
print ('[NeII] 12.8m')
print ('[NeIII] 15.5m')
print ('[SIII] 18.7m')
print ('[SIII] 33.7m')
print ('[OIII] 52m')
print ('[NIII] 57m')
print ('[OIII] 88m') 
print ('[NII] 122m')
print ('[NII] 205m')
print ('')
print ('with their corresponding labels and errors in adjacent columns')
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



if sed==1:
   grid1 = np.loadtxt('C17_popstar_ir_v2.1.dat')
   grid2 = np.loadtxt('C17_popstar_logU_adapted_emp_ir_v2.1.dat')
   grid3 = np.loadtxt('C17_popstar_logU-NO_adapted_emp_ir_v2.1.dat')
   if inter == 0:
      sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation'
      print ('No interpolation for the POPSTAR models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      print ('')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. interpolation'
      print ('Interpolation for the POPSTAR models is going to be used.')
      print ('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
      print ('')
      res_NO = 0.125
elif sed==2:
   grid1 = np.loadtxt('C17_bpass_ir_v2.1.dat')
   grid2 = np.loadtxt('C17_bpass_logU_adapted_emp_ir_v2.1.dat')
   grid3 = np.loadtxt('C17_bpass_logU-NO_adapted_emp_ir_v2.1.dat')
   if inter == 0:
      sed_type = 'BPASS a_IMF = 1.35, M_up = 300, age = 1Myr. No interpolation'
      print ('No interpolation for theBPASS models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      print ('')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'BPASS v.2.1, a_IMF = 1.35, M_up = 300, age = 1Myr. Interpolation'
      print ('Interpolation for theBPASS  models is going to be used.')
      print ('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
      print ('')
      res_NO = 0.125


# Input file reading



grids = []
OHffs = []
eOHffs = []
NOffs = []
eNOffs = []
logUffs = []
elogUffs = []

Label_ID = False
Label_HI_4m = False
Label_eHI_4m = False
Label_HI_7m = False
Label_eHI_7m = False
Label_SIV = False
Label_eSIV = False
Label_HI_12m = False
Label_eHI_12m = False
Label_NeII = False
Label_eNeII = False
Label_NeIII = False
Label_eNeIII = False
Label_SIII_18m = False
Label_eSIII_18m = False
Label_SIII_33m = False
Label_eSIII_33m = False
Label_OIII_52m = False
Label_eOIII_52m = False
Label_NIII = False
Label_eNIII = False
Label_OIII_88m = False
Label_eOIII_88m = False
Label_NII_122m = False
Label_eNII_122m = False
Label_NII_205m = False
Label_eNII_205m = False

for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == 'HI_4m':
      Label_HI_4m = True
   if input1.dtype.names[col] == 'eHI_4m':
      Label_eHI_4m = True
   if input1.dtype.names[col] == 'HI_7m':
      Label_HI_7m = True
   if input1.dtype.names[col] == 'eHI_7m':
      Label_eHI_7m = True
   if input1.dtype.names[col] == 'SIV_10m':
      Label_SIV = True
   if input1.dtype.names[col] == 'eSIV_10m':
      Label_eSIV = True
   if input1.dtype.names[col] == 'HI_12m':
      Label_HI_12m = True
   if input1.dtype.names[col] == 'eHI_12m':
      Label_eHI_12m = True
   if input1.dtype.names[col] == 'NeII_12m':
      Label_NeII = True
   if input1.dtype.names[col] == 'eNeII_12m':
      Label_eNeII = True
   if input1.dtype.names[col] == 'NeIII_15m':
      Label_NeIII = True
   if input1.dtype.names[col] == 'eNeIII_15m':
      Label_eNeIII = True
   if input1.dtype.names[col] == 'SIII_18m':
      Label_SIII_18m = True
   if input1.dtype.names[col] == 'eSIII_18m':
      Label_eSIII_18m = True
   if input1.dtype.names[col] == 'SIII_33m':
      Label_SIII_33m = True
   if input1.dtype.names[col] == 'eSIII_33m':
      Label_eSIII_33m = True
   if input1.dtype.names[col] == 'OIII_52m':
      Label_OIII_52m = True
   if input1.dtype.names[col] == 'eOIII_52m':
      Label_eOIII_52m = True
   if input1.dtype.names[col] == 'NIII_57m':
      Label_NIII = True
   if input1.dtype.names[col] == 'eNIII_57m':
      Label_eNIII = True
   if input1.dtype.names[col] == 'OIII_88m':
      Label_OIII_88m = True
   if input1.dtype.names[col] == 'eOIII_88m':
      Label_eOIII_88m = True
   if input1.dtype.names[col] == 'NII_122m':
      Label_NII_122m = True
   if input1.dtype.names[col] == 'eNII_122m':
      Label_eNII_122m = True
   if input1.dtype.names[col] == 'NII_205m':
      Label_NII_205m = True
   if input1.dtype.names[col] == 'eNII_205m':
      Label_eNII_205m = True


if Label_ID == False:
   Names = np.arange(1,input1.size+1,1)
else:
   Names = input1['ID']
if Label_HI_4m == False:
   HI_4m = np.zeros(input1.size)
else:
   HI_4m = input1['HI_4m']
if Label_eHI_4m == False:
   eHI_4m = np.zeros(input1.size)
else:
   eHI_4m = input1['eHI_4m']
if Label_HI_7m == False:
   HI_7m = np.zeros(input1.size)
else:
   HI_7m = input1['HI_7m']
if Label_eHI_7m == False:
   eHI_7m = np.zeros(input1.size)
else:
   eHI_7m = input1['eHI_7m']
if Label_SIV == False:
   SIV_10m = np.zeros(input1.size)
else:
   SIV_10m = input1['SIV_10m']
if Label_eSIV == False:
   eSIV_10m = np.zeros(input1.size)
else:
   eSIV_10m = input1['eSIV_10m']
if Label_HI_12m == False:
   HI_12m = np.zeros(input1.size)
else:
   HI_12m = input1['HI_12m']
if Label_eHI_12m == False:
   eHI_12m = np.zeros(input1.size)
else:
   eHI_12m = input1['eHI_12m']
if Label_NeII == False:
   NeII_12m = np.zeros(input1.size)
else:
   NeII_12m = input1['NeII_12m']
if Label_eNeII == False:
   eNeII_12m = np.zeros(input1.size)
else:
   eNeII_12m = input1['eNeII_12m']
if Label_NeIII == False:
   NeIII_15m = np.zeros(input1.size)
else:
   NeIII_15m = input1['NeIII_15m']
if Label_eNeIII == False:
   eNeIII_15m = np.zeros(input1.size)
else:
   eNeIII_15m = input1['eNeIII_15m']
if Label_SIII_18m == False:
   SIII_18m = np.zeros(input1.size)
else:
   SIII_18m = input1['SIII_18m']
if Label_eSIII_18m == False:
   eSIII_18m = np.zeros(input1.size)
else:
   eSIII_18m = input1['eSIII_18m']
if Label_SIII_33m == False:
   SIII_33m = np.zeros(input1.size)
else:
   SIII_33m = input1['SIII_33m']
if Label_eSIII_33m == False:
   eSIII_33m = np.zeros(input1.size)
else:
   eSIII_33m = input1['eSIII_33m']
if Label_OIII_52m == False:
   OIII_52m = np.zeros(input1.size)
else:
   OIII_52m = input1['OIII_52m']
if Label_eOIII_52m == False:
   eOIII_52m = np.zeros(input1.size)
else:
   eOIII_52m = input1['eOIII_52m']
if Label_NIII == False:
   NIII_57m = np.zeros(input1.size)
else:
   NIII_57m = input1['NIII_57m']
if Label_eNIII == False:
   eNIII_57m = np.zeros(input1.size)
else:
   eNIII_57m = input1['eNIII_57m']
if Label_OIII_88m == False:
   OIII_88m = np.zeros(input1.size)
else:
   OIII_88m = input1['OIII_88m']
if Label_eOIII_88m == False:
   eOIII_88m = np.zeros(input1.size)
else:
   eOIII_88m = input1['eOIII_88m']
if Label_NII_122m == False:
   NII_122m = np.zeros(input1.size)
else:
   NII_122m = input1['NII_122m']
if Label_eNII_122m == False:
   eNII_122m = np.zeros(input1.size)
else:
   eNII_122m = input1['eNII_122m']
if Label_NII_205m == False:
   NII_205m = np.zeros(input1.size)
else:
   NII_205m = input1['NII_205m']
if Label_eNII_205m == False:
   eNII_205m = np.zeros(input1.size)
else:
   eNII_205m = input1['eNII_205m']

output = np.zeros(input1.size, dtype=[('ID', 'U12'), ('HI_4m', float),('eHI_4m', float),('HI_7m', float),('eHI_7m', float),('SIV_10m', float),('eSIV_10m', float),('HI_12m', float),('eHI_12m', float),('NeII_12m', float),('eNeII_12m', float),('NeIII_15m', float),('eNeIII_15m', float),('SIII_18m', float),('eSIII_18m', float),('SIII_33m', float),('eSIII_33m', float),('OIII_52m', float),('eOIII_52m', float),('NIII_57m', float),('eNIII_57m', float),('OIII_88m', float),('eOIII_88m', float),('NII_122m', float),('eNII_122m', float),('NII_205m', float),('eNII_205m', float),('grid', int),('OH', float),('eOH', float),('NO', float),('eNO', float),('logU', float),('elogU', float)] )

output['ID'] = Names
output['HI_4m'] = HI_4m
output['eHI_4m'] = eHI_4m
output['HI_7m'] = HI_7m
output['eHI_7m'] = eHI_7m
output['SIV_10m'] = SIV_10m
output['eSIV_10m'] = eSIV_10m
output['HI_12m'] = HI_12m
output['eHI_12m'] = eHI_12m
output['NeII_12m'] = NeII_12m
output['eNeII_12m'] = eNeII_12m
output['NeIII_15m'] = NeIII_15m
output['eNeIII_15m'] = eNeIII_15m
output['SIII_18m'] = SIII_18m
output['eSIII_18m'] = eSIII_18m
output['SIII_33m'] = SIII_33m
output['eSIII_33m'] = eSIII_33m
output['OIII_52m'] = OIII_52m
output['eOIII_52m'] = eOIII_52m
output['NIII_57m'] = NIII_57m
output['eNIII_57m'] = eNIII_57m
output['OIII_88m'] = OIII_88m
output['eOIII_88m'] = eOIII_88m
output['NII_122m'] = NII_122m
output['eNII_122m'] = eNII_122m
output['NII_205m'] = NII_205m
output['eNII_205m'] = eNII_205m


print ('Reading grids ....')
print ('')
print ('')
print ('----------------------------------------------------------------')
print ('(%)   ID    Grid  12+log(O/H)  log(N/O)    log(U)')
print ('-----------------------------------------------------------------')




# Beginning of loop of calculation

count = 0
for tab in range(0,len(input1),1):

   count = count + 1


   OH_mc = []
   NO_mc = []
   logU_mc = []
   OHe_mc = []
   NOe_mc = []
   logUe_mc = []  


# Selection of grid
   
   if NIII_57m[tab] > 0 and (OIII_52m[tab] > 0 or OIII_88m[tab] > 0):
      grid = grid2
      grid_type = 2
      grids.append(2)
   else:
      grid = grid3
      grid_type = 3
      grids.append(3)      

# Calculation of N/O

   if NIII_57m[tab] == 0 or (OIII_52m[tab] == 0 and OIII_88m[tab] == 0):
      NOff = -10
      eNOff = 0
   else:
      for monte in range(0,n,1):
         NO_p = 0
         den_NO = 0
         NO_e = 0
         den_NO_e = 0
         tol_max = 1e2

         OIII_52_obs = 0
         if OIII_52m[tab] > 0:
            while OIII_52_obs <= 0:
               OIII_52_obs = np.random.normal(OIII_52m[tab],eOIII_52m[tab]+1e-5)
         NIII_57_obs = 0
         if NIII_57m[tab] > 0:
            while NIII_57_obs <= 0:
               NIII_57_obs = np.random.normal(NIII_57m[tab],eNIII_57m[tab]+1e-5)
         OIII_88_obs = 0
         if OIII_88m[tab] > 0:
            while OIII_88_obs <= 0:
               OIII_88_obs = np.random.normal(OIII_88m[tab],eOIII_88m[tab]+1e-5)
         N3O3a_obs = -10
         if OIII_52_obs > 0:
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
            igrid = interpolate(igrid,1,NO-eNO-0.125,NO+eNO+0.125,10)


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
                  CHI_N3O3b = (np.log10(index[9]/index[10])- N3O3b_obs)**2/np.log10(index[9]/index[10]+1e-5)


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
            grid_mac.append(index[15])

         grid_c = np.reshape(grid_mac,(int(len(grid_mac)/16),16))



# Calculation of O/H and logU


   for monte in range(0,n,1):

      OH_p = 0
      logU_p = 0
      den_OH = 0
      OH_e = 0
      logU_e = 0
      den_OH_e = 0
      tol_max = 1e2

      HI_4m_obs = 0
      if HI_4m[tab] > 0:
         while HI_4m_obs <= 0:
            HI_4m_obs = np.random.normal(HI_4m[tab],eHI_4m[tab]+1e-5)
      HI_7m_obs = 0
      if HI_7m[tab] > 0:
         while HI_7m_obs <= 0:
            HI_7m_obs = np.random.normal(HI_7m[tab],eHI_7m[tab]+1e-5)
      SIV_10m_obs = 0
      if SIV_10m[tab] > 0:
         while SIV_10m_obs <= 0:
            SIV_10m_obs = np.random.normal(SIV_10m[tab],eSIV_10m[tab]+1e-5)
      HI_12m_obs = 0
      if HI_12m[tab] > 0:
         while HI_12m_obs <= 0:
            HI_12m_obs = np.random.normal(HI_12m[tab],eHI_12m[tab]+1e-5)
      NeII_12m_obs = 0
      if NeII_12m[tab] > 0:
         while NeII_12m_obs <= 0:
            NeII_12m_obs = np.random.normal(NeII_12m[tab],eNeII_12m[tab]+1e-3)
      NeIII_15m_obs = 0
      if NeIII_15m[tab] > 0:
         while NeIII_15m_obs <= 0:
            NeIII_15m_obs = np.random.normal(NeIII_15m[tab],eNeIII_15m[tab]+1e-3)
      SIII_18m_obs = 0
      if SIII_18m[tab] > 0:
         while SIII_18m_obs <= 0:
            SIII_18m_obs = np.random.normal(SIII_18m[tab],eSIII_18m[tab]+1e-3)
      SIII_33m_obs = 0
      if SIII_33m[tab] > 0:
         while SIII_33m_obs <= 0:
            SIII_33m_obs = np.random.normal(SIII_33m[tab],eSIII_33m[tab]+1e-3)
      OIII_52m_obs = 0
      if OIII_52m[tab] > 0:
         while OIII_52m_obs <= 0:
            OIII_52m_obs = np.random.normal(OIII_52m[tab],eOIII_52m[tab]+1e-5)
      NIII_57m_obs = 0
      if NIII_57m[tab] > 0:
         while NIII_57m_obs <= 0:
            NIII_57m_obs = np.random.normal(NIII_57m[tab],eNIII_57m[tab]+1e-5)
      OIII_88m_obs = 0
      if OIII_88m[tab] > 0:
         while OIII_88m_obs <= 0:
            OIII_88m_obs = np.random.normal(OIII_88m[tab],eOIII_88m[tab]+1e-5)
      NII_122m_obs = 0
      if NII_122m[tab] > 0:
         while NII_122m_obs <= 0:
            NII_122m_obs = np.random.normal(NII_122m[tab],eNII_122m[tab]+1e-3)
      NII_205m_obs = 0
      if NII_205m[tab] > 0:
         while NII_205m_obs <= 0:
            NII_205m_obs = np.random.normal(NII_205m[tab],eNII_205m[tab]+1e-3)
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
         N2N3a_obs = -10
      else:
         N2N3a_obs = np.log10((NII_122m_obs / NIII_57m_obs))
      if NII_122m_obs == 0 or OIII_52m_obs == 0:
         O3N2a_obs = -10
      else:
         O3N2a_obs = np.log10((OIII_52m_obs/NII_122m_obs))
      if NII_122m_obs == 0 or OIII_88m_obs == 0:
         O3N2b_obs = -10
      else:
         O3N2b_obs = np.log10((OIII_88m_obs/NII_122m_obs))
      if HI_4m_obs == 0 or NIII_57m_obs == 0 or NII_205m_obs == 0:
         N23d_obs = -10
      else:
         N23d_obs = np.log10((NII_205m_obs + NIII_57m_obs)/HI_4m_obs)
      if HI_7m_obs == 0 or NIII_57m_obs == 0 or NII_205m_obs == 0:
         N23e_obs = -10
      else:
         N23e_obs = np.log10((NII_205m_obs + NIII_57m_obs)/HI_7m_obs)
      if HI_12m_obs == 0 or NIII_57m_obs == 0 or NII_205m_obs == 0:
         N23f_obs = -10
      else:
         N23f_obs = np.log10((NII_205m_obs + NIII_57m_obs)/HI_12m_obs)
      if NII_205m_obs == 0 or NIII_57m_obs == 0:
         N2N3b_obs = -10
      else:
         N2N3b_obs = np.log10((NII_205m_obs / NIII_57m_obs))
      if NII_205m_obs == 0 or OIII_52m_obs == 0:
         O3N2c_obs = -10
      else:
         O3N2c_obs = np.log10((OIII_52m_obs/NII_205m_obs))
      if NII_205m_obs == 0 or OIII_88m_obs == 0:
         O3N2d_obs = -10
      else:
         O3N2d_obs = np.log10((OIII_88m_obs/NII_205m_obs))

      if  Ne23a_obs == -10 and Ne23b_obs == -10 and Ne23c_obs == -10 and S34a_obs == -10 and S34b_obs == -10 and S34c_obs == -10 and S34d_obs == -10 and S34e_obs == -10 and S34f_obs == -10 and Ne2Ne3_obs == -10 and S3S4a_obs == -10 and S3S4b_obs == -10 and N23a_obs == -10 and N23b_obs == -10 and N23c_obs == -10 and N2N3a_obs == -10 and O3N2a_obs == -10 and O3N2b_obs == -10 and N23d_obs == -10 and N23e_obs == -10 and N23f_obs == -10 and N2N3b_obs == -10 and O3N2c_obs == -10 and O3N2d_obs == -10: 
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
         CHI_N2N3a = 0
         CHI_O3N2a = 0
         CHI_O3N2b = 0
         CHI_N23d = 0
         CHI_N23e = 0
         CHI_N23f = 0
         CHI_N2N3b = 0
         CHI_O3N2c = 0
         CHI_O3N2d = 0
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
            if N2N3a_obs == -10: 
               CHI_N2N3a = 0
            elif index[12] == 0 or index[14] == 0: 
               CHI_N2N3a = tol_max
            else:   
               CHI_N2N3a = (np.log10(index[14]/index[12])- N2N3a_obs)**2/np.log10((index[14]/index[12]))
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
            if N23d_obs == -10: 
               CHI_N23d = 0
            elif index[3] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23d = tol_max
            else:   
               CHI_N23d = (np.log10((index[12]+index[15])/index[3])- N23d_obs)**2/np.log10((index[12]+index[15])/index[3])
            if N23e_obs == -10: 
               CHI_N23e = 0
            elif index[4] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23e = tol_max
            else:   
               CHI_N23e = (np.log10((index[12]+index[15])/index[4])- N23e_obs)**2/np.log10((index[12]+index[15])/index[4])
            if N23f_obs == -10: 
               CHI_N23f = 0
            elif index[6] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23f = tol_max
            else:   
               CHI_N23f = (np.log10((index[12]+index[15])/index[6])- N23f_obs)**2/np.log10((index[12]+index[15])/index[6])
            if N2N3b_obs == -10: 
               CHI_N2N3b = 0
            elif index[12] == 0 or index[15] == 0: 
               CHI_N2N3b = tol_max
            else:   
               CHI_N2N3b = (np.log10(index[15]/index[12])- N2N3b_obs)**2/np.log10((index[15]/index[12]))
            if O3N2c_obs == -10: 
               CHI_O3N2c = 0
            elif index[11] == 0 or index[15] == 0: 
               CHI_O3N2c = tol_max
            else:   
               CHI_O3N2c = (np.log10(index[11]/index[15])- O3N2c_obs)**2/np.log10((index[11]/index[15]))
            if O3N2d_obs == -10: 
               CHI_O3N2d = 0
            elif index[13] == 0 or index[15] == 0: 
               CHI_O3N2d = tol_max
            else:   
               CHI_O3N2d = (np.log10(index[13]/index[15])- O3N2d_obs)**2/np.log10((index[13]/index[15]))

            CHI_OH = (CHI_Ne23a**2 + CHI_Ne23b**2 + CHI_Ne23c**2 + CHI_Ne2Ne3**2 + CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S3S4a**2 + CHI_S3S4b**2+CHI_N23a**2+CHI_N23b**2+CHI_N23c**2+CHI_N2N3a**2+CHI_O3N2a**2+CHI_O3N2b**2+CHI_N23d**2+CHI_N23e**2+CHI_N23f**2+CHI_N2N3b**2+CHI_O3N2c**2+CHI_O3N2d**2)**0.5

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

      if  Ne23a_obs == -10 and Ne23b_obs == -10 and Ne23c_obs == -10 and S34a_obs == -10 and S34b_obs == -10 and S34c_obs == -10 and S34d_obs == -10 and S34e_obs == -10 and S34f_obs == -10 and Ne2Ne3_obs == -10 and S3S4a_obs == -10 and S3S4b_obs == -10 and N23a_obs == -10 and N23b_obs == -10 and N23c_obs == -10 and N2N3a_obs == -10 and O3N2a_obs == -10 and O3N2b_obs == -10 and N23d_obs == -10 and N23e_obs == -10 and N23f_obs == -10 and N2N3b_obs == -10 and O3N2b_obs == -10 and O3N2d_obs == -10: 
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
         CHI_N2N3a = 0
         CHI_O3N2a = 0
         CHI_O3N2b = 0
         CHI_N23d = 0
         CHI_N23e = 0
         CHI_N23f = 0
         CHI_N2N3b = 0
         CHI_O3N2c = 0
         CHI_O3N2d = 0
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
            if N2N3a_obs == -10: 
               CHI_N2N3a = 0
            elif index[12] == 0 or index[14] == 0: 
               CHI_N2N3a = tol_max
            else:   
               CHI_N2N3a = (np.log10(index[14]/index[12])- N2N3a_obs)**2/np.log10((index[14]/index[12]))
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
            if N23d_obs == -10: 
               CHI_N23d = 0
            elif index[3] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23d = tol_max
            else:   
               CHI_N23d = (np.log10((index[12]+index[15])/index[3])- N23d_obs)**2/np.log10((index[12]+index[15])/index[3])
            if N23e_obs == -10: 
               CHI_N23e = 0
            elif index[4] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23e = tol_max
            else:   
               CHI_N23e = (np.log10((index[12]+index[15])/index[4])- N23e_obs)**2/np.log10((index[12]+index[15])/index[4])
            if N23f_obs == -10: 
               CHI_N23f = 0
            elif index[6] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23f = tol_max
            else:   
               CHI_N23f = (np.log10((index[12]+index[15])/index[6])- N23f_obs)**2/np.log10((index[12]+index[15])/index[6])
            if N2N3b_obs == -10: 
               CHI_N2N3b = 0
            elif index[12] == 0 or index[15] == 0: 
               CHI_N2N3b = tol_max
            else:   
               CHI_N2N3b = (np.log10(index[15]/index[12])- N2N3b_obs)**2/np.log10((index[15]/index[12]))
            if O3N2c_obs == -10: 
               CHI_O3N2c = 0
            elif index[11] == 0 or index[15] == 0: 
               CHI_O3N2c = tol_max
            else:   
               CHI_O3N2c = (np.log10(index[11]/index[15])- O3N2c_obs)**2/np.log10((index[11]/index[15]))
            if O3N2d_obs == -10: 
               CHI_O3N2d = 0
            elif index[13] == 0 or index[15] == 0: 
               CHI_O3N2d = tol_max
            else:   
               CHI_O3N2d = (np.log10(index[13]/index[15])- O3N2d_obs)**2/np.log10((index[13]/index[15]))

            CHI_OH = (CHI_Ne23a**2 + CHI_Ne23b**2 + CHI_Ne23c**2 + CHI_Ne2Ne3**2 + CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S3S4a**2 + CHI_S3S4b**2+CHI_N23a**2+CHI_N23b**2+CHI_N23c**2+CHI_N2N3a**2+CHI_O3N2a**2+CHI_O3N2b**2+CHI_N23d**2+CHI_N23e**2+CHI_N23f**2+CHI_N2N3b**2+CHI_O3N2c**2+CHI_O3N2d**2)**0.5


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
         igrid = interpolate(grid_c,2,logU-elogU-0.25,logU+elogU+0.25,10)
         igrid = igrid[np.lexsort((igrid[:,1],igrid[:,2]))]
         igrid = interpolate(igrid,0,OH-eOH-0.1,OH+eOH+0.1,10)
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
         CHI_N2N3a = 0
         CHI_O3N2a = 0
         CHI_O3N2b = 0
         CHI_N23d = 0
         CHI_N23e = 0
         CHI_N23f = 0
         CHI_N2N3b = 0
         CHI_O3N2c = 0
         CHI_O3N2d = 0
         CHI_OH = 0
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
            if N2N3a_obs == -10: 
               CHI_N2N3a = 0
            elif index[12] == 0 or index[14] == 0: 
               CHI_N2N3a = tol_max
            else:   
               CHI_N2N3a = (np.log10(index[14]/index[12])- N2N3a_obs)**2/np.log10((index[14]/index[12]))
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
            if N23d_obs == -10: 
               CHI_N23d = 0
            elif index[3] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23d = tol_max
            else:   
               CHI_N23d = (np.log10((index[12]+index[15])/index[3])- N23d_obs)**2/np.log10((index[12]+index[15])/index[3])
            if N23e_obs == -10: 
               CHI_N23e = 0
            elif index[4] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23e = tol_max
            else:   
               CHI_N23e = (np.log10((index[12]+index[15])/index[4])- N23e_obs)**2/np.log10((index[12]+index[15])/index[4])
            if N23f_obs == -10: 
               CHI_N23f = 0
            elif index[6] == 0 or index[12] == 0 or index[15] == 0: 
               CHI_N23f = tol_max
            else:   
               CHI_N23f = (np.log10((index[12]+index[15])/index[6])- N23f_obs)**2/np.log10((index[12]+index[15])/index[6])
            if N2N3b_obs == -10: 
               CHI_N2N3b = 0
            elif index[12] == 0 or index[15] == 0: 
               CHI_N2N3b = tol_max
            else:   
               CHI_N2N3b = (np.log10(index[15]/index[12])- N2N3b_obs)**2/np.log10((index[15]/index[12]))
            if O3N2c_obs == -10: 
               CHI_O3N2c = 0
            elif index[11] == 0 or index[15] == 0: 
               CHI_O3N2c = tol_max
            else:   
               CHI_O3N2c = (np.log10(index[11]/index[15])- O3N2c_obs)**2/np.log10((index[11]/index[15]))
            if O3N2d_obs == -10: 
               CHI_O3N2d = 0
            elif index[13] == 0 or index[15] == 0: 
               CHI_O3N2d = tol_max
            else:   
               CHI_O3N2d = (np.log10(index[13]/index[15])- O3N2d_obs)**2/np.log10((index[13]/index[15]))

            CHI_OH = (CHI_Ne23a**2 + CHI_Ne23b**2 + CHI_Ne23c**2 + CHI_Ne2Ne3**2 + CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S3S4a**2 + CHI_S3S4b**2+CHI_N23a**2+CHI_N23b**2+CHI_N23c**2+CHI_N2N3a**2+CHI_O3N2a**2+CHI_O3N2b**2+CHI_N23d**2+CHI_N23e**2+CHI_N23f**2+CHI_N2N3b**2+CHI_O3N2c**2+CHI_O3N2d**2)**0.5



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

   




   
   OHffs.append(OHff)
   eOHffs.append(eOHff)
   NOffs.append(NOff)
   eNOffs.append(eNOff)
   logUffs.append(logUff)
   elogUffs.append(elogUff)
         

   if input0.size == 1 and tab==0: continue
   print (round(100*(count)/float(len(input1)),1),'%',Names[tab],grid_type,'', round(OHff,2), round(eOHff,2),'',round(NOff,2), round(eNOff,2), '',round(logUff,2), round(elogUff,2))

output['grid'] = grids
output['OH'] = OHffs
output['eOH'] = eOHffs
output['NO'] = NOffs
output['eNO'] = eNOffs
output['logU'] = logUffs
output['elogU'] = elogUffs





lineas_header = [' HII-CHI-mistry-IR v.2.1 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','Bra  eBra   Pfa   ePfa  S4     eS4   Hua   eHua Ne2-12  eNe2-12  Ne3-15  eNe3-15  S3-18  eS3-18  S3-33  eS3-33  O3-52 eO3-52  N3-57  eN3-57  O3-88  eO3-88   N2-122  eN2-122 N2-205  eN2-205  i O/H     eO/H  N/O    eN/O  logU   elogU']

header = '\n'.join(lineas_header)

np.savetxt(input00+'_hcm-ir-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*26+['%i']+['%.2f']*6),header=header)
print ('________________________________')
print ('Results are stored in ' + input00 + '_hcm-ir-output.dat')
                           
      

      
