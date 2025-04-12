# Filename: HCm_v 5.0.py



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
print ('This is HII-CHI-mistry v. 5.0')
print (' See Perez-Montero, E. (2014) for details')
print  (' Insert the name of your input text file with all or some of the following columns:')
print (' 3727 [OII], 3868 [NeIII], 4363 [OIII], 5007 [OIII], 6584 [NII], 6725 [SII]')
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
ver_np = np.fromstring(np.__version__,sep=' ')
try:
   if ver_np[0] < 1.13 or ver_np[0] > 1.3:
      input0 = np.genfromtxt(input00,dtype=None,names=True)
   else:
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
   print ('(2) AGN, double component, a(OX) = -0.8, a(UV) = -1.0')
   print ('(3) AGN, double component, a(OX) = -1.2, a(UV) = -1.0')
   print('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      sed = raw_input('Choose SED of the models:')
   else:
      sed = input('Choose SED of the models:')
   if sed == '1' or sed == '2' or sed == '3' : question = False 
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
   grid1 = np.loadtxt('C17_popstar_v5.0.dat')
   grid2 = np.loadtxt('C17_popstar_logU_adapted_emp_v5.0.dat')
   grid3 = np.loadtxt('C17_popstar_logU-NO_adapted_emp_v5.0.dat')
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
      res_NO = 0.0125
elif sed==2:
   grid1 = np.loadtxt('C17_agn_a08_v5.0.dat')
   grid2 = np.loadtxt('C17_agn_a08_v5.0.dat')
   grid3 = np.loadtxt('C17_agn_a08_NO_adapted_emp_v5.0.dat')
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -0.8. No interpolation'
      print ('No interpolation for the AGN a(ox) = -0.8 models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      print ('')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -0.8. Interpolation'
      print ('Interpolation for the AGN a(ox) = -0.8 models is going to be used.')
      print ('The grid has a resolution of 0.01dex for O/H and 0.0125 dex for N/O')
      print ('')
      res_NO = 0.0125
elif sed==3:
   grid1 = np.loadtxt('C17_agn_a12_v5.0.dat')
   grid2 = np.loadtxt('C17_agn_a12_v5.0.dat')
   grid3 = np.loadtxt('C17_agn_a12_NO_adapted_emp_v5.0.dat')
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.2. No interpolation'
      print ('No interpolation for the AGN a(ox) = -1.2 models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      print ('')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.2. Interpolation'
      print ('Interpolation for the AGN a(ox) = -1.2 models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
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
Label_OII = False
Label_eOII = False
Label_NeIII = False
Label_eNeIII = False
Label_OIIIa = False
Label_eOIIIa = False
Label_OIII_4959 = False
Label_eOIII_4959 = False
Label_OIII_5007 = False
Label_eOIII_5007 = False
Label_NII = False
Label_eNII = False
Label_SII = False
Label_eSII = False
Label_SII_6716 = False
Label_eSII_6716 = False
Label_SII_6731 = False
Label_eSII_6731 = False

for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == 'OII_3727':
      Label_OII = True
   if input1.dtype.names[col] == 'eOII_3727':
      Label_eOII = True
   if input1.dtype.names[col] == 'NeIII_3868':
      Label_NeIII = True
   if input1.dtype.names[col] == 'eNeIII_3868':
      Label_eNeIII = True
   if input1.dtype.names[col] == 'OIII_4363':
      Label_OIIIa = True
   if input1.dtype.names[col] == 'eOIIIa_4363':
      Label_eOIIIa = True
   if input1.dtype.names[col] == 'OIII_4959':
      Label_OIII_4959 = True
   if input1.dtype.names[col] == 'eOIII_4959':
      Label_eOIII_4959 = True
   if input1.dtype.names[col] == 'OIII_5007':
      Label_OIII_5007 = True
   if input1.dtype.names[col] == 'eOIII_5007':
      Label_eOIII_5007 = True
   if input1.dtype.names[col] == 'NII_6584':
      Label_NII = True
   if input1.dtype.names[col] == 'eNII_6584':
      Label_eNII = True
   if input1.dtype.names[col] == 'SII_6725':
      Label_SII = True
   if input1.dtype.names[col] == 'eSII_6725':
      Label_eSII = True
   if input1.dtype.names[col] == 'SII_6716':
      Label_SII_6716 = True
   if input1.dtype.names[col] == 'SII_6731':
      Label_SII_6731 = True
   if input1.dtype.names[col] == 'eSII_6716':
      Label_eSII_6716 = True
   if input1.dtype.names[col] == 'eSII_6731':
      Label_eSII_6731 = True

if Label_ID == False:
   Names = np.arange(1,input1.size+1,1)
else:
   Names = input1['ID']
if Label_OII == False:
   OII_3727 = np.zeros(input1.size)
else:
   OII_3727 = input1['OII_3727']
if Label_eOII == False:
   eOII_3727 = np.zeros(input1.size)
else:
   eOII_3727 = input1['eOII_3727']
if Label_NeIII == False:
   NeIII_3868 = np.zeros(input1.size)
else:
   NeIII_3868 = input1['NeIII_3868']
if Label_eNeIII == False:
   eNeIII_3868 = np.zeros(input1.size)
else:
   eNeIII_3868 = input1['eNeIII_3868']
if Label_OIIIa == False:
   OIII_4363 = np.zeros(input1.size)
else:
   OIII_4363 = input1['OIII_4363']
if Label_eOIIIa == False:
   eOIII_4363 = np.zeros(input1.size)
else:
   eOIII_4363 = input1['eOIII_4363']
if Label_OIII_4959 == False and Label_OIII_5007 == False:
   OIII_5007 = np.zeros(input1.size)
elif Label_OIII_4959 == False and Label_OIII_5007 == True:
   OIII_5007 = input1['OIII_5007']
elif Label_OIII_4959 == True and Label_OIII_5007 == False:
   OIII_5007 = 3*input1['OIII_4959']
else:
   OIII_5007 = (input1['OIII_5007']+input1['OIII_4959'])/1.33
if Label_eOIII_4959 == False and Label_eOIII_5007 == False:
   eOIII_5007 = np.zeros(input1.size)
elif Label_eOIII_4959 == False and Label_eOIII_5007 == True:
   eOIII_5007 = input1['eOIII_5007']
elif Label_eOIII_4959 == True and Label_eOIII_5007 == False:
   eOIII_5007 = 4*input1['eOIII_4959']
else:
   eOIII_5007 = (input1['eOIII_5007']+input1['eOIII_4959'])/1.33
if Label_NII == False:
   NII_6584 = np.zeros(input1.size)
else:
   NII_6584 = input1['NII_6584']
if Label_eNII == False:
   eNII_6584 = np.zeros(input1.size)
else:
   eNII_6584 = input1['eNII_6584']
if Label_SII == False and (Label_SII_6716 == False or Label_SII_6731 == False):
   SII_6725 = np.zeros(input1.size)
elif Label_SII == True:
   SII_6725 = input1['SII_6725']
else:
   SII_6725 = input1['SII_6716']+input1['SII_6731']
if Label_eSII == False and (Label_eSII_6716 == False or Label_eSII_6731 == False):
   eSII_6725 = np.zeros(input1.size)
elif Label_eSII == True:
   eSII_6725 = input1['eSII_6725']
else:
   eSII_6725 = input1['eSII_6716']+input1['eSII_6731']

output = np.zeros(input1.size, dtype=[('ID', 'U12'), ('OII_3727', float),('eOII_3727', float),('NeIII_3868', float),('eNeIII_3868', float),('OIII_4363', float),('eOIII_4363', float),('OIII_5007', float),('eOIII_5007', float),('NII_6584', float),('eNII_6584', float),('SII_6725', float),('eSII_6725', float),('grid', int),('OH', float),('eOH', float),('NO', float),('eNO', float),('logU', float),('elogU', float)] )

output['ID'] = Names
output['OII_3727'] = OII_3727
output['eOII_3727'] = eOII_3727
output['NeIII_3868'] = NeIII_3868
output['eNeIII_3868'] = eNeIII_3868
output['OIII_4363'] = OIII_4363
output['eOIII_4363'] = eOIII_4363
output['OIII_5007'] = OIII_5007
output['eOIII_5007'] = eOIII_5007
output['NII_6584'] = NII_6584
output['eNII_6584'] = eNII_6584
output['SII_6725'] = SII_6725
output['eSII_6725'] = eSII_6725




print ('Reading grids ....')
print ('')
print ('')
print ('----------------------------------------------------------------')
print ('(%)   ID    Grid  12+log(O/H)  log(N/O)    log(U)')
print ('-----------------------------------------------------------------')


# Beginning of loop of calculation

count = 0
for tab in range(0,input1.size,1):
   count = count + 1


   OH_mc = []
   NO_mc = []
   logU_mc = []
   OHe_mc = []
   NOe_mc = []
   logUe_mc = []  


# Selection of grid
   
   if OIII_4363[tab] > 0 and OIII_5007[tab] > 0:
      grid = grid1
      grid_type = 1
      grids.append(1)
   elif NII_6584[tab] > 0 and (OII_3727[tab] > 0 or SII_6725[tab] > 0):
      grid = grid2
      grid_type = 2
      grids.append(2)
   else:
      grid = grid3
      grid_type = 3
      grids.append(3)      

# Calculation of N/O

   if NII_6584[tab] == 0 or (OII_3727[tab] == 0 and SII_6725[tab] == 0):
      NOff = -10
      eNOff = 0
   else:
      for monte in range(0,n,1):
         NO_p = 0
         den_NO = 0
         NO_e = 0
         den_NO_e = 0
         tol_max = 1e2

         OII_3727_obs = 0
         if OII_3727[tab] == 0:
            OII_3727_obs = 0
         else:
            while OII_3727_obs <= 0:
               OII_3727_obs = np.random.normal(OII_3727[tab],eOII_3727[tab]+1e-3)
         OIII_4363_obs = 0
         if OIII_4363[tab] == 0:
            OIII_4363_obs = 0
         else:
            while OIII_4363_obs <= 0:
               OIII_4363_obs = np.random.normal(OIII_4363[tab],eOIII_4363[tab]+1e-3)
         OIII_5007_obs = 0
         if OIII_5007[tab] == 0:
            OIII_5007_obs = 0
         else:
            while OIII_5007_obs <= 0:
               OIII_5007_obs = np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-3)
         if OIII_4363_obs == 0  or OIII_5007_obs == 0:
            ROIII_obs = 0
         else:
            ROIII_obs = OIII_5007_obs / OIII_4363_obs
         NII_6584_obs = 0
         if NII_6584[tab] == 0:
            NII_6584_obs = 0
         else:
            while NII_6584_obs <= 0:
               NII_6584_obs = np.random.normal(NII_6584[tab],eNII_6584[tab]+1e-3)
         SII_6725_obs = 0
         if SII_6725[tab] == 0:
               SII_6725_obs = 0
         else:
            while SII_6725_obs <= 0:
               SII_6725_obs = np.random.normal(SII_6725[tab],eSII_6725[tab]+1e-3)
            if SII_6725_obs <= 0: SII_6725_obs = 0
         if NII_6584_obs == 0 or OII_3727_obs == 0:
            N2O2_obs = -10
         else:
            N2O2_obs = np.log10(NII_6584_obs / OII_3727_obs)
         if NII_6584_obs == 0 or SII_6725_obs == 0:
            N2S2_obs = -10
         else:   
            N2S2_obs = np.log10(NII_6584_obs / SII_6725_obs)

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
               CHI_N2O2 =(np.log10(index[7]/index[3]) - N2O2_obs)**2/(abs(np.log10(index[7]/index[3])+1e-3))
            if N2S2_obs == -10: 
               CHI_N2S2 = 0
            elif index[7] == 0 or index[8] == 0:
               CHI_N2S2 = tol_max
            else:
               CHI_N2S2 =(np.log10(index[7]/index[8]) - N2S2_obs)**2/(abs(np.log10(index[7]/index[8])+1e-3))
            CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2)**0.5


            if CHI_NO == 0:
               NO_p = NO_p
               den_NO = den_NO
            else:
               NO_p = index[1] / (CHI_NO) + NO_p
               den_NO = 1 / (CHI_NO) + den_NO


         NO = NO_p / den_NO 


# Calculation of N/O error

   
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
               CHI_N2O2 =(np.log10(index[7]/index[3]) - N2O2_obs)**2/(abs(np.log10(index[7]/index[3])+1e-3))
            if N2S2_obs == -10: 
               CHI_N2S2 = 0
            elif index[7] == 0 or index[8] == 0:
               CHI_N2S2 = tol_max
            else:
               CHI_N2S2 =(np.log10(index[7]/index[8]) - N2S2_obs)**2/(abs(np.log10(index[7]/index[8])+1e-3))

            CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2)**0.5


            if CHI_NO == 0:
               NO_e = NO_e
               den_NO_e = den_NO_e  
            else:
               NO_e = (index[1] - NO)**2 / (CHI_NO) + NO_e
               den_NO_e = 1 / (CHI_NO) + den_NO_e  


         eNO = NO_e / den_NO_e 


#Iterations for the interpolation mode

         if inter == 0 or NO == -10:
            NOf = NO
         elif inter == 1:
            igrid = grid[np.lexsort((grid[:,0],grid[:,2]))]
            igrid = interpolate(igrid,1,NO-eNO-0.125,NO+eNO+0.125,10)

            CHI_ROIII = 0
            CHI_N2O2 = 0
            CHI_N2S2 = 0
            CHI_NO = 0
            NO_p = 0
            den_NO = 0

            for index in igrid:
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
                  CHI_N2O2 =(np.log10(index[7]/index[3]) - N2O2_obs)**2/(abs(np.log10(index[7]/index[3])+1e-3))
               if N2S2_obs == -10: 
                  CHI_N2S2 = 0
               elif index[7] == 0 or index[8] == 0:
                  CHI_N2S2 = tol_max
               else:
                  CHI_N2S2 =(np.log10(index[7]/index[8]) - N2S2_obs)**2/(abs(np.log10(index[7]/index[8])+1e-3))


               CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2)**0.5
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
         grid_c = np.reshape(grid_mac,(int(len(grid_mac)/9),9))

# Calculation of O/H and logU


   for monte in range(0,n,1):

      OH_p = 0
      logU_p = 0
      den_OH = 0
      OH_e = 0
      logU_e = 0
      den_OH_e = 0
      tol_max = 1e2

      OII_3727_obs = 0
      if OII_3727[tab] == 0:
         OII_3727_obs = 0
      else:
         while OII_3727_obs <= 0:
            OII_3727_obs = np.random.normal(OII_3727[tab],eOII_3727[tab]+1e-3)
      NeIII_3868_obs = 0
      if NeIII_3868[tab] == 0:
         NeIII_3868_obs = 0
      else:
         while NeIII_3868_obs <= 0:
            NeIII_3868_obs = np.random.normal(NeIII_3868[tab],eNeIII_3868[tab]+1e-3)
      OIII_4363_obs = 0
      if OIII_4363[tab] == 0:
         OIII_4363_obs = 0
      else:
         while OIII_4363_obs <= 0:
            OIII_4363_obs = np.random.normal(OIII_4363[tab],eOIII_4363[tab]+1e-3)
      OIII_5007_obs = 0
      if OIII_5007[tab] == 0:
         OIII_5007_obs = 0
      else:
         while OIII_5007_obs <= 0:
            OIII_5007_obs = np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-3)
      if OIII_4363_obs == 0  or OIII_5007_obs == 0:
         ROIII_obs = 0
      else:
         ROIII_obs = OIII_5007_obs / OIII_4363_obs
      NII_6584_obs = 0
      if NII_6584[tab] == 0:
         NII_6584_obs = 0
      else:
         while NII_6584_obs <= 0:
            NII_6584_obs = np.random.normal(NII_6584[tab],eNII_6584[tab]+1e-3)
      SII_6725_obs = 0
      if SII_6725[tab] == 0:
            SII_6725_obs = 0
      else:
         while SII_6725_obs <= 0:
            SII_6725_obs = np.random.normal(SII_6725[tab],eSII_6725[tab]+1e-3)
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
      if OIII_5007_obs == 0 or SII_6725_obs == 0:
         O3S2_obs = -10
      else:
         O3S2_obs = np.log10( OIII_5007_obs / SII_6725_obs )

      if R23_obs == -10 and NII_6584_obs == 0 and ROIII_obs == 0 and R2Ne3_obs == -10 and O3S2_obs == -10:
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
         CHI_O3S2 = 0
         CHI_OH = 0
         for index in grid_c:
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
               CHI_R23 = (np.log10(index[3]+index[6])-R23_obs)**2/   (np.abs(np.log10(index[3]+index[6]+1e-3)))
            if OII_3727_obs == 0 or NeIII_3868_obs == 0:
               CHI_O2Ne3 = 0
               CHI_R2Ne3 = 0
            elif index[3] == 0 or index[4] == 0:
               CHI_O2Ne3 = tol_max
               CHI_R2Ne3 = tol_max
            else:
               CHI_O2Ne3 = (index[3]/index[4] - O2Ne3_obs)**2/(index[3]/index[4])
               CHI_R2Ne3 = (np.log10(index[3]+index[4])-R2Ne3_obs)**2/   (np.abs(np.log10(index[3]+index[4]+1e-3)))
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif index[6] == 0 or index[7] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index[6]/index[7]) - O3N2_obs)**2/(np.abs(np.log10(index[6]/index[7]+1e-3)))
            if OIII_5007_obs == 0 or SII_6725_obs == 0:
               CHI_O3S2 = 0
            elif index[6] == 0 or index[8] == 0:
               CHI_O3S2 = tol_max
            else:
               CHI_O3S2 = (np.log10(index[6]/index[8]) - O3S2_obs)**2/(np.abs(np.log10(index[6]/index[8]+1e-3)))


            if ROIII_obs > 0:
               CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROIII_obs == 0 and NII_6584_obs > 0:
               CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O3N2**2 + CHI_O3S2**2 )**0.5 
            elif ROIII_obs == 0 and NII_6584_obs == 0 and OIII_5007_obs > 0:
               CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O3S2**2)**0.5
            elif ROIII_obs == 0 and OIII_5007_obs == 0:
               CHI_OH = (CHI_O2Ne3**2 + CHI_R2Ne3**2 )**0.5

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

      if R23_obs == -10 and NII_6584_obs == 0 and ROIII_obs == 0 and R2Ne3_obs == -10 and O3S2_obs == -10:
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
         CHI_O3S2 = 0
         CHI_OH = 0
         for index in grid_c:
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
               CHI_R23 = (np.log10(index[3]+index[6])-R23_obs)**2/   (np.abs(np.log10(index[3]+index[6]+1e-3)))
            if OII_3727_obs == 0 or NeIII_3868_obs == 0:
               CHI_O2Ne3 = 0
               CHI_R2Ne3 = 0
            elif index[3] == 0 or index[4] == 0:
               CHI_O2Ne3 = tol_max
               CHI_R2Ne3 = tol_max
            else:
               CHI_O2Ne3 = (index[3]/index[4] - O2Ne3_obs)**2/(index[3]/index[4])
               CHI_R2Ne3 = (np.log10(index[3]+index[4])-R2Ne3_obs)**2/   (np.abs(np.log10(index[3]+index[4]+1e-3)))
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif index[6] == 0 or index[7] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index[6]/index[7]) - O3N2_obs)**2/(np.abs(np.log10(index[6]/index[7]+1e-3)))
            if OIII_5007_obs == 0 or SII_6725_obs == 0:
               CHI_O3S2 = 0
            elif index[6] == 0 or index[8] == 0:
               CHI_O3S2 = tol_max
            else:
               CHI_O3S2 = (np.log10(index[6]/index[8]) - O3S2_obs)**2/(np.abs(np.log10(index[6]/index[8]+1e-3)))


            if ROIII_obs > 0:
               CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2  )**0.5
            elif ROIII_obs == 0 and NII_6584_obs > 0:
               CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O3N2**2 + CHI_O3S2**2)**0.5 
            elif ROIII_obs == 0 and NII_6584_obs == 0 and OIII_5007_obs > 0:
               CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O3S2**2)**0.5
            else:
               CHI_OH = (CHI_O2Ne3**2 + CHI_R2Ne3**2 )**0.5

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

         CHI_ROIII = 0
         CHI_OIII = 0
         CHI_OII = 0
         CHI_NII = 0
         CHI_O2O3 = 0
         CHI_R23 = 0
         CHI_O3N2 = 0
         CHI_O2Ne3 = 0
         CHI_R2Ne3 = 0
         CHI_O3S2 = 0
         CHI_OH = 0
         OH_p = 0
         logU_p = 0
         den_OH = 0
      

         for index in igrid:
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
               CHI_R23 = (np.log10(index[3]+index[6])-R23_obs)**2/(np.abs(np.log10(index[3]+index[6]+1e-3)))
            if OII_3727_obs == 0 or NeIII_3868_obs == 0:
               CHI_O2Ne3 = 0
               CHI_R2Ne3 = 0
            elif index[3] == 0 or index[4] == 0:
               CHI_O2Ne3 = tol_max
               CHI_R2Ne3 = tol_max
            else:
               CHI_O2Ne3 = (index[3]/index[4] - O2Ne3_obs)**2/(index[3]/index[4])
               CHI_R2Ne3 = (np.log10(index[3]+index[4])-R2Ne3_obs)**2/(np.abs(np.log10(index[3]+index[4]+1e-3)))
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif index[6] == 0 or index[7] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index[6]/index[7]) - O3N2_obs)**2/(np.abs(np.log10(index[6]/index[7]+1e-3)))
            if OIII_5007_obs == 0 or SII_6725_obs == 0:
               CHI_O3S2 = 0
            elif index[6] == 0 or index[8] == 0:
               CHI_O3S2 = tol_max
            else:
               CHI_O3S2 = (np.log10(index[6]/index[8]) - O3S2_obs)**2/(np.abs(np.log10(index[6]/index[8]+1e-3)))




            if ROIII_obs > 0:
               CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2)**0.5
            elif NII_6584_obs > 0 and OII_3727_obs > 0:
               CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O2Ne3**2 + CHI_R2Ne3**2)**0.5 
            elif NII_6584_obs > 0 and OII_3727_obs == 0:
               CHI_OH = (CHI_NII**2 + CHI_O3N2**2 + CHI_O3S2**2)**0.5
            elif NII_6584_obs == 0:
               CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O2Ne3**2 + CHI_R2Ne3**2 + CHI_O3S2**2 )**0.5

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
   elogUff = (np.std(logU_mc)**2+np.std(logUe_mc)**2)**0.5




   
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

if input0.size == 1:  output = np.delete(output,obj=1,axis=0)


lineas_header = [' HII-CHI-mistry v.5.0 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','ID   O2Hb eO2Hb  Ne3Hb  eNeHb O3aHb  eO3aHb O3nHb  eO3nHb N2Hb   eN2Hb  S2Hb   eS2Hb  i O/H     eO/H  N/O    eN/O  logU   elogU']

header = '\n'.join(lineas_header)


np.savetxt(input00+'_hcm-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*12+['%i']+['%.2f']*6),header=header)
print ('________________________________')
print ('Results are stored in ' + input00 + '_hcm-output.dat')


