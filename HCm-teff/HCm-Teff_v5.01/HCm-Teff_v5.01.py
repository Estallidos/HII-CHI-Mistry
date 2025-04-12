# Filename: HCm-Teff_v5.0.py

import string
import numpy as np
import sys
#sys.stderr = open('errorlog.txt', 'w')


#Function for interpolation of grids

def interpolate(grid,z,zmin,zmax,n):
   ncol = 10
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
print ('This is HII-CHI-mistry-Teff v. 5.0')
print (' See Perez-Montero et al (2019) for details')
print ( ' Insert the name of your input text file with all or some of the following columns:')
print ('12+log(O/H), 3727 [OII], 5007 [OIII], 6725 [SII], 9069 [SIII], 4471 HeI, 5876 HeI, and 4686 HeII')
print ('with their corresponding labels and errors in adjacent columns')
print ('---------------------------------------------------------------------')
print ('')



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
   sys.exit()
print ('')

if input0.size == 1:
   input1 = np.stack((input0,input0))
else:
   input1 = input0


# Iterations for Montecarlo error derivation
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
   print ('---------------------------------------------------------------------')
   print ('(1) Effective temperature, WM-Basic stellar atmospheres (30-60 kK)')
   print ('(2) Effective temperature, Black body (30-90 kK)')
   print ('(3) Photon escape fraction, BPASS cluster atmospheres, age = 4 Myr, Mup = 300, x = 1.35, w/binaries')
   print ('---------------------------------------------------------------------')

   if int(sys.version[0]) < 3:
      sed = raw_input('Choose parameter and models:')
   else:
      sed = input('Choose parameter and models:')
   if sed == '1' or sed == '2' or sed == '3': question = False

print ('')
question = True
while question:
   print ('---------------------------------------------------------------------')
   if sed == '3': 
      geo = '2'
      question = False
   else:
      print ('(1) Plane-parallel geometry')
      print ('(2) Spherical geometry')
      print ('---------------------------------------------------------------------')
      if int(sys.version[0]) < 3:
         geo = raw_input('Choose geometry of the models:')
      else:
         geo = input('Choose geometry of the models:')
      if geo == '1' or geo == '2': question = False


if geo == '3':
   inter = '0'
else:
   question = True
   while question:
      if int(sys.version[0]) < 3:
         inter = raw_input('Choose models [0] No interpolated [1] Interpolated: ')
      else:
         inter = input('Choose models [0] No interpolated [1] Interpolated: ')
      if inter == '0' or inter == '1': question = False
   print ('')


sed = int(sed)
geo = int(geo)
inter = int(inter)


if geo==1 and sed == 1:
   bin = 99
   grid = np.loadtxt('C17_WMb_Teff_30-60_pp.dat')
   if inter == 0:
      sed_type = 'WM-Basic stellar atmosphere. Plane-parallel geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic models with plane-paralell geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic stellar atmosphere. Plane-parallel geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic models with plane-paralell geometry and interpolation')
elif geo==2 and sed == 1:
   bin = 99
   grid = np.loadtxt('C17_WMb_Teff_30-60_sph.dat')
   if inter == 0:
      sed_type = 'WM-Basic stellar atmosphereSpherical geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic stellar atmosphere. Spherical geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic models with spherical geometry and interpolation')

elif geo==1 and sed == 2:
   bin = 132
   grid = np.loadtxt('C17_bb_Teff_30-90_pp.dat')
   if inter == 0:
      sed_type = 'Black body. Plane-parallel geometry. Not interpolated'
      print ('Teff and U calculation using black body models with plane-parallel geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'Black body. Plane-parallel geometry. Interpolated'
      print ('Teff and U calculation using black body models with plane-parallel geometry and interpolation')
elif geo==2 and sed == 2:
   bin = 132
   grid = np.loadtxt('C17_bb_Teff_30-90_sph.dat')
   if inter == 0:
      sed_type = 'Black body. spherical geometry. Not interpolated'
      print ('Teff and U calculation using black body models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'Black body. spherical geometry. Interpolated'
      print ('Teff and U calculation using black body models with spehrical geometry and interpolation')

elif geo==1 and sed == 3:
   sed_type = 'BPASS cluster atmospheres, Mup = 300. Plane-parallel geometry'
   bin = 63
   grid = np.loadtxt('C17_bpass21_imf135_300_pp_esc.dat')
   print ('f_abs and U calculation using BPASS models with plane-parallel geometry')
elif geo==2 and sed == 3:
   bin = 77
   grid = np.loadtxt('C17_bpass21_imf135_300_sph_esc.dat')
   if inter == 0:
      sed_type ='BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. Spherical geometry. Not interpolated'
      print ('F_abs and U calculation usingBPASS models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type ='BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. Spherical geometry. Interpolated'
      print ('F_abs and U calculation usingBPASS models with spherical geometry and interpolation')


OHffs = []
eOHffs = []
Teffs = []
eTeffs = []
logUffs = []
elogUffs = []

Label_ID = False
Label_OH = False
Label_eOH = False
Label_OII = False
Label_eOII = False
Label_OIII_4959 = False
Label_eOIII_4959 = False
Label_OIII_5007 = False
Label_eOIII_5007 = False
Label_SII = False
Label_eSII = False
Label_SII_6716 = False
Label_eSII_6716 = False
Label_SII_6731 = False
Label_eSII_6731 = False
Label_SIII_9069 = False
Label_eSIII_9069 = False
Label_SIII_9532 = False
Label_eSIII_9532 = False
Label_HeI_4471 = False
Label_eHeI_4471 = False
Label_HeI_5876 = False
Label_eHeI_5876 = False
Label_HeII_4686 = False
Label_eHeII_4686 = False


for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == '12logOH':
      Label_OH = True
   if input1.dtype.names[col] == 'e12logOH':
      Label_eOH = True
   if input1.dtype.names[col] == 'OII_3727':
      Label_OII = True
   if input1.dtype.names[col] == 'eOII_3727':
      Label_eOII = True
   if input1.dtype.names[col] == 'OIII_4959':
      Label_OIII_4959 = True
   if input1.dtype.names[col] == 'eOIII_4959':
      Label_eOIII_4959 = True
   if input1.dtype.names[col] == 'OIII_5007':
      Label_OIII_5007 = True
   if input1.dtype.names[col] == 'eOIII_5007':
      Label_eOIII_5007 = True
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
   if input1.dtype.names[col] == 'SIII_9069':
      Label_SIII_9069 = True
   if input1.dtype.names[col] == 'eSIII_9069':
      Label_eSIII_9069 = True
   if input1.dtype.names[col] == 'SIII_9532':
      Label_SIII_9532 = True
   if input1.dtype.names[col] == 'eSIII_9532':
      Label_eSIII_9532 = True
   if input1.dtype.names[col] == 'HeI_4471':
      Label_HeI_4471 = True
   if input1.dtype.names[col] == 'eHeI_4471':
      Label_eHeI_4471 = True
   if input1.dtype.names[col] == 'HeI_5876':
      Label_HeI_5876 = True
   if input1.dtype.names[col] == 'eHeI_5876':
      Label_eHeI_5876 = True
   if input1.dtype.names[col] == 'HeII_4686':
      Label_HeII_4686 = True
   if input1.dtype.names[col] == 'eHeII_4686':
      Label_eHeII_4686 = True


if Label_ID == False:
   Names = np.arange(1,input1.size+1,1)
else:
   Names = input1['ID']
if Label_OH == False:
   logOH = np.zeros(input1.size)
else:
   logOH = input1['12logOH']
if Label_eOH == False:
   elogOH = np.zeros(input1.size)
else:
   elogOH = input1['e12logOH']
if Label_OII == False:
   OII_3727 = np.zeros(input1.size)
else:
   OII_3727 = input1['OII_3727']
if Label_eOII == False:
   eOII_3727 = np.zeros(input1.size)
else:
   eOII_3727 = input1['eOII_3727']
if Label_OIII_4959 == False and Label_OIII_5007 == False:
   OIII_5007 = np.zeros(input1.size)
elif Label_OIII_4959 == False and Label_OIII_5007 == True:
   OIII_5007 = input1['OIII_5007']
elif Label_OIII_4959 == True and Label_OIII_5007 == False:
   OIII_5007 = 4*input1['OIII_4959']
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
if Label_SIII_9069 == False and Label_SIII_9532 == False:
   SIII_9069 = np.zeros(input1.size)
elif Label_SIII_9069 == False and Label_SIII_9532 == True:
   SIII_9069 = input1['SIII_9069']/2.44
elif Label_SIII_9069 == True and Label_SIII_9532 == False:
   SIII_9069 = input1['SIII_9069']
else:
   SIII_9069 = (input1['SIII_9069']+input1['SIII_9532'])/3.44
if Label_eSIII_9069 == False and Label_eSIII_9532 == False:
   eSIII_9069 = np.zeros(input1.size)
elif Label_eSIII_9069 == False and Label_eSIII_9532 == True:
   eSIII_9069 = input1['eSIII_9069']/2.44
elif Label_eSIII_9069 == True and Label_eSIII_9532 == False:
   eSIII_9069 = input1['eSIII_9069']
else:
   eSIII_9069 = (input1['eSIII_9069']+input1['eSIII_9532'])/3.44
if Label_HeI_4471 == False:
   HeI_4471 = np.zeros(input1.size)
else:
   HeI_4471 = input1['HeI_4471']
if Label_eHeI_4471 == False:
   eHeI_4471 = np.zeros(input1.size)
else:
   eHeI_4471 = input1['eHeI_4471']
if Label_HeI_5876 == False:
   HeI_5876 = np.zeros(input1.size)
else:
   HeI_5876 = input1['HeI_5876']
if Label_eHeI_5876 == False:
   eHeI_5876 = np.zeros(input1.size)
else:
   eHeI_5876 = input1['eHeI_5876']
if Label_HeII_4686 == False:
   HeII_4686 = np.zeros(input1.size)
else:
   HeII_4686 = input1['HeII_4686']
if Label_eHeII_4686 == False:
   eHeII_4686 = np.zeros(input1.size)
else:
   eHeII_4686 = input1['eHeII_4686']


output = np.zeros(input1.size, dtype=[('ID', 'U12'), ('OII_3727', float),('eOII_3727', float),('OIII_5007', float),('eOIII_5007', float),('HeI_4471', float),('eHeI_4471', float),('HeI_5876', float),('eHeI_5876', float),('HeII_4686', float),('eHeII_4686', float),('SII_6725', float),('eSII_6725', float),('SIII_9069', float),('eSIII_9069', float),('OH', float),('eOH', float),('Teff', float),('eTeff', float),('logU', float),('elogU', float)] )

output['ID'] = Names
output['OII_3727'] = OII_3727
output['eOII_3727'] = eOII_3727
output['OIII_5007'] = OIII_5007
output['eOIII_5007'] = eOIII_5007
output['HeI_4471'] = HeI_4471
output['eHeI_4471'] = eHeI_4471
output['HeI_5876'] = HeI_5876
output['eHeI_5876'] = eHeI_5876
output['HeII_4686'] = HeII_4686
output['eHeII_4686'] = eHeII_4686
output['SII_6725'] = SII_6725
output['eSII_6725'] = eSII_6725
output['SIII_9069'] = SIII_9069
output['eSIII_9069'] = eSIII_9069


print ('Reading grids ....')
print ('')
print ('')
print ('---------------------------------------')
if sed < 3:
   print( '(%)   ID.   12+log(O/H)  T_eff(K)    log(U)')
else:
   print( '(%)   ID     12+log(O/H)  f_abs       log(U)')

print ('---------------------------------------')


# Beginning of loop of calculation

count = 0
for tab in range(0,len(input1),1):

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


      OII_3727_obs = 0
      if OII_3727[tab] > 0:
         while OII_3727_obs <= 0:
            OII_3727_obs = np.random.normal(OII_3727[tab],eOII_3727[tab]+1e-3)
      OIII_5007_obs = 0
      if OIII_5007[tab] > 0:
         while OIII_5007_obs <= 0:
            OIII_5007_obs = np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-3)
      SII_6725_obs = 0
      if SII_6725[tab] > 0:
         while SII_6725_obs <= 0:
            SII_6725_obs = np.random.normal(SII_6725[tab],eSII_6725[tab]+1e-3)
      SIII_9069_obs = 0
      if SIII_9069[tab] > 0:
         while SIII_9069_obs <= 0:
            SIII_9069_obs = np.random.normal(SIII_9069[tab],eSIII_9069[tab]+1e-3)
      HeI_4471_obs = 0
      if HeI_4471[tab] > 0:
         while HeI_4471_obs <= 0:
            HeI_4471_obs = np.random.normal(HeI_4471[tab],eHeI_4471[tab]+1e-3)
      HeI_5876_obs = 0
      if HeI_5876[tab] > 0:
         while HeI_5876_obs <= 0:
            HeI_5876_obs = np.random.normal(HeI_5876[tab],eHeI_5876[tab]+1e-3)
      HeII_4686_obs = 0
      if HeII_4686[tab] > 0:
         while HeII_4686_obs <= 0:
            HeII_4686_obs = np.random.normal(HeII_4686[tab],eHeII_4686[tab]+1e-3)
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
      if SII_6725_obs == 0 or OIII_5007_obs == 0:
         S2O3_obs = -10
      else:
         S2O3_obs = np.log10(SII_6725_obs / OIII_5007_obs )
      if HeI_4471_obs == 0 or HeII_4686_obs == 0:
         He12a_obs = -10
      else:
         He12a_obs = np.log10(HeI_4471_obs / HeII_4686_obs )
      if HeI_5876_obs == 0 or HeII_4686_obs == 0:
         He12b_obs = -10
      else:
         He12b_obs = np.log10(HeI_5876_obs / HeII_4686_obs )

      

# Interpolation of grid at specific O/H


      if logOH[tab] > 0:
         OH = np.random.normal(logOH[tab],elogOH[tab]+1e-3)
         OH_mc.append(OH)


         grid_T0 = []
         if OH <= 7.1:
            OH = 7.1
            i0 = 0
            i1 = bin
         elif OH >= 7.1 and OH < 7.4:
            i0 = 0
            i1 = bin
         elif OH >= 7.4 and OH < 7.7:
            i0 = bin
            i1 = 2*bin
         elif OH >= 7.7 and OH < 8.0:
            i0 = 2*bin
            i1 = 3*bin
         elif OH >= 8.0 and OH < 8.3:
            i0 = 3*bin
            i1 = 4*bin
         elif OH >= 8.3 and OH < 8.6:
            i0 = 4*bin
            i1 = 5*bin
         elif OH >= 8.6 and OH < 8.9:
            i0 = 5*bin
            i1 = 6*bin
         elif OH >= 8.9:
            OH = 8.9
            i0 = 5*bin
            i1 = 6*bin
       
         for x in range(0,bin):
            for y in range(0,10):
               grid_T0.append(grid[i0+x,y]*np.abs(0.3-OH+grid[i0,0])/0.3+grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)
            
      #         grid_T0.append(grid[i0+x,y]*np.abs(0.3-grid[i0,0]+OH)/0.3 + grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)


         grid_T = np.reshape(grid_T0,(bin,10))

      else:
         OH = 0
         OH_mc.append(OH)
         grid_T = grid
   

      np.savetxt('int_models.dat',grid_T,fmt='%.2f')


# Calculation of T and log U

      if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10:
         Teff = 0
         logU = 0
      else:
         CHI_O2O3 = 0
         CHI_S2S3 = 0
         CHI_S23 = 0
         CHI_S2O3 = 0
         CHI_He12a = 0
         CHI_He12b = 0
         CHI_He12 = 0

         for index in grid_T:
            if index[9] == 0 and HeII_4686_obs > 0: continue
            if S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index[5] == 0 or index[6] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            else:
               CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/S2S3_obs
               CHI_S23 = (index[5]+index[6]-S23_obs)**2/S23_obs
            if S2O3_obs == -10:
               CHI_S2O3 = 0
            elif index[5] == 0 or index[4] == 0:
               CHI_S2O3 = tol_max
            else:
               CHI_S2O3 = (np.log10(index[5]/index[4]) - S2O3_obs)**2/S2O3_obs
            if O2O3_obs == -10:
               CHI_O2O3 = 0
            elif index[3] == 0 or index[4] == 0:
               CHI_O2O3 = tol_max
            else:
               CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/O2O3_obs
            if He12a_obs == -10:
               CHI_He12a = 0
            elif index[7] == 0 or index[9] == 0:
               CHI_He12a = tol_max
            else:
               CHI_He12a = (np.log10(index[7]/index[9]) - He12a_obs)**2/He12a_obs
            if He12b_obs == -10:
               CHI_He12b = 0
            elif index[8] == 0 or index[9] == 0:
               CHI_He12b = tol_max
            else:
               CHI_He12b = (np.log10(index[8]/index[9]) - He12b_obs)**2/He12b_obs
            if CHI_He12a == 0 and CHI_He12b == 0:
               CHI_HE12 = 0
            elif CHI_He12a == 0:
               CHI_He12 = CHI_He12b
            elif CHI_He12b == 0:
               CHI_He12 = CHI_He12a
            else:
               CHI_He12 = (CHI_He12a + CHI_He12b)/2


            if OII_3727_obs == 0:
               CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2)**0.5
            elif SIII_9069_obs == 0 and HeII_4686_obs == 0:
               CHI_Teff = (CHI_O2O3**2 + CHI_S2O3**2)**0.5
            elif SIII_9069_obs == 0 and HeII_4686_obs > 0:
               CHI_Teff = (CHI_He12**2 + CHI_O2O3**2 )**0.5
            else:
               CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2)**0.5


         
            Teff_p = index[1]*(1/CHI_Teff)**2 + Teff_p
            logU_p = index[2] *(1/CHI_Teff)**2 + logU_p         
            den_Teff = (1/CHI_Teff)**2 + den_Teff
         Teff = Teff_p / den_Teff 
         logU = logU_p / den_Teff


# Calculation of T and log U errors


      if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10:
         eTeff = 0
         elogU = 0
      else:
         CHI_O2O3 = 0
         CHI_S2S3 = 0
         CHI_S23 = 0
         CHI_S2O3 = 0
         CHI_He12a = 0
         CHI_He12b = 0
         CHI_He12 = 0

         for index in grid_T:
            if index[9] == 0 and HeII_4686_obs > 0: continue
            if S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index[5] == 0 or index[6] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            else:
               CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/S2S3_obs
               CHI_S23 = (index[5]+index[6]-S23_obs)**2/S23_obs
            if S2O3_obs == -10:
               CHI_S2O3 = 0
            elif index[5] == 0 or index[4] == 0:
               CHI_S2O3 = tol_max
            else:
               CHI_S2O3 = (np.log10(index[5]/index[4]) - S2O3_obs)**2/S2O3_obs
            if O2O3_obs == -10:
               CHI_O2O3 = 0
            elif index[3] == 0 or index[4] == 0:
               CHI_O2O3 = tol_max
            else:
               CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/O2O3_obs
            if He12a_obs == -10:
               CHI_He12a = 0
            elif index[7] == 0 or index[9] == 0:
               CHI_He12a = tol_max
            else:
               CHI_He12a = (np.log10(index[7]/index[9]) - He12a_obs)**2/He12a_obs
            if He12b_obs == -10:
               CHI_He12b = 0
            elif index[8] == 0 or index[9] == 0:
               CHI_He12b = tol_max
            else:
               CHI_He12b = (np.log10(index[8]/index[9]) - He12b_obs)**2/He12b_obs
            if CHI_He12a == 0 and CHI_He12b == 0:
               CHI_HE12 = 0
            elif CHI_He12a == 0:
               CHI_He12 = CHI_He12b
            elif CHI_He12b == 0:
               CHI_He12 = CHI_He12a
            else:
               CHI_He12 = (CHI_He12a + CHI_He12b)/2


            if OII_3727_obs == 0:
               CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2)**0.5
            elif SIII_9069_obs == 0 and HeII_4686_obs == 0:
               CHI_Teff = (CHI_O2O3**2 + CHI_S2O3**2)**0.5
            elif SIII_9069_obs == 0 and HeII_4686_obs > 0:
               CHI_Teff = (CHI_He12**2 + CHI_O2O3**2 )**0.5
            else:
               CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2)**0.5




         
            if sed < 3:
               Teff_e = np.abs(index[1] - Teff) * (1/CHI_Teff) **2+ Teff_e
            else:
               Teff_e = np.abs(np.log10(index[1]+1e-5) - np.log10(Teff)) * (1/CHI_Teff)**2 + Teff_e
            logU_e = np.abs(index[2] - logU) * (1/CHI_Teff)**2 + logU_e         
            den_Teff_e = 1 * (1/CHI_Teff)**2 + den_Teff_e


         eTeff = Teff_e / den_Teff_e 
         if sed == 3:
            eTeff = Teff*np.log10(eTeff)
         elogU = logU_e / den_Teff_e



#Iterations for the interpolation mode

         if inter == 0:
            Teff = Teff
            logU = logU
         elif inter == 1:
            igrid = grid_T[np.lexsort((grid_T[:,0],grid_T[:,2]))]
            igrid = interpolate(igrid,1,Teff-eTeff-1e4,Teff+eTeff+1e4,10)
            igrid = igrid[np.lexsort((igrid[:,0],igrid[:,1]))]
            igrid = interpolate(igrid,2,logU-elogU-0.25,logU+elogU+0.25,10)

#            np.savetxt('int_models.dat',igrid,fmt='%.2f')


            if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10:
               Teff = 0
               logU = 0
            else:
               CHI_O2O3 = 0
               CHI_S2S3 = 0
               CHI_S23 = 0
               CHI_S2O3 = 0
               CHI_He12a = 0
               CHI_He12b = 0
               CHI_He12 = 0


               for index in igrid:
                  if index[9] == 0 and HeII_4686_obs > 0: continue
                  if S2S3_obs == -10:
                     CHI_S2S3 = 0
                     CHI_S23 = 0
                  elif index[5] == 0 or index[6] == 0:
                     CHI_S2S3 = tol_max
                     CHI_S23 = tol_max    
                  else:
                     CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/S2S3_obs
                     CHI_S23 = (index[5]+index[6]-S23_obs)**2/S23_obs
                  if S2O3_obs == -10:
                     CHI_S2O3 = 0
                  elif index[5] == 0 or index[4] == 0:
                     CHI_S2O3 = tol_max
                  else:
                     CHI_S2O3 = (np.log10(index[5]/index[4]) - S2O3_obs)**2/S2O3_obs
                  if O2O3_obs == -10:
                     CHI_O2O3 = 0
                  elif index[3] == 0 or index[4] == 0:
                     CHI_O2O3 = tol_max
                  else:
                     CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/O2O3_obs
                  if He12a_obs == -10:
                     CHI_He12a = 0
                  elif index[7] == 0 or index[9] == 0:
                     CHI_He12a = tol_max
                  else:
                     CHI_He12a = (np.log10(index[7]/index[9]) - He12a_obs)**2/He12a_obs
                  if He12b_obs == -10:
                     CHI_He12b = 0
                  elif index[8] == 0 or index[9] == 0:
                     CHI_He12b = tol_max
                  else:
                     CHI_He12b = (np.log10(index[8]/index[9]) - He12b_obs)**2/He12b_obs
                  if CHI_He12a == 0 and CHI_He12b == 0:
                     CHI_HE12 = 0
                  elif CHI_He12a == 0:
                     CHI_He12 = CHI_He12b
                  elif CHI_He12b == 0:
                     CHI_He12 = CHI_He12a
                  else:
                     CHI_He12 = (CHI_He12a + CHI_He12b)/2

   

                  if OII_3727_obs == 0:
                     CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2)**0.5
                  elif SIII_9069_obs == 0 and HeII_4686_obs == 0:
                     CHI_Teff = (CHI_O2O3**2 + CHI_S2O3**2)**0.5
                  elif SIII_9069_obs == 0 and HeII_4686_obs > 0:
                     CHI_Teff = (CHI_He12**2 + CHI_O2O3**2 )**0.5
                  else:
                     CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2)**0.5
      
   
         
                  if sed < 3:
                     Teff_p = index[1]*(1/CHI_Teff) + Teff_p
                  else:
                     Teff_p = index[1]*(1/CHI_Teff) + Teff_p
                  logU_p = index[2] *(1/CHI_Teff) + logU_p         
                  den_Teff = (1/CHI_Teff) + den_Teff
      
               Teff = Teff_p / den_Teff 
               if sed == 3:
                  Teff = Teff
               logU = logU_p / den_Teff
   


      Teff_mc.append(Teff)
      logU_mc.append(logU)
      eTeff_mc.append(eTeff)
      elogU_mc.append(elogU)


   

   if logOH[tab] > 0:
      OHf = logOH[tab]
      eOHf = elogOH[tab]
   else: 
      OHf = 0
      eOHf = 0
   Tefff = np.mean(Teff_mc)
   eTefff = (np.std(Teff_mc)**2 + np.mean(eTeff_mc)**2)**0.5
   logUf = np.mean(logU_mc)
   elogUf = (np.std(logU_mc)**2+np.mean(elogU_mc)**2)**0.5

   

   OHffs.append(OHf)
   eOHffs.append(eOHf)
   Teffs.append(Tefff)
   eTeffs.append(eTefff)
   logUffs.append(logUf)
   elogUffs.append(elogUf)  

   if input0.size == 1 and tab==0: continue

   if sed == 3:
      print (round(100*(count)/float(len(input1)),1),'%','',Names[tab], round(OHf,2), round(eOHf,2), round(Tefff,2),round(eTefff,2),round(logUf,2),round(elogUf,2)) 
   else:
      print (round(100*(count)/float(len(input1)),1),'%','',Names[tab], round(OHf,2), round(eOHf,2), 100*int(Tefff/100),100*int(eTefff/100),round(logUf,2),round(elogUf,2)) 


output['OH'] = OHffs
output['eOH'] = eOHffs
output['Teff'] = Teffs
output['eTeff'] = eTeffs
output['logU'] = logUffs
output['elogU'] = elogUffs

if input0.size == 1:  output = np.delete(output,obj=1,axis=0)



if sed == 3:
   lineas_header = [' HII-CHI-mistry v.5.0 output file',' Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','ID. O2Hb eO2Hb O3Hb  eO3Hb  eO3Hb S2Hb  eS2Hb S3Hb  eS3Hb 4471Hb e4471Hb 5678Hb e5648Hb He2Hb eHe2Hb O/H   eO/H  fabs  efabs logU elogU']
else:
   lineas_header = [' HII-CHI-mistry-Teff v.5.0 output file',' Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','ID.   O2Hb eO2Hb O3Hb  eO3Hb 4471Hb e4471Hb 5876Hb e5876Hb He2Hb eHe2Hb  S2Hb  eS2Hb S3Hb  eS3Hb  O/H   eO/H  Teff  eTeff logU elogU']


header = '\n'.join(lineas_header)


if sed == 3:
   np.savetxt(input00+'_hcm-teff-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*14+['%.2f']*2+['%.2f']*2+['%.2f']*2),header=header)
else:
   np.savetxt(input00+'_hcm-teff-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*14+['%.2f']*2+['%.0f']*2+['%.2f']*2),header=header)
print ('________________________________')
print ('Results are stored in '+input00+'_hcm-teff-output.dat')


