# Filename: HCm_Teff_v4.1.py

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
print ('This is HII-CHI-mistry_Teff v. 4.1')
print (' See Perez-Montero et al (2019) for details')
print ( ' Insert the name of your input text file with the following columns:')
print ('12+log(O/H), 3727 [OII], 5007 [OIII], 6725 [SII], 9069 [SIII], 4471 HeI, 5876 HeI, and 4686 HeII')
print ('with their corresponding errors in adjacent columns')
print ('relative to Hbeta or 0 for missing information.')
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
try:
   input0 = np.loadtxt(input00)
   if (input0.ndim == 1 and input0.shape[0] != 16) or (input0.ndim > 1 and input0.shape[1] != 16):
      print ('The input file does not have 16 columns. Please check')
      sys.exit()
   print ('The input file is:'+input00)
except:
   print ('Input file error: It does not exist or has wrong format')
   sys.exit()
print ('')

output = []


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

if input0.shape == (16,):
   input1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,input0[0],input0[1],input0[2],input0[3],input0[4],input0[5],input0[6],input0[7],input0[8],input0[9],input0[10],input0[11],input0[12],input0[13],input0[14],input0[15]]
   input = np.reshape(input1,(2,14))
else:
   input = input0


print ('Reading grids ....')
print ('')
print ('')
print ('---------------------------------------')
if sed < 3:
   print( '(%)   12+log(O/H)  T_eff(K)    log(U)')
else:
   print( '(%)   12+log(O/H)  f_abs       log(U)')

print ('---------------------------------------')


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
      HeI_4471_obs = np.random.normal(tab[10],tab[11]+1e-3)
      if HeI_4471_obs <= 0: heI_4471_obs = 0
      HeI_5876_obs = np.random.normal(tab[12],tab[13]+1e-3)
      if HeI_5876_obs <= 0: heI_5876_obs = 0
      HeII_4686_obs = np.random.normal(tab[14],tab[15]+1e-3)
      if HeII_4686_obs <= 0: HeII_4686_obs = 0
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
      if tab[4] == 0 or tab[6] == 0 or SII_6725_obs == 0 or OIII_5007_obs == 0:
         S2O3_obs = -10
      else:
         S2O3_obs = np.log10(SII_6725_obs / OIII_5007_obs )
      if tab[10] == 0 or tab[14] == 0 or HeI_4471_obs == 0 or HeII_4686_obs == 0:
         He12a_obs = -10
      else:
         He12a_obs = np.log10(HeI_4471_obs / HeII_4686_obs )
      if tab[12] == 0 or tab[14] == 0 or HeI_5876_obs == 0 or HeII_4686_obs == 0:
         He12b_obs = -10
      else:
         He12b_obs = np.log10(HeI_5876_obs / HeII_4686_obs )

      

# Interpolation of grid at specific O/H


      if tab[0] > 0:
         OH = np.random.normal(tab[0],tab[1]+1e-3)
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
            if index[9] > 0 and HeII_4686_obs == -10: continue
            if index[5] == 0 or index[6] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            elif S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            else:
               CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/S2S3_obs
               CHI_S23 = (index[5]+index[6]-S23_obs)**2/S23_obs
            if index[5] == 0 or index[4] == 0:
               CHI_S2O3 = tol_max
            elif S2O3_obs == -10:
               CHI_S2O3 = 0
            else:
               CHI_S2O3 = (np.log10(index[5]/index[4]) - S2O3_obs)**2/S2O3_obs
            if index[3] == 0 or index[4] == 0:
               CHI_O2O3 = tol_max
            elif O2O3_obs == -10:
               CHI_O2O3 = 0
            else:
               CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/O2O3_obs
            if index[7] == 0 or index[9] == 0:
               CHI_He12a = tol_max
            elif He12a_obs == -10:
               CHI_He12a = 0
            else:
               CHI_He12a = (np.log10(index[7]/index[9]) - He12a_obs)**2/He12a_obs
            if index[8] == 0 or index[9] == 0:
               CHI_He12b = tol_max
            elif He12b_obs == -10:
               CHI_He12b = 0
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
            if index[9] > 0 and HeII_4686_obs == -10: continue
            if index[5] == 0 or index[6] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            elif S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            else:
               CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/S2S3_obs
               CHI_S23 = (index[5]+index[6]-S23_obs)**2/S23_obs
            if index[5] == 0 or index[4] == 0:
               CHI_S2O3 = tol_max
            elif S2O3_obs == -10:
               CHI_S2O3 = 0
            else:
               CHI_S2O3 = (np.log10(index[5]/index[4]) - S2O3_obs)**2/S2O3_obs
            if index[3] == 0 or index[4] == 0:
               CHI_O2O3 = tol_max
            elif O2O3_obs == -10:
               CHI_O2O3 = 0
            else:
               CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/O2O3_obs
            if index[7] == 0 or index[9] == 0:
               CHI_He12a = tol_max
            elif He12a_obs == -10:
               CHI_He12a = 0
            else:
               CHI_He12a = (np.log10(index[7]/index[9]) - He12a_obs)**2/He12a_obs
            if index[8] == 0 or index[9] == 0:
               CHI_He12b = tol_max
            elif He12b_obs == -10:
               CHI_He12b = 0
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
               Teff_e = np.abs(index[1] - Teff) * (1/CHI_Teff) + Teff_e
            else:
               Teff_e = np.abs(np.log10(index[1]+1e-5) - np.log10(Teff)) * (1/CHI_Teff) + Teff_e
            logU_e = np.abs(index[2] - logU) * (1/CHI_Teff) + logU_e         
            den_Teff_e = 1 * (1/CHI_Teff) + den_Teff_e


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

            np.savetxt('int_models.dat',igrid,fmt='%.2f')


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
                  if index[9] > 0 and HeII_4686_obs == -10: continue
                  if index[5] == 0 or index[6] == 0:
                     CHI_S2S3 = tol_max
                     CHI_S23 = tol_max    
                  elif S2S3_obs == -10:
                     CHI_S2S3 = 0
                     CHI_S23 = 0
                  else:
                     CHI_S2S3 = (np.log10(index[5]/index[6]) - S2S3_obs)**2/S2S3_obs
                     CHI_S23 = (index[5]+index[6]-S23_obs)**2/S23_obs
                  if index[5] == 0 or index[4] == 0:
                     CHI_S2O3 = tol_max
                  elif S2O3_obs == -10:
                     CHI_S2O3 = 0
                  else:
                     CHI_S2O3 = (np.log10(index[5]/index[4]) - S2O3_obs)**2/S2O3_obs
                  if index[3] == 0 or index[4] == 0:
                     CHI_O2O3 = tol_max
                  elif O2O3_obs == -10:
                     CHI_O2O3 = 0
                  else:
                     CHI_O2O3 = (np.log10(index[3]/index[4]) - O2O3_obs)**2/O2O3_obs
                  if index[7] == 0 or index[9] == 0:
                     CHI_He12a = tol_max
                  elif He12a_obs == -10:
                     CHI_He12a = 0
                  else:
                     CHI_He12a = (np.log10(index[7]/index[9]) - He12a_obs)**2/He12a_obs
                  if index[8] == 0 or index[9] == 0:
                     CHI_He12b = tol_max
                  elif He12b_obs == -10:
                     CHI_He12b = 0
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
   output.append(tab[10])
   output.append(tab[11])
   output.append(tab[12])
   output.append(tab[13])
   output.append(tab[14])
   output.append(tab[15])

   output.append(OHf)
   output.append(eOHf)
   output.append(Tefff)
   output.append(eTefff)
   output.append(logUf)
   output.append(elogUf)  

   if input0.shape >= (16,) and count == 1: continue

   if sed == 3:
      print round(100*(count)/float(len(input)),1),'%','', round(OHf,2), round(eOHf,2), round(Tefff,2),round(eTefff,2),round(logUf,2),round(elogUf,2) 
   else:
      print round(100*(count)/float(len(input)),1),'%','', round(OHf,2), round(eOHf,2), 100*int(Tefff/100),100*int(eTefff/100),round(logUf,2),round(elogUf,2) 


out = np.reshape(output,(len(input),20))
if input0.shape == (16,): out = np.delete(out,obj=0,axis=0)

if sed == 3:
   lineas_header = [' HII-CHI-mistry v.4.1 output file',' Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','O2Hb eO2Hb O3Hb  eO3Hb  eO3Hb S2Hb  eS2Hb S3Hb  eS3Hb 4471Hb e4471Hb 5678Hb e5648Hb He2Hb eHe2Hb O/H   eO/H  fabs  efabs logU elogU']
else:
   lineas_header = [' HII-CHI-mistry-Teff v.4.1 output file',' Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'','O2Hb eO2Hb O3Hb  eO3Hb  eO3Hb S2Hb  eS2Hb S3Hb  eS3Hb 4471Hb e4471Hb 5678Hb e5678Hb He2Hb eHe2Hb O/H   eO/H  Teff  eTeff logU elogU']


header = '\n'.join(lineas_header)


if sed == 3:
   np.savetxt(input00+'_hcm-output.dat',out,fmt=' '.join(['%.3f']*14+['%.2f']*2+['%.3f']*2+['%.2f']*2),header=header)
else:
   np.savetxt(input00+'_hcm-output.dat',out,fmt=' '.join(['%.3f']*14+['%.2f']*2+['%.0f']*2+['%.2f']*2),header=header)
print ('________________________________')
print ('Results are stored in '+input00+'_hcm-output.dat')


