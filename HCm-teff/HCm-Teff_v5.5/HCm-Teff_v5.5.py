# Filename: HCm-Teff_v5.5.py

#####################
###### IMPORTS ######
#####################

import string
import numpy as np
import sys
sys.stderr = open('errorlog.txt', 'w')


#######################
###### FUNCTIONS ######
#######################

#Function for interpolation of grids
def interpolate(grid,z,zmin,zmax,n,param):
   if param == 1:
      label_t = 'Teff'
      file_1 = 'C17_WMb_Teff_30-60_pp.dat'
   else:
      label_t = 'f_abs'
      file_1 = 'C17_bpass21_imf135_300_sph_esc_Zg.dat'
   n_comments = 0
   with open('Libraries_teff/'+file_1, 'r') as file1:
      for line in file1:
         if line[0] == '#':
            n_comments += 1
   auxiliar_labels = np.genfromtxt('Libraries_teff/'+file_1, dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
   ncol = len(auxiliar_labels)
   vec = []
   type_list_names = []
   for col in auxiliar_labels:
      inter = 0
      no_inter = 0
      type_list_names.append((col, float))
      for row in range(0,len(grid)):
         if grid[z][row] < zmin or grid[z][row] > zmax: continue
         if z == 'logUsp': x = '12logOH'; y = label_t
         if z == label_t: x = '12logOH'; y = 'logUsp'
         if z == '12logOH': x = label_t; y = 'logUsp'
         if row == (len(grid)-1):
            vec.append(grid[col][row])
            no_inter = no_inter + 1
         elif grid[x][row] < grid[x][row+1] or grid[y][row] < grid[y][row+1] :
            vec.append(grid[col][row])
            no_inter = no_inter + 1
         else:
            inter = inter + 1
            for index in range(0,n):
               i = grid[col][row]+(index)*(grid[col][row+1]-grid[col][row])/n
               vec.append(i)
   out_aux = np.transpose(np.reshape(vec,(-1,n*inter+no_inter)))
   out = np.zeros(out_aux.shape[0], dtype=type_list_names)
   for col_n in range(0, len(auxiliar_labels)):
      out[auxiliar_labels[col_n]] = out_aux[:, col_n]
   return out

#################
##### LOGO ######
#################

def print_logo():
    logo_lines = [
        " ",
        " ",
        "    =================================",
        "    |                               |",
        "    |   H    H   CCCCC  m       m   |",
        "    |   H    H  C       m m   m m   |",
        "    |   HHHHHH  C       m   m   m   |",
        "    |   H    H  C       m       m   |",
        "    |   H    H   CCCCC  m       m   |",
        "    |                               |",
        "    =================================",
         
    ]

    full_name = "                       HII-CHI-Mistry"

    # Print the logo in huge letters
    for line in logo_lines:
        print(line)

    # Print the full name in smaller font
    print(full_name)
    print(" ")

if __name__ == "__main__":
    print_logo()

################################
###### INITIAL ITERATIONS ######
################################

#Description of the code
print (' ---------------------------------------------------------------------')
print ('This is HII-CHI-mistry-Teff v. 5.5')
print (' See Perez-Montero et al (2019) for details')
print ( ' Insert the name of your input text file with all or some of the following columns:')
print ('12+log(O/H)')
print (' 1549 CIV]')
print (' 1909 CIII]')
print (' 3727 [OII]')
print ('4471 Hei')
print ('4686 HeII')
print ('4740 [ArIV]')
print ('4959,5007 [OIII]')
print ('5876 HeI')
print ('6584 [NII]')
print ('6678 HeI')
print ('6717+31 [SII]')
print ('7135 [ArIII]')
print ('9069, 9532 [SIII]') 
print ('')
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
      #Counting comments:
      n_comments = 0
      with open(input00, 'r') as file:
         for line in file:
            if line[0] == '#':
               n_comments += 1
      input0 = np.genfromtxt(input00,dtype=None,names=True, skip_header=n_comments)
   else:
      #Counting comments:
      n_comments = 0
      with open(input00, 'r') as file:
         for line in file:
            if line[0] == '#':
               n_comments += 1
      input0 = np.genfromtxt(input00,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
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
if len(sys.argv) < 3:
   n = 25
else:
   n = int(sys.argv[2])
print ('The number of iterations for MonteCarlo simulation is: ',n)
print ('')


#################################################
###### NON INTERACTIVE OR INTERACTIVE CODE ######
#################################################

interactive = True #Change this value to False to run the code in non-interactive mode

#Questions (inputs fromt terminal)
question1 = interactive #Election of the quantities to be estimated
question2 = interactive #Question for the grids of models (only if "param" is set to 1)
question3 = interactive #Question for the geometry of photoionization models (only if "param" is set to 1)
question4 = interactive #Question for the grid of models (only if "param" is set to 2)
question5 = interactive #Question to use or not interpolation for the grids

#Set values to operate
if question1 == False:
   param = 2 #Choose quantity: (1) Effective temperature and ionization parameter; (2) photon absorption fraction and ionization parameter
if question2 == False and param == 1:
   sed = 2 #Choose grid of models: (1) WM-Basic  (30-60 kK); (2) WM-Basic  (30-60 kK) and Rauch (80-120 kK) stellar atmospheres; (3) Black body (30-100 kK)
if question3 == False and param == 3:
   geo = 1 #Choose value for the geometry: (1) Plane-parallel geometry; (2) Spherical geometry
if question4 == False and param == 2:
   sed_fabs = 1 #Choose grid of models: (1) BPASS cluster atmospheres, age = 4 Myr, Mup = 300, x = 1.35, w/binaries, Z* = Zg; (2) BPASS cluster atmospheres, age = 4 Myr, Mup = 300, x = 1.35, w/binaries, Z* = 1e-5'; (3) Black Body, T = 1e5 K
if question5 == False:
   inter = 0 #Choose value to perform interpolation: (0) no interpolation; (1) interpolation.  Replace value for the given option

   
#############################################
###### SELECTION OF THE GRID OF MODELS ######
#############################################

# Interface with the user
print ('')
while question1:
   print ('-------------------------------------------------')
   print ('------------')
   print ('(1) Effective temperature and ionization parameter')
   print ('(2) Photon absorption fraction and ionization parameter')
   print ('')
   #print('Other SED')
   #print('---------')
   #print ('(3) Different library')
   print('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      param = raw_input('Choose derived parameters: ')
   else:
      param = input('Choose derived parameters: ')
   if param == '1' or param == '2': question1 = False

print ('')
if param == '1': 
   sed_fabs = 0
   while question2:
      print ('---------------------------------------------------------------------')
      print ('(1) WM-Basic  (30-60 kK) ')
      print ('(2) WM-Basic  (30-60 kK) and Rauch (80-120 kK) stellar atmospheres')
      print ('(3) Black body (30-100 kK)')
      print ('---------------------------------------------------------------------')

      if int(sys.version[0]) < 3:
         sed = raw_input('Choose models:')
      else:
         sed = input('Choose and models:')
      if sed == '1' or sed == '2' or sed == '3': question2 = False
      print ('')
   while question3:
      print ('(1) Plane-parallel geometry')
      print ('(2) Spherical geometry')
      print ('---------------------------------------------------------------------')
      if int(sys.version[0]) < 3:
         geo = raw_input('Choose geometry of the models: ')
      else:
         geo = input('Choose geometry of the models: ')
      if geo == '1' or geo == '2': question3 = False
   print('')


if param == '2':
   print ('')
   sed = 0
   geo = 0
   while question4:
      print ('---------------------------------------------------------------------')
      print ('(1) BPASS cluster atmospheres, age = 4 Myr, Mup = 300, x = 1.35, w/binaries, Z* = Zg')
      print ('(2) BPASS cluster atmospheres, age = 4 Myr, Mup = 300, x = 1.35, w/binaries, Z* = 1e-5')
      print ('(3) Black Body, T = 1e5 K')
      print ('---------------------------------------------------------------------')

      if int(sys.version[0]) < 3:
         sed_fabs = raw_input('Choose models:')
      else:
         sed_fabs = input('Choose models:')
      if sed_fabs == '1' or sed_fabs == '2' or sed_fabs == '3': question4 = False


while question5:
   if int(sys.version[0]) < 3:
      inter = raw_input('Choose models [0] No interpolated [1] Interpolated: ')
   else:
      inter = input('Choose models [0] No interpolated [1] Interpolated: ')
   if inter == '0' or inter == '1': question5 = False
print ('')


param = int(param)
sed = int(sed)
sed_fabs = int(sed_fabs)
geo = int(geo)
inter = int(inter)


if geo==1 and sed == 1:
   bin = 99
   file_lib = 'C17_WMb_Teff_30-60_pp.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'WM-Basic stellar atmosphere. Plane-parallel geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic models with plane-paralell geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic stellar atmosphere. Plane-parallel geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic models with plane-paralell geometry and interpolation')
elif geo==2 and sed == 1:
   bin = 99
   file_lib = 'C17_WMb_Teff_30-60_sph.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'WM-Basic stellar atmosphereSpherical geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic stellar atmosphere. Spherical geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic models with spherical geometry and interpolation')

elif geo==1 and sed == 2:
   bin = 132
   file_lib = 'C17_WMb+Rauch_Teff_30-120_pp.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'WM-Basic and Rauch stellar atmosphere. Plane-parallel geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic and Rauch models with plane-parallel geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic and Rauch stellar atmosphere. Plane-parallel geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic and Rauch models with plane-parallel geometry and interpolation')
elif geo==2 and sed == 2:
   bin = 132
   file_lib = 'C17_WMb+Rauch_Teff_30-120_sph.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'WM-Basic  and Rauch stellar atmosphereSpherical geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic and Rauch models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic and Rauch stellar atmosphere. Spherical geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic and Rauch models with spherical geometry and interpolation')

elif geo==1 and sed == 3:
   bin = 143
   file_lib = 'C17_bb_Teff_30-100_pp.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Black body. Plane-parallel geometry. Not interpolated'
      print ('Teff and U calculation using black body models with plane-parallel geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'Black body. Plane-parallel geometry. Interpolated'
      print ('Teff and U calculation using black body models with plane-parallel geometry and interpolation')
elif geo==2 and sed == 3:
   bin = 143
   file_lib = 'C17_bb_Teff_30-100_sph.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Black body. spherical geometry. Not interpolated'
      print ('Teff and U calculation using black body models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'Black body. spherical geometry. Interpolated'
      print ('Teff and U calculation using black body models with spehrical geometry and interpolation')

elif sed_fabs == 1:
   sed_type = 'BPASS cluster atmospheres, Mup = 300 with binaries, Z = Zg. Spherical geometry'
   bin = 88
   file_lib = 'C17_bpass21_imf135_300_sph_esc_Zg.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type ='BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. and Z* = Zg. Not interpolated'
      print ('F_abs and U calculation using BPASS models withz* = Zg and non-interpolation')
   elif inter == 1:
      sed_type ='BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. and Z* = Zg. Interpolated'
      print ('F_abs and U calculation using BPASS models with Z* = Zg and interpolation')

elif sed_fabs == 2:
   bin = 88
   file_lib = 'C17_bpass21_imf135_300_sph_esc_Z5.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type ='BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. Z* = 1e-5. Not interpolated'
      print ('F_abs and U calculation using BPASS models with Z* = 1e-5 and non-interpolation')
   elif inter == 1:
      sed_type ='BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. Z* = 1e-5. Interpolated'
      print ('F_abs and U calculation using BPASS models with Z* = 1e-5 and interpolation')

elif sed_fabs == 3:
   bin = 88
   file_lib = 'C17_bb_T10.0_sph_esc.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type ='Black-body, T* = 1e5 K. Not interpolated'
      print ('F_abs and U calculation using Black-body 1e5 K models  and non-interpolation')
   elif inter == 1:
      sed_type ='Black-body, T* = 1e5 K.  Interpolated'
      print ('F_abs and U calculation using Black-body T* 1e5 K and interpolation')


#Different library
elif sed==4:
   file_lib = new_library
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+new_library, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1  
   grid_aux = np.genfromtxt('Libraries_teff/'+new_library,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'User file ' + new_library + ' used as library for the models no interpolated'
      print ('No interpolation for the library '+new_library)
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'User file ' + new_library + ' used as library for the models interpolated'
      print ('Interpolation for the library '+new_library)

###########################################
###### SUMMARY OF THE GRID OF MODELS ######
###########################################

#Final grid: no limitations
grid = grid_aux
print ('-------------------------------------------------')
print ('Summary of the models')
print ('---------------------')
print('The following grid is going to be used:')
print('- Full library: '+file_lib)
print('       Total number of models: ' + str(len(grid)))
print ('-------------------------------------------------')
print ('')


#################################################
###### CREATING ARRAY TO STORE ESTIMATIONS ######
#################################################

OHffs = []
eOHffs = []
Teffs = []
eTeffs = []
logUffs = []
elogUffs = []

#Labels to check information provided in the input file
Label_ID = False
Label_OH = False
Label_eOH = False
Label_CIV = False
Label_eCIV = False
Label_CIII = False
Label_eCIII = False
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
Label_HeI_6678 = False
Label_eHeI_6678 = False
Label_HeII_4686 = False
Label_eHeII_4686 = False
Label_ArIV_4740 = False
Label_eArIV_4740 = False
Label_ArIII_7135 = False
Label_eArIII_7135 = False
Label_NII_6584 = False
Label_eNII_6584 = False


for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == '12logOH':
      Label_OH = True
   if input1.dtype.names[col] == 'e12logOH':
      Label_eOH = True
   if input1.dtype.names[col] == 'CIV_1549':
      Label_CIV = True
   if input1.dtype.names[col] == 'eCIV_1549':
      Label_eCIV = True
   if input1.dtype.names[col] == 'CIII_1909':
      Label_CIII = True
   if input1.dtype.names[col] == 'eCIII_1909':
      Label_eCIII = True
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
   if input1.dtype.names[col] == 'HeI_6678':
      Label_HeI_6678 = True
   if input1.dtype.names[col] == 'eHeI_6678':
      Label_eHeI_6678 = True
   if input1.dtype.names[col] == 'HeII_4686':
      Label_HeII_4686 = True
   if input1.dtype.names[col] == 'eHeII_4686':
      Label_eHeII_4686 = True
   if input1.dtype.names[col] == 'ArIV_4740':
      Label_ArIV_4740 = True
   if input1.dtype.names[col] == 'eArIV_4740':
      Label_eArIV_4740 = True
   if input1.dtype.names[col] == 'ArIII_7135':
      Label_ArIII_7135 = True
   if input1.dtype.names[col] == 'eArIII_7135':
      Label_eArIII_7135 = True
   if input1.dtype.names[col] == 'NII_6584':
      Label_NII_6584 = True
   if input1.dtype.names[col] == 'eNII_6584':
      Label_eNII_6584 = True


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
if Label_CIV == False:
   CIV_1549 = np.zeros(input1.size)
else:
   CIV_1549 = input1['CIV_1549']
if Label_eCIV == False:
   eCIV_1549 = np.zeros(input1.size)
else:
   eCIV_1549 = input1['eCIV_1549']
if Label_CIII == False:
   CIII_1909 = np.zeros(input1.size)
else:
   CIII_1909 = input1['CIII_1909']
if Label_eCIII == False:
   eCIII_1909 = np.zeros(input1.size)
else:
   eCIII_1909 = input1['eCIII_1909']
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
if Label_OIII_4959 == False:
   OIII_4959 = np.zeros(input1.size)
else:
   OIII_4959 = input1['OIII_4959']
if Label_eOIII_4959 == False:
   eOIII_4959 = np.zeros(input1.size)
else:
   eOIII_4959 = input1['eOIII_4959']
if Label_SII_6716 == False:
   SII_6716 = np.zeros(input1.size)
else:
   SII_6716 = input1['SII_6716']
if Label_eSII_6716 == False:
   eSII_6716 = np.zeros(input1.size)
else:
   eSII_6716 = input1['eSII_6716']
if Label_SII_6731 == False:
   SII_6731 = np.zeros(input1.size)
else:
   SII_6731 = input1['SII_6731']
if Label_eSII_6731 == False:
   eSII_6731 = np.zeros(input1.size)
else:
   eSII_6731 = input1['eSII_6731']
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
   SIII_9069 = input1['SIII_9532']/2.44
elif Label_SIII_9069 == True and Label_SIII_9532 == False:
   SIII_9069 = input1['SIII_9069']
else:
   SIII_9069 = (input1['SIII_9069']+input1['SIII_9532'])/3.44
if Label_eSIII_9069 == False and Label_eSIII_9532 == False:
   eSIII_9069 = np.zeros(input1.size)
elif Label_eSIII_9069 == False and Label_eSIII_9532 == True:
   eSIII_9069 = input1['eSIII_9532']/2.44
elif Label_eSIII_9069 == True and Label_eSIII_9532 == False:
   eSIII_9069 = input1['eSIII_9069']
else:
   eSIII_9069 = (input1['eSIII_9069']+input1['eSIII_9532'])/3.44
if Label_SIII_9532 == False:
   SIII_9532 = np.zeros(input1.size)
else:
   SIII_9532 = input1['SIII_9532']
if Label_eSIII_9532 == False:
   eSIII_9532 = np.zeros(input1.size)
else:
   eSIII_9532 = input1['eSIII_9532']
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
if Label_HeI_6678 == False:
   HeI_6678 = np.zeros(input1.size)
else:
   HeI_6678 = input1['HeI_6678']
if Label_eHeI_6678 == False:
   eHeI_6678 = np.zeros(input1.size)
else:
   eHeI_6678 = input1['eHeI_6678']
if Label_HeII_4686 == False:
   HeII_4686 = np.zeros(input1.size)
else:
   HeII_4686 = input1['HeII_4686']
if Label_eHeII_4686 == False:
   eHeII_4686 = np.zeros(input1.size)
else:
   eHeII_4686 = input1['eHeII_4686']
if Label_ArIV_4740 == False:
   ArIV_4740 = np.zeros(input1.size)
else:
   ArIV_4740 = input1['ArIV_4740']
if Label_eArIV_4740 == False:
   eArIV_4740 = np.zeros(input1.size)
else:
   eArIV_4740 = input1['eArIV_4740']
if Label_ArIII_7135 == False:
   ArIII_7135 = np.zeros(input1.size)
else:
   ArIII_7135 = input1['ArIII_7135']
if Label_eArIII_7135 == False:
   eArIII_7135 = np.zeros(input1.size)
else:
   eArIII_7135 = input1['eArIII_7135']
if Label_NII_6584 == False:
   NII_6584 = np.zeros(input1.size)
else:
   NII_6584 = input1['NII_6584']
if Label_eNII_6584 == False:
   eNII_6584 = np.zeros(input1.size)
else:
   eNII_6584 = input1['eNII_6584']

################################################################
###### OUTPUT FORMAT AND INFORMATION: ONLY EMISSION LINES ######
################################################################

#Creation of output only with information from inputs
aux_list = []
aux_list.append(('ID','U12'))
if Label_CIV == True:
   aux_list.append(('CIV_1549', float))
if Label_eCIV == True:
   aux_list.append(('eCIV_1549', float))
if Label_CIII == True:
   aux_list.append(('CIII_1909', float))
if Label_eCIII == True:
   aux_list.append(('eCIII_1909', float))
if Label_OII == True:
   aux_list.append(('OII_3727', float))
if Label_eOII == True:
   aux_list.append(('eOII_3727', float))
if Label_OIII_4959 == True:
   aux_list.append(('OIII_4959', float))
if Label_eOIII_4959 == True:
   aux_list.append(('eOIII_4959', float))
if Label_OIII_5007 == True:
   aux_list.append(('OIII_5007', float))
if Label_eOIII_5007 == True:
   aux_list.append(('eOIII_5007', float))
if Label_SII == True:
   aux_list.append(('SII_6725', float))
if Label_eSII == True:
   aux_list.append(('eSII_6725', float))
if Label_SII_6716 == True:
   aux_list.append(('SII_6716', float))
if Label_eSII_6716 == True:
   aux_list.append(('eSII_6716', float))
if Label_SII_6731 == True:
   aux_list.append(('SII_6731', float))
if Label_eSII_6731 == True:
   aux_list.append(('eSII_6731', float))
if Label_SIII_9069 == True:
   aux_list.append(('SIII_9069', float))
if Label_eSIII_9069 == True:
   aux_list.append(('eSIII_9069', float))
if Label_SIII_9532 == True:
   aux_list.append(('SIII_9532', float))
if Label_eSIII_9532 == True:
   aux_list.append(('eSIII_9532', float))
if Label_HeI_4471 == True:
   aux_list.append(('HeI_4471', float))
if Label_eHeI_4471 == True:
   aux_list.append(('eHeI_4471', float))
if Label_HeI_5876 == True:
   aux_list.append(('HeI_5876', float))
if Label_eHeI_5876 == True:
   aux_list.append(('eHeI_5876', float))
if Label_HeI_6678 == True:
   aux_list.append(('HeI_6678', float))
if Label_eHeI_6678 == True:
   aux_list.append(('eHeI_6678', float))
if Label_HeII_4686 == True:
   aux_list.append(('HeII_4686', float))
if Label_eHeII_4686 == True:
   aux_list.append(('eHeII_4686', float))
if Label_ArIV_4740 == True:
   aux_list.append(('ArIV_4740', float))
if Label_eArIV_4740 == True:
   aux_list.append(('eArIV_4740', float))
if Label_ArIII_7135 == True:
   aux_list.append(('ArIII_7135', float))
if Label_eArIII_7135 == True:
   aux_list.append(('eArIII_7135', float))
if Label_NII_6584 == True:
   aux_list.append(('NII_6584', float))
if Label_eNII_6584 == True:
   aux_list.append(('eNII_6584', float))

aux_list.append(('OH', float))
aux_list.append(('eOH', float))
if param == 2:
   aux_list.append(('f_abs', float))
   aux_list.append(('ef_abs', float))
else:
   aux_list.append(('Teff', float))
   aux_list.append(('eTeff', float))
aux_list.append(('logU', float))
aux_list.append(('elogU', float))
output = np.zeros(input1.size, dtype=aux_list)

output['ID'] = Names
if Label_CIV == True:
   output['CIV_1549'] = CIV_1549
if Label_eCIV == True:
   output['eCIV_1549'] = eCIV_1549
if Label_CIII == True:
   output['CIII_1909'] = CIII_1909
if Label_eCIII == True:
   output['eCIII_1909'] = eCIII_1909
if Label_OII == True:
   output['OII_3727'] = OII_3727
if Label_eOII == True:
   output['eOII_3727'] = eOII_3727
if Label_OIII_4959 == True:
   output['OIII_4959'] = OIII_4959
if Label_eOIII_4959 == True:
   output['eOIII_4959'] = eOIII_4959
if Label_OIII_5007 == True:
   output['OIII_5007'] = OIII_5007
if Label_eOIII_5007 == True:
   output['eOIII_5007'] = eOIII_5007
if Label_HeI_4471 == True:
   output['HeI_4471'] = HeI_4471
if Label_eHeI_4471 == True:
   output['eHeI_4471'] = eHeI_4471
if Label_HeI_5876 == True:
   output['HeI_5876'] = HeI_5876
if Label_eHeI_5876 == True:
   output['eHeI_5876'] = eHeI_5876
if Label_HeI_6678 == True:
   output['HeI_6678'] = HeI_6678
if Label_eHeI_6678 == True:
   output['eHeI_6678'] = eHeI_6678
if Label_HeII_4686 == True:
   output['HeII_4686'] = HeII_4686
if Label_eHeII_4686 == True:
   output['eHeII_4686'] = eHeII_4686
if Label_SII == True:
   output['SII_6725'] = SII_6725
if Label_eSII == True:
   output['eSII_6725'] = eSII_6725
if Label_SII_6716 == True:
   output['SII_6716'] = SII_6716
if Label_eSII_6716 == True:
   output['eSII_6716'] = eSII_6716
if Label_SII_6731 == True:
   output['SII_6731'] = SII_6731
if Label_eSII_6731 == True:
   output['eSII_6731'] = eSII_6731
if Label_SIII_9069 == True:
   output['SIII_9069'] = SIII_9069
if Label_eSIII_9069 == True:
   output['eSIII_9069'] = eSIII_9069
if Label_SIII_9532 == True:
   output['SIII_9532'] = SIII_9532
if Label_eSIII_9532 == True:
   output['eSIII_9532'] = eSIII_9532
if Label_ArIV_4740 == True:
   output['ArIV_4740'] = ArIV_4740
if Label_eArIV_4740 == True:
   output['eArIV_4740'] = eArIV_4740
if Label_ArIII_7135 == True:
   output['ArIII_7135'] = ArIII_7135
if Label_eArIII_7135 == True:
   output['eArIII_7135'] = eArIII_7135
if Label_NII_6584 == True:
   output['NII_6584'] = NII_6584
if Label_eNII_6584 == True:
   output['eNII_6584'] = eNII_6584


################################################
###### ESTIMATIONS OF PHYSICAL PARAMETERS ######
################################################

print ('Reading grids ....')
print ('')
print ('')
print ('---------------------------------------')
if param == 1:
   print( '(%)   ID.   12+log(O/H)  T_eff(K)    log(U)')
else:
   print( '(%)   ID     12+log(O/H)  f_abs       log(U)')

print ('---------------------------------------')


# Beginning of loop of calculation
count = 0
for tab in range(0,len(input1),1):
   try:
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

         CIV_1549_obs = 0
         if CIV_1549[tab] > 0:
            while CIV_1549_obs <= 0:
               CIV_1549_obs = np.random.normal(CIV_1549[tab],eCIV_1549[tab]+1e-3)
         CIII_1909_obs = 0
         if CIII_1909[tab] > 0:
            while CIII_1909_obs <= 0:
               CIII_1909_obs = np.random.normal(CIII_1909[tab],eCIII_1909[tab]+1e-3)
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
         HeI_6678_obs = 0
         if HeI_6678[tab] > 0:
            while HeI_6678_obs <= 0:
               HeI_6678_obs = np.random.normal(HeI_6678[tab],eHeI_6678[tab]+1e-3)
         HeII_4686_obs = 0
         if HeII_4686[tab] > 0:
            while HeII_4686_obs <= 0:
               HeII_4686_obs = np.random.normal(HeII_4686[tab],eHeII_4686[tab]+1e-3)
         ArIV_4740_obs = 0
         if ArIV_4740[tab] > 0:
            while ArIV_4740_obs <= 0:
               ArIV_4740_obs = np.random.normal(ArIV_4740[tab],eArIV_4740[tab]+1e-3)
         ArIII_7135_obs = 0
         if ArIII_7135[tab] > 0:
            while ArIII_7135_obs <= 0:
               ArIII_7135_obs = np.random.normal(ArIII_7135[tab],eArIII_7135[tab]+1e-3)
         NII_6584_obs = 0
         if NII_6584[tab] > 0:
            while NII_6584_obs <= 0:
               NII_6584_obs = np.random.normal(NII_6584[tab],eNII_6584[tab]+1e-3)

         if CIV_1549_obs == 0 or CIII_1909_obs == 0:
            C3C4_obs = -10
         else:
            C3C4_obs = np.log10(CIII_1909_obs / CIV_1549_obs) 
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
         if HeI_6678_obs == 0 or HeII_4686_obs == 0:
            He12c_obs = -10
         else:
            He12c_obs = np.log10(HeI_6678_obs / HeII_4686_obs )
         if SII_6725_obs == 0 or ArIII_7135_obs == 0:
            S2Ar3_obs = -10
         else:
            S2Ar3_obs = np.log10(SII_6725_obs / ArIII_7135_obs )
         if ArIV_4740_obs == 0 or ArIII_7135_obs == 0:
            Ar3Ar4_obs = -10
         else:
            Ar3Ar4_obs = np.log10( ArIII_7135_obs/ArIV_4740_obs )
         if NII_6584_obs == 0 or OIII_5007_obs == 0:
            N2O3_obs = -10
         else:
            N2O3_obs = np.log10(NII_6584_obs / OIII_5007_obs )
         if NII_6584_obs == 0 or SIII_9069_obs == 0:
            N2S3_obs = -10
         else:
            N2S3_obs = np.log10(NII_6584_obs / SIII_9069_obs )
         if NII_6584_obs == 0 or ArIII_7135_obs == 0:
            N2Ar3_obs = -10
         else:
            N2Ar3_obs = np.log10(NII_6584_obs / ArIII_7135_obs )

         

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
               for y in grid.dtype.names:
                  grid_T0.append(grid[y][i0+x]*np.abs(0.3-OH+grid['12logOH'][i0])/0.3+grid[y][i1+x]*np.abs(0.3-grid['12logOH'][i1]+OH)/0.3)
                  #grid_T0.append(grid[i0+x,y]*np.abs(0.3-OH+grid[i0,0])/0.3+grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)
               
         #         grid_T0.append(grid[i0+x,y]*np.abs(0.3-grid[i0,0]+OH)/0.3 + grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)


            grid_T_aux = np.reshape(grid_T0,(bin,16))
            grid_T = np.zeros(grid_T_aux.shape[0], dtype=grid.dtype)
            for col_n in range(0, len(grid.dtype.names)):
               grid_T[grid.dtype.names[col_n]] = grid_T_aux[:, col_n]

         else:
            OH = 0
            OH_mc.append(OH)
            grid_T = grid
      

         #np.savetxt('int_models.dat',grid_T,fmt='%.2f')


   # Calculation of T and log U

         if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10 and He12c_obs == -10 and S2Ar3_obs == -10 and Ar3Ar4_obs == -10 and N2O3_obs == -10 and N2S3_obs == -10 and N2Ar3_obs == -10 and S2O3_obs == -10 and C3C4_obs == -10:
            if param == 1:
               Teff = 0
            else:
               Teff = -10
            logU = 0
         else:
            CHI_C3C4 = 0
            CHI_O2O3 = 0
            CHI_S2S3 = 0
            CHI_S23 = 0
            CHI_S2O3 = 0
            CHI_He12a = 0
            CHI_He12b = 0
            CHI_He12c = 0
            CHI_He12 = 0
            CHI_S2Ar3 = 0
            CHI_Ar3Ar4 = 0
            CHI_N2O3 = 0
            CHI_N2Ar3 = 0
            CHI_N2S3 = 0
            for index in grid_T:
               if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
               if C3C4_obs == -10:
                  CHI_C3C4 = 0
               elif index['CIV_1549'] == 0 or index['CIII_1909'] == 0:
                  CHI_C3C4 = tol_max
               else:
                  CHI_C3C4 = (np.log10(index['CIII_1909']/index['CIV_1549']) - C3C4_obs)**2/np.abs(C3C4_obs)
               if S2S3_obs == -10:
                  CHI_S2S3 = 0
                  CHI_S23 = 0
               elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
                  CHI_S2S3 = tol_max
                  CHI_S23 = tol_max    
               else:
                  CHI_S2S3 = (np.log10(index['SII_671731']/index['SIII_9069']) - S2S3_obs)**2/np.abs(S2S3_obs)
                  CHI_S23 = (index['SII_671731']+index['SIII_9069']-S23_obs)**2/np.abs(S23_obs)
               if S2O3_obs == -10:
                  CHI_S2O3 = 0
               elif index['SII_671731'] == 0 or index['OIII_5007'] == 0:
                  CHI_S2O3 = tol_max
               else:
                  CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/np.abs(S2O3_obs)
               if O2O3_obs == -10:
                  CHI_O2O3 = 0
               elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                  CHI_O2O3 = tol_max
               else:
                  CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/np.abs(O2O3_obs)
               if He12a_obs == -10:
                  CHI_He12a = 0
               elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
                  CHI_He12a = tol_max
               else:
                  CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/np.abs(He12a_obs)
               if He12b_obs == -10:
                  CHI_He12b = 0
               elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
                  CHI_He12b = tol_max
               else:
                  CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/np.abs(He12b_obs)
               if He12c_obs == -10:
                  CHI_He12c = 0
               elif index['HeI_6678'] == 0 or index['HeII_4686'] == 0:
                  CHI_He12c = tol_max
               else:
                  CHI_He12c = (np.log10(index['HeI_6678']/index['HeII_4686']) - He12c_obs)**2/np.abs(He12c_obs)
               if CHI_He12a == 0 and CHI_He12b == 0 and CHI_He12c == 0:
                  CHI_HE12 = 0
               elif CHI_He12b >= 0:
                  CHI_He12 = CHI_He12b
               elif CHI_He12c > 0:
                  CHI_He12 = CHI_He12c
               else:
                  CHI_He12 = CHI_He12a 
               if S2Ar3_obs == -10:
                  CHI_S2Ar3 = 0
               elif index['SII_671731'] == 0 or index['ArIII_7135'] == 0:
                  CHI_S2Ar3 = tol_max
               else:
                  CHI_S2Ar3 = (np.log10(index['SII_671731']/index['ArIII_7135']) - S2Ar3_obs)**2/np.abs(S2Ar3_obs)
               if Ar3Ar4_obs == -10:
                  CHI_Ar3Ar4 = 0
               elif index['ArIII_7135'] == 0 or index['ArIV_4740'] == 0:
                  CHI_Ar3Ar4 = tol_max
               else:
                  CHI_Ar3Ar4 = (np.log10(index['ArIII_7135']/index['ArIV_4740']) - Ar3Ar4_obs)**2/np.abs(Ar3Ar4_obs)
               if N2O3_obs == -10:
                  CHI_N2O3 = 0
               elif index['NII_6584'] == 0 or index['OIII_5007'] == 0:
                  CHI_N2O3 = tol_max
               else:
                  CHI_N2O3 = (np.log10(index['NII_6584']/index['OIII_5007']) - N2O3_obs)**2/np.abs(N2O3_obs)
               if N2S3_obs == -10:
                  CHI_N2S3 = 0
               elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
                  CHI_N2S3 = tol_max
               else:
                  CHI_N2S3 = (np.log10(index['NII_6584']/index['SIII_9069']) - N2S3_obs)**2/np.abs(N2S3_obs)
               if N2Ar3_obs == -10:
                  CHI_N2Ar3 = 0
               elif index['NII_6584'] == 0 or index['ArIII_7135'] == 0:
                  CHI_N2Ar3 = tol_max
               else:
                  CHI_N2Ar3 = (np.log10(index['NII_6584']/index['ArIII_7135']) - N2Ar3_obs)**2/np.abs(N2Ar3_obs)



               if OII_3727_obs == 0 and SIII_9069_obs >0:
                  CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
               elif SIII_9069_obs == 0 and OII_3727_obs > 0:
                  CHI_Teff = (CHI_O2O3**2 + CHI_S2Ar3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
               elif SIII_9069_obs == 0 and OII_3727_obs == 0:
                  CHI_Teff = (CHI_He12**2 + CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
               else:
                  CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5


               if CHI_Teff == 0:
                  Teff_p = Teff_p
                  logU_p = logU_p
                  den_Teff = den_Teff
               else:
                  if param == 1:
                     Teff_p = index['Teff']*(1/CHI_Teff)**2 + Teff_p
                  else:
                     Teff_p = index['f_abs']*(1/CHI_Teff)**2 + Teff_p 
                  logU_p = index['logUsp'] *(1/CHI_Teff)**2 + logU_p         
                  den_Teff = (1/CHI_Teff)**2 + den_Teff

            if Teff_p == 0:
               if param == 1:
                  Teff = 0
               else:
                  Teff = -10
               logU = 0
            else:
               Teff = Teff_p / den_Teff 
               logU = logU_p / den_Teff


   # Calculation of T and log U errors


         if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10 and He12c_obs == -10 and S2Ar3_obs == -10 and Ar3Ar4_obs == -10 and N2O3_obs == -10 and N2S3_obs == -10 and N2Ar3_obs == -10 and S2O3_obs == -10 and C3C4_obs == -10:
            if param == 1:
               eTeff = 0
            else:
               eTeff = -10
            elogU = 0
         else:
            CHI_C3C4 = 0
            CHI_O2O3 = 0
            CHI_S2S3 = 0
            CHI_S23 = 0
            CHI_S2O3 = 0
            CHI_He12a = 0
            CHI_He12b = 0
            CHI_He12c = 0
            CHI_He12 = 0
            CHI_S2Ar3 = 0
            CHI_Ar3Ar4 = 0
            CHI_N2O3 = 0
            CHI_N2S3 = 0
            CHI_N2Ar3 = 0
            for index in grid_T:
               if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
               if index['ArIV_4740'] > 0 and ArIV_4740_obs == -10: continue
               if C3C4_obs == -10:
                  CHI_C3C4 = 0
               elif index['CIV_1549'] == 0 or index['CIII_1909'] == 0:
                  CHI_C3C4 = tol_max
               else:
                  CHI_C3C4 = (np.log10(index['CIII_1909']/index['CIV_1549']) - C3C4_obs)**2/np.abs(C3C4_obs)
               if S2S3_obs == -10:
                  CHI_S2S3 = 0
                  CHI_S23 = 0
               elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
                  CHI_S2S3 = tol_max
                  CHI_S23 = tol_max    
               else:
                  CHI_S2S3 = (np.log10(index['SII_671731']/index['SIII_9069']) - S2S3_obs)**2/np.abs(S2S3_obs)
                  CHI_S23 = (index['SII_671731']+index['SIII_9069']-S23_obs)**2/np.abs(S23_obs)
               if S2O3_obs == -10:
                  CHI_S2O3 = 0
               elif index['SII_671731'] == 0 or index['OIII_5007'] == 0:
                  CHI_S2O3 = tol_max
               else:
                  CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/np.abs(S2O3_obs)
               if O2O3_obs == -10:
                  CHI_O2O3 = 0
               elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                  CHI_O2O3 = tol_max
               else:
                  CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/np.abs(O2O3_obs)
               if He12a_obs == -10:
                  CHI_He12a = 0
               elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
                  CHI_He12a = tol_max
               else:
                  CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/np.abs(He12a_obs)
               if He12b_obs == -10:
                  CHI_He12b = 0
               elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
                  CHI_He12b = tol_max
               else:
                  CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/np.abs(He12b_obs)
               if He12c_obs == -10:
                  CHI_He12c = 0
               elif index['HeI_6678'] == 0 or index['HeII_4686'] == 0:
                  CHI_He12c = tol_max
               else:
                  CHI_He12c = (np.log10(index['HeI_6678']/index['HeII_4686']) - He12c_obs)**2/np.abs(He12c_obs)
               if CHI_He12a == 0 and CHI_He12b == 0 and CHI_He12c == 0:
                  CHI_HE12 = 0
               elif CHI_He12b >= 0:
                  CHI_He12 = CHI_He12b
               elif CHI_He12c > 0:
                  CHI_He12 = CHI_He12c
               else:
                  CHI_He12 = CHI_He12a 
               if S2Ar3_obs == -10:
                  CHI_S2Ar3 = 0
               elif index['SII_671731'] == 0 or index['ArIII_7135'] == 0:
                  CHI_S2Ar3 = tol_max
               else:
                  CHI_S2Ar3 = (np.log10(index['SII_671731']/index['ArIII_7135']) - S2Ar3_obs)**2/np.abs(S2Ar3_obs)
               if Ar3Ar4_obs == -10:
                  CHI_Ar3Ar4 = 0
               elif index['ArIII_7135'] == 0 or index['ArIV_4740'] == 0:
                  CHI_Ar3Ar4 = tol_max
               else:
                  CHI_Ar3Ar4 = (np.log10(index['ArIII_7135']/index['ArIV_4740']) - Ar3Ar4_obs)**2/np.abs(Ar3Ar4_obs)
               if N2O3_obs == -10:
                  CHI_N2O3 = 0
               elif index['NII_6584'] == 0 or index['OIII_5007'] == 0:
                  CHI_N2O3 = tol_max
               else:
                  CHI_N2O3 = (np.log10(index['NII_6584']/index['OIII_5007']) - N2O3_obs)**2/np.abs(N2O3_obs)
               if N2S3_obs == -10:
                  CHI_N2S3 = 0
               elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
                  CHI_N2S3 = tol_max
               else:
                  CHI_N2S3 = (np.log10(index['NII_6584']/index['SIII_9069']) - N2S3_obs)**2/np.abs(N2S3_obs)
               if N2Ar3_obs == -10:
                  CHI_N2Ar3 = 0
               elif index['NII_6584'] == 0 or index['ArIII_7135'] == 0:
                  CHI_N2Ar3 = tol_max
               else:
                  CHI_N2Ar3 = (np.log10(index['NII_6584']/index['ArIII_7135']) - N2Ar3_obs)**2/np.abs(N2Ar3_obs)



               if OII_3727_obs == 0 and SIII_9069_obs >0:
                  CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
               elif SIII_9069_obs == 0 and OII_3727_obs > 0:
                  CHI_Teff = (CHI_O2O3**2 + CHI_S2Ar3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
               elif SIII_9069_obs == 0 and OII_3727_obs == 0:
                  CHI_Teff = (CHI_He12**2 + CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
               else:
                  CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5





            
               if CHI_Teff == 0:
                  Teff_e = Teff_e
                  logU_e = logU_e
                  den_Teff_e = den_Teff_e
               else:
                  if param == 1:
                     Teff_e = np.abs(index['Teff'] - Teff) * (1/CHI_Teff) **2+ Teff_e
                  else:
                     Teff_e = np.abs((index['f_abs']+1e-5) - (Teff)) * (1/CHI_Teff)**2 + Teff_e
                  logU_e = np.abs(index['logUsp'] - logU) * (1/CHI_Teff)**2 + logU_e         
                  den_Teff_e = 1 * (1/CHI_Teff)**2 + den_Teff_e


            if Teff_e == 0:
               if param == 1:
                  eTeff = 0
               else:
                  eTeff = -10
               elogU = 0
            else:
               if param == 1:
                  eTeff = Teff_e / den_Teff_e 
               else: 
                  eTeff = Teff_e / den_Teff_e 
               elogU = logU_e / den_Teff_e



   #Iterations for the interpolation mode

            if inter == 0:
               Teff = Teff
               logU = logU
            elif inter == 1:
               igrid = grid_T[np.lexsort((grid_T['12logOH'],grid_T['logUsp']))]
               if param == 1:
                  igrid = interpolate(igrid,'Teff',Teff-eTeff-eTeff/10,Teff+eTeff+eTeff/10,10, param)
                  igrid = igrid[np.lexsort((igrid['12logOH'],igrid['Teff']))]
               else:
                  igrid = interpolate(igrid,'f_abs',Teff-eTeff-eTeff/20,Teff+eTeff+eTeff/20,20, param)
                  igrid = igrid[np.lexsort((igrid['12logOH'],igrid['f_abs']))]
               igrid = interpolate(igrid,'logUsp',logU-elogU-0.25,logU+elogU+0.25,10, param)

   #            np.savetxt('int_models.dat',igrid,fmt='%.2f')


               if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10 and He12c_obs == -10 and S2Ar3_obs == -10 and Ar3Ar4_obs == -10 and N2O3_obs == -10 and N2S3_obs == -10 and N2Ar3_obs == -10 and S2O3_obs == -10 and C3C4_obs == -10:
                  Teff = 0
                  logU = 0
               else:
                  CHI_C3C4 = 0
                  CHI_O2O3 = 0
                  CHI_S2S3 = 0   count = count + 1
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


      CIV_1549_obs = 0
      if CIV_1549[tab] > 0:
         while CIV_1549_obs <= 0:
            CIV_1549_obs = np.random.normal(CIV_1549[tab],eCIV_1549[tab]+1e-3)
      CIII_1909_obs = 0
      if CIII_1909[tab] > 0:
         while CIII_1909_obs <= 0:
            CIII_1909_obs = np.random.normal(CIII_1909[tab],eCIII_1909[tab]+1e-3)
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
      HeI_6678_obs = 0
      if HeI_6678[tab] > 0:
         while HeI_6678_obs <= 0:
            HeI_6678_obs = np.random.normal(HeI_6678[tab],eHeI_6678[tab]+1e-3)
      HeII_4686_obs = 0
      if HeII_4686[tab] > 0:
         while HeII_4686_obs <= 0:
            HeII_4686_obs = np.random.normal(HeII_4686[tab],eHeII_4686[tab]+1e-3)
      ArIV_4740_obs = 0
      if ArIV_4740[tab] > 0:
         while ArIV_4740_obs <= 0:
            ArIV_4740_obs = np.random.normal(ArIV_4740[tab],eArIV_4740[tab]+1e-3)
      ArIII_7135_obs = 0
      if ArIII_7135[tab] > 0:
         while ArIII_7135_obs <= 0:
            ArIII_7135_obs = np.random.normal(ArIII_7135[tab],eArIII_7135[tab]+1e-3)
      NII_6584_obs = 0
      if NII_6584[tab] > 0:
         while NII_6584_obs <= 0:
            NII_6584_obs = np.random.normal(NII_6584[tab],eNII_6584[tab]+1e-3)

      if CIV_1549_obs == 0 or CIII_1909_obs == 0:
         C3C4_obs = -10
      else:
         C3C4_obs = np.log10(CIII_1909_obs / CIV_1549_obs) 
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
      if HeI_6678_obs == 0 or HeII_4686_obs == 0:
         He12c_obs = -10
      else:
         He12c_obs = np.log10(HeI_6678_obs / HeII_4686_obs )
      if SII_6725_obs == 0 or ArIII_7135_obs == 0:
         S2Ar3_obs = -10
      else:
         S2Ar3_obs = np.log10(SII_6725_obs / ArIII_7135_obs )
      if ArIV_4740_obs == 0 or ArIII_7135_obs == 0:
         Ar3Ar4_obs = -10
      else:
         Ar3Ar4_obs = np.log10( ArIII_7135_obs/ArIV_4740_obs )
      if NII_6584_obs == 0 or OIII_5007_obs == 0:
         N2O3_obs = -10
      else:
         N2O3_obs = np.log10(NII_6584_obs / OIII_5007_obs )
      if NII_6584_obs == 0 or SIII_9069_obs == 0:
         N2S3_obs = -10
      else:
         N2S3_obs = np.log10(NII_6584_obs / SIII_9069_obs )
      if NII_6584_obs == 0 or ArIII_7135_obs == 0:
         N2Ar3_obs = -10
      else:
         N2Ar3_obs = np.log10(NII_6584_obs / ArIII_7135_obs )

      

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
            for y in grid.dtype.names:
               grid_T0.append(grid[y][i0+x]*np.abs(0.3-OH+grid['12logOH'][i0])/0.3+grid[y][i1+x]*np.abs(0.3-grid['12logOH'][i1]+OH)/0.3)
               #grid_T0.append(grid[i0+x,y]*np.abs(0.3-OH+grid[i0,0])/0.3+grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)
            
      #         grid_T0.append(grid[i0+x,y]*np.abs(0.3-grid[i0,0]+OH)/0.3 + grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)


         grid_T_aux = np.reshape(grid_T0,(bin,16))
         grid_T = np.zeros(grid_T_aux.shape[0], dtype=grid.dtype)
         for col_n in range(0, len(grid.dtype.names)):
             grid_T[grid.dtype.names[col_n]] = grid_T_aux[:, col_n]

      else:
         OH = 0
         OH_mc.append(OH)
         grid_T = grid
   

      #np.savetxt('int_models.dat',grid_T,fmt='%.2f')


# Calculation of T and log U

      if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10 and He12c_obs == -10 and S2Ar3_obs == -10 and Ar3Ar4_obs == -10 and N2O3_obs == -10 and N2S3_obs == -10 and N2Ar3_obs == -10 and S2O3_obs == -10 and C3C4_obs == -10:
         if param == 1:
            Teff = 0
         else:
            Teff = -10
         logU = 0
      else:
         CHI_C3C4 = 0
         CHI_O2O3 = 0
         CHI_S2S3 = 0
         CHI_S23 = 0
         CHI_S2O3 = 0
         CHI_He12a = 0
         CHI_He12b = 0
         CHI_He12c = 0
         CHI_He12 = 0
         CHI_S2Ar3 = 0
         CHI_Ar3Ar4 = 0
         CHI_N2O3 = 0
         CHI_N2Ar3 = 0
         CHI_N2S3 = 0
         for index in grid_T:
            if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
            if C3C4_obs == -10:
               CHI_C3C4 = 0
            elif index['CIV_1549'] == 0 or index['CIII_1909'] == 0:
               CHI_C3C4 = tol_max
            else:
               CHI_C3C4 = (np.log10(index['CIII_1909']/index['CIV_1549']) - C3C4_obs)**2/np.abs(C3C4_obs)
            if S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            else:
               CHI_S2S3 = (np.log10(index['SII_671731']/index['SIII_9069']) - S2S3_obs)**2/np.abs(S2S3_obs)
               CHI_S23 = (index['SII_671731']+index['SIII_9069']-S23_obs)**2/np.abs(S23_obs)
            if S2O3_obs == -10:
               CHI_S2O3 = 0
            elif index['SII_671731'] == 0 or index['OIII_5007'] == 0:
               CHI_S2O3 = tol_max
            else:
               CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/np.abs(S2O3_obs)
            if O2O3_obs == -10:
               CHI_O2O3 = 0
            elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
               CHI_O2O3 = tol_max
            else:
               CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/np.abs(O2O3_obs)
            if He12a_obs == -10:
               CHI_He12a = 0
            elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
               CHI_He12a = tol_max
            else:
               CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/np.abs(He12a_obs)
            if He12b_obs == -10:
               CHI_He12b = 0
            elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
               CHI_He12b = tol_max
            else:
               CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/np.abs(He12b_obs)
            if He12c_obs == -10:
               CHI_He12c = 0
            elif index['HeI_6678'] == 0 or index['HeII_4686'] == 0:
               CHI_He12c = tol_max
            else:
               CHI_He12c = (np.log10(index['HeI_6678']/index['HeII_4686']) - He12c_obs)**2/np.abs(He12c_obs)
            if CHI_He12a == 0 and CHI_He12b == 0 and CHI_He12c == 0:
               CHI_HE12 = 0
            elif CHI_He12b >= 0:
               CHI_He12 = CHI_He12b
            elif CHI_He12c > 0:
               CHI_He12 = CHI_He12c
            else:
               CHI_He12 = CHI_He12a 
            if S2Ar3_obs == -10:
               CHI_S2Ar3 = 0
            elif index['SII_671731'] == 0 or index['ArIII_7135'] == 0:
               CHI_S2Ar3 = tol_max
            else:
               CHI_S2Ar3 = (np.log10(index['SII_671731']/index['ArIII_7135']) - S2Ar3_obs)**2/np.abs(S2Ar3_obs)
            if Ar3Ar4_obs == -10:
               CHI_Ar3Ar4 = 0
            elif index['ArIII_7135'] == 0 or index['ArIV_4740'] == 0:
               CHI_Ar3Ar4 = tol_max
            else:
               CHI_Ar3Ar4 = (np.log10(index['ArIII_7135']/index['ArIV_4740']) - Ar3Ar4_obs)**2/np.abs(Ar3Ar4_obs)
            if N2O3_obs == -10:
               CHI_N2O3 = 0
            elif index['NII_6584'] == 0 or index['OIII_5007'] == 0:
               CHI_N2O3 = tol_max
            else:
               CHI_N2O3 = (np.log10(index['NII_6584']/index['OIII_5007']) - N2O3_obs)**2/np.abs(N2O3_obs)
            if N2S3_obs == -10:
               CHI_N2S3 = 0
            elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
               CHI_N2S3 = tol_max
            else:
               CHI_N2S3 = (np.log10(index['NII_6584']/index['SIII_9069']) - N2S3_obs)**2/np.abs(N2S3_obs)
            if N2Ar3_obs == -10:
               CHI_N2Ar3 = 0
            elif index['NII_6584'] == 0 or index['ArIII_7135'] == 0:
               CHI_N2Ar3 = tol_max
            else:
               CHI_N2Ar3 = (np.log10(index['NII_6584']/index['ArIII_7135']) - N2Ar3_obs)**2/np.abs(N2Ar3_obs)



            if OII_3727_obs == 0 and SIII_9069_obs >0:
               CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
            elif SIII_9069_obs == 0 and OII_3727_obs > 0:
               CHI_Teff = (CHI_O2O3**2 + CHI_S2Ar3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
            elif SIII_9069_obs == 0 and OII_3727_obs == 0:
               CHI_Teff = (CHI_He12**2 + CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
            else:
               CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5


            if CHI_Teff == 0:
               Teff_p = Teff_p
               logU_p = logU_p
               den_Teff = den_Teff
            else:
               if param == 1:
                  Teff_p = index['Teff']*(1/CHI_Teff)**2 + Teff_p
               else:
                  Teff_p = index['f_abs']*(1/CHI_Teff)**2 + Teff_p 
               logU_p = index['logUsp'] *(1/CHI_Teff)**2 + logU_p         
               den_Teff = (1/CHI_Teff)**2 + den_Teff

         if Teff_p == 0:
            if param == 1:
               Teff = 0
            else:
               Teff = -10
            logU = 0
         else:
            Teff = Teff_p / den_Teff 
            logU = logU_p / den_Teff


# Calculation of T and log U errors


      if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10 and He12c_obs == -10 and S2Ar3_obs == -10 and Ar3Ar4_obs == -10 and N2O3_obs == -10 and N2S3_obs == -10 and N2Ar3_obs == -10 and S2O3_obs == -10 and C3C4_obs == -10:
         if param == 1:
            eTeff = 0
         else:
            eTeff = -10
         elogU = 0
      else:
         CHI_C3C4 = 0
         CHI_O2O3 = 0
         CHI_S2S3 = 0
         CHI_S23 = 0
         CHI_S2O3 = 0
         CHI_He12a = 0
         CHI_He12b = 0
         CHI_He12c = 0
         CHI_He12 = 0
         CHI_S2Ar3 = 0
         CHI_Ar3Ar4 = 0
         CHI_N2O3 = 0
         CHI_N2S3 = 0
         CHI_N2Ar3 = 0
         for index in grid_T:
            if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
            if index['ArIV_4740'] > 0 and ArIV_4740_obs == -10: continue
            if C3C4_obs == -10:
               CHI_C3C4 = 0
            elif index['CIV_1549'] == 0 or index['CIII_1909'] == 0:
               CHI_C3C4 = tol_max
            else:
               CHI_C3C4 = (np.log10(index['CIII_1909']/index['CIV_1549']) - C3C4_obs)**2/np.abs(C3C4_obs)
            if S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            else:
               CHI_S2S3 = (np.log10(index['SII_671731']/index['SIII_9069']) - S2S3_obs)**2/np.abs(S2S3_obs)
               CHI_S23 = (index['SII_671731']+index['SIII_9069']-S23_obs)**2/np.abs(S23_obs)
            if S2O3_obs == -10:
               CHI_S2O3 = 0
            elif index['SII_671731'] == 0 or index['OIII_5007'] == 0:
               CHI_S2O3 = tol_max
            else:
               CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/np.abs(S2O3_obs)
            if O2O3_obs == -10:
               CHI_O2O3 = 0
            elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
               CHI_O2O3 = tol_max
            else:
               CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/np.abs(O2O3_obs)
            if He12a_obs == -10:
               CHI_He12a = 0
            elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
               CHI_He12a = tol_max
            else:
               CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/np.abs(He12a_obs)
            if He12b_obs == -10:
               CHI_He12b = 0
            elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
               CHI_He12b = tol_max
            else:
               CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/np.abs(He12b_obs)
            if He12c_obs == -10:
               CHI_He12c = 0
            elif index['HeI_6678'] == 0 or index['HeII_4686'] == 0:
               CHI_He12c = tol_max
            else:
               CHI_He12c = (np.log10(index['HeI_6678']/index['HeII_4686']) - He12c_obs)**2/np.abs(He12c_obs)
            if CHI_He12a == 0 and CHI_He12b == 0 and CHI_He12c == 0:
               CHI_HE12 = 0
            elif CHI_He12b >= 0:
               CHI_He12 = CHI_He12b
            elif CHI_He12c > 0:
               CHI_He12 = CHI_He12c
            else:
               CHI_He12 = CHI_He12a 
            if S2Ar3_obs == -10:
               CHI_S2Ar3 = 0
            elif index['SII_671731'] == 0 or index['ArIII_7135'] == 0:
               CHI_S2Ar3 = tol_max
            else:
               CHI_S2Ar3 = (np.log10(index['SII_671731']/index['ArIII_7135']) - S2Ar3_obs)**2/np.abs(S2Ar3_obs)
            if Ar3Ar4_obs == -10:
               CHI_Ar3Ar4 = 0
            elif index['ArIII_7135'] == 0 or index['ArIV_4740'] == 0:
               CHI_Ar3Ar4 = tol_max
            else:
               CHI_Ar3Ar4 = (np.log10(index['ArIII_7135']/index['ArIV_4740']) - Ar3Ar4_obs)**2/np.abs(Ar3Ar4_obs)
            if N2O3_obs == -10:
               CHI_N2O3 = 0
            elif index['NII_6584'] == 0 or index['OIII_5007'] == 0:
               CHI_N2O3 = tol_max
            else:
               CHI_N2O3 = (np.log10(index['NII_6584']/index['OIII_5007']) - N2O3_obs)**2/np.abs(N2O3_obs)
            if N2S3_obs == -10:
               CHI_N2S3 = 0
            elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
               CHI_N2S3 = tol_max
            else:
               CHI_N2S3 = (np.log10(index['NII_6584']/index['SIII_9069']) - N2S3_obs)**2/np.abs(N2S3_obs)
            if N2Ar3_obs == -10:
               CHI_N2Ar3 = 0
            elif index['NII_6584'] == 0 or index['ArIII_7135'] == 0:
               CHI_N2Ar3 = tol_max
            else:
               CHI_N2Ar3 = (np.log10(index['NII_6584']/index['ArIII_7135']) - N2Ar3_obs)**2/np.abs(N2Ar3_obs)



            if OII_3727_obs == 0 and SIII_9069_obs >0:
               CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
            elif SIII_9069_obs == 0 and OII_3727_obs > 0:
               CHI_Teff = (CHI_O2O3**2 + CHI_S2Ar3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
            elif SIII_9069_obs == 0 and OII_3727_obs == 0:
               CHI_Teff = (CHI_He12**2 + CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
            else:
               CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5





         
            if CHI_Teff == 0:
               Teff_e = Teff_e
               logU_e = logU_e
               den_Teff_e = den_Teff_e
            else:
               if param == 1:
                  Teff_e = np.abs(index['Teff'] - Teff) * (1/CHI_Teff) **2+ Teff_e
               else:
                  Teff_e = np.abs((index['f_abs']+1e-5) - (Teff)) * (1/CHI_Teff)**2 + Teff_e
               logU_e = np.abs(index['logUsp'] - logU) * (1/CHI_Teff)**2 + logU_e         
               den_Teff_e = 1 * (1/CHI_Teff)**2 + den_Teff_e


         if Teff_e == 0:
            if param == 1:
               eTeff = 0
            else:
               eTeff = -10
            elogU = 0
         else:
            if param == 1:
               eTeff = Teff_e / den_Teff_e 
            else: 
               eTeff = Teff_e / den_Teff_e 
            elogU = logU_e / den_Teff_e



#Iterations for the interpolation mode

         if inter == 0:
            Teff = Teff
            logU = logU
         elif inter == 1:
            igrid = grid_T[np.lexsort((grid_T['12logOH'],grid_T['logUsp']))]
            if param == 1:
               igrid = interpolate(igrid,'Teff',Teff-eTeff-eTeff/10,Teff+eTeff+eTeff/10,10, param)
               igrid = igrid[np.lexsort((igrid['12logOH'],igrid['Teff']))]
            else:
               igrid = interpolate(igrid,'f_abs',Teff-eTeff-eTeff/20,Teff+eTeff+eTeff/20,20, param)
               igrid = igrid[np.lexsort((igrid['12logOH'],igrid['f_abs']))]
            igrid = interpolate(igrid,'logUsp',logU-elogU-0.25,logU+elogU+0.25,10, param)

#            np.savetxt('int_models.dat',igrid,fmt='%.2f')


            if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10 and He12c_obs == -10 and S2Ar3_obs == -10 and Ar3Ar4_obs == -10 and N2O3_obs == -10 and N2S3_obs == -10 and N2Ar3_obs == -10 and S2O3_obs == -10 and C3C4_obs == -10:
               Teff = 0
               logU = 0
            else:
               CHI_C3C4 = 0
               CHI_O2O3 = 0
               CHI_S2S3 = 0
               CHI_S23 = 0
               CHI_S2O3 = 0
               CHI_He12a = 0
               CHI_He12b = 0
               CHI_He12c = 0
               CHI_He12 = 0
               CHI_S2Ar3 = 0
               CHI_Ar3Ar4 = 0
               CHI_N2O3 = 0
               CHI_N2S3 = 0
               CHI_N2Ar3 = 0
               for index in igrid:
                  if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
                  if C3C4_obs == -10:
                     CHI_C3C4 = 0
                  elif index['CIV_1549'] == 0 or index['CIII_1909'] == 0:
                     CHI_C3C4 = tol_max
                  else:
                     CHI_C3C4 = (np.log10(index['CIII_1909']/index['CIV_1549']) - C3C4_obs)**2/C3C4_obs
                  if S2S3_obs == -10:
                     CHI_S2S3 = 0
                     CHI_S23 = 0
                  elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
                     CHI_S2S3 = tol_max
                     CHI_S23 = tol_max    
                  else:
                     CHI_S2S3 = (np.log10(index['SII_671731']/index['SIII_9069']) - S2S3_obs)**2/S2S3_obs
                     CHI_S23 = (index['SII_671731']+index[6]-S23_obs)**2/S23_obs
                  if S2O3_obs == -10:
                     CHI_S2O3 = 0
                  elif index['SII_671731'] == 0 or index['OIII_5007'] == 0:
                     CHI_S2O3 = tol_max
                  else:
                     CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/np.abs(S2O3_obs)
                  if O2O3_obs == -10:
                     CHI_O2O3 = 0
                  elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                     CHI_O2O3 = tol_max
                  else:
                     CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/np.abs(O2O3_obs)
                  if He12a_obs == -10:
                     CHI_He12a = 0
                  elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
                     CHI_He12a = tol_max
                  else:
                     CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/np.abs(He12a_obs)
                  if He12b_obs == -10:
                     CHI_He12b = 0
                  elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
                     CHI_He12b = tol_max
                  else:
                     CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/np.abs(He12b_obs)
                  if He12c_obs == -10:
                     CHI_He12c = 0
                  elif index['HeI_6678'] == 0 or index['HeII_4686'] == 0:
                     CHI_He12c = tol_max
                  else:
                     CHI_He12c = (np.log10(index['HeI_6678']/index['HeII_4686']) - He12c_obs)**2/np.abs(He12c_obs)
                  if CHI_He12a == 0 and CHI_He12b == 0 and CHI_He12c == 0:
                     CHI_HE12 = 0
                  elif CHI_He12b >= 0:
                     CHI_He12 = CHI_He12b
                  elif CHI_He12c > 0:
                     CHI_He12 = CHI_He12c
                  else:
                     CHI_He12 = CHI_He12a 
                  if S2Ar3_obs == -10:
                     CHI_S2Ar3 = 0
                  elif index['SII_671731'] == 0 or index['ArIII_7135'] == 0:
                     CHI_S2Ar3 = tol_max
                  else:
                     CHI_S2Ar3 = (np.log10(index['SII_671731']/index['ArIII_7135']) - S2Ar3_obs)**2/np.abs(S2Ar3_obs)
                  if Ar3Ar4_obs == -10:
                     CHI_Ar3Ar4 = 0
                  elif index['ArIII_7135'] == 0 or index['ArIV_4740'] == 0:
                     CHI_Ar3Ar4 = tol_max
                  else:
                     CHI_Ar3Ar4 = (np.log10(index['ArIII_7135']/index['ArIV_4740']) - Ar3Ar4_obs)**2/np.abs(Ar3Ar4_obs)
                  if N2O3_obs == -10:
                     CHI_N2O3 = 0
                  elif index['NII_6584'] == 0 or index['OIII_5007'] == 0:
                     CHI_N2O3 = tol_max
                  else:
                     CHI_N2O3 = (np.log10(index['NII_6584']/index['OIII_5007']) - N2O3_obs)**2/np.abs(N2O3_obs)
                  if N2S3_obs == -10:
                     CHI_N2S3 = 0
                  elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
                     CHI_N2S3 = tol_max
                  else:
                     CHI_N2S3 = (np.log10(index['NII_6584']/index['SIII_9069']) - N2S3_obs)**2/np.abs(N2S3_obs)
                  if N2Ar3_obs == -10:
                     CHI_N2Ar3 = 0
                  elif index['NII_6584'] == 0 or index['ArIII_7135'] == 0:
                     CHI_N2Ar3 = tol_max
                  else:
                     CHI_N2Ar3 = (np.log10(index['NII_6584']/index['ArIII_7135']) - N2Ar3_obs)**2/np.abs(N2Ar3_obs)


                  if OII_3727_obs == 0 and SIII_9069_obs >0:
                     CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
                  elif SIII_9069_obs == 0 and OII_3727_obs > 0:
                     CHI_Teff = (CHI_O2O3**2 + CHI_S2Ar3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
                  elif SIII_9069_obs == 0 and OII_3727_obs == 0:
                     CHI_Teff = (CHI_He12**2 +  CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
                  else:
                     CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5

                  if CHI_Teff == 0:
                     Teff_p = Teff_p
                     logU_p = logU_p
                  else:
                     if param == 1:
                        Teff_p = index['Teff']*(1/CHI_Teff) + Teff_p
                     else:
                        Teff_p = index['f_abs']*(1/CHI_Teff) + Teff_p
                     logU_p = index['logUsp'] *(1/CHI_Teff) + logU_p         
                     den_Teff = (1/CHI_Teff) + den_Teff
      
               if Teff_p == 0:
                  if param == 1:
                     Teff = 0
                  else:
                     Teff = -10
                  logU = 0
               else:
                  Teff = Teff_p / den_Teff 
                  if param == 2:
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
                  CHI_S23 = 0
                  CHI_S2O3 = 0
                  CHI_He12a = 0
                  CHI_He12b = 0
                  CHI_He12c = 0
                  CHI_He12 = 0
                  CHI_S2Ar3 = 0
                  CHI_Ar3Ar4 = 0
                  CHI_N2O3 = 0
                  CHI_N2S3 = 0
                  CHI_N2Ar3 = 0
                  for index in igrid:
                     if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
                     if C3C4_obs == -10:
                        CHI_C3C4 = 0
                     elif index['CIV_1549'] == 0 or index['CIII_1909'] == 0:
                        CHI_C3C4 = tol_max
                     else:
                        CHI_C3C4 = (np.log10(index['CIII_1909']/index['CIV_1549']) - C3C4_obs)**2/C3C4_obs
                     if S2S3_obs == -10:
                        CHI_S2S3 = 0
                        CHI_S23 = 0
                     elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
                        CHI_S2S3 = tol_max
                        CHI_S23 = tol_max    
                     else:
                        CHI_S2S3 = (np.log10(index['SII_671731']/index['SIII_9069']) - S2S3_obs)**2/S2S3_obs
                        CHI_S23 = (index['SII_671731']+index[6]-S23_obs)**2/S23_obs
                     if S2O3_obs == -10:
                        CHI_S2O3 = 0
                     elif index['SII_671731'] == 0 or index['OIII_5007'] == 0:
                        CHI_S2O3 = tol_max
                     else:
                        CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/np.abs(S2O3_obs)
                     if O2O3_obs == -10:
                        CHI_O2O3 = 0
                     elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                        CHI_O2O3 = tol_max
                     else:
                        CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/np.abs(O2O3_obs)
                     if He12a_obs == -10:
                        CHI_He12a = 0
                     elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
                        CHI_He12a = tol_max
                     else:
                        CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/np.abs(He12a_obs)
                     if He12b_obs == -10:
                        CHI_He12b = 0
                     elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
                        CHI_He12b = tol_max
                     else:
                        CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/np.abs(He12b_obs)
                     if He12c_obs == -10:
                        CHI_He12c = 0
                     elif index['HeI_6678'] == 0 or index['HeII_4686'] == 0:
                        CHI_He12c = tol_max
                     else:
                        CHI_He12c = (np.log10(index['HeI_6678']/index['HeII_4686']) - He12c_obs)**2/np.abs(He12c_obs)
                     if CHI_He12a == 0 and CHI_He12b == 0 and CHI_He12c == 0:
                        CHI_HE12 = 0
                     elif CHI_He12b >= 0:
                        CHI_He12 = CHI_He12b
                     elif CHI_He12c > 0:
                        CHI_He12 = CHI_He12c
                     else:
                        CHI_He12 = CHI_He12a 
                     if S2Ar3_obs == -10:
                        CHI_S2Ar3 = 0
                     elif index['SII_671731'] == 0 or index['ArIII_7135'] == 0:
                        CHI_S2Ar3 = tol_max
                     else:
                        CHI_S2Ar3 = (np.log10(index['SII_671731']/index['ArIII_7135']) - S2Ar3_obs)**2/np.abs(S2Ar3_obs)
                     if Ar3Ar4_obs == -10:
                        CHI_Ar3Ar4 = 0
                     elif index['ArIII_7135'] == 0 or index['ArIV_4740'] == 0:
                        CHI_Ar3Ar4 = tol_max
                     else:
                        CHI_Ar3Ar4 = (np.log10(index['ArIII_7135']/index['ArIV_4740']) - Ar3Ar4_obs)**2/np.abs(Ar3Ar4_obs)
                     if N2O3_obs == -10:
                        CHI_N2O3 = 0
                     elif index['NII_6584'] == 0 or index['OIII_5007'] == 0:
                        CHI_N2O3 = tol_max
                     else:
                        CHI_N2O3 = (np.log10(index['NII_6584']/index['OIII_5007']) - N2O3_obs)**2/np.abs(N2O3_obs)
                     if N2S3_obs == -10:
                        CHI_N2S3 = 0
                     elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
                        CHI_N2S3 = tol_max
                     else:
                        CHI_N2S3 = (np.log10(index['NII_6584']/index['SIII_9069']) - N2S3_obs)**2/np.abs(N2S3_obs)
                     if N2Ar3_obs == -10:
                        CHI_N2Ar3 = 0
                     elif index['NII_6584'] == 0 or index['ArIII_7135'] == 0:
                        CHI_N2Ar3 = tol_max
                     else:
                        CHI_N2Ar3 = (np.log10(index['NII_6584']/index['ArIII_7135']) - N2Ar3_obs)**2/np.abs(N2Ar3_obs)


                     if OII_3727_obs == 0 and SIII_9069_obs >0:
                        CHI_Teff = (CHI_S2S3**2 + CHI_He12**2 +  CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
                     elif SIII_9069_obs == 0 and OII_3727_obs > 0:
                        CHI_Teff = (CHI_O2O3**2 + CHI_S2Ar3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
                     elif SIII_9069_obs == 0 and OII_3727_obs == 0:
                        CHI_Teff = (CHI_He12**2 +  CHI_S2O3**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5
                     else:
                        CHI_Teff = (CHI_S2S3**2 + CHI_O2O3**2 + CHI_He12**2 + CHI_Ar3Ar4**2 + CHI_N2O3**2 + CHI_N2S3**2 + CHI_N2Ar3**2 + CHI_C3C4**2)**0.5

                     if CHI_Teff == 0:
                        Teff_p = Teff_p
                        logU_p = logU_p
                     else:
                        if param == 1:
                           Teff_p = index['Teff']*(1/CHI_Teff) + Teff_p
                        else:
                           Teff_p = index['f_abs']*(1/CHI_Teff) + Teff_p
                        logU_p = index['logUsp'] *(1/CHI_Teff) + logU_p         
                        den_Teff = (1/CHI_Teff) + den_Teff
         
                  if Teff_p == 0:
                     if param == 1:
                        Teff = 0
                     else:
                        Teff = -10
                     logU = 0
                  else:
                     Teff = Teff_p / den_Teff 
                     if param == 2:
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
      
   except:
      OHf = 9999
      eOHf = 9999
      Tefff = -9999
      eTefff = -9999
      logUf = 9999
      elogUf = 9999
   

   OHffs.append(OHf)
   eOHffs.append(eOHf)
   Teffs.append(Tefff)
   eTeffs.append(eTefff)
   logUffs.append(logUf)
   elogUffs.append(elogUf)  

   ##################################
   # Displaying results in terminal #
   ##################################
   
   if input0.size == 1 and tab==0: continue

   if param == 2:
      print (round(100*(count)/float(len(input1)),1),'%','',Names[tab], round(OHf,2), round(eOHf,2), round(Tefff,2),round(eTefff,2),round(logUf,2),round(elogUf,2)) 
   else:
      print (round(100*(count)/float(len(input1)),1),'%','',Names[tab], round(OHf,2), round(eOHf,2), 100*int(Tefff/100),100*int(eTefff/100),round(logUf,2),round(elogUf,2)) 


####################################################
###### OUTPUT FORMAT AND INFORMATION: RESULTS ######
####################################################

output['OH'] = OHffs
output['eOH'] = eOHffs
if param == 1:
   output['Teff'] = Teffs
   output['eTeff'] = eTeffs
else:
   output['f_abs'] = Teffs
   output['ef_abs'] = eTeffs
output['logU'] = logUffs
output['elogU'] = elogUffs


if input0.size == 1:  output = np.delete(output,obj=1,axis=0)

lineas_header = [' HII-CHI-mistry-Teff v.5.5 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type, 'Library file used : '+file_lib, '']

line_label = '{:10}  '.format(output.dtype.names[0])
for ind2 in range(1, len(output.dtype.names)-6):
   line_label += '{:10}  '.format(output.dtype.names[ind2])

if param == 1:
   line_label += '{:10}   {:10}   {:10}   {:10}   {:10}   {:10}'.format('O/H',    'eO/H',  'Teff' ,   'eTeff' , 'logU',   'elogU')
else:
   line_label += '{:10}   {:10}   {:10}   {:10}   {:10}   {:10}'.format('O/H',    'eO/H',  'f_abs' ,   'ef_abs' , 'logU',   'elogU')

lineas_header.append(line_label)
header = '\n'.join(lineas_header)
if param == 1:
   np.savetxt('.'.join(input00.split('.')[:-1])+'_hcm-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*(len(output.dtype.names)-7)+['%.2f']*2+['%i']*2+['%.2f']*2), header=header)
else:
   np.savetxt('.'.join(input00.split('.')[:-1])+'_hcm-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*(len(output.dtype.names)-7)+['%.2f']*2+['%.2f']*2+['%.2f']*2), header=header)

lines_stor = []
with open('.'.join(input00.split('.')[:-1])+'_hcm-output.dat', 'r+') as output_file:
   for line in output_file:
      lines_stor.append(line)

file_overwrite = open('.'.join(input00.split('.')[:-1])+'_hcm-output.dat', 'r+')
file_overwrite.seek(0)
for line_n in lines_stor:  
   if line_n[0] == '#' and line_n[2:4] == 'ID':
      file_overwrite.write(line_n[2:])
   else:
      file_overwrite.write(line_n)
file_overwrite.truncate()   
file_overwrite.close()

print ('-------------------------------------------------')
print ('Results are stored in ' + '.'.join(input00.split('.')[:-1])+'_hcm-output.dat')
print ('-------------------------------------------------')


#############################################
###### INFORMATION AND CONTACT DETAILS ######
#############################################

# Enrique Perez-Montero, epm@iaa.es
# Borja Perez-Diaz, bperez@iaa.es


#################
###### END ######
#################
