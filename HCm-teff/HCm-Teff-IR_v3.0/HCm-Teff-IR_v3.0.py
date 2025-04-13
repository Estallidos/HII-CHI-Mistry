# Filename: HCm-Teff-IR_v3.0.py

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
def interpolate(grid,z,zmin,zmax,n):
   label_t = 'T_eff'
   file_1 = 'C17_WMb-IR_Teff_30-60_pp.dat'
   n_comments = 0
   with open('Libraries_teff-IR/'+file_1, 'r') as file1:
      for line in file1:
         if line[0] == '#':
            n_comments += 1
   auxiliar_labels = np.genfromtxt('Libraries_teff-IR/'+file_1, dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
   ncol = len(auxiliar_labels)
   vec = []
   type_list_names = []
   for col in auxiliar_labels:
      inter = 0
      no_inter = 0
      type_list_names.append((col, float))
      for row in range(0,len(grid)):
         if grid[z][row] < zmin or grid[z][row] > zmax: continue
         if z == 'logU': x = '12logOH'; y = label_t
         if z == label_t: x = '12logOH'; y = 'logU'
         if z == '12logOH': x = label_t; y = 'logU'
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
print ('This is HII-CHI-mistry-Teff-IR v. 3.0')
print (' See Perez-Montero et al (2024) for details')
print ( ' Insert the name of your input text file with all or some of the following columns:')
print ('')
print ('12+log(O/H)')
print ('[ArII] 7.0m')
print ('[ArV] 7.9m')
print ('[ArII] 9.0m')
print ('[SIV] 10.5m')
print ('[NeII] 12.8m')
print ('[ArV] 13.1m')
print ('[NeV] 14.3m')
print ('[NeIII] 15.5m')
print ('[SIII] 18.7m')
print ('[NeV] 24.3m')
print ('[OIV] 25.9m')
print ('[SIII] 33.7m')
print ('[OIII] 52m')
print ('[NIII] 57m')
print ('[OIII] 88m') 
print ('[NII] 122m')
print ('[NII] 205m')
print ('')
print ('in arbitrary units or 0 for missing information')
print ('with their corresponding labels and errors.')
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

#############################################
###### NON INTERACTIVE OR INTERACTIVE CODE ######
#############################################

interactive = True #Change this value to False to run the code in non-interactive mode

#Questions (inputs fromt terminal)
question1 = interactive #Question for Teff or alpha_ox calculation
question2 = interactive #For Teff question for the grids of models
question3 = interactive #For Teff question for the geometry of photoionization models
question4 = interactive #For alpha question for the fraction of free electrons
question5 = interactive #For alpha question for the presence of dust
question6 = interactive #Question to use or not interpolation for the grids

#Set values to operate
if question1 == False:
   param = 2 #Choose (1) Teff calculation or (2) alpha_OX calculation in AGN
if question2 == False:
   sed = 2 #Choose grid of models: (1) WM-Basic  (30-60 kK); (2) WM-Basic  (30-60 kK) and Rauch (80-120 kK) stellar atmospheres; (3) Black body (30-100 kK)
if question3 == False:
   geo = 1 #Choose value for the geometry: (1) Plane-parallel geometry; (2) Spherical geometry
if question4 == False:
   efrac = 2 #Choose fraction of free electrons as stopping criterion: (1) 2%, (2) 98%, or (3) 99.9%
if question5 == False:
   grains = 1 #Choose (1) presence or (2) absence of dust
if question6 == False:
   inter = 0 #Choose value to perform interpolation: (0) no interpolation; (1) interpolation.  Replace value for the given option

   
#############################################
###### SELECTION OF THE GRID OF MODELS ######
#############################################

# Interface with the user
print ('')
print ('')
while question1:
   print ('-------------------------------------------------')
   print ('------------')
   print ('(1) Effective temperature and ionization parameter')
   print ('(2) alpha(OX) AGN power law, T =1.5e5 K, alpha(UV) = -0.5')


   print ('')
   print('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      param = raw_input('Choose derived parameters: ')
   else:
      param = input('Choose derived parameters: ')
   if param == '1' or param == '2': question1 = False

print ('')
if param == '1': 
   while question2:
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
   sed = 0
   geo = 0
   #FRACTION OF FREE ELECTRONS
   while question4:
      if int(sys.version[0]) < 3:
         efrac = raw_input('Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons  [3] 99.9% free electrons: ')
      else:
         efrac = input('Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons  [3]  99.9% free electrons:  ')
      if efrac == '1' or efrac == '2' or efrac == '3': question4 = False
      efrac = int(efrac)
   #Presence or absence of dust in the models

   while question5:
      if int(sys.version[0]) < 3:
         grains = raw_input('Choose AGN models with [1] or without [2] dust grains: ')
      else:
         grains = input('Choose AGN models with [1] or without [2] dust grains:  ')
      if grains == '1' or grains == '2': question5 = False
      grains = int(grains)
   print ('')



while question6:
   if int(sys.version[0]) < 3:
      inter = raw_input('Choose models [0] No interpolated [1] Interpolated: ')
   else:
      inter = input('Choose models [0] No interpolated [1] Interpolated: ')
   if inter == '0' or inter == '1': question6 = False
print ('')


param = int(param)
sed = int(sed)
geo = int(geo)
inter = int(inter)



if param == 1 and geo==1 and sed == 1:
   bin = 99
   file_lib = 'C17_WMb-IR_Teff_30-60_pp.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'WM-Basic stellar atmosphere. Plane-parallel geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic models with plane-paralell geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic stellar atmosphere. Plane-parallel geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic models with plane-paralell geometry and interpolation')
elif param == 1 and geo==2 and sed == 1:
   bin = 99
   file_lib = 'C17_WMb-IR_Teff_30-60_sph.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'WM-Basic stellar atmosphereSpherical geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic stellar atmosphere. Spherical geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic models with spherical geometry and interpolation')

elif param == 1 and geo==1 and sed == 2:
   bin = 132
   file_lib = 'C17_WMb+Rauch-IR_Teff_30-120_pp.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'WM-Basic and Rauch stellar atmosphere. Plane-parallel geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic and Rauch models with plane-parallel geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic and Rauch stellar atmosphere. Plane-parallel geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic and Rauch models with plane-parallel geometry and interpolation')
elif param == 1 and geo==2 and sed == 2:
   bin = 132
   file_lib = 'C17_WMb+Rauch-IR_Teff_30-120_sph.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'WM-Basic  and Rauch stellar atmosphereSpherical geometry. Not interpolated'
      print ('Teff and U calculation using WM-Basic and Rauch models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'WM-Basic and Rauch stellar atmosphere. Spherical geometry. Interpolated'
      print ('Teff and U calculation using WM-Basic and Rauch models with spherical geometry and interpolation')

elif param == 1 and geo==1 and sed == 3:
   bin = 143
   file_lib = 'C17_bb-IR_Teff_30-100_pp.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Black body. Plane-parallel geometry. Not interpolated'
      print ('Teff and U calculation using black body models with plane-parallel geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'Black body. Plane-parallel geometry. Interpolated'
      print ('Teff and U calculation using black body models with plane-parallel geometry and interpolation')
elif param == 1 and geo==2 and sed == 3:
   bin = 143
   file_lib = 'C17_bb-IR_Teff_30-100_sph.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Black body. spherical geometry. Not interpolated'
      print ('Teff and U calculation using black body models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type = 'Black body. spherical geometry. Interpolated'
      print ('Teff and U calculation using black body models with spehrical geometry and interpolation')

elif param  == 2 and efrac == 2 and grains == 1:
   sed_type = 'Double composite AGN and free electron fraction = 98% with dust grains. No interpolation.'
   bin = 78
   file_lib = 'C17_agn_IR_efrac098.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   print ('a(OX) and U calculation using AGN models with 98% free electrons with dust grains')

elif param  == 2 and efrac == 1 and grains == 1:
   sed_type = 'Double composite AGN and free electron fraction = 2% with dust grains. No interpolation.'
   bin = 78
   file_lib = 'C17_agn_IR_efrac002.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   print ('a(OX) and U calculation using AGN models with 2% free electrons with dust grains')


elif param  == 2 and efrac == 3 and grains == 1:
   sed_type = 'Double composite AGN and free electron fraction = 99.9% with dust grains. No interpolation.'
   bin = 78
   file_lib = 'C17_agn_IR_efrac0999.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   print ('a(OX) and U calculation using AGN models with 99.9% free electrons with dust grains')

elif param  == 2 and efrac == 2 and grains == 2:
   sed_type = 'Double composite AGN and free electron fraction = 98% without dust grains. No interpolation.'
   bin = 78
   file_lib = 'C17_agn_IR_efrac098_nograins.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   print ('a(OX) and U calculation using AGN models with 98% free electrons without dust grains')


elif param  == 2 and efrac == 1 and grains == 2:
   sed_type = 'Double composite AGN and free electron fraction = 2% with dust grains. No interpolation.'
   bin = 78
   file_lib = 'C17_agn_IR_efrac002_nograins.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   print ('a(OX) and U calculation using AGN models with 2% free electrons without dust grains')

elif param  == 2 and efrac == 3 and grains == 2:
   sed_type = 'Double composite AGN and free electron fraction = 99.9% without dust grains. No interpolation.'
   bin = 78
   file_lib = 'C17_agn_IR_efrac0999_nograins.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   print ('a(OX) and U calculation using AGN models with 99.9% free electrons without dust grains')
elif param  == 2 and efrac == 2 and grains == 2:
   sed_type = 'Double composite AGN and free electron fraction = 98% without dust grains. No interpolation.'
   bin = 78
   file_lib = 'C17_agn_IR_efrac098_nograins.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff-IR/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff-IR/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   print ('a(OX) and U calculation using AGN models with 99.9% free electrons without dust grains')


############################################
###### SUMMARY OF THE GRID OF MODELS ######
###########################################

#Final grid: no limitations
grid = grid_aux
print ('-------------------------------------------------')
print ('Summary of the models')
print ('---------------------')
print(' ')
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
Label_SIV = False
Label_eSIV = False
Label_NeII = False
Label_eNeII = False
Label_NeV_14m = False
Label_eNeV_14m = False
Label_NeV_24m = False
Label_eNeV_24m = False
Label_NeIII = False
Label_eNeIII = False
Label_SIII_18m = False
Label_eSIII_18m = False
Label_OIV = False
Label_eOIV = False
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
Label_ArII = False
Label_eArII = False
Label_ArIII = False
Label_eArIII = False
Label_ArV_8m = False
Label_eArV_8m = False
Label_ArV_13m = False
Label_eArV_13m = False


for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == '12logOH':
      Label_OH = True
   if input1.dtype.names[col] == 'e12logOH':
      Label_eOH = True
   if input1.dtype.names[col] == 'SIV_10m':
      Label_SIV = True
   if input1.dtype.names[col] == 'eSIV_10m':
      Label_eSIV = True
   if input1.dtype.names[col] == 'NeII_12m':
      Label_NeII = True
   if input1.dtype.names[col] == 'eNeII_12m':
      Label_eNeII = True
   if input1.dtype.names[col] == 'NeV_14m':
      Label_NeV_14m = True
   if input1.dtype.names[col] == 'eNeV_14m':
      Label_eNeV_14m = True
   if input1.dtype.names[col] == 'NeV_24m':
      Label_NeV_24m = True
   if input1.dtype.names[col] == 'eNeV_24m':
      Label_eNeV_24m = True
   if input1.dtype.names[col] == 'NeIII_15m':
      Label_NeIII = True
   if input1.dtype.names[col] == 'eNeIII_15m':
      Label_eNeIII = True
   if input1.dtype.names[col] == 'SIII_18m':
      Label_SIII_18m = True
   if input1.dtype.names[col] == 'eSIII_18m':
      Label_eSIII_18m = True
   if input1.dtype.names[col] == 'OIV_26m':
      Label_OIV = True
   if input1.dtype.names[col] == 'eOIV_26m':
      Label_eOIV = True
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
   if input1.dtype.names[col] == 'ArII_7m':
      Label_ArII = True
   if input1.dtype.names[col] == 'eArII_7m':
      Label_eArII = True
   if input1.dtype.names[col] == 'ArIII_9m':
      Label_ArIII = True
   if input1.dtype.names[col] == 'eArIII_9m':
      Label_eArIII = True
   if input1.dtype.names[col] == 'ArV_8m':
      Label_ArV_8m = True
   if input1.dtype.names[col] == 'eArV_8m':
      Label_eArV_8m = True
   if input1.dtype.names[col] == 'ArV_13m':
      Label_ArV_13m = True
   if input1.dtype.names[col] == 'eArV_13m':
      Label_eArV_13m = True




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
if Label_SIV == False:
   SIV_10m = np.zeros(input1.size)
else:
   SIV_10m = input1['SIV_10m']
if Label_eSIV == False:
   eSIV_10m = np.zeros(input1.size)
else:
   eSIV_10m = input1['eSIV_10m']
if Label_NeII == False:
   NeII_12m = np.zeros(input1.size)
else:
   NeII_12m = input1['NeII_12m']
if Label_eNeII == False:
   eNeII_12m = np.zeros(input1.size)
else:
   eNeII_12m = input1['eNeII_12m']
if Label_NeV_14m == False:
   NeV_14m = np.zeros(input1.size)
else:
   NeV_14m = input1['NeV_14m']
if Label_eNeV_14m == False:
   eNeV_14m = np.zeros(input1.size)
else:
   eNeV_14m = input1['eNeV_14m']
if Label_NeV_24m == False:
   NeV_24m = np.zeros(input1.size)
else:
   NeV_24m = input1['NeV_24m']
if Label_eNeV_24m == False:
   eNeV_24m = np.zeros(input1.size)
else:
   eNeV_24m = input1['eNeV_24m']
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
if Label_OIV == False:
   OIV_26m = np.zeros(input1.size)
else:
   OIV_26m = input1['OIV_26m']
if Label_eOIV == False:
   eOIV_26m = np.zeros(input1.size)
else:
   eOIV_26m = input1['eOIV_26m']
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
if Label_ArII == False:
   ArII_7m = np.zeros(input1.size)
else:
   ArII_7m = input1['ArII_7m']
if Label_eArII == False:
   eArII_7m = np.zeros(input1.size)
else:
   eArII_7m = input1['eArII_7m']
if Label_ArIII == False:
   ArIII_9m = np.zeros(input1.size)
else:
   ArIII_9m = input1['ArIII_9m']
if Label_eArIII == False:
   eArIII_9m = np.zeros(input1.size)
else:
   eArIII_9m = input1['eArIII_9m']
if Label_ArV_8m == False:
   ArV_8m = np.zeros(input1.size)
else:
   ArV_8m = input1['ArV_8m']
if Label_eArV_8m == False:
   eArV_8m = np.zeros(input1.size)
else:
   eArV_8m = input1['eArV_8m']
if Label_ArV_13m == False:
   ArV_13m = np.zeros(input1.size)
else:
   ArV_13m = input1['ArV_13m']
if Label_eArV_13m == False:
   eArV_13m = np.zeros(input1.size)
else:
   eArV_13m = input1['eArV_13m']

################################################################
###### OUTPUT FORMAT AND INFORMATION: ONLY EMISSION LINES ######
################################################################

#Creation of output only with information from inputs
aux_list = []
aux_list.append(('ID','U12'))
if Label_ArII == True:
   aux_list.append(('ArII_7m', float))
if Label_eArII == True:
   aux_list.append(('eArII_7m', float))
if Label_ArV_8m == True:
   aux_list.append(('ArV_8m', float))
if Label_eArV_8m == True:
   aux_list.append(('eArV_8m', float))
if Label_ArIII == True:
   aux_list.append(('ArIII_9m', float))
if Label_eArIII == True:
   aux_list.append(('eArIII_9m', float))
if Label_SIV == True:
   aux_list.append(('SIV_10m', float))
if Label_eSIV == True:
   aux_list.append(('eSIV_10m', float))
if Label_NeII == True:
   aux_list.append(('NeII_12m', float))
if Label_eNeII == True:
   aux_list.append(('eNeII_12m', float))
if Label_ArV_13m == True:
   aux_list.append(('ArV_13m', float))
if Label_eArV_13m == True:
   aux_list.append(('eArV_13m', float))
if Label_NeV_14m == True:
   aux_list.append(('NeV_14m', float))
if Label_eNeV_14m == True:
   aux_list.append(('eNeV_14m', float))
if Label_NeIII == True:
   aux_list.append(('NeIII_15m', float))
if Label_eNeIII == True:
   aux_list.append(('eNeIII_15m', float))
if Label_SIII_18m == True:
   aux_list.append(('SIII_18m', float))
if Label_eSIII_18m == True:
   aux_list.append(('eSIII_18m', float))
if Label_NeV_24m == True:
   aux_list.append(('NeV_24m', float))
if Label_eNeV_24m == True:
   aux_list.append(('eNeV_24m', float))
if Label_OIV == True:
   aux_list.append(('OIV_26m', float))
if Label_eOIV == True:
   aux_list.append(('eOIV_26m', float))
if Label_SIII_33m == True:
   aux_list.append(('SIII_33m', float))
if Label_eSIII_33m == True:
   aux_list.append(('eSIII_33m', float))
if Label_OIII_52m == True:
   aux_list.append(('OIII_52m', float))
if Label_eOIII_52m == True:
   aux_list.append(('eOIII_52m', float))
if Label_NIII == True:
   aux_list.append(('NIII_57m', float))
if Label_eNIII == True:
   aux_list.append(('eNIII_57m', float))
if Label_OIII_88m == True:
   aux_list.append(('OIII_88m', float))
if Label_eOIII_88m == True:
   aux_list.append(('eOIII_88m', float))
if Label_NII_122m == True:
   aux_list.append(('NII_122m', float))
if Label_eNII_122m == True:
   aux_list.append(('eNII_122m', float))
if Label_NII_205m == True:
   aux_list.append(('NII_205m', float))
if Label_eNII_205m == True:
   aux_list.append(('eNII_205m', float))


aux_list.append(('OH', float))
aux_list.append(('eOH', float))
aux_list.append(('Teff', float))
aux_list.append(('eTeff', float))
aux_list.append(('logU', float))
aux_list.append(('elogU', float))
output = np.zeros(input1.size, dtype=aux_list)



output['ID'] = Names
if Label_ArII == True:
   output['ArII_7m'] = ArII_7m
if Label_eArII == True:
   output['eArII_7m'] = eArII_7m
if Label_ArV_8m == True:
   output['ArV_8m'] = ArV_8m
if Label_eArV_8m == True:
   output['eArV_8m'] = eArV_8m
if Label_ArIII == True:
   output['ArIII_9m'] = ArIII_9m
if Label_eArIII == True:
   output['eArIII_9m'] = eArIII_9m
if Label_SIV == True:
   output['SIV_10m'] = SIV_10m
if Label_eSIV == True:
   output['eSIV_10m'] = eSIV_10m
if Label_NeII == True:
   output['NeII_12m'] = NeII_12m
if Label_eNeII == True:
   output['eNeII_12m'] = eNeII_12m
if Label_ArV_13m == True:
   output['ArV_13m'] = ArV_13m
if Label_eArV_13m == True:
   output['eArV_13m'] = eArV_13m
if Label_NeV_14m == True:
   output['NeV_14m'] = NeV_14m
if Label_eNeV_14m == True:
   output['eNeV_14m'] = eNeV_14m
if Label_NeIII == True:
   output['NeIII_15m'] = NeIII_15m
if Label_eNeIII == True:
   output['eNeIII_15m'] = eNeIII_15m
if Label_SIII_18m == True:
   output['SIII_18m'] = SIII_18m
if Label_eSIII_18m == True:
   output['eSIII_18m'] = eSIII_18m
if Label_NeV_24m == True:
   output['NeV_24m'] = NeV_24m
if Label_eNeV_24m == True:
   output['eNeV_24m'] = eNeV_24m
if Label_OIV == True:
   output['OIV_26m'] = OIV_26m
if Label_eOIV == True:
   output['eOIV_26m'] = eOIV_26m
if Label_SIII_33m == True:
   output['SIII_33m'] = SIII_33m
if Label_eSIII_33m == True:
   output['eSIII_33m'] = eSIII_33m
if Label_OIII_52m == True:
   output['OIII_52m'] = OIII_52m
if Label_eOIII_52m == True:
   output['eOIII_52m'] = eOIII_52m
if Label_NIII == True:
   output['NIII_57m'] = NIII_57m
if Label_eNIII == True:
   output['eNIII_57m'] = eNIII_57m
if Label_OIII_88m == True:
   output['OIII_88m'] = OIII_88m
if Label_eOIII_88m == True:
   output['eOIII_88m'] = eOIII_88m
if Label_NII_122m == True:
   output['NII_122m'] = NII_122m
if Label_eNII_122m == True:
   output['eNII_122m'] = eNII_122m
if Label_NII_205m == True:
   output['NII_205m'] = NII_205m
if Label_eNII_205m == True:
   output['eNII_205m'] = eNII_205m

################################################
###### ESTIMATIONS OF PHYSICAL PARAMETERS ######
################################################

print ('Reading grids ....')
print ('')
print ('')
print ('---------------------------------------')
if param == 1:
   print( '(%)   ID.   12+log(O/H)  T_eff(K)    log(U)')
elif param == 2:
   print( '(%)   ID     12+log(O/H)  a_ox       log(U)')
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

         ArII_7m_obs = 0
         if ArII_7m[tab] > 0:
            while ArII_7m_obs <= 0:
               ArII_7m_obs = np.random.normal(ArII_7m[tab],eArII_7m[tab]+1e-5)
         ArV_8m_obs = 0
         if ArV_8m[tab] > 0:
            while ArV_8m_obs <= 0:
               ArV_8m_obs = np.random.normal(ArV_8m[tab],eArV_8m[tab]+1e-5)
         ArIII_9m_obs = 0
         if ArIII_9m[tab] > 0:
            while ArIII_9m_obs <= 0:
               ArIII_9m_obs = np.random.normal(ArIII_9m[tab],eArIII_9m[tab]+1e-5)
         SIV_10m_obs = 0
         if SIV_10m[tab] > 0:
            while SIV_10m_obs <= 0:
               SIV_10m_obs = np.random.normal(SIV_10m[tab],eSIV_10m[tab]+1e-5)
         NeII_12m_obs = 0
         if NeII_12m[tab] > 0:
            while NeII_12m_obs <= 0:
               NeII_12m_obs = np.random.normal(NeII_12m[tab],eNeII_12m[tab]+1e-3)
         ArV_13m_obs = 0
         if ArV_13m[tab] > 0:
            while ArV_13m_obs <= 0:
               ArV_13m_obs = np.random.normal(ArV_13m[tab],eArV_13m[tab]+1e-5)
         NeV_14m_obs = 0
         if NeV_14m[tab] > 0:
            while NeV_14m_obs <= 0:
               NeV_14m_obs = np.random.normal(NeV_14m[tab],eNeV_14m[tab]+1e-3)
         NeV_24m_obs = 0
         if NeV_24m[tab] > 0:
            while NeV_24m_obs <= 0:
               NeV_24m_obs = np.random.normal(NeV_24m[tab],eNeV_24m[tab]+1e-3)
         NeIII_15m_obs = 0
         if NeIII_15m[tab] > 0:
            while NeIII_15m_obs <= 0:
               NeIII_15m_obs = np.random.normal(NeIII_15m[tab],eNeIII_15m[tab]+1e-3)
         SIII_18m_obs = 0
         if SIII_18m[tab] > 0:
            while SIII_18m_obs <= 0:
               SIII_18m_obs = np.random.normal(SIII_18m[tab],eSIII_18m[tab]+1e-3)
         OIV_26m_obs = 0
         if OIV_26m[tab] > 0:
            while OIV_26m_obs <= 0:
               OIV_26m_obs = np.random.normal(OIV_26m[tab],eOIV_26m[tab]+1e-3)
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
         if ArII_7m_obs == 0 or ArIII_9m_obs == 0:
            Ar2Ar3_obs = -10
         else:
            Ar2Ar3_obs = np.log10((ArII_7m_obs / ArIII_9m_obs))
         if ArII_7m_obs == 0 or ArIII_9m_obs == 0 or ArV_8m_obs == 0:
            Ar23Ar5a_obs = -10
         else:
            Ar23Ar5a_obs = np.log10((ArII_7m_obs + ArIII_9m_obs)/ArV_8m_obs)
         if ArII_7m_obs == 0 or ArIII_9m_obs == 0 or ArV_13m_obs == 0:
            Ar23Ar5b_obs = -10
         else:
            Ar23Ar5b_obs = np.log10((ArII_7m_obs + ArIII_9m_obs)/ArV_13m_obs)
         if Ar23Ar5a_obs > -10 or Ar23Ar5b_obs > -10:
            Ar23Ar5_obs = 0
         else:
            Ar23Ar5_obs = -10
         if ArIII_9m_obs == 0 or ArV_8m_obs == 0:
            Ar3Ar5a_obs = -10
         else:
            Ar3Ar5a_obs = np.log10(ArIII_9m_obs/ArV_8m_obs)
         if ArIII_9m_obs == 0 or ArV_13m_obs == 0:
            Ar3Ar5b_obs = -10
         else:
            Ar3Ar5b_obs = np.log10(ArIII_9m_obs/ArV_13m_obs)
         if Ar3Ar5a_obs > -10 or Ar3Ar5b_obs > -10:
            Ar3Ar5_obs = 0
         else:
            Ar3Ar5_obs = -10
         if NeII_12m_obs == 0 or NeIII_15m_obs == 0:
            Ne2Ne3_obs = -10
         else:
            Ne2Ne3_obs = np.log10((NeII_12m_obs / NeIII_15m_obs))
         if NeII_12m_obs == 0 or NeIII_15m_obs == 0 or NeV_14m_obs == 0:
            Ne23Ne5a_obs = -10
         else:
            Ne23Ne5a_obs = np.log10((NeII_12m_obs + NeIII_15m_obs)/NeV_14m_obs)
         if NeII_12m_obs == 0 or NeIII_15m_obs == 0 or NeV_24m_obs == 0:
            Ne23Ne5b_obs = -10
         else:
            Ne23Ne5b_obs = np.log10((NeII_12m_obs + NeIII_15m_obs)/NeV_24m_obs)
         if Ne23Ne5a_obs > -10 or Ne23Ne5b_obs > -10:
            Ne23Ne5_obs = 0
         else:
            Ne23Ne5_obs = -10
         if NeIII_15m_obs == 0 or NeV_14m_obs == 0:
            Ne3Ne5a_obs = -10
         else:
            Ne3Ne5a_obs = np.log10(NeIII_15m_obs/NeV_14m_obs)
         if NeIII_15m_obs == 0 or NeV_24m_obs == 0:
            Ne3Ne5b_obs = -10
         else:
            Ne3Ne5b_obs = np.log10(NeIII_15m_obs/NeV_24m_obs)
         if Ne3Ne5a_obs > -10 or Ne3Ne5b_obs > -10:
            Ne3Ne5_obs = 0
         else:
            Ne3Ne5_obs = -10
         if SIV_10m_obs == 0  or SIII_18m_obs == 0:
            S3S4a_obs = -10
         else:
            S3S4a_obs = np.log10((SIII_18m_obs / SIV_10m_obs))
         if SIV_10m_obs == 0  or SIII_33m_obs == 0:
            S3S4b_obs = -10
         else:
            S3S4b_obs = np.log10(SIII_33m_obs / SIV_10m_obs)
         if S3S4a_obs > -10 or S3S4b_obs > -10:
            S3S4_obs = 0
         else:
            S3S4_obs = -10
         if NII_122m_obs == 0 or NIII_57m_obs == 0:
            N2N3a_obs = -10
         else:
            N2N3a_obs = np.log10((NII_122m_obs / NIII_57m_obs))
         if OIV_26m_obs == 0 or OIII_52m_obs == 0:
            O3O4a_obs = -10
         else:
            O3O4a_obs = np.log10((OIII_52m_obs/OIV_26m_obs))
         if OIV_26m_obs == 0 or OIII_88m_obs == 0:
            O3O4b_obs = -10
         else:
            O3O4b_obs = np.log10((OIII_88m_obs/OIV_26m_obs))
         if O3O4a_obs > -10 or O3O4b_obs > -10:
            O3O4_obs = 0
         else:
            O3O4_obs = -10
         if OIV_26m_obs == 0 or SIV_10m_obs == 0 or SIII_18m_obs == 0:
            S34O4a_obs = -10
         else:
            S34O4a_obs = np.log10((SIII_18m_obs+SIV_10m_obs)/OIV_26m_obs)
         if OIV_26m_obs == 0 or SIV_10m_obs == 0 or SIII_33m_obs == 0:
            S34O4b_obs = -10
         else:
            S34O4b_obs = np.log10((SIII_33m_obs+SIV_10m_obs)/OIV_26m_obs)
         if S34O4a_obs > -10 or S34O4b_obs > -10:
            S34O4_obs = 0
         else:
            S34O4_obs = -10
         if NII_205m_obs == 0 or NIII_57m_obs == 0:
            N2N3b_obs = -10
         else:
            N2N3b_obs = np.log10((NII_205m_obs / NIII_57m_obs))
         if N2N3a_obs > -10 or N2N3b_obs > -10:
            N2N3_obs = 0
         else:
            N2N3_obs = -10


         lista_param_sf = [Ar2Ar3_obs, Ne2Ne3_obs, S3S4_obs, N2N3_obs]
         lista_param_agn1 = [Ar3Ar5_obs, Ne3Ne5_obs, O3O4_obs]
         lista_param_agn2 = [Ar3Ar5_obs, Ne3Ne5_obs, S34O4_obs]
         counter_sf = sum(1 for parametro in lista_param_sf if parametro > -10)
         counter_agn1 = sum(1 for parametro in lista_param_agn1 if parametro > -10)
         counter_agn2 = sum(1 for parametro in lista_param_agn2 if parametro > -10)




   # Interpolation of grid at specific O/H


         if logOH[tab] > 0:
            OH = np.random.normal(logOH[tab],elogOH[tab]+1e-3)
            OH_mc.append(OH)


            grid_T0 = []
            if param == 1:
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
            elif param == 2:
               if OH <= 8.1:
                  OH = 8.1
                  i0 = 0
                  i1 = bin
               elif OH >= 8.1 and OH < 8.4:
                  i0 = 0
                  i1 = bin
               elif OH >= 8.4 and OH < 8.7:
                  i0 = bin
                  i1 = 2*bin
               elif OH >= 8.7 and OH < 9.0:
                  i0 = 2*bin
                  i1 = 3*bin
               elif OH >= 9.0:
                  OH = 9.0
                  i0 = 2*bin
                  i1 = 3*bin




            for x in range(0,bin):
               for y in grid.dtype.names:
                  grid_T0.append(grid[y][i0+x]*np.abs(0.3-OH+grid['12logOH'][i0])/0.3+grid[y][i1+x]*np.abs(0.3-grid['12logOH'][i1]+OH)/0.3)
                  #grid_T0.append(grid[i0+x,y]*np.abs(0.3-OH+grid[i0,0])/0.3+grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)
               
         #         grid_T0.append(grid[i0+x,y]*np.abs(0.3-grid[i0,0]+OH)/0.3 + grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)


            grid_T_aux = np.reshape(grid_T0,(bin,20))
            grid_T = np.zeros(grid_T_aux.shape[0], dtype=grid.dtype)
            for col_n in range(0, len(grid.dtype.names)):
               grid_T[grid.dtype.names[col_n]] = grid_T_aux[:, col_n]

         else:
            OH = 0
            OH_mc.append(OH)
            grid_T = grid
      

#         np.savetxt('int_models.dat',grid_T,fmt='%.2f')

   # Calculation of T and log U

         if (param == 1 and counter_sf < 2) or (param == 2 and counter_agn1 < 2 and counter_agn2 < 2):
            Teff = 0
            logU = 0
         else:
            CHI_Ar2Ar3 = 0
            CHI_Ar23Ar5 = 0
            CHI_Ar23Ar5a = 0
            CHI_Ar23Ar5b = 0
            CHI_Ar3Ar5 = 0
            CHI_Ar3Ar5a = 0
            CHI_Ar3Ar5b = 0
            CHI_Ne2Ne3 = 0
            CHI_Ne23Ne5 = 0
            CHI_Ne23Ne5a = 0
            CHI_Ne23Ne5b = 0
            CHI_Ne3Ne5 = 0
            CHI_Ne3Ne5a = 0
            CHI_Ne3Ne5b = 0
            CHI_S3S4 = 0
            CHI_S3S4a = 0
            CHI_S3S4b = 0
            CHI_O3O4 = 0
            CHI_O3O4a = 0
            CHI_O3O4b = 0
            CHI_S34O4 = 0
            CHI_S34O4a = 0
            CHI_S34O4b = 0
            CHI_N2N3 = 0
            CHI_N2N3a = 0
            CHI_N2N3b = 0
            for index in grid_T:
               if S3S4a_obs == -10: 
                  CHI_S3S4a = 0
               elif index['SIV_10m'] == 0 or index['SIII_18m'] == 0 : 
                  CHI_S3S4a = tol_max
               else:   
                  CHI_S3S4a = (np.log10(index['SIII_18m']/index['SIV_10m'])- S3S4a_obs)**2/np.log10((index['SIII_18m']/index['SIV_10m']))
               if S3S4b_obs == -10: 
                  CHI_S3S4b = 0
               elif index['SIV_10m'] == 0 or index['SIII_33m'] == 0 : 
                  CHI_S3S4b = tol_max
               else:   
                  CHI_S3S4b = (np.log10(index['SIII_33m']/index['SIV_10m'])- S3S4b_obs)**2/np.log10((index['SIII_33m']/index['SIV_10m']))
               if CHI_S3S4a == 0 and CHI_S3S4b == 0:
                  CHI_S3S4 = 0
               elif CHI_S3S4a > 0 and CHI_S3S4b == 0:
                  CHI_S3S4 = CHI_S3S4a
               elif CHI_S3S4a == 0 and CHI_S3S4b > 0:
                  CHI_S3S4 = CHI_S3S4b
               else:
                  CHI_S3S4 = (CHI_S3S4a + CHI_S3S4b)/2 
               if Ar2Ar3_obs == -10: 
                  CHI_Ar2Ar3 = 0
               elif index['ArII_7m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar2Ar3 = tol_max
               else:   
                  CHI_Ar2Ar3 = (np.log10(index['ArII_7m']/index['ArIII_9m'])- Ar2Ar3_obs)**2/np.log10((index['ArII_7m']/index['ArIII_9m']))
               if Ar23Ar5a_obs == -10: 
                  CHI_Ar23Ar5a = 0
               elif index['ArII_7m'] == 0 or index['ArV_8m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar23Ar5a = tol_max
               else:   
                  CHI_Ar23Ar5a = (np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_8m'])- Ar23Ar5a_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_8m'])
               if Ar23Ar5b_obs == -10: 
                  CHI_Ar23Ar5b = 0
               elif index['ArII_7m'] == 0 or index['ArV_13m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar23Ar5b = tol_max
               else:   
                  CHI_Ar23Ar5b = (np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_13m'])- Ar23Ar5b_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_13m'])
               if CHI_Ar23Ar5a == 0 and CHI_Ar23Ar5b == 0:
                  CHI_Ar23Ar5 = 0
               elif CHI_Ar23Ar5a > 0 and CHI_Ar23Ar5b == 0:
                  CHI_Ar23Ar5 = CHI_Ar23Ar5a
               elif CHI_Ar23Ar5a == 0 and CHI_Ar23Ar5b > 0:
                  CHI_Ar23Ar5 = CHI_Ar23Ar5b
               else:
                  CHI_Ar23Ar5 = (CHI_Ar23Ar5a + CHI_Ar23Ar5b)/2 
               if Ar3Ar5a_obs == -10: 
                  CHI_Ar3Ar5a = 0
               elif index['ArV_8m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar3Ar5a = tol_max
               else:   
                  CHI_Ar3Ar5a = (np.log10((index['ArIII_9m'])/index['ArV_8m'])- Ar3Ar5a_obs)**2/np.log10((index['ArIII_9m'])/index['ArV_8m'])
               if Ar3Ar5b_obs == -10: 
                  CHI_Ar3Ar5b = 0
               elif index['ArV_13m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar3Ar5b = tol_max
               else:   
                  CHI_Ar3Ar5b = (np.log10((index['ArIII_9m'])/index['ArV_13m'])- Ar3Ar5b_obs)**2/np.log10((index['ArIII_9m'])/index['ArV_13m'])
               if CHI_Ar3Ar5a == 0 and CHI_Ar3Ar5b == 0:
                  CHI_Ar3Ar5 = 0
               elif CHI_Ar3Ar5a > 0 and CHI_Ar3Ar5b == 0:
                  CHI_Ar3Ar5 = CHI_Ar3Ar5a
               elif CHI_Ar3Ar5a == 0 and CHI_Ar3Ar5b > 0:
                  CHI_Ar3Ar5 = CHI_Ar3Ar5b
               else:
                  CHI_Ar3Ar5 = (CHI_Ar3Ar5a + CHI_Ar3Ar5b)/2 
               if Ne2Ne3_obs == -10: 
                  CHI_Ne2Ne3 = 0
               elif index['NeII_12m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne2Ne3 = tol_max
               else:   
                  CHI_Ne2Ne3 = (np.log10(index['NeII_12m']/index['NeIII_15m'])- Ne2Ne3_obs)**2/np.log10((index['NeII_12m']/index['NeIII_15m']))
               if Ne23Ne5a_obs == -10: 
                  CHI_Ne23Ne5a = 0
               elif index['NeII_12m'] == 0 or index['NeV_14m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne23Ne5a = tol_max
               else:   
                  CHI_Ne23Ne5a = (np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_14m'])- Ne23Ne5a_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_14m'])
               if Ne23Ne5b_obs == -10: 
                  CHI_Ne23Ne5b = 0
               elif index['NeII_12m'] == 0 or index['NeV_24m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne23Ne5b = tol_max
               else:   
                  CHI_Ne23Ne5b = (np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_24m'])- Ne23Ne5b_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_24m'])
               if CHI_Ne23Ne5a == 0 and CHI_Ne23Ne5b == 0:
                  CHI_Ne23Ne5 = 0
               elif CHI_Ne23Ne5a > 0 and CHI_Ne23Ne5b == 0:
                  CHI_Ne23Ne5 = CHI_Ne23Ne5a
               elif CHI_Ne23Ne5a == 0 and CHI_Ne23Ne5b > 0:
                  CHI_Ne23Ne5 = CHI_Ne23Ne5b
               else:
                  CHI_Ne23Ne5 = (CHI_Ne23Ne5a + CHI_Ne23Ne5b)/2 
               if Ne3Ne5a_obs == -10: 
                  CHI_Ne3Ne5a = 0
               elif index['NeV_14m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne3Ne5a = tol_max
               else:   
                  CHI_Ne3Ne5a = (np.log10((index['NeIII_15m'])/index['NeV_14m'])- Ne3Ne5a_obs)**2/np.log10((index['NeIII_15m'])/index['NeV_14m'])
               if Ne3Ne5b_obs == -10: 
                  CHI_Ne3Ne5b = 0
               elif index['NeV_24m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne3Ne5b = tol_max
               else:   
                  CHI_Ne3Ne5b = (np.log10((index['NeIII_15m'])/index['NeV_24m'])- Ne3Ne5b_obs)**2/np.log10((index['NeIII_15m'])/index['NeV_24m'])
               if CHI_Ne3Ne5a == 0 and CHI_Ne3Ne5b == 0:
                  CHI_Ne3Ne5 = 0
               elif CHI_Ne3Ne5a > 0 and CHI_Ne3Ne5b == 0:
                  CHI_Ne3Ne5 = CHI_Ne3Ne5a
               elif CHI_Ne3Ne5a == 0 and CHI_Ne3Ne5b > 0:
                  CHI_Ne3Ne5 = CHI_Ne3Ne5b
               else:
                  CHI_Ne3Ne5 = (CHI_Ne3Ne5a + CHI_Ne3Ne5b)/2 
               if O3O4a_obs == -10: 
                  CHI_O3O4a = 0
               elif index['OIV_26m'] == 0 or index['OIII_52m'] == 0: 
                  CHI_O3O4a = tol_max
               else:   
                  CHI_O3O4a = (np.log10(index['OIII_52m']/index['OIV_26m'])- O3O4a_obs)**2/np.log10((index['OIII_52m']/index['OIV_26m']))
               if O3O4b_obs == -10: 
                  CHI_O3O4b = 0
               elif index['OIV_26m'] == 0 or index['OIII_88m'] == 0: 
                  CHI_O3O4b = tol_max
               else:   
                  CHI_O3O4b = (np.log10(index['OIII_88m']/index['OIV_26m'])- O3O4b_obs)**2/np.log10((index['OIII_88m']/index['OIV_26m']))
               if CHI_O3O4a == 0 and CHI_O3O4b == 0:
                  CHI_O3O4 = 0
               elif CHI_O3O4a > 0 and CHI_O3O4b == 0:
                  CHI_O3O4 = CHI_O3O4a
               elif CHI_O3O4a == 0 and CHI_O3O4b > 0:
                  CHI_O3O4 = CHI_O3O4b
               else:
                  CHI_O3O4 = (CHI_O3O4a + CHI_O3O4b)/2 
               if S34O4a_obs == -10: 
                  print('a')
                  CHI_S34O4a = 0
               elif index['OIV_26m'] == 0 or index['SIII_18m'] == 0 or index['SIV_10m'] == 0: 
                  CHI_S34O4a = tol_max
               else:   
                  CHI_S34O4a = (np.log10((index['SIII_18m']+index['SIV_10m'])/index['OIV_26m'])- S34O4a_obs)**2/np.log10(((index['SIII_18m']+index['SIV_10m'])/index['OIV_26m']))
               if S34O4b_obs == -10: 
                  CHI_S34O4b = 0
               elif index['OIV_26m'] == 0 or index['SIII_33m'] == 0 or index['SIV_10m'] == 0: 
                  CHI_S34O4b = tol_max
               else:   
                  CHI_S34O4b = (np.log10((index['SIII_33m']+index['SIV_10m'])/index['OIV_26m'])- S34O4b_obs)**2/np.log10((index['SIII_33m']+index['SIV_10m'])/index['OIV_26m'])
               if CHI_S34O4a == 0 and CHI_S34O4b == 0:
                  CHI_S34O4 = 0
               elif CHI_S34O4a > 0 and CHI_S34O4b == 0:
                  CHI_S34O4 = CHI_S34O4a
               elif CHI_S34O4a == 0 and CHI_S34O4b > 0:
                  CHI_S34O4 = CHI_S34O4b
               else:
                  CHI_S34O4 = (CHI_S34O4a + CHI_S34O4b)/2 
               if N2N3a_obs == -10: 
                  CHI_N2N3a = 0
               elif index['NIII_57m'] == 0 or index['NII_122m'] == 0: 
                  CHI_N2N3a = tol_max
               else:   
                  CHI_N2N3a = (np.log10(index['NII_122m']/index['NIII_57m'])- N2N3a_obs)**2/np.log10((index['NII_122m']/index['NIII_57m']))
               if N2N3b_obs == -10: 
                  CHI_N2N3b = 0
               elif index['NIII_57m'] == 0 or index['NII_205m'] == 0: 
                  CHI_N2N3b = tol_max
               else:   
                  CHI_N2N3b = (np.log10(index['NII_205m']/index['NIII_57m'])- N2N3b_obs)**2/np.log10((index['NII_205m']/index['NIII_57m']))
               if CHI_N2N3a == 0 and CHI_N2N3b == 0:
                  CHI_N2N3 = 0
               elif CHI_N2N3a > 0 and CHI_N2N3b == 0:
                  CHI_N2N3 = CHI_N2N3a
               elif CHI_N2N3a == 0 and CHI_N2N3b > 0:
                  CHI_N2N3 = CHI_N2N3b
               else:
                  CHI_N2N3 = (CHI_N2N3a + CHI_N2N3b)/2 


               if param == 1:
                  CHI_Teff = (CHI_S3S4**2 +  CHI_Ne2Ne3**2 + CHI_N2N3**2 + CHI_Ar2Ar3**2)**0.5
               else:
                  if Ar23Ar5_obs > -10 and Ne23Ne5_obs > -10 and O3O4_obs > -10:
                     CHI_Teff = (CHI_Ne23Ne5**2 + CHI_O3O4**2 + CHI_Ar23Ar5**2)**0.5
                  elif Ar23Ar5_obs > -10 and Ne23Ne5_obs > -10 and O3O4_obs == -10:
                     CHI_Teff = (CHI_Ne23Ne5**2 + CHI_S34O4**2 + CHI_Ar23Ar5**2)**0.5
                  elif Ar23Ar5_obs == -10 and Ne23Ne5_obs > -10 and O3O4_obs > -10:
                     CHI_Teff = (CHI_Ne23Ne5**2 + CHI_O3O4**2 + CHI_Ar3Ar5**2)**0.5
                  elif Ar23Ar5_obs == -10 and Ne23Ne5_obs > -10 and O3O4_obs == -10:
                     CHI_Teff = (CHI_Ne23Ne5**2 + CHI_S34O4**2 + CHI_Ar3Ar5**2)**0.5
                  elif Ar23Ar5_obs == -10 and Ne23Ne5_obs == -10 and O3O4_obs > -10:
                     CHI_Teff = (CHI_Ne3Ne5**2 + CHI_O3O4**2 + CHI_Ar3Ar5**2)**0.5
                  elif Ar23Ar5_obs == -10 and Ne23Ne5_obs == -10 and O3O4_obs == -10:
                     CHI_Teff = (CHI_Ne3Ne5**2 + CHI_S34O4**2 + CHI_Ar3Ar5**2)**0.5

               if param == 1:
                  Teff_p = index['T_eff']*(1/CHI_Teff)**2 + Teff_p
               else:
                  Teff_p = index['a_ox']*(1/CHI_Teff)**2 + Teff_p
               logU_p = index['logU'] *(1/CHI_Teff)**2 + logU_p         
               den_Teff = (1/CHI_Teff)**2 + den_Teff
            Teff = Teff_p / den_Teff 
            logU = logU_p / den_Teff

   # Calculation of T and log U errors


         if (param == 1 and counter_sf < 2) or (param == 2 and counter_agn1 < 2 and counter_agn2 < 2):
            eTeff = 0
            elogU = 0
         else:
            CHI_Ar2Ar3 = 0
            CHI_Ar23Ar5 = 0
            CHI_Ar23Ar5a = 0
            CHI_Ar23Ar5b = 0
            CHI_Ar3Ar5 = 0
            CHI_Ar3Ar5a = 0
            CHI_Ar3Ar5b = 0
            CHI_Ne2Ne3 = 0
            CHI_Ne23Ne5 = 0
            CHI_Ne23Ne5a = 0
            CHI_Ne23Ne5b = 0
            CHI_Ne3Ne5 = 0
            CHI_Ne3Ne5a = 0
            CHI_Ne3Ne5b = 0
            CHI_S3S4 = 0
            CHI_S3S4a = 0
            CHI_S3S4b = 0
            CHI_O3O4 = 0
            CHI_O3O4a = 0
            CHI_O3O4b = 0
            CHI_S34O4 = 0
            CHI_S34O4a = 0
            CHI_S34O4b = 0
            CHI_N2N3 = 0
            CHI_N2N3a = 0
            CHI_N2N3b = 0

            for index in grid_T:
               if S3S4a_obs == -10: 
                  CHI_S3S4a = 0
               elif index['SIV_10m'] == 0 or index['SIII_18m'] == 0 : 
                  CHI_S3S4a = tol_max
               else:   
                  CHI_S3S4a = (np.log10(index['SIII_18m']/index['SIV_10m'])- S3S4a_obs)**2/np.log10((index['SIII_18m']/index['SIV_10m']))
               if S3S4b_obs == -10: 
                  CHI_S3S4b = 0
               elif index['SIV_10m'] == 0 or index['SIII_33m'] == 0 : 
                  CHI_S3S4b = tol_max
               else:   
                  CHI_S3S4b = (np.log10(index['SIII_33m']/index['SIV_10m'])- S3S4b_obs)**2/np.log10((index['SIII_33m']/index['SIV_10m']))
               if CHI_S3S4a == 0 and CHI_S3S4b == 0:
                  CHI_S3S4 = 0
               elif CHI_S3S4a > 0 and CHI_S3S4b == 0:
                  CHI_S3S4 = CHI_S3S4a
               elif CHI_S3S4a == 0 and CHI_S3S4b > 0:
                  CHI_S3S4 = CHI_S3S4b
               else:
                  CHI_S3S4 = (CHI_S3S4a + CHI_S3S4b)/2 
               if Ar2Ar3_obs == -10: 
                  CHI_Ar2Ar3 = 0
               elif index['ArII_7m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar2Ar3 = tol_max
               else:   
                  CHI_Ar2Ar3 = (np.log10(index['ArII_7m']/index['ArIII_9m'])- Ar2Ar3_obs)**2/np.log10((index['ArII_7m']/index['ArIII_9m']))
               if Ar23Ar5a_obs == -10: 
                  CHI_Ar23Ar5a = 0
               elif index['ArII_7m'] == 0 or index['ArV_8m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar23Ar5a = tol_max
               else:   
                  CHI_Ar23Ar5a = (np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_8m'])- Ar23Ar5a_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_8m'])
               if Ar23Ar5b_obs == -10: 
                  CHI_Ar23Ar5b = 0
               elif index['ArII_7m'] == 0 or index['ArV_13m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar23Ar5b = tol_max
               else:   
                  CHI_Ar23Ar5b = (np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_13m'])- Ar23Ar5b_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_13m'])
               if CHI_Ar23Ar5a == 0 and CHI_Ar23Ar5b == 0:
                  CHI_Ar23Ar5 = 0
               elif CHI_Ar23Ar5a > 0 and CHI_Ar23Ar5b == 0:
                  CHI_Ar23Ar5 = CHI_Ar23Ar5a
               elif CHI_Ar23Ar5a == 0 and CHI_Ar23Ar5b > 0:
                  CHI_Ar23Ar5 = CHI_Ar23Ar5b
               else:
                  CHI_Ar23Ar5 = (CHI_Ar23Ar5a + CHI_Ar23Ar5b)/2 
               if Ar3Ar5a_obs == -10: 
                  CHI_Ar3Ar5a = 0
               elif index['ArV_8m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar3Ar5a = tol_max
               else:   
                  CHI_Ar3Ar5a = (np.log10((index['ArIII_9m'])/index['ArV_8m'])- Ar3Ar5a_obs)**2/np.log10((index['ArIII_9m'])/index['ArV_8m'])
               if Ar3Ar5b_obs == -10: 
                  CHI_Ar3Ar5b = 0
               elif index['ArV_13m'] == 0 or index['ArIII_9m'] == 0: 
                  CHI_Ar3Ar5b = tol_max
               else:   
                  CHI_Ar3Ar5b = (np.log10((index['ArIII_9m'])/index['ArV_13m'])- Ar3Ar5b_obs)**2/np.log10((index['ArIII_9m'])/index['ArV_13m'])
               if CHI_Ar3Ar5a == 0 and CHI_Ar3Ar5b == 0:
                  CHI_Ar3Ar5 = 0
               elif CHI_Ar3Ar5a > 0 and CHI_Ar3Ar5b == 0:
                  CHI_Ar3Ar5 = CHI_Ar3Ar5a
               elif CHI_Ar3Ar5a == 0 and CHI_Ar3Ar5b > 0:
                  CHI_Ar3Ar5 = CHI_Ar3Ar5b
               else:
                  CHI_Ar3Ar5 = (CHI_Ar3Ar5a + CHI_Ar3Ar5b)/2 
               if Ne2Ne3_obs == -10: 
                  CHI_Ne2Ne3 = 0
               elif index['NeII_12m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne2Ne3 = tol_max
               else:   
                  CHI_Ne2Ne3 = (np.log10(index['NeII_12m']/index['NeIII_15m'])- Ne2Ne3_obs)**2/np.log10((index['NeII_12m']/index['NeIII_15m']))
               if Ne23Ne5a_obs == -10: 
                  CHI_Ne23Ne5a = 0
               elif index['NeII_12m'] == 0 or index['NeV_14m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne23Ne5a = tol_max
               else:   
                  CHI_Ne23Ne5a = (np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_14m'])- Ne23Ne5a_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_14m'])
               if Ne23Ne5b_obs == -10: 
                  CHI_Ne23Ne5b = 0
               elif index['NeII_12m'] == 0 or index['NeV_24m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne23Ne5b = tol_max
               else:   
                  CHI_Ne23Ne5b = (np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_24m'])- Ne23Ne5b_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_24m'])
               if CHI_Ne23Ne5a == 0 and CHI_Ne23Ne5b == 0:
                  CHI_Ne23Ne5 = 0
               elif CHI_Ne23Ne5a > 0 and CHI_Ne23Ne5b == 0:
                  CHI_Ne23Ne5 = CHI_Ne23Ne5a
               elif CHI_Ne23Ne5a == 0 and CHI_Ne23Ne5b > 0:
                  CHI_Ne23Ne5 = CHI_Ne23Ne5b
               else:
                  CHI_Ne23Ne5 = (CHI_Ne23Ne5a + CHI_Ne23Ne5b)/2 
               if Ne3Ne5a_obs == -10: 
                  CHI_Ne3Ne5a = 0
               elif index['NeV_14m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne3Ne5a = tol_max
               else:   
                  CHI_Ne3Ne5a = (np.log10((index['NeIII_15m'])/index['NeV_14m'])- Ne3Ne5a_obs)**2/np.log10((index['NeIII_15m'])/index['NeV_14m'])
               if Ne3Ne5b_obs == -10: 
                  CHI_Ne3Ne5b = 0
               elif index['NeV_24m'] == 0 or index['NeIII_15m'] == 0: 
                  CHI_Ne3Ne5b = tol_max
               else:   
                  CHI_Ne3Ne5b = (np.log10((index['NeIII_15m'])/index['NeV_24m'])- Ne3Ne5b_obs)**2/np.log10((index['NeIII_15m'])/index['NeV_24m'])
               if CHI_Ne3Ne5a == 0 and CHI_Ne3Ne5b == 0:
                  CHI_Ne3Ne5 = 0
               elif CHI_Ne3Ne5a > 0 and CHI_Ne3Ne5b == 0:
                  CHI_Ne3Ne5 = CHI_Ne3Ne5a
               elif CHI_Ne3Ne5a == 0 and CHI_Ne3Ne5b > 0:
                  CHI_Ne3Ne5 = CHI_Ne3Ne5b
               else:
                  CHI_Ne3Ne5 = (CHI_Ne3Ne5a + CHI_Ne3Ne5b)/2 
               if O3O4a_obs == -10: 
                  CHI_O3O4a = 0
               elif index['OIV_26m'] == 0 or index['OIII_52m'] == 0: 
                  CHI_O3O4a = tol_max
               else:   
                  CHI_O3O4a = (np.log10(index['OIII_52m']/index['OIV_26m'])- O3O4a_obs)**2/np.log10((index['OIII_52m']/index['OIV_26m']))
               if O3O4b_obs == -10: 
                  CHI_O3O4b = 0
               elif index['OIV_26m'] == 0 or index['OIII_88m'] == 0: 
                  CHI_O3O4b = tol_max
               else:   
                  CHI_O3O4b = (np.log10(index['OIII_88m']/index['OIV_26m'])- O3O4b_obs)**2/np.log10((index['OIII_88m']/index['OIV_26m']))
               if CHI_O3O4a == 0 and CHI_O3O4b == 0:
                  CHI_O3O4 = 0
               elif CHI_O3O4a > 0 and CHI_O3O4b == 0:
                  CHI_O3O4 = CHI_O3O4a
               elif CHI_O3O4a == 0 and CHI_O3O4b > 0:
                  CHI_O3O4 = CHI_O3O4b
               else:
                  CHI_O3O4 = (CHI_O3O4a + CHI_O3O4b)/2 
               if S34O4a_obs == -10: 
                  CHI_S34O4a = 0
               elif index['OIV_26m'] == 0 or index['SIII_18m'] == 0 or index['SIV_10m'] == 0: 
                  CHI_S34O4a = tol_max
               else:   
                  CHI_S34O4a = (np.log10((index['SIII_18m']+index['SIV_10m'])/index['OIV_26m'])- S34O4a_obs)**2/np.log10(((index['SIII_18m']+index['SIV_10m'])/index['OIV_26m']))
               if S34O4b_obs == -10: 
                  CHI_S34O4b = 0
               elif index['OIV_26m'] == 0 or index['SIII_33m'] == 0 or index['SIV_10m'] == 0: 
                  CHI_S34O4b = tol_max
               else:   
                  CHI_S34O4b = (np.log10((index['SIII_33m']+index['SIV_10m'])/index['OIV_26m'])- S34O4b_obs)**2/np.log10((index['SIII_33m']+index['SIV_10m'])/index['OIV_26m'])
               if CHI_S34O4a == 0 and CHI_S34O4b == 0:
                  CHI_S34O4 = 0
               elif CHI_S34O4a > 0 and CHI_S34O4b == 0:
                  CHI_S34O4 = CHI_S34O4a
               elif CHI_S34O4a == 0 and CHI_S34O4b > 0:
                  CHI_S34O4 = CHI_S34O4b
               else:
                  CHI_S34O4 = (CHI_S34O4a + CHI_S34O4b)/2 
               if N2N3a_obs == -10: 
                  CHI_N2N3a = 0
               elif index['NIII_57m'] == 0 or index['NII_122m'] == 0: 
                  CHI_N2N3a = tol_max
               else:   
                  CHI_N2N3a = (np.log10(index['NII_122m']/index['NIII_57m'])- N2N3a_obs)**2/np.log10((index['NII_122m']/index['NIII_57m']))
               if N2N3b_obs == -10: 
                  CHI_N2N3b = 0
               elif index['NIII_57m'] == 0 or index['NII_205m'] == 0: 
                  CHI_N2N3b = tol_max
               else:   
                  CHI_N2N3b = (np.log10(index['NII_205m']/index['NIII_57m'])- N2N3b_obs)**2/np.log10((index['NII_205m']/index['NIII_57m']))
               if CHI_N2N3a == 0 and CHI_N2N3b == 0:
                  CHI_N2N3 = 0
               elif CHI_N2N3a > 0 and CHI_N2N3b == 0:
                  CHI_N2N3 = CHI_N2N3a
               elif CHI_N2N3a == 0 and CHI_N2N3b > 0:
                  CHI_N2N3 = CHI_N2N3b
               else:
                  CHI_N2N3 = (CHI_N2N3a + CHI_N2N3b)/2 

               if param == 1:
                  CHI_Teff = (CHI_S3S4**2 +  CHI_Ne2Ne3**2 + CHI_N2N3**2 + CHI_Ar2Ar3**2)**0.5
               else:
                  if Ar23Ar5_obs > -10 and Ne23Ne5_obs > -10 and O3O4_obs > -10:
                     CHI_Teff = (CHI_Ne23Ne5**2 + CHI_O3O4**2 + CHI_Ar23Ar5**2)**0.5
                  elif Ar23Ar5_obs > -10 and Ne23Ne5_obs > -10 and O3O4_obs == -10:
                     CHI_Teff = (CHI_Ne23Ne5**2 + CHI_S34O4**2 + CHI_Ar23Ar5**2)**0.5
                  elif Ar23Ar5_obs == -10 and Ne23Ne5_obs > -10 and O3O4_obs > -10:
                     CHI_Teff = (CHI_Ne23Ne5**2 + CHI_O3O4**2 + CHI_Ar3Ar5**2)**0.5
                  elif Ar23Ar5_obs == -10 and Ne23Ne5_obs > -10 and O3O4_obs == -10:
                     CHI_Teff = (CHI_Ne23Ne5**2 + CHI_S34O4**2 + CHI_Ar3Ar5**2)**0.5
                  elif Ar23Ar5_obs == -10 and Ne23Ne5_obs == -10 and O3O4_obs > -10:
                     CHI_Teff = (CHI_Ne3Ne5**2 + CHI_O3O4**2 + CHI_Ar3Ar5**2)**0.5
                  elif Ar23Ar5_obs == -10 and Ne23Ne5_obs == -10 and O3O4_obs == -10:
                     CHI_Teff = (CHI_Ne3Ne5**2 + CHI_S34O4**2 + CHI_Ar3Ar5**2)**0.5


               if param == 1:
                  Teff_e = np.abs(index['T_eff'] - Teff) * (1/CHI_Teff) **2+ Teff_e
               else:
                  Teff_e = np.abs(index['a_ox'] - Teff) * (1/CHI_Teff)**2 + Teff_e
               logU_e = np.abs(index['logU'] - logU) * (1/CHI_Teff)**2 + logU_e         
               den_Teff_e = 1 * (1/CHI_Teff)**2 + den_Teff_e


            eTeff = Teff_e / den_Teff_e 
            elogU = logU_e / den_Teff_e



   #Iterations for the interpolation mode

            if inter == 0:
               Teff = Teff
               logU = logU
            elif inter == 1:
               igrid = grid_T[np.lexsort((grid_T['12logOH'],grid_T['logU']))]
               if param == 2:
                  igrid = interpolate(igrid, 'a_ox', Teff-eTeff-0.02, Teff+eTeff+0.02, 10, param)
                  igrid = igrid[np.lexsort((igrid['12logOH'],igrid['a_ox']))]
               else:
                  igrid = interpolate(igrid,'T_eff',Teff-eTeff-eTeff/10,Teff+eTeff+eTeff/10,10, param)
                  igrid = igrid[np.lexsort((igrid['12logOH'],igrid['T_eff']))]
               igrid = interpolate(igrid,'logU',logU-elogU-0.25,logU+elogU+0.25,10, param)


   #            np.savetxt('int_models.dat',igrid,fmt='%.2f')



               if (param == 1 and counter_sf < 2) or (param == 2 and counter_agn1 < 2 and counter_agn2 < 2):
                  Teff = 0
                  logU = 0
               else:
                  CHI_Ar2Ar3 = 0
                  CHI_Ar23Ar5 = 0
                  CHI_Ar23Ar5a = 0
                  CHI_Ar23Ar5b = 0
                  CHI_Ar3Ar5 = 0
                  CHI_Ar3Ar5a = 0
                  CHI_Ar3Ar5b = 0
                  CHI_Ne2Ne3a = 0
                  CHI_Ne23Ne5a = 0
                  CHI_Ne23Ne5b = 0
                  CHI_Ne23Ne5 = 0
                  CHI_Ne3Ne5 = 0
                  CHI_Ne3Ne5a = 0
                  CHI_Ne3Ne5b = 0
                  CHI_S3S4a = 0
                  CHI_S3S4b = 0
                  CHI_S3S4 = 0
                  CHI_N2N3a = 0
                  CHI_O3O4a = 0
                  CHI_O3O4b = 0
                  CHI_O3O4 = 0
                  CHI_S34O4 = 0
                  CHI_S34O4a = 0
                  CHI_S34O4b = 0
                  CHI_N2N3b = 0
                  CHI_n2N3 = 0

                  for index in igrid:
                     if S3S4a_obs == -10: 
                        CHI_S3S4a = 0
                     elif index['SIV_10m'] == 0 or index['SIII_18m'] == 0 : 
                        CHI_S3S4a = tol_max
                     else:   
                        CHI_S3S4a = (np.log10(index['SIII_18m']/index['SIV_10m'])- S3S4a_obs)**2/np.log10((index['SIII_18m']/index['SIV_10m']))
                     if S3S4b_obs == -10: 
                        CHI_S3S4b = 0
                     elif index['SIV_10m'] == 0 or index['SIII_33m'] == 0 : 
                        CHI_S3S4b = tol_max
                     else:   
                        CHI_S3S4b = (np.log10(index['SIII_33m']/index['SIV_10m'])- S3S4b_obs)**2/np.log10((index['SIII_33m']/index['SIV_10m']))
                     if CHI_S3S4a == 0 and CHI_S3S4b == 0:
                        CHI_S3S4 = 0
                     elif CHI_S3S4a > 0 and CHI_S3S4b == 0:
                        CHI_S3S4 = CHI_S3S4a
                     elif CHI_S3S4a == 0 and CHI_S3S4b > 0:
                        CHI_S3S4 = CHI_S3S4b
                     else:
                        CHI_S3S4 = (CHI_S3S4a + CHI_S3S4b)/2 
                     if Ar2Ar3_obs == -10: 
                        CHI_Ar2Ar3 = 0
                     elif index['ArII_7m'] == 0 or index['ArIII_9m'] == 0: 
                        CHI_Ar2Ar3 = tol_max
                     else:   
                        CHI_Ar2Ar3 = (np.log10(index['ArII_7m']/index['ArIII_9m'])-                                                 Ar2Ar3_obs)**2/np.log10((index['ArII_7m']/index['ArIII_9m']))
                     if Ar23Ar5a_obs == -10: 
                        CHI_Ar23Ar5a = 0
                     elif index['ArII_7m'] == 0 or index['ArV_8m'] == 0 or index['ArIII_9m'] == 0: 
                        CHI_Ar23Ar5a = tol_max
                     else:   
                        CHI_Ar23Ar5a = (np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_8m'])- Ar23Ar5a_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_8m'])
                     if Ar23Ar5b_obs == -10: 
                        CHI_Ar23Ar5b = 0
                     elif index['ArII_7m'] == 0 or index['ArV_13m'] == 0 or index['ArIII_9m'] == 0: 
                        CHI_Ar23Ar5b = tol_max
                     else:   
                        CHI_Ar23Ar5b = (np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_13m'])- Ar23Ar5b_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m'])/index['ArV_13m'])
                     if CHI_Ar23Ar5a == 0 and CHI_Ar23Ar5b == 0:
                        CHI_Ar23Ar5 = 0
                     elif CHI_Ar23Ar5a > 0 and CHI_Ar23Ar5b == 0:
                        CHI_Ar23Ar5 = CHI_Ar23Ar5a
                     elif CHI_Ar23Ar5a == 0 and CHI_Ar23Ar5b > 0:
                        CHI_Ar23Ar5 = CHI_Ar23Ar5b
                     else:
                        CHI_Ar23Ar5 = (CHI_Ar23Ar5a + CHI_Ar23Ar5b)/2 
                     if Ar3Ar5a_obs == -10: 
                        CHI_Ar3Ar5a = 0
                     elif index['ArV_8m'] == 0 or index['ArIII_9m'] == 0:
                        CHI_Ar3Ar5a = tol_max
                     else:   
                        CHI_Ar3Ar5a = (np.log10((index['ArIII_9m'])/index['ArV_8m'])- Ar3Ar5a_obs)**2/np.log10((index['ArIII_9m'])/index['ArV_8m'])
                     if Ar3Ar5b_obs == -10: 
                        CHI_Ar3Ar5b = 0
                     elif index['ArV_13m'] == 0 or index['ArIII_9m'] == 0: 
                        CHI_Ar3Ar5b = tol_max
                     else:   
                        CHI_Ar3Ar5b = (np.log10((index['ArIII_9m'])/index['ArV_13m'])- Ar3Ar5b_obs)**2/np.log10((index['ArIII_9m'])/index['ArV_13m'])
                     if CHI_Ar3Ar5a == 0 and CHI_Ar3Ar5b == 0:
                        CHI_Ar3Ar5 = 0
                     elif CHI_Ar3Ar5a > 0 and CHI_Ar3Ar5b == 0:
                        CHI_Ar3Ar5 = CHI_Ar3Ar5a
                     elif CHI_Ar3Ar5a == 0 and CHI_Ar3Ar5b > 0:
                        CHI_Ar3Ar5 = CHI_Ar3Ar5b
                     else:
                        CHI_Ar3Ar5 = (CHI_Ar3Ar5a + CHI_Ar3Ar5b)/2 
                     if Ne2Ne3_obs == -10: 
                        CHI_Ne2Ne3 = 0
                     elif index['NeII_12m'] == 0 or index['NeIII_15m'] == 0: 
                        CHI_Ne2Ne3 = tol_max
                     else:   
                        CHI_Ne2Ne3 = (np.log10(index['NeII_12m']/index['NeIII_15m'])- Ne2Ne3_obs)**2/np.log10((index['NeII_12m']/index['NeIII_15m']))
                     if Ne23Ne5a_obs == -10: 
                        CHI_Ne23Ne5a = 0
                     elif index['NeII_12m'] == 0 or index['NeV_14m'] == 0 or index['NeIII_15m'] == 0: 
                        CHI_Ne23Ne5a = tol_max
                     else:   
                        CHI_Ne23Ne5a = (np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_14m'])- Ne23Ne5a_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_14m'])
                     if Ne23Ne5b_obs == -10: 
                        CHI_Ne23Ne5b = 0
                     elif index['NeII_12m'] == 0 or index['NeV_24m'] == 0 or index['NeIII_15m'] == 0: 
                        CHI_Ne23Ne5b = tol_max
                     else:   
                        CHI_Ne23Ne5b = (np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_24m'])- Ne23Ne5b_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m'])/index['NeV_24m'])
                     if CHI_Ne23Ne5a == 0 and CHI_Ne23Ne5b == 0:
                        CHI_Ne23Ne5 = 0
                     elif CHI_Ne23Ne5a > 0 and CHI_Ne23Ne5b == 0:
                        CHI_Ne23Ne5 = CHI_Ne23Ne5a
                     elif CHI_Ne23Ne5a == 0 and CHI_Ne23Ne5b > 0:
                        CHI_Ne23Ne5 = CHI_Ne23Ne5b
                     else:
                        CHI_Ne23Ne5 = (CHI_Ne23Ne5a + CHI_Ne23Ne5b)/2 
                     if Ne3Ne5a_obs == -10: 
                        CHI_Ne3Ne5a = 0
                     elif index['NeV_14m'] == 0 or index['NeIII_15m'] == 0: 
                        CHI_Ne3Ne5a = tol_max
                     else:   
                        CHI_Ne3Ne5a = (np.log10((index['NeIII_15m'])/index['NeV_14m'])- Ne3Ne5a_obs)**2/np.log10((index['NeIII_15m'])/index['NeV_14m'])
                     if Ne3Ne5b_obs == -10: 
                        CHI_Ne3Ne5b = 0
                     elif index['NeV_24m'] == 0 or index['NeIII_15m'] == 0: 
                        CHI_Ne3Ne5b = tol_max
                     else:   
                        CHI_Ne3Ne5b = (np.log10((index['NeIII_15m'])/index['NeV_24m'])- Ne3Ne5b_obs)**2/np.log10((index['NeIII_15m'])/index['NeV_24m'])
                     if CHI_Ne3Ne5a == 0 and CHI_Ne3Ne5b == 0:
                        CHI_Ne3Ne5 = 0
                     elif CHI_Ne3Ne5a > 0 and CHI_Ne3Ne5b == 0:
                        CHI_Ne3Ne5 = CHI_Ne3Ne5a
                     elif CHI_Ne3Ne5a == 0 and CHI_Ne3Ne5b > 0:
                        CHI_Ne3Ne5 = CHI_Ne3Ne5b
                     else:
                        CHI_Ne3Ne5 = (CHI_Ne3Ne5a + CHI_Ne3Ne5b)/2 
                     if O3O4a_obs == -10: 
                        CHI_O3O4a = 0
                     elif index['OIV_26m'] == 0 or index['OIII_52m'] == 0: 
                        CHI_O3O4a = tol_max
                     else:   
                        CHI_O3O4a = (np.log10(index['OIII_52m']/index['OIV_26m'])- O3O4a_obs)**2/np.log10((index['OIII_52m']/index['OIV_26m']))
                     if O3O4b_obs == -10: 
                        CHI_O3O4b = 0
                     elif index['OIV_26m'] == 0 or index['OIII_88m'] == 0: 
                        CHI_O3O4b = tol_max
                     else:   
                        CHI_O3O4b = (np.log10(index['OIII_88m']/index['OIV_26m'])- O3O4b_obs)**2/np.log10((index['OIII_88m']/index['OIV_26m']))
                     if CHI_O3O4a == 0 and CHI_O3O4b == 0:
                        CHI_O3O4 = 0
                     elif CHI_O3O4a > 0 and CHI_O3O4b == 0:
                        CHI_O3O4 = CHI_O3O4a
                     elif CHI_O3O4a == 0 and CHI_O3O4b > 0:
                        CHI_O3O4 = CHI_O3O4b
                     else:
                        CHI_O3O4 = (CHI_O3O4a + CHI_O3O4b)/2 
                     if S34O4a_obs == -10: 
                        CHI_S34O4a = 0
                     elif index['OIV_26m'] == 0 or index['SIII_18m'] == 0 or index['SIV_10m'] == 0: 
                        CHI_S34O4a = tol_max
                     else:   
                        CHI_S34O4a = (np.log10((index['SIII_18m']+index['SIV_10m'])/index['OIV_26m'])- S34O4a_obs)**2/np.log10(((index['SIII_18m']+index['SIV_10m'])/index['OIV_26m']))
                     if S34O4b_obs == -10: 
                        CHI_S34O4b = 0
                     elif index['OIV_26m'] == 0 or index['SIII_33m'] == 0 or index['SIV_10m'] == 0: 
                        CHI_S34O4b = tol_max
                     else:   
                        CHI_S34O4b = (np.log10((index['SIII_33m']+index['SIV_10m'])/index['OIV_26m'])- S34O4b_obs)**2/np.log10((index['SIII_33m']+index['SIV_10m'])/index['OIV_26m'])
                     if CHI_S34O4a == 0 and CHI_S34O4b == 0:
                        CHI_S34O4 = 0
                     elif CHI_S34O4a > 0 and CHI_S34O4b == 0:
                        CHI_S34O4 = CHI_S34O4a
                     elif CHI_S34O4a == 0 and CHI_S34O4b > 0:
                        CHI_S34O4 = CHI_S34O4b
                     else:
                        CHI_S34O4 = (CHI_S34O4a + CHI_S34O4b)/2 
                     if N2N3a_obs == -10: 
                        CHI_N2N3a = 0
                     elif index['NIII_57m'] == 0 or index['NII_122m'] == 0: 
                        CHI_N2N3a = tol_max
                     else:   
                        CHI_N2N3a = (np.log10(index['NII_122m']/index['NIII_57m'])- N2N3a_obs)**2/np.log10((index['NII_122m']/index['NIII_57m']))
                     if N2N3b_obs == -10: 
                        CHI_N2N3b = 0
                     elif index['NIII_57m'] == 0 or index['NII_205m'] == 0: 
                        CHI_N2N3b = tol_max
                     else:   
                        CHI_N2N3b = (np.log10(index['NII_205m']/index['NIII_57m'])- N2N3b_obs)**2/np.log10((index['NII_205m']/index['NIII_57m']))
                     if CHI_N2N3a == 0 and CHI_N2N3b == 0:
                        CHI_N2N3 = 0
                     elif CHI_N2N3a > 0 and CHI_N2N3b == 0:
                        CHI_N2N3 = CHI_N2N3a
                     elif CHI_N2N3a == 0 and CHI_N2N3b > 0:
                        CHI_N2N3 = CHI_N2N3b
                     else:
                        CHI_N2N3 = (CHI_N2N3a + CHI_N2N3b)/2 

                     if param == 1:
                        CHI_Teff = (CHI_S3S4**2 +  CHI_Ne2Ne3**2 + CHI_N2N3**2 + CHI_Ar2Ar3**2)**0.5
                     else:
                        if Ar23Ar5_obs > -10 and Ne23Ne5_obs > -10 and O3O4_obs > -10:
                           CHI_Teff = (CHI_Ne23Ne5**2 + CHI_O3O4**2 + CHI_Ar23Ar5**2)**0.5
                        elif Ar23Ar5_obs > -10 and Ne23Ne5_obs > -10 and O3O4_obs == -10:
                           CHI_Teff = (CHI_Ne23Ne5**2 + CHI_S34O4**2 + CHI_Ar23Ar5**2)**0.5
                        elif Ar23Ar5_obs == -10 and Ne23Ne5_obs > -10 and O3O4_obs > -10:
                           CHI_Teff = (CHI_Ne23Ne5**2 + CHI_O3O4**2 + CHI_Ar3Ar5**2)**0.5
                        elif Ar23Ar5_obs == -10 and Ne23Ne5_obs > -10 and O3O4_obs == -10:
                           CHI_Teff = (CHI_Ne23Ne5**2 + CHI_S34O4**2 + CHI_Ar3Ar5**2)**0.5
                        elif Ar23Ar5_obs == -10 and Ne23Ne5_obs == -10 and O3O4_obs > -10:
                           CHI_Teff = (CHI_Ne3Ne5**2 + CHI_O3O4**2 + CHI_Ar3Ar5**2)**0.5
                        elif Ar23Ar5_obs == -10 and Ne23Ne5_obs == -10 and O3O4_obs == -10:
                           CHI_Teff = (CHI_Ne3Ne5**2 + CHI_S34O4**2 + CHI_Ar3Ar5**2)**0.5


                     if param == 1:
                        Teff_p = index['T_eff']*(1/CHI_Teff)**2 + Teff_p
                     else:
                        Teff_p = index['a_ox']*(1/CHI_Teff)**2 + Teff_p
                     logU_p = index['logU'] *(1/CHI_Teff)**2 + logU_p         
                     den_Teff = (1/CHI_Teff)**2 + den_Teff
      
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


output['OH'] = OHffs
output['eOH'] = eOHffs
output['Teff'] = Teffs
output['eTeff'] = eTeffs
output['logU'] = logUffs
output['elogU'] = elogUffs

if input0.size == 1:  output = np.delete(output,obj=1,axis=0)


#np.savetxt('qq',output,fmt='%s')

lineas_header = [' HII-CHI-mistry-Teff-IR v.2.4 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type, 'Library file used : '+file_lib, '']

line_label = '{:10}  '.format(output.dtype.names[0])
for ind2 in range(1, len(output.dtype.names)-6):
   line_label += '{:10}  '.format(output.dtype.names[ind2])

if param == 1:
   line_label += '{:10}   {:10}   {:10}   {:10}   {:10}   {:10}'.format( 'O/H',    'eO/H',  'Teff' ,   'eTeff' , 'logU',   'elogU')
   lineas_header.append(line_label)
   header = '\n'.join(lineas_header)
   np.savetxt('.'.join(input00.split('.')[:-1])+'_hcm-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*(len(output.dtype.names)-7)+['%.2f']*2+['%i']*2+['%.2f']*2), header=header)
else:
   line_label += '{:10}   {:10}   {:10}   {:10}   {:10}   {:10}'.format( 'O/H',    'eO/H',  'a_OX' ,   'ea_OX' , 'logU',   'elogU')
   lineas_header.append(line_label)
   header = '\n'.join(lineas_header)
   np.savetxt('.'.join(input00.split('.')[:-1])+'_hcm-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*(len(output.dtype.names)-7)+['%.2f']*6), header=header)


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

print ('________________________________')
print ('Results are stored in ' + input00 + '_hcm-output.dat')


#############################################
###### INFORMATION AND CONTACT DETAILS ######
#############################################

# Enrique Perez-Montero, epm@iaa.es
# Borja Perez-Diaz, bperez@iaa.es


#################
###### END ######
#################
