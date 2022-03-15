# Filename: HCm-Teff_v5.1.py

import string
import numpy as np
import sys
#sys.stderr = open('errorlog.txt', 'w')


#Function for interpolation of grids

def interpolate(grid,z,zmin,zmax,n,sed_t):
   if sed_t != 3:
      label_t = 'Teff'
      file_1 = 'C17_WMb_Teff_30-60_pp.dat'
   else:
      label_t = 'f_abs'
      file_1 = 'C17_bpass21_imf135_300_pp_esc.dat'
   n_comments = 0
   with open('Libraries_teff/'+file_1, 'r') as file1:
      for line in file1:
         if line[0] == '#':
            n_comments += 1
   auxiliar_labels = np.genfromtxt('Libraries_teff/'+file_1, dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
   ncol = len(auxiliar_labels)
   vec = []
   if z == 2:
      label_z = 'logUsp'
   if z == 1:
      label_z = label_t
   if z == 0:
      label_z = '12logOH'
   type_list_names = []
   for col in auxiliar_labels:
      inter = 0
      no_inter = 0
      type_list_names.append((col, float))
      for row in range(0,len(grid)):
         if grid[label_z][row] < zmin or grid[label_z][row] > zmax: continue
         if z == 2: x = '12logOH'; y = label_t
         if z == 1: x = '12logOH'; y = 'logUsp'
         if z == 0: x = label_t; y = 'logUsp'
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






print (' ---------------------------------------------------------------------')
print ('This is HII-CHI-mistry-Teff v. 5.1')
print (' See Perez-Montero et al (2019) for details')
print ( ' Insert the name of your input text file with all or some of the following columns:')
print ('12+log(O/H)')
print ('[OII] 3727')
print ('[OIII] 5007')
print ('[SII] 6725')
print ('[SIII] 9069')
print ('HeI 4471')
print ('HeI 5876')
print ('HeII 4686')
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
   print ('-------------------------------------------------')
   print ('Default SEDs')
   print ('------------')
   print ('(1) Effective temperature, WM-Basic stellar atmospheres (30-60 kK)')
   print ('(2) Effective temperature, Black body (30-90 kK)')
   print ('(3) Photon escape fraction, BPASS cluster atmospheres, age = 4 Myr, Mup = 300, x = 1.35, w/binaries')
   print ('')
   #print('Other SED')
   #print('---------')
   #print ('(4) Different library')
   print('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      sed = raw_input('Choose parameter and models: ')
   else:
      sed = input('Choose parameter and models: ')
   if sed == '1' or sed == '2' or sed == '3': question = False

print ('')
question = True
while question:
   print ('---------------------------------------------------------------------')
   if sed == '3': #or sed == '4': 
      geo = '2'
      question = False
   else:
      print ('(1) Plane-parallel geometry')
      print ('(2) Spherical geometry')
      print ('---------------------------------------------------------------------')
      if int(sys.version[0]) < 3:
         geo = raw_input('Choose geometry of the models: ')
      else:
         geo = input('Choose geometry of the models: ')
      if geo == '1' or geo == '2': question = False
print('')


#Particular file introduced by the user (not tested)
if sed == '4':
   question = True
   while question:
      print('Introduce name of the file containing the models. It must be located in the folder "Libraries_teff".')
      print(' ')
      if int(sys.version[0]) < 3:
         new_library = raw_input('Name of file: ')
      else:
         new_library = input('Name of file: ')
 
      #Searching for the file
      try:
         #Counting comments:
         n_comments = 0
         with open('Libraries_teff/'+new_library, 'r') as file:
            for line in file:
               if line[0] == '#':
                  n_comments += 1
         library_user = np.genfromtxt('Libraries_teff/'+new_library, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         print(' ')
         print('Loading library '+new_library+'. Checking correct format of the file.')
         question = False
      except:
         print(' ')
         print('Library was not found in folder "Libraries_teff" or file does not exist.')
   question = True
   while question:
      try:
         #Counting comments:
         n_comments = 0
         with open('Libraries_teff/'+new_library, 'r') as file:
            for line in file:
               if line[0] == '#':
                  n_comments += 1
         library_user = np.genfromtxt('Libraries_teff/'+new_library, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         #Checking correct format:
         #Counting comments:
         n_comments = 0
         with open('Libraries_teff/C17_bb_Teff_30-90_pp.dat', 'r') as file:
            for line in file:
               if line[0] == '#':
                  n_comments += 1
         auxiliar_labels = np.genfromtxt('Libraries_teff/C17_bb_Teff_30-90_pp.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
         missing_labels = []
         for label in auxiliar_labels:
            if label in library_user.dtype.names:
               continue
            else:
               missing_labels.append(label)
         #Displaying message for the user:
         print('Succesfully reading of the file')
         if len(missing_labels) == 0:
            print('File presents the correct format')
            question = False
         else:
            print('File does not present the correct format. The following columns are missing:')
            for need_label in missing_labels:
               print('- '+need_label)
            print('More details on the correct format for the library are found in readme file.')
            print(' ')
            print('Reintroduce name of the file with fixed format:')
            print(' ')
            if int(sys.version[0]) < 3:
               new_library = raw_input('Name of file: ')
            else:
               new_library = input('Name of file: ')
      except:
         print('Something went wrong while reading file. Please, reintroduce name of the file:')
         print(' ')
         if int(sys.version[0]) < 3:
            new_library = raw_input('Name of file: ')
         else:
            new_library = input('Name of file: ')

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
   file_lib = 'C17_bb_Teff_30-90_pp.dat'
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
elif geo==2 and sed == 2:
   bin = 132
   file_lib = 'C17_bb_Teff_30-90_sph.dat'
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

elif geo==1 and sed == 3:
   sed_type = 'BPASS cluster atmospheres, Mup = 300. Plane-parallel geometry'
   bin = 63
   file_lib = 'C17_bpass21_imf135_300_pp_esc.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   print ('f_abs and U calculation using BPASS models with plane-parallel geometry')
elif geo==2 and sed == 3:
   bin = 77
   file_lib = 'C17_bpass21_imf135_300_sph_esc.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_teff/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_teff/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type ='BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. Spherical geometry. Not interpolated'
      print ('F_abs and U calculation usingBPASS models with spherical geometry and non-interpolation')
   elif inter == 1:
      sed_type ='BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. Spherical geometry. Interpolated'
      print ('F_abs and U calculation usingBPASS models with spherical geometry and interpolation')


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

#Final grid: no limitations
grid = grid_aux
print(' ')
print('The following grid is going to be used:')
print('- Full library: '+file_lib)
print('       Total number of models: ' + str(len(grid)))
print(' ')

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
   if input1.dtype.names[col] == 'SII_6717':
      Label_SII_6716 = True
   if input1.dtype.names[col] == 'SII_6731':
      Label_SII_6731 = True
   if input1.dtype.names[col] == 'eSII_6717':
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
   SII_6725 = input1['SII_6717']+input1['SII_6731']
if Label_eSII == False and (Label_eSII_6716 == False or Label_eSII_6731 == False):
   eSII_6725 = np.zeros(input1.size)
elif Label_eSII == True:
   eSII_6725 = input1['eSII_6725']
else:
   eSII_6725 = input1['eSII_6717']+input1['eSII_6731']
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


#Creation of output only with information from inputs
aux_list = []
aux_list.append(('ID','U12'))
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
   aux_list.append(('SII_6717', float))
if Label_eSII_6716 == True:
   aux_list.append(('eSII_6717', float))
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
if Label_HeII_4686 == True:
   aux_list.append(('HeII_4686', float))
if Label_eHeII_4686 == True:
   aux_list.append(('eHeII_4686', float))

aux_list.append(('OH', float))
aux_list.append(('eOH', float))
if sed == 3:
   aux_list.append(('f_abs', float))
   aux_list.append(('ef_abs', float))
else:
   aux_list.append(('Teff', float))
   aux_list.append(('eTeff', float))
aux_list.append(('logU', float))
aux_list.append(('elogU', float))
output = np.zeros(input1.size, dtype=aux_list)

output['ID'] = Names
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
if Label_HeII_4686 == True:
   output['HeII_4686'] = HeII_4686
if Label_eHeII_4686 == True:
   output['eHeII_4686'] = eHeII_4686
if Label_SII == True:
   output['SII_6725'] = SII_6725
if Label_eSII == True:
   output['eSII_6725'] = eSII_6725
if Label_SII_6716 == True:
   output['SII_6717'] = SII_6716
if Label_eSII_6716 == True:
   output['eSII_6717'] = eSII_6716
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


print ('Reading grids ....')
print ('')
print ('')
print ('---------------------------------------')
if sed != 3:
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
            for y in grid.dtype.names:
               grid_T0.append(grid[y][i0+x]*np.abs(0.3-OH+grid['12logOH'][i0])/0.3+grid[y][i1+x]*np.abs(0.3-grid['12logOH'][i1]+OH)/0.3)
               #grid_T0.append(grid[i0+x,y]*np.abs(0.3-OH+grid[i0,0])/0.3+grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)
            
      #         grid_T0.append(grid[i0+x,y]*np.abs(0.3-grid[i0,0]+OH)/0.3 + grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)


         grid_T_aux = np.reshape(grid_T0,(bin,10))
         grid_T = np.zeros(grid_T_aux.shape[0], dtype=grid.dtype)
         for col_n in range(0, len(grid.dtype.names)):
             grid_T[grid.dtype.names[col_n]] = grid_T_aux[:, col_n]

      else:
         OH = 0
         OH_mc.append(OH)
         grid_T = grid
   

      #np.savetxt('int_models.dat',grid_T,fmt='%.2f')


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
            if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
            if S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            else:
               CHI_S2S3 = (np.log10(index['SII_671731']/index['SIII_9069']) - S2S3_obs)**2/S2S3_obs
               CHI_S23 = (index['SII_671731']+index['SIII_9069']-S23_obs)**2/S23_obs
            if S2O3_obs == -10:
               CHI_S2O3 = 0
            elif index['SII_671731'] == 0 or index['OIII_5007'] == 0:
               CHI_S2O3 = tol_max
            else:
               CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/S2O3_obs
            if O2O3_obs == -10:
               CHI_O2O3 = 0
            elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
               CHI_O2O3 = tol_max
            else:
               CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/O2O3_obs
            if He12a_obs == -10:
               CHI_He12a = 0
            elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
               CHI_He12a = tol_max
            else:
               CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/He12a_obs
            if He12b_obs == -10:
               CHI_He12b = 0
            elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
               CHI_He12b = tol_max
            else:
               CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/He12b_obs
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


            if sed != 3:
               Teff_p = index['Teff']*(1/CHI_Teff)**2 + Teff_p
            else:
               Teff_p = index['f_abs']*(1/CHI_Teff)**2 + Teff_p 
            logU_p = index['logUsp'] *(1/CHI_Teff)**2 + logU_p         
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
            if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
            if S2S3_obs == -10:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max    
            else:
               CHI_S2S3 = (np.log10(index['SII_671731']/index['SIII_9069']) - S2S3_obs)**2/S2S3_obs
               CHI_S23 = (index['SII_671731']+index['SIII_9069']-S23_obs)**2/S23_obs
            if S2O3_obs == -10:
               CHI_S2O3 = 0
            elif index['SII_671731'] == 0 or index['OIII_5007'] == 0:
               CHI_S2O3 = tol_max
            else:
               CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/S2O3_obs
            if O2O3_obs == -10:
               CHI_O2O3 = 0
            elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
               CHI_O2O3 = tol_max
            else:
               CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/O2O3_obs
            if He12a_obs == -10:
               CHI_He12a = 0
            elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
               CHI_He12a = tol_max
            else:
               CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/He12a_obs
            if He12b_obs == -10:
               CHI_He12b = 0
            elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
               CHI_He12b = tol_max
            else:
               CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/He12b_obs
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




         
            if sed != 3:
               Teff_e = np.abs(index['Teff'] - Teff) * (1/CHI_Teff) **2+ Teff_e
            else:
               Teff_e = np.abs(np.log10(index['f_abs']+1e-5) - np.log10(Teff)) * (1/CHI_Teff)**2 + Teff_e
            logU_e = np.abs(index['logUsp'] - logU) * (1/CHI_Teff)**2 + logU_e         
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
            igrid = grid_T[np.lexsort((grid_T['12logOH'],grid_T['logUsp']))]
            igrid = interpolate(igrid,1,Teff-eTeff-1e4,Teff+eTeff+1e4,10, sed)
            if sed != 3:
               igrid = igrid[np.lexsort((igrid['12logOH'],igrid['Teff']))]
            else:
               igrid = igrid[np.lexsort((igrid['12logOH'],igrid['f_abs']))]
            igrid = interpolate(igrid,2,logU-elogU-0.25,logU+elogU+0.25,10, sed)

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
                  if index['HeII_4686'] == 0 and HeII_4686_obs > 0: continue
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
                     CHI_S2O3 = (np.log10(index['SII_671731']/index['OIII_5007']) - S2O3_obs)**2/S2O3_obs
                  if O2O3_obs == -10:
                     CHI_O2O3 = 0
                  elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                     CHI_O2O3 = tol_max
                  else:
                     CHI_O2O3 = (np.log10(index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/O2O3_obs
                  if He12a_obs == -10:
                     CHI_He12a = 0
                  elif index['HeI_4471'] == 0 or index['HeII_4686'] == 0:
                     CHI_He12a = tol_max
                  else:
                     CHI_He12a = (np.log10(index['HeI_4471']/index['HeII_4686']) - He12a_obs)**2/He12a_obs
                  if He12b_obs == -10:
                     CHI_He12b = 0
                  elif index['HeI_5876'] == 0 or index['HeII_4686'] == 0:
                     CHI_He12b = tol_max
                  else:
                     CHI_He12b = (np.log10(index['HeI_5876']/index['HeII_4686']) - He12b_obs)**2/He12b_obs
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
      
   
         
                  if sed != 3:
                     Teff_p = index['Teff']*(1/CHI_Teff) + Teff_p
                  else:
                     Teff_p = index['f_abs']*(1/CHI_Teff) + Teff_p
                  logU_p = index['logUsp'] *(1/CHI_Teff) + logU_p         
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
if sed != 3:
   output['Teff'] = Teffs
   output['eTeff'] = eTeffs
else:
   output['f_abs'] = Teffs
   output['ef_abs'] = eTeffs
output['logU'] = logUffs
output['elogU'] = elogUffs


if input0.size == 1:  output = np.delete(output,obj=1,axis=0)


lineas_header = [' HII-CHI-mistry-Teff v.5.1 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type, 'Library file used : '+file_lib, '']

line_label = '{:10}  '.format(output.dtype.names[0])
for ind2 in range(1, len(output.dtype.names)-7):
   line_label += '{:10}  '.format(output.dtype.names[ind2])

if sed != 3:
   line_label += '{:10}   {:10}   {:10}   {:10}   {:10}   {:10}   {:10}'.format('i', 'O/H',    'eO/H',  'Teff' ,   'eTeff' , 'logU',   'elogU')
else:
   line_label += '{:10}   {:10}   {:10}   {:10}   {:10}   {:10}   {:10}'.format('i', 'O/H',    'eO/H',  'f_abs' ,   'ef_abs' , 'logU',   'elogU')

lineas_header.append(line_label)
header = '\n'.join(lineas_header)
np.savetxt('.'.join(input00.split('.')[:-1])+'_hcm-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*(len(output.dtype.names)-8)+['%i']+['%.2f']*6), header=header)

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
