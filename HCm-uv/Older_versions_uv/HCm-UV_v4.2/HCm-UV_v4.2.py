# Filename: HCm_UV_v4.2.py

import string
import numpy as np
import sys
#sys.stderr = open('errorlog.txt', 'w')
import warnings
warnings.filterwarnings("ignore")

#Function for interpolation of grids

def interpolate(grid,z,zmin,zmax,n):
   #Columns of the library
   n_comments = 0
   with open('Libraries_uv/C17_POPSTAR_1myr.dat', 'r') as file1:
      for line in file1:
         if line[0] == '#':
            n_comments += 1
   auxiliar_labels = np.genfromtxt('Libraries_uv/C17_POPSTAR_1myr.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
   ncol = len(auxiliar_labels)
   vec = []
   if z == 2:
      label_z = 'logU'
   if z == 1:
      label_z = 'logCO'
   if z == 0:
      label_z = '12logOH'
   type_list_names = []
   for col in auxiliar_labels:
      inter = 0
      no_inter = 0
      type_list_names.append((col, float))
      for row in range(0,len(grid)):
         if grid[label_z][row] < zmin or grid[label_z][row] > zmax: continue
         if z == 2: x = '12logOH'; y = 'logCO'
         if z == 1: x = '12logOH'; y = 'logU'
         if z == 0: x = 'logCO'; y = 'logU'
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
print (' This is HII-CHI-mistry for UV version 4.2')
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
   #Counting comments:
   n_comments = 0
   with open(input00, 'r') as file2:
      for line in file2:
         if line[0] == '#':
            n_comments += 1		
   input0 = np.genfromtxt(input00,dtype=None,names=True, encoding = 'ascii', skip_header = n_comments)
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
   print('Default SEDs')
   print('------------')
   print ('(1) POPSTAR with Chabrier IMF, age = 1 Myr')
   print ('(2) BPASS v.2.1 a_IMF = 1.35, Mup = 300, age = 1Myr')
   print ('')
   print('Other SED')
   print('---------')
   print ('(3) Different library')
   print('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      sed = raw_input('Choose SED of the models: ')
   else:
      sed = input('Choose SED of the models: ')
   if sed == '1' or sed == '2' or sed == '3': question = False 
print ('')


#Particular file introduced by the user
if sed == '3':
   question = True
   while question:
      print('Introduce name of the file containing the models. It must be located in the folder "Libraries_uv".')
      print(' ')
      if int(sys.version[0]) < 3:
         new_library = raw_input('Name of file: ')
      else:
         new_library = input('Name of file: ')
 
      #Searching for the file
      try:
         #Counting comments:
         n_comments = 0
         with open('Libraries_uv/'+new_library, 'r') as file3:
            for line in file3:
               if line[0] == '#':
                  n_comments += 1
         library_user = np.genfromtxt('Libraries_uv/'+new_library, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         print(' ')
         print('Loading library '+new_library+'. Checking correct format of the file.')
         question = False
      except:
         print(' ')
         print('Library was not found in folder "Libraries_uv" or file does not exist.')
   question = True
   while question:
      try:
         #Counting comments:
         n_comments = 0
         with open('Libraries_uv/'+new_library, 'r') as file4:
            for line in file4:
               if line[0] == '#':
                  n_comments += 1
         library_user = np.genfromtxt('Libraries_uv/'+new_library, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         #Checking correct format:
         #Counting comments:
         n_comments = 0
         with open('Libraries_uv/C17_POPSTAR_1myr.dat', 'r') as file5:
            for line in file5:
               if line[0] == '#':
                  n_comments += 1
         auxiliar_labels = np.genfromtxt('Libraries_uv/C17_POPSTAR_1myr.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
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
print(' ')
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

#POPSTAR MODEL
if sed==1:
   file_lib = 'C17_POPSTAR_1myr.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file6:
      for line in file6:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation'
      print ('No interpolation for the POPSTAR models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      print ('')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. Interpolation'
      print ('Interpolation for the POPSTAR models is going to be used.')
      print ('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
      print ('')
      res_NO = 0.125
        
#BPASS MODEL
elif sed==2:
   file_lib = 'C17_BPASS_IMF135_mup300_1myr.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file7:
      for line in file7:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)  
   if inter == 0:
      sed_type = 'BPASS a_IMF = 1.35, M_up = 300, age = 1Myr. No interpolation'
      print ('No interpolation for the BPASS models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      print ('')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'BPASS v.2.1, a_IMF = 1.35, M_up = 300, age = 1Myr. Interpolation'
      print ('Interpolation for the BPASS  models is going to be used.')
      print ('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
      print ('')
      res_NO = 0.125

#Different library
elif sed==3:
   file_lib = new_library
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+new_library, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1  
   grid_aux = np.genfromtxt('Libraries_uv/'+new_library,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'User file ' + new_library + ' used as library for the models no interpolated'
      print ('No interpolation for the library '+new_library)
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'User file ' + new_library + ' used as library for the models interpolated'
      print ('Interpolation for the library '+new_library)




#Reading constraints and creating library with constraints
print(' ')
print('Select a file with the constraint laws to be used to limit the grid of models when the measurement of a quantity is impossible without any relation.')
print(' ')

print ('')
question = True
while question:
   print('-------------------------------------------------')
   print('Default constraints')
   print('------------')
   print ('(1) Constraints for Star-Forming Galaxies')
   print ('(2) Constraints for Extreme Emission Line Galaxies')
   print ('')
   print('Other constraints')
   print('---------')
   print ('(3) Different constraint file')
   print('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      const = raw_input('Choose constraint for the grids: ')
   else:
      const = input('Choose constraint for the grids: ')
   if const == '1' or const == '2' or const == '3': question = False 
print ('')

#Particular file introduced by the user
if const == '3':
   question = True
   while question:
      print('Introduce name of the file containing the constraints for the grids. It must be located in the folder "Constraints".')
      print(' ')
      if int(sys.version[0]) < 3:
         new_const = raw_input('Name of file: ')
      else:
         new_const = input('Name of file: ')
 
      #Searching for the file
      try:
         #Counting comments:
         n_comments = 0
         with open('Constraints/'+new_const, 'r') as file9:
            for line in file9:
               if line[0] == '#':
                  n_comments += 1
         const_user = np.genfromtxt('Constraints/'+new_const, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         print(' ')
         print('Loading constraint file '+new_const+'. Checking correct format of the file.')
         question = False
      except:
         print(' ')
         print('File was not found in folder "Constraints" or file does not exist.')
   question = True
   while question:
      try:
         #Counting comments:
         n_comments = 0
         with open('Constraints/'+new_const, 'r') as file10:
            for line in file10:
               if line[0] == '#':
                  n_comments += 1
         const_user = np.genfromtxt('Constraints/'+new_const, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         #Checking correct format:
         #Counting comments:
         n_comments = 0
         with open('Constraints/template_OH.dat', 'r') as file11:
            for line in file11:
               if line[0] == '#':
                  n_comments += 1
         auxiliar_labels = np.genfromtxt('Constraints/template_OH.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
         missing_labels = []
         for label in auxiliar_labels:
            if label in const_user.dtype.names:
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
               new_const = raw_input('Name of file: ')
            else:
               new_const = input('Name of file: ')
      except:
         print('Something went wrong while reading file. Please, reintroduce name of the file:')
         print(' ')
         if int(sys.version[0]) < 3:
            new_const = raw_input('Name of file: ')
         else:
            new_const = input('Name of file: ')

#Valuable columns of the files
uv_lin = ['12logOH', 'logCO', 'logU', 'Lya_1216', 'CIV_1549', 'HeII_1640', 'OIII_1665', 'CIII_1909', 'OIII_5007']
lin_uv_label = ['12+log(O/H)', 'log(C/O)', 'log(U)', 'Lya_1216', 'CIV_1549', 'HeII_1640', 'OIII_1665', 'CIII_1909', 'OIII_5007']

#Generation of grids with constraints laws:
if const == '1' or const == '2' or const == '3':
   #First grid does not change
   grid1 = grid_aux
   file_lib_2 = file_lib

#Generating libraries for the constraints in the files
if const == '1': #Star-Forming Galaxies
   const_file = 'template_OH.dat'
   name_const = 'Constraints/template_OH.dat'
   n_comments = 0
   with open(name_const, 'r') as file12:
      for line in file12:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == '2':
   const_file = 'template_OH_eelg.dat'
   name_const = 'Constraints/template_OH_eelg.dat'
   n_comments = 0
   with open(name_const, 'r') as file13:
      for line in file13:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == '3':
   const_file = new_const
   name_const = 'Constraints/'+new_const
   n_comments = 0
   with open(name_const, 'r') as file14:
      for line in file14:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

#Limiting the grids:
lin_uv_val = []
#The initial grid need to be constrained in the ioniztion parameter

#Retrieving each column of the grid
for label in uv_lin:
   aux_col = grid1[label].tolist()
   lin_uv_val.append(aux_col)
   #Creation of the grids
   name_OH_U = '.'.join(file_lib_2.split('.')[0:-1])+'_OH_U_constrained.'+file_lib.split('.')[-1]
   name_OH_U_CO = '.'.join(file_lib_2.split('.')[0:-1])+'_OH_U_CO_constrained.'+file_lib.split('.')[-1]
   file_open = open('Libraries_uv/'+ name_OH_U, 'w') #OH and U relation
   file_open_2 = open('Libraries_uv/'+name_OH_U_CO, 'w') #OH, NO and U relation
   file_open.write('#Constrained by relation between 12+log(O/H) and log(U)\n')
   file_open_2.write('#Constrained by relation between 12+log(O/H), log(U) and log(C/O)\n')
   #Header line
   label_line = '{:15} '.format(lin_uv_label[0].replace(' ',''))
   for ind in range(1, len(lin_uv_label)-1):
      label_line += '\t {:15} '.format(lin_uv_label[ind].replace(' ',''))
   label_line += '\t {:15}\n'.format(lin_uv_label[-1].replace(' ','')) 
   file_open.write(label_line)
   file_open_2.write(label_line)
#Values:
for ind_val in range(0, len(lin_uv_val[0])):
   index_desired = np.where(const_data['12logOH'] == lin_uv_val[0][ind_val])[0][0] #Searching for constrain in given value of O/H
   if lin_uv_val[2][ind_val] <= const_data['logU_max'][index_desired] and lin_uv_val[2][ind_val] >= const_data['logU_min'][index_desired]:
      val_line = '{:7.7f} '.format(lin_uv_val[0][ind_val])
      for ind2 in range(1, len(lin_uv_label)-1):
         val_line += '\t {:7.7f} '.format(lin_uv_val[ind2][ind_val])
      val_line += '\t {:7.7f}\n'.format(lin_uv_val[-1][ind_val])
      file_open.write(val_line)
   if lin_uv_val[2][ind_val] <= const_data['logU_max'][index_desired] and lin_uv_val[2][ind_val] >= const_data['logU_min'][index_desired] and lin_uv_val[1][ind_val] <= const_data['logCO_max'][index_desired] and lin_uv_val[1][ind_val] >= const_data['logCO_min'][index_desired]:
      val_line = '{:7.7f} '.format(lin_uv_val[0][ind_val])
      for ind2 in range(1, len(lin_uv_label)-1):
         val_line += '\t {:7.7f} '.format(lin_uv_val[ind2][ind_val])
      val_line += '\t {:7.7f}\n'.format(lin_uv_val[-1][ind_val])
      file_open_2.write(val_line)
file_open.close()
file_open_2.close()
#Counting comments:
n_comments = 0
with open('Libraries_uv/'+name_OH_U, 'r') as file15:
   for line in file15:
      if line[0] == '#':
         n_comments += 1
grid2 = np.genfromtxt('Libraries_uv/'+name_OH_U,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
n_comments = 0
with open('Libraries_uv/'+name_OH_U_CO, 'r') as file:
   for line in file:
      if line[0] == '#':
         n_comments += 1
grid3 = np.genfromtxt('Libraries_uv/'+name_OH_U_CO,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

#Residual in CO
if inter==0:
   res_CO = np.max([sorted(set(grid1['logCO']))[ind+1]-sorted(set(grid1['logCO']))[ind] for ind in range(0, len(set(grid1['logCO']))-1)])
if inter==1:
   res_CO = np.max([sorted(set(grid1['logCO']))[ind+1]-sorted(set(grid1['logCO']))[ind] for ind in range(0, len(set(grid1['logCO']))-1)])/10
   
print('Libraries generated with the constraints. The following grids are going to be used:')
print('- Full library (Grid#1): '+file_lib_2)
print('       Total number of models: ' + str(len(grid1)))
print('- Library constrained by 12+log(O/H) - log(U) relation (Grid#2): '+name_OH_U)
print('       Total number of models: ' + str(len(grid2)))
print('- Library constrained by 12+log(O/H) - log(U) - log(C/O) relation (Grid#3): '+name_OH_U_CO)
print('       Total number of models: ' + str(len(grid3)))
print(' ')

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


#Creation of output only with information from inputs
aux_list = []
aux_list.append(('ID','U12'))
if Label_Lya == True:
   aux_list.append(('Lya_1216', float))
if Label_eLya == True:
   aux_list.append(('eLya_1216', float))
if Label_CIV == True:
   aux_list.append(('CIV_1549', float))
if Label_eCIV == True:
   aux_list.append(('eCIV_1549', float))
if Label_HeII == True:
   aux_list.append(('HeII_1640', float))
if Label_eHeII == True:
   aux_list.append(('eHeII_1640', float))
if Label_OIII_1665 == True:
   aux_list.append(('OIII_1665', float))
if Label_eOIII_1665 == True:
   aux_list.append(('eOIII_1665', float))
if Label_CIII == True:
   aux_list.append(('CIII_1909', float))
if Label_eCIII == True:
   aux_list.append(('eCIII_1909', float))
if Label_Hbeta == True:
   aux_list.append(('Hb_4861', float))
if Label_eHbeta == True:
   aux_list.append(('eHb_4861', float))
if Label_OIII_5007 == True:
   aux_list.append(('OIII_5007', float))
if Label_eOIII_5007 == True:
   aux_list.append(('eOIII_5007', float))



aux_list.append(('grid', int))
aux_list.append(('OH', float))
aux_list.append(('eOH', float))
aux_list.append(('CO', float))
aux_list.append(('eCO', float))
aux_list.append(('logU', float))
aux_list.append(('elogU', float))
output = np.zeros(input1.size, dtype=aux_list)



#output = np.zeros(input1.size, dtype=[('ID', 'U12'), ('Lya_1216', float),('eLya_1216', float),('CIV_1549', float),('eCIV_1549', float),('HeII_1640', float),('eHeII_1640', float),('OIII_1665', float),('eOIII_1665', float),('CIII_1909', float),('eCIII_1909', float),('Hb_4861', float),('eHb_4861', float),('OIII_5007', float),('eOIII_5007', float),('grid', int),('OH', float),('eOH', float),('CO', float),('eCO', float),('logU', float),('elogU', float)] )

output['ID'] = Names
if Label_Lya == True:
   output['Lya_1216'] = Lya_1216
if Label_eLya == True:
   output['eLya_1216'] = eLya_1216
if Label_CIV == True:
   output['CIV_1549'] = CIV_1549
if Label_eCIV == True:
   output['eCIV_1549'] = eCIV_1549
if Label_HeII == True:
   output['HeII_1640'] = HeII_1640
if Label_eHeII == True:
   output['eHeII_1640'] = eHeII_1640
if Label_OIII_1665 == True:
   output['OIII_1665'] = OIII_1665
if Label_eOIII_1665 == True:
   output['eOIII_1665'] = eOIII_1665
if Label_CIII == True:
   output['CIII_1909'] = CIII_1909
if Label_eCIII == True:
   output['CIII_1909'] = CIII_1909
if Label_Hbeta == True:
   output['Hb_4861'] = Hb_4861
if Label_eHbeta == True:
   output['eHb_4861'] = eHb_4861
if Label_OIII_5007 == True:
   output['OIII_5007'] = OIII_5007
if Label_eOIII_5007 == True:
   output['eOIII_5007'] = eOIII_5007

print ('Reading grids ....')
print ('')
print ('')
print ('----------------------------------------------------------------')
print ('(%)   ID    Grid  12+log(O/H)  log(C/O)    log(U)')
print ('----------------------------------------------------------------')

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
      if Lya_1216[tab] <= 0:
         Lya_1216_obs = 0
      else:
         while Lya_1216_obs <= 0:
            Lya_1216_obs = np.random.normal(Lya_1216[tab],eLya_1216[tab]+1e-5)
      CIV_1549_obs = 0
      if CIV_1549[tab] <= 0:
         CIV_1549_obs = 0
      else:
         while CIV_1549_obs <= 0:
            CIV_1549_obs = np.random.normal(CIV_1549[tab],eCIV_1549[tab]+1e-5)
      HeII_1640_obs = 0
      if HeII_1640[tab] <= 0:
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
      if CIII_1909[tab] <= 0:
         CIII_1909_obs = 0
      else:
         while CIII_1909_obs <= 0:
            CIII_1909_obs = np.random.normal(CIII_1909[tab],eCIII_1909[tab]+1e-5)
      Hb_4861_obs = 0
      if Hb_4861[tab] <= 0:
         Hb_4861_obs = 0
      else:
         while Hb_4861_obs <= 0:
            Hb_4861_obs = np.random.normal(Hb_4861[tab],eHb_4861[tab]+1e-5)
      OIII_5007_obs = 0
      if OIII_5007[tab] <= 0:
         OIII_5007_obs = 0
      else:
         while OIII_5007_obs <= 0:
            OIII_5007_obs = np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-5)
      if OIII_1665_obs <= 0 or OIII_5007_obs <= 0:
         ROIII_obs = 0
      else:
         ROIII_obs = OIII_5007_obs/OIII_1665_obs
      if Lya_1216_obs <= 0 or CIII_1909_obs <= 0:
         C34_obs = 0
      else:
         C34_obs = (CIII_1909_obs + CIV_1549_obs) / (Lya_1216_obs)
      if HeII_1640_obs <= 0 or CIII_1909_obs <= 0:
         C34He2_obs = 0
      else:
         C34He2_obs = (CIII_1909_obs + CIV_1549_obs) / (HeII_1640_obs)
      if CIII_1909_obs <= 0 or OIII_1665_obs <= 0:
         C3O3_obs = -10
      else:   
         C3O3_obs = np.log10((CIII_1909_obs) / (OIII_1665_obs))
      if CIII_1909_obs <= 0 or CIV_1549_obs <= 0:
         C3C4_obs = 0
      else:
         C3C4_obs = (CIII_1909_obs/CIV_1549_obs)
      if CIII_1909_obs <= 0 or Hb_4861_obs <= 0:
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
            elif index['OIII_1665'] == 0 or index['OIII_5007'] == 0:
               CHI_ROIII = tol_max
            else:
               CHI_ROIII = (index['OIII_5007']/index['OIII_1665'] - ROIII_obs)**2/(index['OIII_5007']/index['OIII_1665'])
            if C3O3_obs == -10:
               CHI_C3O3 = 0          
            elif index['CIII_1909'] == 0 or index['OIII_1665'] == 0:
               CHI_C3O3 = tol_max
            else:
               CHI_C3O3 =(np.log10((index['CIII_1909'])/index['OIII_1665']) - C3O3_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_1665']+1e-5))

            CHI_CO = (CHI_ROIII**2 + CHI_C3O3**2 )**0.5


            if CHI_CO == 0:
               CO_p = CO_p
               den_CO = den_CO
            else:
               CO_p = index['logCO'] /np.exp(CHI_CO) + CO_p
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
            elif index['OIII_1665'] == 0 or index['OIII_5007'] == 0:
               CHI_ROIII = tol_max
            else:
               CHI_ROIII = (index['OIII_5007']/index['OIII_1665'] - ROIII_obs)**2/(index['OIII_5007']/index['OIII_1665'])
            if C3O3_obs == -10:
               CHI_C3O3 = 0          
            elif index['CIII_1909'] == 0 or index['OIII_1665'] == 0:
               CHI_C3O3 = tol_max
            else:
               CHI_C3O3 =(np.log10((index['CIII_1909'])/index['OIII_1665']) - C3O3_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_1665']+1e-5))

            CHI_CO = (CHI_ROIII**2 + CHI_C3O3**2 )**0.5


            if CHI_CO == 0:
               CO_e = CO_e
               den_CO_e = den_CO_e  
            else:
               CO_e = (index['logCO'] - CO)**2 / np.exp(CHI_CO) + CO_e
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
            if CO > -10 and np.abs(index['logCO'] - CO) > np.abs(eCO+0.125):
               continue
            if CIV_1549_obs > 0 and index['CIV_1549'] == 0:
               continue
            if HeII_1640_obs > 0 and index['HeII_1640'] == 0:
               continue
            else:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index['OIII_1665'] == 0 or index['OIII_5007'] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index['OIII_5007']/index['OIII_1665'] - ROIII_obs)**2/(index['OIII_5007']/index['OIII_1665'])
               if C34_obs == 0:
                  CHI_C34 = 0
               elif index['Lya_1216'] == 0 or index['CIII_1909'] == 0:
                  CHI_C34 = tol_max
               else:
                  CHI_C34 = ((index['CIII_1909']+index['CIV_1549'])/index['Lya_1216'] - C34_obs)**2/((index['CIII_1909']+index['CIV_1549'])/index['Lya_1216'])
               if C34He2_obs == 0:
                  CHI_C34He2 = 0
               elif index['HeII_1640'] == 0 or index['CIII_1909'] == 0:
                  CHI_C34He2 = tol_max
               else:
                  CHI_C34He2 = ((index['CIII_1909']+index['CIV_1549'])/index['HeII_1640'] - C34He2_obs)**2/((index['CIII_1909']+index['CIV_1549'])/index['HeII_1640'])
               if C34Hb_obs == 0:
                  CHI_C34Hb = 0
               elif index['CIII_1909'] == 0:
                  CHI_C34Hb = tol_max
               else:
                  CHI_C34Hb = (index['CIII_1909']+index['CIV_1549'] - C34Hb_obs)**2/(index['CIII_1909']+index['CIV_1549'])
               if C3C4_obs == 0:
                  CHI_C3C4 = 0
               elif index['CIV_1549'] == 0 or index['CIII_1909'] == 0:
                  CHI_C3C4 = tol_max
               else:
                  CHI_C3C4 = (index['CIII_1909']/index['CIV_1549'] - C3C4_obs)**2/(index['CIII_1909']/index['CIV_1549'])

               if C34Hb_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34Hb**2  + CHI_C3C4**2)**0.5
               else:	
                  CHI_OH = (CHI_ROIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_C3C4**2 )**0.5



            if CHI_OH == 0:
               OH_p = OH_p
               logU_p = logU_p
               den_OH = den_OH
            else:
               OH_p = index['12logOH'] / np.exp(CHI_OH) + OH_p
               logU_p = index['logU'] / np.exp(CHI_OH) + logU_p
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
            if CO > -10 and np.abs(index['logCO'] - CO) > np.abs(eCO+res_CO):
               continue
            if CIV_1549_obs > 0 and index['CIV_1549'] == 0:
               continue
            if HeII_1640_obs > 0 and index['HeII_1640'] == 0:
               continue
            else:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index['OIII_1665'] == 0 or index['OIII_5007'] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index['OIII_5007']/index['OIII_1665'] - ROIII_obs)**2/(index['OIII_5007']/index['OIII_1665'])
               if C34_obs == 0:
                  CHI_C34 = 0
               elif index['Lya_1216'] == 0 or index['CIII_1909'] == 0:
                  CHI_C34 = tol_max
               else:
                  CHI_C34 = ((index['CIII_1909']+index['CIV_1549'])/index['Lya_1216'] - C34_obs)**2/((index['CIII_1909']+index['CIV_1549'])/index['Lya_1216'])
               if C34He2_obs == 0:
                  CHI_C34He2 = 0
               elif index['HeII_1640'] == 0 or index['CIII_1909'] == 0:
                  CHI_C34He2 = tol_max
               else:
                  CHI_C34He2 = ((index['CIII_1909']+index['CIV_1549'])/index['HeII_1640'] - C34He2_obs)**2/((index['CIII_1909']+index['CIV_1549'])/index['HeII_1640'])
               if C34Hb_obs == 0:
                  CHI_C34Hb = 0
               elif index['CIII_1909'] == 0:
                  CHI_C34Hb = tol_max
               else:
                  CHI_C34Hb = (index['CIII_1909']+index['CIV_1549'] - C34Hb_obs)**2/(index['CIII_1909']+index['CIV_1549'])
               if C3C4_obs == 0:
                  CHI_C3C4 = 0
               elif index['CIV_1549'] == 0 or index['CIII_1909'] == 0:
                  CHI_C3C4 = tol_max
               else:
                  CHI_C3C4 = (index['CIII_1909']/index['CIV_1549'] - C3C4_obs)**2/(index['CIII_1909']/index['CIV_1549'])


               if C34Hb_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34Hb**2 + CHI_C3C4**2)**0.5
               else:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_C3C4**2  )**0.5

            if CHI_OH == 0:
               OH_e = OH_e
               logU_e = logU_e
               den_OH_e = den_OH_e
            else:
               OH_e = (index['12logOH'] - OH)**2 /np.exp(CHI_OH) + OH_e
               logU_e = (index['logU'] - logU)**2 /np.exp(CHI_OH) + logU_e
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
            igrid = igrid[np.lexsort((igrid['logCO'],igrid['logU']))]
            igrid = interpolate(igrid,0,OH-eOH-0.1,OH+eOH+0.1,10)
         if CO == -10:
            igrid = igrid
         else:
            igrid = igrid[np.lexsort((igrid['12logOH'],igrid['logU']))]
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
            elif index['OIII_1665'] == 0 or index['OIII_5007'] == 0:
               CHI_ROIII = tol_max
            else:
               CHI_ROIII = (index['OIII_5007']/index['OIII_1665'] - ROIII_obs)**2/(index['OIII_5007']/index['OIII_1665'])
            if C3O3_obs == -10:
               CHI_C3O3 = 0          
            elif index['CIII_1909'] == 0 or index['OIII_1665'] == 0:
               CHI_C3O3 = tol_max
            else:
               CHI_C3O3 =(np.log10((index['CIII_1909'])/index['OIII_1665']) - C3O3_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_1665']+1e-5))
            if C34_obs == 0:
               CHI_C34 = 0
            elif index['Lya_1216'] == 0:
               CHI_C34 = tol_max
            else:
               CHI_C34 = ((index['CIV_1549']+index['CIII_1909'])/index['Lya_1216'] - C34_obs)**2/((index['CIV_1549']+index['CIII_1909'])/index['Lya_1216'])
            if C34Hb_obs == 0:
               CHI_C34Hb = 0
            elif index['CIV_1549'] == 0:
               CHI_C34Hb = tol_max
            else:
               CHI_C34Hb = (index['CIV_1549']+index['CIII_1909'] - C34_obs)**2/(index['CIV_1549']+index['CIII_1909'])
            if C3C4_obs == 0:
               CHI_C3C4 = 0
            elif index['CIII_1909'] == 0 or index['CIV_1549'] == 0:
               CHI_C3C4 = tol_max
            else:
               CHI_C3C4 = (index['CIV_1549']/index['CIII_1909'] - C3C4_obs)**2/(index['CIV_1549']/index['CIII_1909'])


               if C34Hb_obs > 0:
                  CHI_OH = (CHI_ROIII**2 + CHI_C34Hb**2 + CHI_C3C4**2)**0.5
               else:	
                  CHI_OH = (CHI_ROIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_C3C4**2 )**0.5



            if CHI_OH == 0:
               OH_p = OH_p
               logU_p = logU_p
               den_OH = den_OH
            else:
               OH_p = index['12logOH'] /np.exp(CHI_OH) + OH_p
               logU_p = index['logU'] /np.exp(CHI_OH) + logU_p
               den_OH = 1 /np.exp(CHI_OH) + den_OH

            CHI_CO = (CHI_ROIII**2 + CHI_C3O3**2 )**0.5

            if CHI_CO == 0:
               CO_p = CO_p
               den_CO = den_CO
            else:
               CO_p = index['logCO'] /np.exp(CHI_CO)**2 + CO_p
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


lineas_header = [' HII-CHI-mistry_UV v.4.2 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'Library file used : '+file_lib_2, 'Template used to constraint grid of models: '+const_file,'']


line_label = '{:10}  '.format(output.dtype.names[0])
for ind2 in range(1, len(output.dtype.names)-7):
   line_label += '{:10}  '.format(output.dtype.names[ind2])


line_label += '{:10}   {:10}   {:10}   {:10}   {:10}   {:10}   {:10}'.format('i', 'O/H',    'eO/H',  'C/O' ,   'eC/O' , 'logU',   'elogU')
lineas_header.append(line_label)
header = '\n'.join(lineas_header)
np.savetxt('.'.join(input00.split('.')[:-1])+'_hcm-uv-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*(len(output.dtype.names)-8)+['%i']+['%.2f']*6), header=header)

lines_stor = []
with open('.'.join(input00.split('.')[:-1])+'_hcm-uv-output.dat', 'r+') as output_file:
   for line in output_file:
      lines_stor.append(line)

file_overwrite = open('.'.join(input00.split('.')[:-1])+'_hcm-uv-output.dat', 'r+')
file_overwrite.seek(0)
for line_n in lines_stor:  
   if line_n[0] == '#' and line_n[2:4] == 'ID':
      file_overwrite.write(line_n[2:])
   else:
      file_overwrite.write(line_n)
file_overwrite.truncate()   
file_overwrite.close()
print ('________________________________')
print ('Results are stored in ' + '.'.join(input00.split('.')[:-1]) + '_hcm-uv-output.dat')


