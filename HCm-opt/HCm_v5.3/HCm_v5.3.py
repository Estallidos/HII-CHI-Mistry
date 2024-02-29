# Filename: HCm_v 5.3.py



import string
import numpy as np
import sys
sys.stderr = open('errorlog.txt', 'w')
import warnings
warnings.filterwarnings("ignore")

def interpolate(grid,z,zmin,zmax,n):
   #Columns of the library
   n_comments = 0
   with open('Libraries_opt/C17_POPSTAR_1myr.dat', 'r') as file1:
      for line in file1:
         if line[0] == '#':
            n_comments += 1
   auxiliar_labels = np.genfromtxt('Libraries_opt/C17_POPSTAR_1myr.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
   ncol = len(auxiliar_labels)
   vec = []
   if z == 2:
      label_z = 'logU'
   if z == 1:
      label_z = 'logNO'
   if z == 0:
      label_z = '12logOH'
   type_list_names = []
   for col in auxiliar_labels:
      inter = 0
      no_inter = 0
      type_list_names.append((col, float))
      for row in range(0,len(grid)):
         if grid[label_z][row] < zmin or grid[label_z][row] > zmax: continue
         if z == 2: x = '12logOH'; y = 'logNO'
         if z == 1: x = '12logOH'; y = 'logU'
         if z == 0: x = 'logNO'; y = 'logU'
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

################################
###### INITIAL ITERATIONS ######
################################

#Description of the code

print (' ---------------------------------------------------------------------')
print ('This is HII-CHI-mistry v. 5.3')
print (' See Perez-Montero, E. (2014) for details')
print  (' Insert the name of your input text file with all or some of the following columns:')
print (' 3727 [OII]')
print (' 3868 [NeIII]')
print (' 4363 [OIII]')
print (' 4959 [OIII]')
print (' 5007 [OIII]')
print (' 5755 [NII]')
print (' 6312 [SIII]')
print (' 6584 [NII]')
print (' 6717,31 [SII]')
print (' 7319,30 [OII]')

print (' 9069 [SIII]')
print (' 9532 [SIII]')
print ('with their corresponding labels and errors in adjacent columns')
print ('---------------------------------------------------------------------')



#Input file reading
if len(sys.argv) == 1:
   if int(sys.version[0]) < 3:
      input00 = raw_input('Insert input file name: ')
   else:
      input00 = input('Insert input file name: ')
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
   print ('The input file is: '+input00)
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


# Interface with the user
print ('')
question = True
while question:
   print ('-------------------------------------------------')
   print ('Default SEDs')
   print ('------------')
   print ('(1) POPSTAR with Chabrier IMF, age = 1 Myr')
   print ('(2) BPASS v.2.1 a_IMF = 1.35, Mup = 300, age = 1Myr')
   print ('(3) AGN, double component,  a(UV) = -1.0')
   print (' ')
   print ('Other SED')
   print ('---------')
   print ('(4) Different library')
   print ('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      sed = raw_input('Choose SED of the models: ')
   else:
      sed = input('Choose SED of the models: ')
   if sed == '1' or sed == '2' or sed == '3' or sed == '4': question = False 
print ('')

# Further questions on the AGN models
if sed == '3':
   #SLOPE ALPHA
   question = True
   while question:
      if int(sys.version[0]) < 3:
         alpha = raw_input('Choose value for alpha(OX) in the AGN models: [1] -0.8 [2] -1.2: ')
      else:
         alpha = input('Choose value for alpha(OX) in the AGN models: [1] -0.8 [2] -1.2:  ')
      if alpha == '1' or alpha == '2': question = False
      alpha = int(alpha)
   print ('')
   #Fraction of free electrons (stopping criteria in the models)
   question = True
   while question:
      if int(sys.version[0]) < 3:
         efrac = raw_input('Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons: ')
      else:
         efrac = input('Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons:  ')
      if efrac == '1' or efrac == '2': question = False
      efrac = int(efrac)
   print ('')
 
#Particular file introduced by the user
if sed == '4':
   question = True
   while question:
      print ('Introduce name of the file containing the models. It must be located in the folder "Libraries_opt".')
      print ('')
      if int(sys.version[0]) < 3:
         new_library = raw_input('Name of file: ')
      else:
         new_library = input('Name of file: ')
 
      #Searching for the file
      try:
         #Counting comments:
         n_comments = 0
         with open('Libraries_opt/'+new_library, 'r') as file3:
            for line in file3:
               if line[0] == '#':
                  n_comments += 1
         library_user = np.genfromtxt('Libraries_opt/'+new_library, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         print ('')
         print ('Loading library '+new_library+'. Checking correct format of the file.')
         question = False
      except:
         print ('')
         print ('Library was not found in folder "Libraries_opt" or file does not exist.')
   question = True
   while question:
      try:
         #Counting comments:
         n_comments = 0
         with open('Libraries_opt/'+new_library, 'r') as file4:
            for line in file4:
               if line[0] == '#':
                  n_comments += 1
         library_user = np.genfromtxt('Libraries_opt/'+new_library, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         #Checking correct format:
         #Counting comments:
         n_comments = 0
         with open('Libraries_opt/C17_POPSTAR_1myr.dat', 'r') as file5:
            for line in file5:
               if line[0] == '#':
                  n_comments += 1
         auxiliar_labels = np.genfromtxt('Libraries_opt/C17_POPSTAR_1myr.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
         missing_labels = []
         for label in auxiliar_labels:
            if label in library_user.dtype.names:
               continue
            else:
               missing_labels.append(label)
         #Displaying message for the user:
         print ('Succesfully reading of the file')
         print ('')
         if len(missing_labels) == 0:
            print ('File presents the correct format')
            question = False
         else:
            print ('File does not present the correct format. The following columns are missing:')
            for need_label in missing_labels:
               print ('- '+need_label)
            print ('More details on the correct format for the library are found in readme file.')
            print ('')
            print ('Reintroduce name of the file with fixed format.')
            print ('')
            if int(sys.version[0]) < 3:
               new_library = raw_input('Name of file: ')
            else:
               new_library = input('Name of file: ')
      except:
         print ('Something went wrong while reading file. Please, reintroduce name of the file.')
         print ('')
         if int(sys.version[0]) < 3:
            new_library = raw_input('Name of file: ')
         else:
            new_library = input('Name of file: ')
            
#Interpolation in the grid of models
question = True
print ('')
while question:
   if int(sys.version[0]) < 3:
      inter = raw_input('Choose models: [0] No interpolated [1] Interpolated: ')
   else:
      inter = input('Choose models: [0] No interpolated [1] Interpolated: ')
   if inter == '0' or inter == '1': question = False
print ('')

sed = int(sed)
inter = int(inter)

#POPSTAR MODEL
if sed==1:
   file_lib = 'C17_POPSTAR_1myr.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_opt/'+file_lib, 'r') as file6:
      for line in file6:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_opt/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
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
   with open('Libraries_opt/'+file_lib, 'r') as file7:
      for line in file7:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_opt/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)  
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
      

#AGN MODEL FOR alpha_OX = -0.8, efrac = 2% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==1 and efrac == 1:
   file_lib = 'C17_AGN_alpha08_efrac02_CNfix.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_opt/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_opt/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 2%. No interpolation. No constraint in ionization parameter.'
      print ('No interpolation for the AGN a(ox) = -0.8 with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O. ')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 2% interpolated. No constraint in ionization parameter.'
      print ('Interpolation for the AGN a(ox) = -0.8 with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
      res_NO = 0.125
    
#AGN MODEL FOR alpha_OX = -0.8, efrac = 98% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==1 and efrac == 2:
   file_lib = 'C17_AGN_alpha08_efrac98_CNfix.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_opt/'+file_lib, 'r') as file9:
      for line in file9:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_opt/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 98%. No interpolation. No constraint in ionization parameter.'
      print ('No interpolation for the AGN a(ox) = -0.8 with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 98% interpolated. No constraint in ionization parameter.'
      print ('Interpolation for the AGN a(ox) = -0.8 with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
      res_NO = 0.125
    
#AGN MODEL FOR alpha_OX = -1.2, efrac = 2% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==2 and efrac == 1:
   file_lib = 'C17_AGN_alpha12_efrac02_CNfix.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_opt/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_opt/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 2%. No interpolation. No constraint in ionization parameter.'
      print ('No interpolation for the AGN a(ox) = -1.2 with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 2% interpolated. No constraint in ionization parameter.'
      print ('Interpolation for the AGN a(ox) = -1.2 with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
      res_NO = 0.125
 
#AGN MODEL FOR alpha_OX = -1.2, efrac = 98% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==2 and efrac == 2:
   file_lib = 'C17_AGN_alpha12_efrac98_CNfix.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_opt/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_opt/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 98%. No interpolation. No constraint in ionization parameter.'
      print ('No interpolation for the AGN a(ox) = -1.2 with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 98% interpolated. No constraint in ionization parameter.'
      print ('Interpolation for the AGN a(ox) = -1.2 with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
        
#Different library
elif sed==4:
   file_lib = new_library
   #Counting comments:
   n_comments = 0
   with open('Libraries_opt/'+file_lib, 'r') as file12:
      for line in file12:
         if line[0] == '#':
            n_comments += 1  
   grid_aux = np.genfromtxt('Libraries_opt/'+file_lib, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'User file ' + file_lib + ' used as library for the models no interpolated'
      print ('No interpolation for the library '+file_lib)
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'User file ' + file_lib + ' used as library for the models interpolated'
      print ('Interpolation for the library '+file_lib)

#Valuable columns of the files
opt_lin = ['12logOH', 'logNO', 'logU', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007', 'NII_5755', 'SIII_6312', 'NII_6584', 'SII_671731', 'OII_7325', 'SIII_9069']
lin_opt_label = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007', 'NII_5755', 'SIII_6312', 'NII_6584', 'SII_6717,31', 'OII_7325', 'SIII_9069']

########################################
###### SORTING THE GRID OF MODELS ######
########################################

print (' ')
print ('Sorting the grid of models')
print (' ')

index_OH_NO_U_sorted = [] #storing the correct order of the indexes

#Sorting abundances 12+log(O/H)
OH_values = grid_aux['12logOH'] #Oxygen abundances
if len(OH_values) != 1:
   sorted_list_OH = sorted(range(len(OH_values)),key=OH_values.__getitem__)
if len(OH_values) == 1:
   sorted_list_OH = [0]

#Sorting abundance ratios log(N/O)
OH_values_diff = list(set(OH_values[sorted_list_OH]))
OH_values_diff.sort() #It is necessary to sort again the list of different elements
for OH_num in OH_values_diff:
   index_OH_fix = np.where(OH_values == OH_num)[0] #Index(es) for a particular abundance 12+log(O/H)
   NO_values = grid_aux['logNO'][index_OH_fix]
   if len(NO_values) != 1:
      sorted_list_NO = sorted(range(len(NO_values)), key=NO_values.__getitem__)
   if len(NO_values) == 1:
      sorted_list_NO = [0]
   NO_values_diff = list(set(NO_values[sorted_list_NO]))
   NO_values_diff.sort() #It s necessary to sort again the list of different elements
   for NO_num in NO_values_diff:
      index_OH_NO_fix = np.where(NO_values == NO_num)[0] #Index(es) for particular abundances 12+log(O/H) and log(N/O)
      #Sorting ionization parameters
      U_values = grid_aux['logU'][index_OH_fix[index_OH_NO_fix]]
      if len(U_values) != 1:
         sorted_list_U = sorted(range(len(U_values)), key=U_values.__getitem__)
      if len(U_values) == 1:
         sorted_list_U = [0]
      index_OH_NO_U = index_OH_fix[index_OH_NO_fix[sorted_list_U]] #Sorted index(es) for U at fixed O/H and N/O
      for index_sort in index_OH_NO_U:
         index_OH_NO_U_sorted.append(index_sort) #Adding index in the correct order

#Generating new library file
list_comments = [] #Storing comments in the file:
with open('Libraries_opt/'+file_lib, 'r') as file_aux:
   for line in file_aux:
      if line[0] == '#':
         list_comments.append(line)

#Storing columns:
lin_opt_col = []
#Retrieving each column of the grid
for label in opt_lin:
   aux_col = grid_aux[label].tolist()
   lin_opt_col.append(aux_col)

#Comments
grid_to_write = open('Libraries_opt/'+file_lib, 'w')
for line_com in list_comments:
   grid_to_write.write(line_com)
#Header line
label_line = '{:15} '.format(lin_opt_label[0].replace(' ',''))
for ind in range(1, len(lin_opt_label)-1):
   label_line += '\t {:15} '.format(lin_opt_label[ind].replace(' ',''))
label_line += '\t {:15}\n'.format(lin_opt_label[-1].replace(' ','')) 
grid_to_write.write(label_line)
#Values:
for ind_val in index_OH_NO_U_sorted:
   val_line = '{:7.7f} '.format(lin_opt_col[0][ind_val])
   for ind2 in range(1, len(lin_opt_label)-1):
      val_line += '\t {:7.7f} '.format(lin_opt_col[ind2][ind_val])
   val_line += '\t {:7.7f}\n'.format(lin_opt_col[-1][ind_val])
   grid_to_write.write(val_line)        
grid_to_write.close()

#Opening sorted grid of models
n_comments = 0
with open('Libraries_opt/'+file_lib, 'r') as file12:
   for line in file12:
      if line[0] == '#':
         n_comments += 1  
grid_aux = np.genfromtxt('Libraries_opt/'+file_lib, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

################################################
###### CONSTRAINTS FOR THE GRID OF MODELS ######
################################################
    
#Reading constraints and creating library with constraints
print ('')
print ('Select a file with the constraints to be used to limit the grid of models when the measurement of a quantity is impossible without any relation.')
print ('')

#Displayig options for the user
print ('')
question = True
while question:
   print ('-------------------------------------------------')
   print ('Default constraints')
   print ('-------------------')
   print ('(1) Constraints for Star-Forming Galaxies')
   print ('(2) Constraints for Extreme Emission Line Galaxies')
   print ('(3) Constraints for AGNs (no restriction in the ionization parameter)')
   print ('(4) Constraints for high ionization AGNs (log(U) > -2.5)')
   print ('(5) Constraints for low ionization AGNs (log(U) < -2.5)')
   print ('')
   print ('Other constraints')
   print ('-----------------')
   print ('(6) Different constraint file')
   print ('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      const = raw_input('Choose constraint for the grids: ')
   else:
      const = input('Choose constraint for the grids: ')
   if const == '1' or const == '2' or const == '3' or const == '4' or const == '5' or const == '6': question = False 
print ('')

#Particular file introduced by the user
if const == '6':
   question = True
   while question:
      print ('Introduce name of the file containing the constraints for the grids. It must be located in the folder "Constraints".')
      print ('')
      if int(sys.version[0]) < 3:
         new_const = raw_input('Name of file: ')
      else:
         new_const = input('Name of file: ')
 
      #Searching for the file
      try:
         #Counting comments:
         n_comments = 0
         with open('Constraints/'+new_const, 'r') as file13:
            for line in file13:
               if line[0] == '#':
                  n_comments += 1
         const_user = np.genfromtxt('Constraints/'+new_const, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         print ('')
         print ('Loading constraint file '+new_const+'. Checking correct format of the file.')
         question = False
      except:
         print ('')
         print ('File was not found in folder "Constraints" or file does not exist.')
   question = True
   while question:
      try:
         #Counting comments:
         n_comments = 0
         with open('Constraints/'+new_const, 'r') as file14:
            for line in file14:
               if line[0] == '#':
                  n_comments += 1
         const_user = np.genfromtxt('Constraints/'+new_const, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         #Checking correct format:
         #Counting comments:
         n_comments = 0
         with open('Constraints/template_OH.dat', 'r') as file15:
            for line in file15:
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
         print ('Succesfully reading of the file')
         if len(missing_labels) == 0:
            print ('File presents the correct format')
            question = False
         else:
            print ('File does not present the correct format. The following columns are missing:')
            for need_label in missing_labels:
               print ('- '+need_label)
            print ('More details on the correct format for the library are found in readme file.')
            print ('')
            print ('Reintroduce name of the file with fixed format.')
            print ('')
            if int(sys.version[0]) < 3:
               new_const = raw_input('Name of file: ')
            else:
               new_const = input('Name of file: ')
      except:
         print('Something went wrong while reading file. Please, reintroduce name of the file.')
         print('')
         if int(sys.version[0]) < 3:
            new_const = raw_input('Name of file: ')
         else:
            new_const = input('Name of file: ')

#Generation of grids with constrained laws:
if const == '1' or const == '2' or const == '3' or const == '6':
   #First grid does not change
   grid1 = grid_aux
   file_lib_2 = file_lib
elif const == '4' or const == '5':
   lin_opt_agn = []
   #The initial grid need to be constrained in the ionization parameter
   if const == '4':
      U_max = 0.0
      U_min = -2.5
      tag = 'high'
   if const == '5':
      U_max = -2.5
      U_min = -4.0
      tag = 'low'
   #Retrieving each column of the grid
   for label in opt_lin:
      aux_col = grid_aux[label].tolist()
      lin_opt_agn.append(aux_col)
   #Creation of the grid
   file_lib_2 = '.'.join(file_lib.split('.')[0:-1])+'_'+tag+'.'+file_lib.split('.')[-1]
   file_open = open('Libraries_opt/'+file_lib_2, 'w')
   file_open.write('#Library constrained for '+tag+' ionization AGNs\n')
   #Header line
   label_line = '{:15} '.format(lin_opt_label[0].replace(' ',''))
   for ind in range(1, len(lin_opt_label)-1):
      label_line += '\t {:15} '.format(lin_opt_label[ind].replace(' ',''))
   label_line += '\t {:15}\n'.format(lin_opt_label[-1].replace(' ','')) 
   file_open.write(label_line)
   #Values:
   for ind_val in range(0, len(lin_opt_agn[0])):
      if lin_opt_agn[2][ind_val] <= U_max and lin_opt_agn[2][ind_val] >= U_min:
         val_line = '{:7.7f} '.format(lin_opt_agn[0][ind_val])
         for ind2 in range(1, len(lin_opt_label)-1):
            val_line += '\t {:7.7f} '.format(lin_opt_agn[ind2][ind_val])
         val_line += '\t {:7.7f}\n'.format(lin_opt_agn[-1][ind_val])
         file_open.write(val_line)        
   file_open.close()
   #Counting comments:
   n_comments = 0
   with open('Libraries_opt/'+file_lib_2, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid1 = np.genfromtxt('Libraries_opt/'+file_lib_2,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

#Generating libraries for the constraints in the files
if const == '1': #Star-Forming Galaxies
   name_const = 'Constraints/template_OH.dat'
   const_file = 'template_OH.dat'
   n_comments = 0
   with open(name_const, 'r') as file16:
      for line in file16:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == '2':
   name_const = 'Constraints/template_OH_eelg.dat'
   const_file = 'template_OH_eelg.dat'
   n_comments = 0
   with open(name_const, 'r') as file17:
      for line in file17:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == '3':
   name_const = 'Constraints/template_OH_agn.dat'
   const_file = 'template_OH_agn.dat'
   n_comments = 0
   with open(name_const, 'r') as file18:
      for line in file18:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == '4':
   name_const = 'Constraints/template_OH_agn_high.dat'
   const_file = 'template_OH_agn_high.dat'
   n_comments = 0
   with open(name_const, 'r') as file19:
      for line in file19:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == '5':
   name_const = 'Constraints/template_OH_agn_low.dat'
   const_file = 'template_OH_agn_low.dat'
   n_comments = 0
   with open(name_const, 'r') as file20:
      for line in file20:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == '6':
   const_file = new_const
   name_const = 'Constraints/'+new_const
   n_comments = 0
   with open(name_const, 'r') as file21:
      for line in file21:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

#Limiting the grids:
lin_opt_val = [] #The initial grid need to be constrained in the ionization parameter
#Retrieving each column of the grid
for label in opt_lin:
   aux_col = grid1[label].tolist()
   lin_opt_val.append(aux_col)
   #Creation of the grids
   name_OH_U = '.'.join(file_lib_2.split('.')[0:-1])+'_OH_U_constrained.'+file_lib.split('.')[-1]
   name_OH_U_NO = '.'.join(file_lib_2.split('.')[0:-1])+'_OH_U_NO_constrained.'+file_lib.split('.')[-1]
   file_open = open('Libraries_opt/'+ name_OH_U, 'w') #OH and U relation
   file_open_2 = open('Libraries_opt/'+name_OH_U_NO, 'w') #OH, NO and U relation
   file_open.write('#Constrained by relation between 12+log(O/H) and log(U)\n')
   file_open_2.write('#Constrained by relation between 12+log(O/H), log(U) and log(N/O)\n')
   #Header line
   label_line = '{:15} '.format(lin_opt_label[0].replace(' ',''))
   for ind in range(1, len(lin_opt_label)-1):
      label_line += '\t {:15} '.format(lin_opt_label[ind].replace(' ',''))
   label_line += '\t {:15}\n'.format(lin_opt_label[-1].replace(' ','')) 
   file_open.write(label_line)
   file_open_2.write(label_line)
#Values:
for ind_val in range(0, len(lin_opt_val[0])):
   index_desired = np.where(const_data['12logOH'] == lin_opt_val[0][ind_val])[0][0] #Searching for constrain in given value of O/H
   if lin_opt_val[2][ind_val] <= const_data['logU_max'][index_desired] and lin_opt_val[2][ind_val] >= const_data['logU_min'][index_desired]:
      val_line = '{:7.7f} '.format(lin_opt_val[0][ind_val])
      for ind2 in range(1, len(lin_opt_label)-1):
         val_line += '\t {:7.7f} '.format(lin_opt_val[ind2][ind_val])
      val_line += '\t {:7.7f}\n'.format(lin_opt_val[-1][ind_val])
      file_open.write(val_line)
   if lin_opt_val[2][ind_val] <= const_data['logU_max'][index_desired] and lin_opt_val[2][ind_val] >= const_data['logU_min'][index_desired] and lin_opt_val[1][ind_val] <= const_data['logNO_max'][index_desired] and lin_opt_val[1][ind_val] >= const_data['logNO_min'][index_desired]:
      val_line = '{:7.7f} '.format(lin_opt_val[0][ind_val])
      for ind2 in range(1, len(lin_opt_label)-1):
         val_line += '\t {:7.7f} '.format(lin_opt_val[ind2][ind_val])
      val_line += '\t {:7.7f}\n'.format(lin_opt_val[-1][ind_val])
      file_open_2.write(val_line)
file_open.close()
file_open_2.close()
#Counting comments:
n_comments = 0
with open('Libraries_opt/'+name_OH_U, 'r') as file22:
   for line in file22:
      if line[0] == '#':
         n_comments += 1
grid2 = np.genfromtxt('Libraries_opt/'+name_OH_U,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
n_comments = 0
with open('Libraries_opt/'+name_OH_U_NO, 'r') as file23:
   for line in file23:
      if line[0] == '#':
         n_comments += 1
grid3 = np.genfromtxt('Libraries_opt/'+name_OH_U_NO,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

#Residual in NO
if inter==0:
   res_NO = np.max([sorted(set(grid1['logNO']))[ind+1]-sorted(set(grid1['logNO']))[ind] for ind in range(0, len(set(grid1['logNO']))-1)])
if inter==1:
   res_NO = np.max([sorted(set(grid1['logNO']))[ind+1]-sorted(set(grid1['logNO']))[ind] for ind in range(0, len(set(grid1['logNO']))-1)])/10

###########################################
###### SUMMARY OF THE GRID OF MODELS ######
###########################################

print ('-------------------------------------------------')
print ('Summary of the models')
print ('---------------------')
print ('Libraries generated with the constraint file: '+const_file+'. The following grids are going to be used:')
print ('- Full library (Grid#1): '+file_lib_2)
print ('       Total number of models: ' + str(len(grid1)))
print ('- Library constrained by 12+log(O/H) - log(U) relation (Grid#2): '+name_OH_U)
print ('       Total number of models: ' + str(len(grid2)))
print ('- Library constrained by 12+log(O/H) - log(U) - log(N/O) relation (Grid#3): '+name_OH_U_NO)
print ('       Total number of models: ' + str(len(grid3)))
print ('-------------------------------------------------')
print ('')

#################################################
###### CREATING ARRAY TO STORE ESTIMATIONS ######
#################################################

grids = []
OHffs = []
eOHffs = []
NOffs = []
eNOffs = []
logUffs = []
elogUffs = []

#Labels to check information provided in the input file

Label_ID = False
Label_OII = False
Label_eOII = False
Label_OIIa = False
Label_eOIIa = False
Label_NeIII = False
Label_eNeIII = False
Label_OIIIa = False
Label_eOIIIa = False
Label_OIII_4959 = False
Label_eOIII_4959 = False
Label_OIII_5007 = False
Label_eOIII_5007 = False
Label_NIIa = False
Label_eNIIa = False
Label_SIIIa = False
Label_eSIIIa = False
Label_NII = False
Label_eNII = False
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

for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == 'OII_3727':
      Label_OII = True
   if input1.dtype.names[col] == 'eOII_3727':
      Label_eOII = True
   if input1.dtype.names[col] == 'OII_7325':
      Label_OIIa = True
   if input1.dtype.names[col] == 'eOII_7325':
      Label_eOIIa = True
   if input1.dtype.names[col] == 'NeIII_3868':
      Label_NeIII = True
   if input1.dtype.names[col] == 'eNeIII_3868':
      Label_eNeIII = True
   if input1.dtype.names[col] == 'OIII_4363':
      Label_OIIIa = True
   if input1.dtype.names[col] == 'eOIII_4363':
      Label_eOIIIa = True
   if input1.dtype.names[col] == 'OIII_4959':
      Label_OIII_4959 = True
   if input1.dtype.names[col] == 'eOIII_4959':
      Label_eOIII_4959 = True
   if input1.dtype.names[col] == 'OIII_5007':
      Label_OIII_5007 = True
   if input1.dtype.names[col] == 'eOIII_5007':
      Label_eOIII_5007 = True
   if input1.dtype.names[col] == 'NII_5755':
      Label_NIIa = True
   if input1.dtype.names[col] == 'eNII_5755':
      Label_eNIIa = True
   if input1.dtype.names[col] == 'SIII_6312':
      Label_SIIIa = True
   if input1.dtype.names[col] == 'eSIII_6312':
      Label_eSIIIa = True
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
   if input1.dtype.names[col] == 'SIII_9069':
      Label_SIII_9069 = True
   if input1.dtype.names[col] == 'eSIII_9069':
      Label_eSIII_9069 = True
   if input1.dtype.names[col] == 'SIII_9532':
      Label_SIII_9532 = True
   if input1.dtype.names[col] == 'eSIII_9532':
      Label_eSIII_9532 = True


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
if Label_OIIa == False:
   OII_7325 = np.zeros(input1.size)
else:
   OII_7325 = input1['OII_7325']
if Label_eOIIa == False:
   eOII_7325 = np.zeros(input1.size)
else:
   eOII_7325 = input1['eOII_7325']
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
if Label_OIII_4959 == False:
    OIII_4959 = np.zeros(input1.size)
else:
    OIII_4959 = input1['OIII_4959']
if Label_eOIII_4959 == False:
    eOIII_4959 = np.zeros(input1.size)
else:
    eOIII_4959 = input1['eOIII_4959']
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
   eOIII_5007 = 3*input1['eOIII_4959']
else:
   eOIII_5007 = (input1['eOIII_5007']+input1['eOIII_4959'])/1.33
if Label_NIIa == False:
   NII_5755 = np.zeros(input1.size)
else:
   NII_5755 = input1['NII_5755']
if Label_eNIIa == False:
   eNII_5755 = np.zeros(input1.size)
else:
   eNII_5755 = input1['eNII_5755']
if Label_SIIIa == False:
   SIII_6312 = np.zeros(input1.size)
else:
   SIII_6312 = input1['SIII_6312']
if Label_eSIIIa == False:
   eSIII_6312 = np.zeros(input1.size)
else:
   eSIII_6312 = input1['eSIII_6312']
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



#Creation of output only with information from inputs
aux_list = []
aux_list.append(('ID','U12'))
if Label_OII == True:
   aux_list.append(('OII_3727', float))
if Label_eOII == True:
   aux_list.append(('eOII_3727', float))
if Label_NeIII == True:
   aux_list.append(('NeIII_3868', float))
if Label_eNeIII == True:
   aux_list.append(('eNeIII_3868', float))
if Label_OIIIa == True:
   aux_list.append(('OIII_4363', float))
if Label_eOIIIa == True:
   aux_list.append(('eOIII_4363', float))
if Label_OIII_4959 == True:
   aux_list.append(('OIII_4959', float))
if Label_eOIII_4959 == True:
   aux_list.append(('eOIII_4959', float))
if Label_OIII_5007 == True:
   aux_list.append(('OIII_5007', float))
if Label_eOIII_5007 == True:
   aux_list.append(('eOIII_5007', float))
if Label_NIIa == True:
   aux_list.append(('NII_5755', float))
if Label_eNIIa == True:
   aux_list.append(('eNII_5755', float))
if Label_SIIIa == True:
   aux_list.append(('SIII_6312', float))
if Label_eSIIIa == True:
   aux_list.append(('eSIII_6312', float))
if Label_NII == True:
   aux_list.append(('NII_6584', float))
if Label_eNII == True:
   aux_list.append(('eNII_6584', float))
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
if Label_OIIa == True:
   aux_list.append(('OII_7325', float))
if Label_eOIIa == True:
   aux_list.append(('eOII_7325', float))
if Label_SIII_9069 == True:
   aux_list.append(('SIII_9069', float))
if Label_eSIII_9069 == True:
   aux_list.append(('eSIII_9069', float))
if Label_SIII_9532 == True:
   aux_list.append(('SIII_9532', float))
if Label_eSIII_9532 == True:
   aux_list.append(('eSIII_9532', float))

aux_list.append(('grid', int))
aux_list.append(('OH', float))
aux_list.append(('eOH', float))
aux_list.append(('NO', float))
aux_list.append(('eNO', float))
aux_list.append(('logU', float))
aux_list.append(('elogU', float))
output = np.zeros(input1.size, dtype=aux_list)

output['ID'] = Names
if Label_OII == True:
   output['OII_3727'] = OII_3727
if Label_eOII == True:
   output['eOII_3727'] = eOII_3727
if Label_OIIa == True:
   output['OII_7325'] = OII_7325
if Label_eOIIa == True:
   output['eOII_7325'] = eOII_7325
if Label_NeIII == True:
   output['NeIII_3868'] = NeIII_3868
if Label_eNeIII == True:
   output['eNeIII_3868'] = eNeIII_3868
if Label_OIIIa == True:
   output['OIII_4363'] = OIII_4363
if Label_eOIIIa == True:
   output['eOIII_4363'] = eOIII_4363
if Label_OIII_4959 == True:
   output['OIII_4959'] = OIII_4959
if Label_eOIII_4959 == True:
   output['eOIII_4959'] = eOIII_4959
if Label_OIII_5007 == True:
   output['OIII_5007'] = OIII_5007
if Label_eOIII_5007 == True:
   output['eOIII_5007'] = eOIII_5007
if Label_NIIa == True:
   output['NII_5755'] = NII_5755
if Label_eNIIa == True:
   output['eNII_5755'] = eNII_5755
if Label_SIIIa == True:
   output['SIII_6312'] = SIII_6312
if Label_eSIIIa == True:
   output['eSIII_6312'] = eSIII_6312
if Label_NII == True:
   output['NII_6584'] = NII_6584
if Label_eNII == True:
   output['eNII_6584'] = eNII_6584
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




print ('Reading grids ....')
print ('')
print ('')
print ('----------------------------------------------------------------')
print ('(%)   ID    Grid  12+log(O/H)   log(N/O)    log(U)')
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
   
   if (OIII_4363[tab] > 0 and OIII_5007[tab] > 0) or (NII_6584[tab] > 0 and NII_5755[tab] > 0) or (SIII_9069[tab] > 0 and SIII_6312[tab] > 0) or (OII_3727[tab] > 0 and OII_7325[tab] > 0):

      grid = grid1
      grid_type = 1
      grids.append(1)
   elif NII_6584[tab] > 0 and (OII_3727[tab] > 0 or SII_6725[tab] > 0 or OII_7325[tab] > 0 or OIII_5007[tab] > 0):
      grid = grid2
      grid_type = 2
      grids.append(2)
   else:
      grid = grid3
      grid_type = 3
      grids.append(3)

# Calculation of N/O

   if NII_6584[tab] <= 0 or (OII_3727[tab] <= 0 and SII_6725[tab] <= 0 and OII_7325[tab] <= 0 and OIII_5007[tab] <=0):
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
         if OII_3727[tab] <= 0:
            OII_3727_obs = 0
         else:
            while OII_3727_obs <= 0:
               OII_3727_obs = np.random.normal(OII_3727[tab],eOII_3727[tab]+1e-5)
         OII_7325_obs = 0
         if OII_7325[tab] <= 0:
            OII_7325_obs = 0
         else:
            while OII_7325_obs <= 0:
               OII_7325_obs = np.random.normal(OII_7325[tab],eOII_7325[tab]+1e-5)
         if OII_3727_obs <= 0  or OII_7325_obs <= 0:
            ROII_obs = 0
         else:
            ROII_obs = OII_3727_obs / OII_7325_obs
         OIII_4363_obs = 0
         if OIII_4363[tab] <= 0:
            OIII_4363_obs = 0
         else:
            while OIII_4363_obs <= 0:
               OIII_4363_obs = np.random.normal(OIII_4363[tab],eOIII_4363[tab]+1e-5)
         OIII_5007_obs = 0
         if OIII_5007[tab] <= 0:
            OIII_5007_obs = 0
         else:
            while OIII_5007_obs <= 0:
               OIII_5007_obs = np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-5)
         if OIII_4363_obs <= 0  or OIII_5007_obs <= 0:
            ROIII_obs = 0
         else:
            ROIII_obs = OIII_5007_obs / OIII_4363_obs
         NII_6584_obs = 0
         if NII_6584[tab] <= 0:
            NII_6584_obs = 0
         else:
            while NII_6584_obs <= 0:
               NII_6584_obs = np.random.normal(NII_6584[tab],eNII_6584[tab]+1e-5)
         NII_5755_obs = 0
         if NII_5755[tab] <= 0:
            NII_5755_obs = 0
         else:
            while NII_5755_obs <= 0:
               NII_5755_obs = np.random.normal(NII_5755[tab],eNII_5755[tab]+1e-5)
         if NII_5755_obs <= 0  or NII_6584_obs <= 0:
            RNII_obs = 0
         else:
            RNII_obs = NII_6584_obs / NII_5755_obs
         SII_6725_obs = 0
         if SII_6725[tab] <= 0:
               SII_6725_obs = 0
         else:
            while SII_6725_obs <= 0:
               SII_6725_obs = np.random.normal(SII_6725[tab],eSII_6725[tab]+1e-5)
         if NII_6584_obs <= 0 or OII_3727_obs <= 0:
            N2O2_obs = -10
         else:
            N2O2_obs = np.log10(NII_6584_obs / OII_3727_obs)
         if NII_6584_obs <= 0 or OII_7325_obs <= 0:
            N2O2a_obs = -10
         else:
            N2O2a_obs = np.log10(NII_6584_obs / OII_7325_obs)
         if NII_6584_obs <= 0 or SII_6725_obs <= 0:
            N2S2_obs = -10
         else:   
            N2S2_obs = np.log10(NII_6584_obs / SII_6725_obs)
         SIII_6312_obs = 0
         if SIII_6312[tab] <= 0:
               SIII_6312_obs = 0
         else:
            while SIII_6312_obs <= 0:
               SIII_6312_obs = np.random.normal(SIII_6312[tab],eSIII_6312[tab]+1e-5)
         SIII_9069_obs = 0
         if SIII_9069[tab] <= 0:
               SIII_9069_obs = 0
         else:
            while SIII_9069_obs <= 0:
               SIII_9069_obs = np.random.normal(SIII_9069[tab],eSIII_9069[tab]+1e-5)
         if SIII_6312_obs <= 0  or SIII_9069_obs <= 0:
            RSIII_obs = 0
         else:
            RSIII_obs = SIII_9069_obs / SIII_6312_obs
         if OIII_5007_obs <= 0 or NII_6584_obs <= 0:
            O3N2_obs = -10
         else:
            O3N2_obs = np.log10( OIII_5007_obs / NII_6584_obs )


         CHI_ROIII = 0
         CHI_RNII = 0
         CHI_RSIII = 0
         CHI_ROII = 0
         CHI_N2O2 = 0
         CHI_N2O2a = 0
         CHI_N2S2 = 0
         CHI_NO = 0
         for index in grid:
            if ROIII_obs == 0: 
               CHI_ROIII = 0
            elif index['OIII_4363'] == 0 or index['OIII_5007'] == 0:
               CHI_ROIII = tol_max
            else:   
               CHI_ROIII = (index['OIII_5007']/index['OIII_4363']- ROIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
            if ROII_obs == 0: 
               CHI_ROII = 0
            elif index['OII_3727'] == 0 or index['OII_7325'] == 0:
               CHI_ROII = tol_max
            else:   
               CHI_ROII = (index['OII_3727']/index['OII_7325']- ROII_obs)**2/(index['OII_3727']/index['OII_7325'])
            if RNII_obs == 0: 
               CHI_RNII = 0
            elif index['NII_5755'] == 0 or index['NII_6584'] == 0:
               CHI_RNII = tol_max
            else:   
               CHI_RNII = (index['NII_6584']/index['NII_5755']- RNII_obs)**2/(index['NII_6584']/index['NII_5755'])
            if RSIII_obs == 0: 
               CHI_RSIII = 0
            elif index['SIII_6312'] == 0 or index['SIII_9069'] == 0:
               CHI_RSIII = tol_max
            else:   
               CHI_RSIII = (index['SIII_9069']/index['SIII_6312']- RSIII_obs)**2/(index['SIII_9069']/index['SIII_6312'])
            if N2O2_obs == -10:
               CHI_N2O2 = 0
            elif index['OII_3727'] == 0 or index['NII_6584'] == 0:
               CHI_N2O2 = tol_max
            else:
               CHI_N2O2 =(np.log10(index['NII_6584']/index['OII_3727']) - N2O2_obs)**2/(abs(np.log10(index['NII_6584']/index['OII_3727'])+1e-5))
            if N2O2a_obs == -10:
               CHI_N2O2a = 0
            elif N2O2_obs > -10:
               CHI_N2O2a = 0
            elif index['OII_7325'] == 0 or index['NII_6584'] == 0:
               CHI_N2O2a = tol_max
            else:
               CHI_N2O2a =(np.log10(index['NII_6584']/index['OII_7325']) - N2O2a_obs)**2/(abs(np.log10(index['NII_6584']/index['OII_7325'])+1e-5))
            if N2S2_obs == -10: 
               CHI_N2S2 = 0
            elif N2O2_obs > -10:
               CHI_N2S2 = 0
            elif index['NII_6584'] == 0 or index['SII_671731'] == 0:
               CHI_N2S2 = tol_max
            else:
               CHI_N2S2 =(np.log10(index['NII_6584']/index['SII_671731']) - N2S2_obs)**2/(abs(np.log10(index['NII_6584']/index['SII_671731'])+1e-5))
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif N2O2_obs > -10:
               CHI_O3N2 = 0
            elif index['OIII_5007'] == 0 or index['NII_6584'] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index['OIII_5007']/index['NII_6584']) - O3N2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['NII_6584']+1e-5)))

            if ROIII_obs > 0:
               CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2 + CHI_O3N2**2 + CHI_N2O2a**2)**0.5
            elif RSIII_obs > 0:
               CHI_NO = (CHI_RSIII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
            elif RNII_obs > 0:
               CHI_NO = (CHI_RNII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
            elif ROII_obs > 0:
               CHI_NO = (CHI_ROII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
            else:
               CHI_NO = (CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5


            if CHI_NO == 0:
               NO_p = NO_p
               den_NO = den_NO
            else:
               NO_p = index[1] / (CHI_NO)**2 + NO_p
               den_NO = 1 / (CHI_NO)**2 + den_NO


         NO = NO_p / den_NO 

# Calculation of N/O error

   
         CHI_ROIII = 0
         CHI_ROII = 0
         CHI_RNII = 0
         CHI_RSIII = 0
         CHI_N2O2 = 0
         CHI_N2O2a = 0
         CHI_N2S2 = 0
         CHI_NO = 0
      
         for index in grid:
            if ROIII_obs == 0: 
               CHI_ROIII = 0
            elif index['OIII_4363'] == 0 or index['OIII_5007'] == 0:
               CHI_ROIII = tol_max
            else:   
               CHI_ROIII = (index['OIII_5007']/index['OIII_4363']- ROIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
            if ROII_obs == 0: 
               CHI_ROII = 0
            elif index['OII_3727'] == 0 or index['OII_7325'] == 0:
               CHI_ROII = tol_max
            else:   
               CHI_ROII = (index['OII_3727']/index['OII_7325']- ROII_obs)**2/(index['OII_3727']/index['OII_7325'])
            if RNII_obs == 0: 
               CHI_RNII = 0
            elif index['NII_5755'] == 0 or index['NII_6584'] == 0:
               CHI_RNII = tol_max
            else:   
               CHI_RNII = (index['NII_6584']/index['NII_5755']- RNII_obs)**2/(index['NII_6584']/index['NII_5755'])
            if RSIII_obs == 0: 
               CHI_RSIII = 0
            elif index['SIII_6312'] == 0 or index['SIII_9069'] == 0:
               CHI_RSIII = tol_max
            else:   
               CHI_RSIII = (index['SIII_9069']/index['SIII_6312']- RSIII_obs)**2/(index['SIII_9069']/index['SIII_6312'])
            if N2O2_obs == -10:
               CHI_N2O2 = 0
            elif index['OII_3727'] == 0 or index['NII_6584'] == 0:
               CHI_N2O2 = tol_max
            else:
               CHI_N2O2 =(np.log10(index['NII_6584']/index['OII_3727']) - N2O2_obs)**2/(abs(np.log10(index['NII_6584']/index['OII_3727'])+1e-5))
            if N2O2a_obs == -10:
               CHI_N2O2a = 0
            elif N2O2_obs > -10:
               CHI_N2O2a = 0
            elif index['OII_7325'] == 0 or index['NII_6584'] == 0:
               CHI_N2O2a = tol_max
            else:
               CHI_N2O2a =(np.log10(index['NII_6584']/index['OII_7325']) - N2O2a_obs)**2/(abs(np.log10(index['NII_6584']/index['OII_7325'])+1e-5))
            if N2S2_obs == -10: 
               CHI_N2S2 = 0
            elif N2O2_obs > -10:
               CHI_N2S2 = 0
            elif index['NII_6584'] == 0 or index['SII_671731'] == 0:
               CHI_N2S2 = tol_max
            else:
               CHI_N2S2 =(np.log10(index['NII_6584']/index['SII_671731']) - N2S2_obs)**2/(abs(np.log10(index['NII_6584']/index['SII_671731'])+1e-5))
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif N2O2_obs > -10:
               CHI_O3N2 = 0
            elif index['OIII_5007'] == 0 or index['NII_6584'] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index['OIII_5007']/index['NII_6584']) - O3N2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['NII_6584']+1e-5)))

            if ROIII_obs > 0:
               CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2 + CHI_O3N2**2 + CHI_N2O2a**2)**0.5
            elif RSIII_obs > 0:
               CHI_NO = (CHI_RSIII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
            elif RNII_obs > 0:
               CHI_NO = (CHI_RNII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
            elif ROII_obs > 0:
               CHI_NO = (CHI_ROII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
            else:
               CHI_NO = (CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5


            if CHI_NO == 0:
               NO_e = NO_e
               den_NO_e = den_NO_e  
            else:
               NO_e = (index['logNO'] - NO)**2 / (CHI_NO)**2 + NO_e
               den_NO_e = 1 / (CHI_NO)**2 + den_NO_e  


         eNO = NO_e / den_NO_e 


#Iterations for the interpolation mode

         if inter == 0 or NO == -10:
            NOf = NO
         elif inter == 1:
            igrid = grid[np.lexsort((grid['12logOH'],grid['logU']))]
            igrid = interpolate(igrid,1,NO-eNO-0.125,NO+eNO+0.125,10)

            CHI_ROIII = 0
            CHI_ROII = 0
            CHI_RNII = 0
            CHI_RSIII = 0
            CHI_N2O2 = 0
            CHI_N2O2a = 0
            CHI_N2S2 = 0
            CHI_NO = 0
            NO_p = 0
            den_NO = 0

            for index in igrid:
               if ROIII_obs == 0: 
                  CHI_ROIII = 0
               elif index['OIII_4363'] == 0:
                  CHI_ROIII = tol_max
               else:   
                  CHI_ROIII = (index['OIII_5007']/index['OIII_4363']- ROIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
               if ROII_obs == 0: 
                  CHI_ROII = 0
               elif index['OII_7325'] == 0:
                  CHI_ROII = tol_max
               else:   
                  CHI_ROII = (index['OII_3727']/index['OII_3727']- ROII_obs)**2/(index['OII_3727']/index['OII_7325'])
               if RNII_obs == 0: 
                  CHI_RNII = 0
               elif index['NII_5755'] == 0 or index['NII_6584'] == 0:
                  CHI_RNII = tol_max
               else:   
                  CHI_RNII = (index['NII_6584']/index['NII_5755']- RNII_obs)**2/(index['NII_6584']/index['NII_5755'])
               if RSIII_obs == 0: 
                  CHI_RSIII = 0
               elif index['SIII_6312'] == 0 or index['SIII_9069'] == 0:
                  CHI_RSIII = tol_max
               else:   
                  CHI_RSIII = (index['SIII_9069']/index['SIII_6312']- RSIII_obs)**2/(index['SIII_9069']/index['SIII_6312'])
               if N2O2_obs == -10:
                  CHI_N2O2 = 0
               elif index['OII_3727'] == 0 or index['NII_6584'] == 0:
                  CHI_N2O2 = tol_max
               else:
                  CHI_N2O2 =(np.log10(index['NII_6584']/index['OII_3727']) - N2O2_obs)**2/(abs(np.log10(index['NII_6584']/index['OII_3727'])+1e-5))
               if N2O2a_obs == -10:
                  CHI_N2O2a = 0
               elif N2O2_obs > -10: 
                  CHI_N2O2a = 0
               elif index['OII_7325'] == 0 or index['NII_6584'] == 0:
                  CHI_N2O2a = tol_max
               else:
                  CHI_N2O2a =(np.log10(index['NII_6584']/index['OII_7325']) - N2O2a_obs)**2/(abs(np.log10(index['NII_6584']/index['OII_7325'])+1e-5))
               if N2S2_obs == -10: 
                  CHI_N2S2 = 0
               elif N2O2_obs > -10: 
                  CHI_N2S2 = 0
               elif index['NII_6584'] == 0 or index['SII_671731'] == 0:
                  CHI_N2S2 = tol_max
               else:
                  CHI_N2S2 =(np.log10(index['NII_6584']/index['SII_671731']) - N2S2_obs)**2/(abs(np.log10(index['NII_6584']/index['SII_671731'])+1e-5))
               if OIII_5007_obs == 0 or NII_6584_obs == 0:
                  CHI_O3N2 = 0
               elif N2O2_obs > -10:
                  CHI_O3N2 = 0
               elif index['OIII_5007'] == 0 or index['NII_6584'] == 0:
                  CHI_O3N2 = tol_max
               else:
                  CHI_O3N2 = (np.log10(index['OIII_5007']/index['NII_6584']) - O3N2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['NII_6584']+1e-5)))

               if RNII_obs > 0:
                  CHI_NO = (CHI_RNII**2 + CHI_N2O2**2 + CHI_N2S2**2 + CHI_O3N2**2 + CHI_N2O2a**2)**0.5
               elif RSIII_obs > 0:
                  CHI_NO = (CHI_RSIII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
               elif ROIII_obs > 0:
                  CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
               elif ROII_obs > 0:
                  CHI_NO = (CHI_ROII**2 + CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5
               else:
                  CHI_NO = (CHI_N2O2**2 + CHI_N2S2**2+ CHI_O3N2**2+ CHI_N2O2a**2)**0.5


               if CHI_NO == 0:
                  NO_p = NO_p
                  den_NO = den_NO
               else:
                  NO_p = index['logNO'] / CHI_NO + NO_p
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
         if np.abs(index['logNO'] - NOff) > np.abs(eNOff+res_NO):
            continue
         else:
            grid_mac.append(index[[name_col for name_col in grid1.dtype.names]])
            format_array = []
            for name_col in grid1.dtype.names:
               format_array.append((name_col, grid1[0][name_col].dtype))
         
            #grid_c_1 = np.reshape(grid_mac,(int(len(grid_mac)/len(grid_aux.dtype.names[:-3])),len(grid_aux.dtype.names[:-3])))
            grid_c = np.array(grid_mac, dtype=format_array)

# Calculation of O/H, and logU


   for monte in range(0,n,1):

      OH_p = 0
      logU_p = 0
      den_OH = 0
      OH_e = 0
      logU_e = 0
      den_OH_e = 0
      tol_max = 1e2
      OII_3727_obs = 0
      if OII_3727[tab] <= 0:
         OII_3727_obs = 0
      else:
         while OII_3727_obs <= 0:
            OII_3727_obs = np.random.normal(OII_3727[tab],eOII_3727[tab]+1e-5)
      OII_7325_obs = 0
      if OII_7325[tab] <= 0:
         OII_7325_obs = 0
      else:
         while OII_7325_obs <= 0:
            OII_7325_obs = np.random.normal(OII_7325[tab],eOII_7325[tab]+1e-5)
      NeIII_3868_obs = 0
      if NeIII_3868[tab] <= 0:
         NeIII_3868_obs = 0
      else:
         while NeIII_3868_obs <= 0:
            NeIII_3868_obs = np.random.normal(NeIII_3868[tab],eNeIII_3868[tab]+1e-5)
      OIII_4363_obs = 0
      if OIII_4363[tab] <= 0:
         OIII_4363_obs = 0
      else:
         while OIII_4363_obs <= 0:
            OIII_4363_obs = np.random.normal(OIII_4363[tab],eOIII_4363[tab]+1e-5)
      OIII_5007_obs = 0
      if OIII_5007[tab] <= 0:
         OIII_5007_obs = 0
      else:
         while OIII_5007_obs <= 0:
            OIII_5007_obs = np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-5)
      if OIII_4363_obs <= 0  or OIII_5007_obs <= 0:
         ROIII_obs = 0
      else:
         ROIII_obs = OIII_5007_obs / OIII_4363_obs
      if OII_7325_obs <= 0  or OII_3727_obs <= 0:
         ROII_obs = 0
      else:
         ROII_obs = OII_3727_obs / OII_7325_obs
      NII_6584_obs = 0
      if NII_6584[tab] <= 0:
         NII_6584_obs = 0
      else:
         while NII_6584_obs <= 0:
            NII_6584_obs = np.random.normal(NII_6584[tab],eNII_6584[tab]+1e-5)
         NII_5755_obs = 0
      if NII_5755[tab] <= 0:
         NII_5755_obs = 0
      else:
         while NII_5755_obs <= 0:
            NII_5755_obs = np.random.normal(NII_5755[tab],eNII_5755[tab]+1e-5)
      if NII_5755_obs <= 0  or NII_6584_obs <= 0:
         RNII_obs = 0
      else:
         RNII_obs = NII_6584_obs / NII_5755_obs
      SII_6725_obs = 0
      if SII_6725[tab] <= 0:
            SII_6725_obs = 0
      else:
         while SII_6725_obs <= 0:
            SII_6725_obs = np.random.normal(SII_6725[tab],eSII_6725[tab]+1e-5)
      if OII_3727_obs <= 0 or OIII_5007_obs<= 0:
         O2O3_obs = 0
         R23_obs = -10
      else:
         R23_obs = np.log10(OII_3727_obs + OIII_5007_obs )
         O2O3_obs = (OII_3727_obs / OIII_5007_obs )
      if OII_3727_obs <= 0 or NeIII_3868_obs<= 0:
         O2Ne3_obs = 0
         R2Ne3_obs = -10
      else:
         O2Ne3_obs = (OII_3727_obs / NeIII_3868_obs )
         R2Ne3_obs = np.log10(OII_3727_obs + NeIII_3868_obs )
      if OIII_5007_obs <= 0 or NII_6584_obs <= 0:
         O3N2_obs = -10
      else:
         O3N2_obs = np.log10( OIII_5007_obs / NII_6584_obs )
      if OIII_5007_obs <= 0 or SII_6725_obs <= 0:
         O3S2_obs = -10
      else:
         O3S2_obs = np.log10( OIII_5007_obs / SII_6725_obs )
      SIII_6312_obs = 0
      if SIII_6312[tab] <= 0:
            SIII_6312_obs = 0
      else:
         while SIII_6312_obs <= 0:
            SIII_6312_obs = np.random.normal(SIII_6312[tab],eSIII_6312[tab]+1e-5)
      SIII_9069_obs = 0
      if SIII_9069[tab] <= 0:
            SIII_9069_obs = 0
      else:
         while SIII_9069_obs <= 0:
            SIII_9069_obs = np.random.normal(SIII_9069[tab],eSIII_9069[tab]+1e-5)
      if SIII_6312_obs <= 0  or SIII_9069_obs <= 0:
         RSIII_obs = 0
      else:
         RSIII_obs = SIII_9069_obs / SIII_6312_obs
      if SII_6725_obs <= 0 or SIII_9069_obs<= 0:
         S2S3_obs = 0
         S23_obs = 0
      else:
         S23_obs = (SII_6725_obs + 3.44*SIII_9069_obs )
         S2S3_obs = (SII_6725_obs / SIII_9069_obs )
      if NII_6584_obs <= 0 or SIII_9069_obs<= 0:
         N2S3_obs = 0
      else:
         N2S3_obs = (NII_6584_obs / SIII_9069_obs )
      if SII_6725_obs <= 0 or SIII_6312_obs<= 0:
         S2S3a_obs = 0
         S23a_obs = 0
      else:
         S23a_obs = (SII_6725_obs + SIII_6312_obs )
         S2S3a_obs = (SII_6725_obs / SIII_6312_obs )


      if R23_obs == -10 and NII_6584_obs == 0 and ROIII_obs == 0 and R2Ne3_obs == -10 and O3S2_obs == -10 and ROII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and S2S3_obs == 0 and S2S3a_obs == 0:
         OH = 0
         logU = 0
      else:
         CHI_ROIII = 0
         CHI_RNII = 0
         CHI_ROII = 0
         CHI_ROIII = 0
         CHI_NII = 0
         CHI_SIII = 0
         CHI_SII = 0
         CHI_OII = 0
         CHI_O2O3 = 0
         CHI_R23 = 0
         CHI_O2Ne3 = 0
         CHI_R2Ne3 = 0
         CHI_O3N2 = 0
         CHI_O3S2 = 0
         CHI_S23 = 0
         CHI_S2S3 = 0
         CHI_N2S3 = 0
         CHI_S23a = 0
         CHI_S2S3a = 0
         CHI_OH = 0
         for index in grid_c:
            if index['OIII_5007'] == 0 and OIII_5007_obs > 0: continue
            if index['OII_3727'] == 0 and OII_3727_obs > 0: continue
            if index['NII_6584'] == 0 and NII_6584_obs > 0: continue
            if index['SIII_9069'] == 0 and SIII_9069_obs > 0: continue
            if index['OIII_4363'] == 0 and OIII_4363_obs > 0: continue
            if index['SIII_6312'] == 0 and SIII_6312_obs > 0: continue
            if ROIII_obs == 0: 
               CHI_ROIII = 0
            elif index['OIII_4363'] == 0:
               CHI_ROIII = tol_max
            else:   
               CHI_ROIII = (index['OIII_5007']/index['OIII_4363']- ROIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
            if ROII_obs == 0: 
               CHI_ROII = 0
            elif index['OII_7325'] == 0:
               CHI_ROII = tol_max
            else:   
               CHI_ROII = (index['OII_3727']/index['OII_3727']- ROII_obs)**2/(index['OII_3727']/index['OII_3727'])
            if RNII_obs == 0: 
               CHI_RNII = 0
            elif index['NII_5755'] == 0 or index['NII_6584'] == 0:
               CHI_RNII = tol_max
            else:   
               CHI_RNII = (index['NII_6584']/index['NII_5755']- RNII_obs)**2/(index['NII_6584']/index['NII_5755'])
            if RSIII_obs == 0: 
               CHI_RSIII = 0
            elif index['SIII_6312'] == 0 or index['SIII_9069'] == 0:
               CHI_RSIII = tol_max
            else:   
               CHI_RSIII = (index['SIII_9069']/index['SIII_6312']- RSIII_obs)**2/(index['SIII_9069']/index['SIII_6312'])
            if OIII_5007_obs == 0:
               CHI_OIII = 0
            elif index['OIII_5007'] == 0:
               CHI_OIII = tol_max
            else:
               CHI_OIII = (index['OIII_5007'] - OIII_5007_obs)**2/index['OIII_5007']
            if SIII_9069_obs == 0:
               CHI_SIII = 0
            elif index['SIII_9069'] == 0:
               CHI_SIII = tol_max
            else:
               CHI_SIII = (index['SIII_9069'] - SIII_9069_obs)**2/index['SIII_9069']
            if SII_6725_obs == 0:
               CHI_SII = 0
            elif index['SII_671731'] == 0:
               CHI_SII = tol_max
            else:
               CHI_SII = (index['SII_671731'] - SII_6725_obs)**2/index['SII_671731']
            if OII_3727_obs == 0:
               CHI_OII = 0
            elif index['OII_3727'] == 0:
               CHI_OII = tol_max
            else:
               CHI_OII = (index['OII_3727'] - OII_3727_obs)**2/index['OII_3727']
            if NII_6584_obs == 0:
               CHI_NII = 0
            elif index['NII_6584'] == 0:
               CHI_NII = tol_max
            else:
               CHI_NII = (index['NII_6584'] - NII_6584_obs)**2/index['NII_6584']
            if OII_3727_obs == 0 or OIII_5007_obs == 0:
               CHI_O2O3 = 0
               CHI_R23 = 0
            elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
               CHI_O2O3 = tol_max
               CHI_R23 = tol_max
            else:
               CHI_O2O3 = (index['OII_3727']/index['OIII_5007'] - O2O3_obs)**2/(index['OII_3727']/index['OIII_5007'])
               CHI_R23 = (np.log10(index['OII_3727']+index['OIII_5007'])-R23_obs)**2/   (np.abs(np.log10(index['OII_3727']+index['OIII_5007']+1e-5)))
            if OII_3727_obs == 0 or NeIII_3868_obs == 0:
               CHI_O2Ne3 = 0
               CHI_R2Ne3 = 0
            elif index['OII_3727'] == 0 or index['NeIII_3868'] == 0:
               CHI_O2Ne3 = tol_max
               CHI_R2Ne3 = tol_max
            else:
               CHI_O2Ne3 = (index['OII_3727']/index['NeIII_3868'] - O2Ne3_obs)**2/(index['OII_3727']/index['NeIII_3868'])
               CHI_R2Ne3 = (np.log10(index['OII_3727']+index['NeIII_3868'])-R2Ne3_obs)**2/   (np.abs(np.log10(index['OII_3727']+index['NeIII_3868']+1e-5)))
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif index['OIII_5007'] == 0 or index['NII_6584'] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index['OIII_5007']/index['NII_6584']) - O3N2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['NII_6584']+1e-5)))
            if OIII_5007_obs == 0 or SII_6725_obs == 0:
               CHI_O3S2 = 0
            elif index['OIII_5007'] == 0 or index['SII_671731'] == 0:
               CHI_O3S2 = tol_max
            else:
               CHI_O3S2 = (np.log10(index['OIII_5007']/index['SII_671731']) - O3S2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['SII_671731']+1e-5)))
            if SII_6725_obs == 0 or SIII_9069_obs == 0:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max
            else:
               CHI_S2S3 = (index['SII_671731']/index['SIII_9069'] - S2S3_obs)**2/(index['SII_671731']/index['SIII_9069'])
               CHI_S23 = ((index['SII_671731']+3.44*index['SIII_9069'])-S23_obs)**2/(np.abs((index['SII_671731']+3.44*index['SIII_9069'])))
            if NII_6584_obs == 0 or SIII_9069_obs == 0:
               CHI_N2S3 = 0
            elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
               CHI_N2S3 = tol_max
            else:
               CHI_N2S3 = (index['NII_6584']/index['SIII_9069'] - N2S3_obs)**2/(index['NII_6584']/index['SIII_9069'])
            if SII_6725_obs == 0 or SIII_6312_obs == 0:
               CHI_S2S3a = 0
               CHI_S23a = 0
            elif SII_6725_obs > 0 and SIII_9069_obs > 0:
               CHI_S2S3a = 0
               CHI_S23a = 0
            elif index['SII_671731'] == 0 or index['SIII_6312'] == 0:
               CHI_S2S3a = tol_max
               CHI_S23a = tol_max
            else:
               CHI_S2S3a = (index['SII_671731']/index['SIII_6312'] - S2S3a_obs)**2/(index['SII_671731']/index['SIII_6312'])
               CHI_S23a = ((index['SII_671731']+index['SIII_6312'])-S23a_obs)**2/   (np.abs((index['SII_671731']+index['SIII_6312']+1e-5)))

            if RSIII_obs > 0 :
               CHI_OH = (CHI_RSIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROIII_obs > 0 :
               CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif RNII_obs > 0 :
               CHI_OH = (CHI_RNII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROII_obs > 0 :
               CHI_OH = (CHI_ROII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROIII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and ROII_obs == 0 and NII_6584_obs > 0:
               CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O3N2**2 + CHI_O3S2**2 + CHI_S2S3**2 + CHI_S2S3a**2 + CHI_O3S2**2)**0.5
            elif ROIII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and NII_6584_obs == 0 and OIII_5007_obs > 0:
               CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O3S2**2 + CHI_S2S3**2 + CHI_S2S3a**2)**0.5
            elif ROIII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and OIII_5007_obs == 0:
               CHI_OH = (CHI_O2Ne3**2 + CHI_R2Ne3**2 + CHI_S2S3**2 + CHI_S2S3a**2)**0.5

            if CHI_OH == 0:
               OH_p = OH_p
               logU_p = logU_p
               den_OH = den_OH
            else:
               OH_p = index['12logOH'] / (CHI_OH)**2 + OH_p
               logU_p = index['logU'] / (CHI_OH)**2 + logU_p
               den_OH = 1 / (CHI_OH)**2 + den_OH


         OH = OH_p / den_OH
         logU = logU_p / den_OH

#Calculation of error of O/H  and logU

      if R23_obs == -10 and NII_6584_obs == 0 and ROIII_obs == 0 and R2Ne3_obs == -10 and O3S2_obs == -10 and RNII_obs == 0 and RSIII_obs == 0 and S2S3_obs == 0 and S2S3a_obs == 0:
         eOH = 0
         elogU = 0
      else:
         CHI_ROIII = 0
         CHI_RNII = 0
         CHI_RSIII = 0
         CHI_NII = 0
         CHI_OIII = 0
         CHI_OII = 0
         CHI_SIII = 0
         CHI_SII = 0
         CHI_O2O3 = 0
         CHI_R23 = 0
         CHI_O2Ne3 = 0
         CHI_R2Ne3 = 0
         CHI_O3N2 = 0
         CHI_O3S2 = 0
         CHI_S23 = 0
         CHI_S2S3 = 0
         CHI_N2S3 = 0
         CHI_S23a = 0
         CHI_S2S3a = 0
         CHI_OH = 0

         for index in grid_c:
            if index['OIII_5007'] == 0 and OIII_5007_obs > 0: continue
            if index['OII_3727'] == 0 and OII_3727_obs > 0: continue
            if index['NII_6584'] == 0 and NII_6584_obs > 0: continue
            if index['SIII_9069'] == 0 and SIII_9069_obs > 0: continue
            if index['OIII_4363'] == 0 and OIII_4363_obs > 0: continue
            if index['SIII_6312'] == 0 and SIII_6312_obs > 0: continue
            if ROIII_obs == 0: 
               CHI_ROIII = 0
            elif index['OIII_4363'] == 0:
               CHI_ROIII = tol_max
            else:   
               CHI_ROIII = (index['OIII_5007']/index['OIII_4363']- ROIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
            if ROII_obs == 0: 
               CHI_ROII = 0
            elif index['OII_7325'] == 0:
               CHI_ROII = tol_max
            else:   
               CHI_ROII = (index['OII_3727']/index['OII_7325']- ROII_obs)**2/(index['OII_3727']/index['OII_7325'])
            if RNII_obs == 0: 
               CHI_RNII = 0
            elif index['NII_5755'] == 0 or index['NII_6584'] == 0:
               CHI_RNII = tol_max
            else:   
               CHI_RNII = (index['NII_6584']/index['NII_5755']- RNII_obs)**2/(index['NII_6584']/index['NII_5755'])
            if RSIII_obs == 0: 
               CHI_RSIII = 0
            elif index['SIII_6312'] == 0 or index['SIII_9069'] == 0:
               CHI_RSIII = tol_max
            else:   
               CHI_RSIII = (index['SIII_9069']/index['SIII_6312']- RSIII_obs)**2/(index['SIII_9069']/index['SIII_6312'])
            if SIII_9069_obs == 0:
               CHI_SIII = 0
            elif index['SIII_9069'] == 0:
               CHI_SIII = tol_max
            else:
               CHI_SIII = (index['SIII_9069'] - SIII_9069_obs)**2/index['SIII_9069']
            if SII_6725_obs == 0:
               CHI_SII = 0
            elif index['SII_671731'] == 0:
               CHI_SII = tol_max
            else:
               CHI_SII = (index['SII_671731'] - SII_6725_obs)**2/index['SII_671731']
            if OIII_5007_obs == 0:
               CHI_OIII = 0
            elif index['OIII_5007'] == 0:
               CHI_OIII = tol_max
            else:
               CHI_OIII = (index['OIII_5007'] - OIII_5007_obs)**2/index['OIII_5007']
            if OII_3727_obs == 0:
               CHI_OII = 0
            elif index['OII_3727'] == 0:
               CHI_OII = tol_max
            else:
               CHI_OII = (index['OII_3727'] - OII_3727_obs)**2/index['OII_3727']
            if NII_6584_obs == 0:
               CHI_NII = 0
            elif index['NII_6584'] == 0:
               CHI_NII = tol_max
            else:
               CHI_NII = (index['NII_6584'] - NII_6584_obs)**2/index['NII_6584']
            if OII_3727_obs == 0 or OIII_5007_obs == 0:
               CHI_O2O3 = 0
               CHI_R23 = 0
            elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
               CHI_O2O3 = tol_max
               CHI_R23 = tol_max
            else:
               CHI_O2O3 = (index['OII_3727']/index['OIII_5007'] - O2O3_obs)**2/(index['OII_3727']/index['OIII_5007'])
               CHI_R23 = (np.log10(index['OII_3727']+index['OIII_5007'])-R23_obs)**2/   (np.abs(np.log10(index['OII_3727']+index['OIII_5007']+1e-5)))
            if OII_3727_obs == 0 or NeIII_3868_obs == 0:
               CHI_O2Ne3 = 0
               CHI_R2Ne3 = 0
            elif index['OII_3727'] == 0 or index['NeIII_3868'] == 0:
               CHI_O2Ne3 = tol_max
               CHI_R2Ne3 = tol_max
            else:
               CHI_O2Ne3 = (index['OII_3727']/index['NeIII_3868'] - O2Ne3_obs)**2/(index['OII_3727']/index['NeIII_3868'])
               CHI_R2Ne3 = (np.log10(index['OII_3727']+index['NeIII_3868'])-R2Ne3_obs)**2/   (np.abs(np.log10(index['OII_3727']+index['NeIII_3868']+1e-5)))
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif index['OIII_5007'] == 0 or index['NII_6584'] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index['OIII_5007']/index['NII_6584']) - O3N2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['NII_6584']+1e-5)))
            if OIII_5007_obs == 0 or SII_6725_obs == 0:
               CHI_O3S2 = 0
            elif OIII_5007_obs > 0 and NII_6584_obs > 0:
               CHI_O3S2 = 0
            elif index['OIII_5007'] == 0 or index['SII_671731'] == 0:
               CHI_O3S2 = tol_max
            else:
               CHI_O3S2 = (np.log10(index['OIII_5007']/index['SII_671731']) - O3S2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['SII_671731']+1e-5)))
            if SII_6725_obs == 0 or SIII_9069_obs == 0:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max
            else:
               CHI_S2S3 = (index['SII_671731']/index['SIII_9069'] - S2S3_obs)**2/(index['SII_671731']/index['SIII_9069'])
               CHI_S23 = ((index['SII_671731']+3.44*index['SIII_9069'])-S23_obs)**2/   (np.abs((index['SII_671731']+3.44*index['SIII_9069']+1e-5)))
            if NII_6584_obs == 0 or SIII_9069_obs == 0:
               CHI_N2S3 = 0
            elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
               CHI_N2S3 = tol_max
            else:
               CHI_N2S3 = (index['NII_6584']/index['SIII_9069'] - N2S3_obs)**2/(index['NII_6584']/index['SIII_9069'])
            if SII_6725_obs == 0 or SIII_6312_obs == 0:
               CHI_S2S3a = 0
               CHI_S23a = 0
            elif SII_6725_obs > 0 and SIII_9069_obs > 0:
               CHI_S2S3a = 0
               CHI_S23a = 0
            elif index['SII_671731'] == 0 or index['SIII_6312'] == 0:
               CHI_S2S3a = tol_max
               CHI_S23a = tol_max
            else:
               CHI_S2S3a = (index['SII_671731']/index['SIII_6312'] - S2S3a_obs)**2/(index['SII_671731']/index['SIII_6312'])
               CHI_S23a = ((index['SII_671731']+index['SIII_6312'])-S23a_obs)**2/   (np.abs((index['SII_671731']+index['SIII_6312']+1e-5)))


            if RSIII_obs > 0 :
               CHI_OH = (CHI_RSIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROIII_obs > 0 :
               CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif RNII_obs > 0 :
               CHI_OH = (CHI_RNII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROII_obs > 0 :
               CHI_OH = (CHI_ROII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROIII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and NII_6584_obs > 0:
               CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O3N2**2 + CHI_O3S2**2 + CHI_S2S3**2 + CHI_S2S3a**2 + CHI_O3S2**2)**0.5
            elif ROIII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and NII_6584_obs == 0 and OIII_5007_obs > 0:
               CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O3S2**2 + CHI_S2S3**2 + CHI_S2S3a**2)**0.5
            elif ROIII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and OIII_5007_obs == 0:
               CHI_OH = (CHI_O2Ne3**2 + CHI_R2Ne3**2 + CHI_S2S3**2 + CHI_S2S3a**2)**0.5


            if CHI_OH == 0:
               OH_e = OH_e
               logU_e = logU_e
               den_OH_e = den_OH_e
            else:
               OH_e = (index['12logOH'] - OH)**2 / (CHI_OH)**2 + OH_e
               logU_e = (index['logU'] - logU)**2 / (CHI_OH)**2 + logU_e
               den_OH_e = 1 / (CHI_OH)**2 + den_OH_e 



         eOH = OH_e / den_OH_e
         elogU = logU_e / den_OH_e 
# Iterations for interpolated models

      if inter == 0 or OH == 0:
         OHf = OH
         logUf = logU
      elif inter == 1:
         igrid = interpolate(grid_c,2,logU-elogU-0.25,logU+elogU+0.25,10)
         igrid = igrid[np.lexsort((igrid['logNO'],igrid['logU']))]
         igrid = interpolate(igrid,0,OH-eOH-0.1,OH+eOH+0.1,10)
         igrid = igrid[np.lexsort((igrid['12logOH'],igrid['logU']))]

         CHI_ROIII = 0
         CHI_ROII = 0
         CHI_RNII = 0
         CHI_RSIII = 0
         CHI_SII = 0
         CHI_SIII = 0
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
         CHI_S23 = 0
         CHI_S2S3 = 0
         CHI_N2S3 = 0
         CHI_S23a = 0
         CHI_S2S3a = 0
         OH_p = 0
         logU_p = 0
         den_OH = 0
      

         for index in igrid:
            if index['OIII_5007'] == 0 and OIII_5007_obs > 0: continue
            if index['OII_3727'] == 0 and OII_3727_obs > 0: continue
            if index['NII_6584'] == 0 and NII_6584_obs > 0: continue
            if index['SIII_9069'] == 0 and SIII_9069_obs > 0: continue
            if index['OIII_4363'] == 0 and OIII_4363_obs > 0: continue
            if index['SIII_6312'] == 0 and SIII_6312_obs > 0: continue
            if ROIII_obs == 0: 
               CHI_ROIII = 0
            elif index['OIII_4363'] == 0:
               CHI_ROIII = tol_max
            else:   
               CHI_ROIII = (index['OIII_5007']/index['OIII_4363']- ROIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
            if ROII_obs == 0: 
               CHI_ROII = 0
            elif index['OII_7325'] == 0:
               CHI_ROII = tol_max
            else:   
               CHI_ROII = (index['OII_3727']/index['OII_7325']- ROII_obs)**2/(index['OII_3727']/index['OII_7325'])
            if RNII_obs == 0: 
               CHI_RNII = 0
            elif index['NII_5755'] == 0 or index['NII_6584'] == 0:
               CHI_RNII = tol_max
            else:   
               CHI_RNII = (index['NII_6584']/index['NII_5755']- RNII_obs)**2/(index['NII_6584']/index['NII_5755'])
            if RSIII_obs == 0: 
               CHI_RSIII = 0
            elif index['SIII_6312'] == 0 or index['SIII_9069'] == 0:
               CHI_RSIII = tol_max
            else:   
               CHI_RSIII = (index['SIII_9069']/index['SIII_6312']- RSIII_obs)**2/(index['SIII_9069']/index['SIII_6312'])
            if SIII_9069_obs == 0:
               CHI_SIII = 0
            elif index['SIII_9069'] == 0:
               CHI_SIII = tol_max
            else:
               CHI_SIII = (index['SIII_9069'] - SIII_9069_obs)**2/index['SIII_9069']
            if SII_6725_obs == 0:
               CHI_SII = 0
            elif index['SII_671731'] == 0:
               CHI_SII = tol_max
            else:
               CHI_SII = (index['SII_671731'] - SII_6725_obs)**2/index['SII_671731']
            if OIII_5007_obs == 0:
               CHI_OIII = 0
            elif index['OIII_5007'] == 0:
               CHI_OIII = tol_max
            else:
               CHI_OIII = (index['OIII_5007'] - OIII_5007_obs)**2/index['OIII_5007']
            if OII_3727_obs == 0:
               CHI_OII = 0
            elif index['OII_3727'] == 0:
               CHI_OII = tol_max
            else:
               CHI_OII = (index['OII_3727'] - OII_3727_obs)**2/index['OII_3727']
            if NII_6584_obs == 0:
               CHI_NII = 0
            elif index['NII_6584'] == 0:
               CHI_NII = tol_max
            else:
               CHI_NII = (index['NII_6584'] - NII_6584_obs)**2/index['NII_6584']
            if OII_3727_obs == 0 or OIII_5007_obs == 0:
               CHI_O2O3 = 0
               CHI_R23 = 0
            elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
               CHI_O2O3 = tol_max
               CHI_R23 = tol_max
            else:
               CHI_O2O3 = (index['OII_3727']/index['OIII_5007'] - O2O3_obs)**2/(index['OII_3727']/index['OIII_5007'])
               CHI_R23 = (np.log10(index['OII_3727']+index['OIII_5007'])-R23_obs)**2/(np.abs(np.log10(index['OII_3727']+index['OIII_5007']+1e-5)))
            if OII_3727_obs == 0 or NeIII_3868_obs == 0:
               CHI_O2Ne3 = 0
               CHI_R2Ne3 = 0
            elif index['OII_3727'] == 0 or index['NeIII_3868'] == 0:
               CHI_O2Ne3 = tol_max
               CHI_R2Ne3 = tol_max
            else:
               CHI_O2Ne3 = (index['OII_3727']/index['NeIII_3868'] - O2Ne3_obs)**2/(index['OII_3727']/index['NeIII_3868'])
               CHI_R2Ne3 = (np.log10(index['OII_3727']+index['NeIII_3868'])-R2Ne3_obs)**2/(np.abs(np.log10(index['OII_3727']+index['NeIII_3868']+1e-5)))
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
               CHI_O3N2 = 0
            elif index['OIII_5007'] == 0 or index['NII_6584'] == 0:
               CHI_O3N2 = tol_max
            else:
               CHI_O3N2 = (np.log10(index['OIII_5007']/index['NII_6584']) - O3N2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['NII_6584']+1e-5)))
            if OIII_5007_obs == 0 or SII_6725_obs == 0:
               CHI_O3S2 = 0
            elif OIII_5007_obs > 0 and NII_6584_obs > 0:
               CHI_O3S2 = 0
            elif index['OIII_5007'] == 0 or index['SII_671731'] == 0:
               CHI_O3S2 = tol_max
            else:
               CHI_O3S2 = (np.log10(index['OIII_5007']/index['SII_671731']) - O3S2_obs)**2/(np.abs(np.log10(index['OIII_5007']/index['SII_671731']+1e-5)))
            if SII_6725_obs == 0 or SIII_9069_obs == 0:
               CHI_S2S3 = 0
               CHI_S23 = 0
            elif index['SII_671731'] == 0 or index['SIII_9069'] == 0:
               CHI_S2S3 = tol_max
               CHI_S23 = tol_max
            else:
               CHI_S2S3 = (index['SII_671731']/index['SIII_9069'] - S2S3_obs)**2/(index['SII_671731']/index['SIII_9069'])
               CHI_S23 = ((index['SII_671731']+3.44*index['SIII_9069'])-S23_obs)**2/   (np.abs((index['SII_671731']+3.44*index['SIII_9069']+1e-5)))
            if NII_6584_obs == 0 or SIII_9069_obs == 0:
               CHI_N2S3 = 0
            elif index['NII_6584'] == 0 or index['SIII_9069'] == 0:
               CHI_N2S3 = tol_max
            else:
               CHI_N2S3 = (index['NII_6584']/index['SIII_9069'] - N2S3_obs)**2/(index['NII_6584']/index['SIII_9069'])
            if SII_6725_obs == 0 or SIII_6312_obs == 0:
               CHI_S2S3a = 0
               CHI_S23a = 0
            elif SII_6725_obs > 0 and SIII_9069_obs > 0:
               CHI_S2S3a = 0
               CHI_S23a = 0
            elif index['SII_671731'] == 0 or index['SIII_6312'] == 0:
               CHI_S2S3a = tol_max
               CHI_S23a = tol_max
            else:
               CHI_S2S3a = (index['SII_671731']/index['SIII_6312'] - S2S3a_obs)**2/(index['SII_671731']/index['SIII_6312'])
               CHI_S23a = ((index['SII_671731']+index['SIII_6312'])-S23a_obs)**2/   (np.abs((index['SII_671731']+index['SIII_6312']+1e-5)))


            if RSIII_obs > 0 :
               CHI_OH = (CHI_RSIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROIII_obs > 0 :
               CHI_OH = (CHI_ROIII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif RNII_obs > 0 :
               CHI_OH = (CHI_RNII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROII_obs > 0 :
               CHI_OH = (CHI_ROII**2 + CHI_NII**2 + CHI_OII**2 + CHI_OIII**2 )**0.5
            elif ROIII_obs == 0 and ROII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and NII_6584_obs > 0:
               CHI_OH = (CHI_NII**2 + CHI_O2O3**2 + CHI_R23**2 + CHI_O3N2**2 + CHI_O3S2**2 + CHI_S2S3**2 + CHI_S2S3a**2 + CHI_O3S2**2)**0.5
            elif ROIII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and NII_6584_obs == 0 and OIII_5007_obs > 0:
               CHI_OH = (CHI_O2O3**2 + CHI_R23**2 + CHI_O3S2**2 + CHI_S2S3**2 + CHI_S2S3a**2)**0.5
            elif ROIII_obs == 0 and RNII_obs == 0 and RSIII_obs == 0 and OIII_5007_obs == 0:
               CHI_OH = (CHI_O2Ne3**2 + CHI_R2Ne3**2 + CHI_S2S3**2 + CHI_S2S3a**2)**0.5



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
   print (round(100*(count)/float(len(input1)),1),'%',Names[tab],grid_type,'', round(OHff,2), round(eOHff,2),'','',round(NOff,2), round(eNOff,2), '',round(logUff,2), round(elogUff,2))


output['grid'] = grids
output['OH'] = OHffs
output['eOH'] = eOHffs
output['NO'] = NOffs
output['eNO'] = eNOffs
output['logU'] = logUffs
output['elogU'] = elogUffs

if input0.size == 1:  output = np.delete(output,obj=1,axis=0)


lineas_header = ['HII-CHI-mistry v.5.3 output file', 'Input file: '+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'Library file used : '+file_lib_2, 'Template used to constraint grid of models: '+const_file,'']

line_label = '{:10}  '.format(output.dtype.names[0])
for ind2 in range(1, len(output.dtype.names)-7):
   line_label += '{:10}  '.format(output.dtype.names[ind2])


line_label += '{:10}   {:10}   {:10}   {:10}   {:10}   {:10}   {:10}'.format('i', 'O/H',    'eO/H',   'N/O' ,   'eN/O' , 'logU',   'elogU')
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


