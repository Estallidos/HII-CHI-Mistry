# Filename: HII-CHCm-IR_v3.2.py

#####################
###### IMPORTS ######
#####################

import string
import numpy as np
import sys
#sys.stderr = open('errorlog.txt', 'w')
import warnings
warnings.filterwarnings("ignore")

#######################
###### FUNCTIONS ######
#######################

#Function for interpolation of grids
def interpolate(grid,z,zmin,zmax,n):
   #Columns of the library
   n_comments = 0
   with open('Libraries_ir/C17_POPSTAR_1myr.dat', 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   auxiliar_labels = np.genfromtxt('Libraries_ir/C17_POPSTAR_1myr.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
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
print ('-------------------------------------------------')
print ('This is HII-CHI-mistry IR v. 3.2')
print (' See Fernandez-Ontiveros et al. (2021) and Perez-Diaz et al. (2022) for details')
print (' Insert the name of your input text file with all or some of the following columns:')
print ('')
#print ('Hb 4861AA')
print ('HI 4.05m')
print ('[ArII] 6.98m')
print ('HI 7.46m')
print ('[ArV] 7.90m')
print ('[ArIII] 8.99m')
print ('[SIV] 10.5m')
print ('HI 12.4m')
print ('[NeII] 12.8m')
print ('[ArV] 13.1m')
print ('[NeV] 14.3m')
print ('[NeIII] 15.5m')
print ('[SIII] 18.7m')
print ('[NeV] 24.2m')
print ('[OIV] 25.9m')
print ('[SIII] 33.7m')
print ('[OIII] 52m')
print ('[NIII] 57m')
print ('[OIII] 88m') 
print ('[NII] 122m')
print ('[NII] 205m')
print ('')
print ('with their corresponding labels and errors in adjacent columns')
print ('-------------------------------------------------')
print (' ')

# Input file reading
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
   with open(input00, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1		
   input0 = np.genfromtxt(input00,dtype=None,names=True, encoding = 'ascii', skip_header = n_comments)
   print ('The input file is: '+input00)
except:
   print ('Input file error: It does not exist or has wrong format')

   sys.exit()

print ('')
print(input00)
if input0.size == 1:
   input1 = np.stack((input0,input0))
else:
   input1 = input0

# Iterations for Montecarlo error derivation
if len(sys.argv) < 3:
   n = 25
else:
   n = int(sys.argv[2])
print ('The number of iterations for MonteCarlo simulation is: '+str(n))
print ('')

#############################################
###### NON INTERACTIVE OR INTERACTIVE CODE ######
#############################################

interactive = True #Change this value to False to run the code in non-interactive mode

#Questions (inputs fromt terminal)
question1 = interactive #Question to select the grids of models
question2 = interactive #Question for the value of alphaOX (only required if "sed" is set to 3)
question3 = interactive #Question for the value of efrac (only required if "sed" is set to 3)
question4 = interactive #Question required to introduce a particular file by the user (only required if "sed" is set 4)
question6 = interactive #Question to use or not interpolation for the grids
question7 = interactive #Question regarding the constrain law to be assumed by the code
question8 = interactive #Question required to tntroduced a particular file by the user (only if "const" is set to 6)

#Set values to operate
if question1 == False:
   sed = 1 #Choose grid of model: (1) for POPSTAR models; (2) for BPASS models; (3) for AGN models; (4) models introduced by the user (will require file). Replace value for the given option
if question2 == False and sed == 3:
   alpha = 1 #Choose value of alpha_OX: (1) for alpha_OX = -0.8; (2) for alpha_OX = -1.2. Replace value with any of these options
if question3 == False and sed == 3:
   efrac = 1 #Choose value for the stopping criteria, i.e., the fraction of free electrons: (1) for efrac = 0.02; (2) for efrac = 0.98
if question4 == False and sed == 4:
   new_library = 'Name_of_the_file' #Introduced name of the file with the grids of models. It must be located under the folder "Libraries_ir"
if question6 == False:
   inter = 0 #Choose value to perform interpolation: (0) no interpolation; (1) interpolation.  Replace value for the given option
if question7 == False:
   const = 1 #Choose value for the constraint laws between O/H, N/O and U that must be assumed when the code does not have enough information: (1) constraints obtained for Star-Forming Galaxies; (2) constraints obtained for Extreme Emission Line Galaxies; (3) constraints between N/O and O/H obtained for AGNs, without restriction in the ionization parameter; (4) constraints between N/O and O/H obtained for AGNs, and log(U) > -2.5; (5) constraints between N/O and O/H obtained for AGNs, and log(U) < -2.5; (6) constraint law introduced by the user (will required a file).  Replace value for the given option
if question8 == False and const == 6:
   new_const = 'Name_of_the_file' #Introduced name of the file with the constraint laws. It must be located under the folder "Constraint"

#############################################
###### SELECTION OF THE GRID OF MODELS ######
#############################################

# Interface with the user
print ('')
while question1:
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
   if sed == '1' or sed == '2' or sed == '3' or sed == '4': question1 = False 
print ('')
sed = int(sed)

# Further questions on the AGN models
if sed == 3:
   #SLOPE ALPHA
   while question2:
      if int(sys.version[0]) < 3:
         alpha = raw_input('Choose value for alpha(OX) in the AGN models: [1] -0.8 [2] -1.2: ')
      else:
         alpha = input('Choose value for alpha(OX) in the AGN models: [1] -0.8 [2] -1.2:  ')
      if alpha == '1' or alpha == '2': question2 = False
   print ('')
   #Fraction of free electrons (stopping criteria in the models)
   while question3:
      if int(sys.version[0]) < 3:
         efrac = raw_input('Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons: ')
      else:
         efrac = input('Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons:  ')
      if efrac == '1' or efrac == '2': question3 = False
   print ('')
   alpha = int(alpha)
   efrac = int(efrac)

#Particular file introduced by the user
if sed == 4:
   while question4:
      print ('Introduce name of the file containing the models. It must be located in the folder "Libraries_ir".')
      print (' ')
      if int(sys.version[0]) < 3:
         new_library = raw_input('Name of file: ')
      else:
         new_library = input('Name of file: ')
 
      #Searching for the file
      try:
         #Counting comments:
         n_comments = 0
         with open('Libraries_ir/'+new_library, 'r') as file:
            for line in file:
               if line[0] == '#':
                  n_comments += 1
         library_user = np.genfromtxt('Libraries_ir/'+new_library, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         print (' ')
         print ('Loading library '+new_library+'. Checking correct format of the file.')
         question4 = False
      except:
         print (' ')
         print ('Library was not found in folder "Libraries_ir" or file does not exist.')
   question5 = True      
   while question5:
      try:
         #Counting comments:
         n_comments = 0
         with open('Libraries_ir/'+new_library, 'r') as file:
            for line in file:
               if line[0] == '#':
                  n_comments += 1
         library_user = np.genfromtxt('Libraries_ir/'+new_library, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         #Checking correct format:
         #Counting comments:
         n_comments = 0
         with open('Libraries_ir/C17_POPSTAR_1myr.dat', 'r') as file:
            for line in file:
               if line[0] == '#':
                  n_comments += 1
         auxiliar_labels = np.genfromtxt('Libraries_ir/C17_POPSTAR_1myr.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
         missing_labels = []
         for label in auxiliar_labels:
            if label in library_user.dtype.names:
               continue
            else:
               missing_labels.append(label)
         #Displaying message for the user:
         print ('Succesfully reading of the file')
         if len(missing_labels) == 0:
            print ('File presents the correct format')
            question5 = False
         else:
            print ('File does not present the correct format. The following columns are missing:')
            for need_label in missing_labels:
               print('- '+need_label)
            print ('More details on the correct format for the library are found in readme file.')
            if interactive == True:
               print (' ')
               print ('Reintroduce name of the file with fixed format.')
               print (' ')
               if int(sys.version[0]) < 3:
                  new_library = raw_input('Name of file: ')
               else:
                  new_library = input('Name of file: ')
            elif interactive == False:
               print (' ')
               print ('Re-run the program. Please correct the format of the library before running the code again')
               sys.exit()
      except:
         if interactive == True:
            print ('Something went wrong while reading file. Please, reintroduce name of the file.')
            print (' ')
            if int(sys.version[0]) < 3:
               new_library = raw_input('Name of file: ')
            else:
               new_library = input('Name of file: ')
         elif interactive == False:
            print ('Something went wrong while reading file. Please, check again the file introduced, before running the code again.')
            sys.exit()

#Interpolation in the grid of models
print ('')
while question6:
   if int(sys.version[0]) < 3:
      inter = raw_input('Choose models [0] No interpolated [1] Interpolated: ')
   else:
      inter = input('Choose models [0] No interpolated [1] Interpolated: ')
   if inter == '0' or inter == '1': question6 = False
print ('')
inter = int(inter)

#POPSTAR MODEL
if sed==1:
   file_lib = 'C17_POPSTAR_1myr.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_ir/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_ir/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
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
   with open('Libraries_ir/'+file_lib, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_ir/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)  
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
   with open('Libraries_ir/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_ir/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
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
   with open('Libraries_ir/'+file_lib, 'r') as file9:
      for line in file9:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_ir/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
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
   with open('Libraries_ir/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_ir/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
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
   with open('Libraries_ir/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_ir/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
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
   with open('Libraries_ir/'+new_library, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1  
   grid_aux = np.genfromtxt('Libraries_ir/'+new_library,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'User file ' + new_library + ' used as library for the models no interpolated'
      print ('No interpolation for the library '+new_library)
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'User file ' + new_library + ' used as library for the models interpolated'
      print ('Interpolation for the library '+new_library)
 
#Valuable columns of the files
ir_lin = ['12logOH', 'logNO', 'logU', 'HI_4m', 'ArII_7m', 'HI_7m', 'ArV_8m', 'ArIII_9m', 'SIV_10m', 'HI_12m', 'NeII_12m', 'ArV_13m', 'NeV_14m', 'NeIII_15m', 'SIII_19m', 'NeV_24m', 'OIV_26m', 'SIII_33m', 'OIII_52m', 'NIII_57m', 'OIII_88m', 'NII_122m', 'NII_205m']
lin_ir_label = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'HI_4m', 'ArII_7m', 'HI_7m', 'ArV_8m', 'ArIII_9m', 'SIV_10m', 'HI_12m', 'NeII_12m', 'ArV_13m', 'NeV_14m', 'NeIII_15m', 'SIII_19m', 'NeV_24m', 'OIV_26m', 'SIII_33m', 'OIII_52m', 'NIII_57m', 'OIII_88m', 'NII_122m', 'NII_205m']

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
with open('Libraries_ir/'+file_lib, 'r') as file_aux:
   for line in file_aux:
      if line[0] == '#':
         list_comments.append(line)

#Storing columns:
lin_ir_col = []
#Retrieving each column of the grid
for label in ir_lin:
   aux_col = grid_aux[label].tolist()
   lin_ir_col.append(aux_col)

#Comments
grid_to_write = open('Libraries_ir/'+file_lib, 'w')
for line_com in list_comments:
   grid_to_write.write(line_com)
#Header line
label_line = '{:15} '.format(lin_ir_label[0].replace(' ',''))
for ind in range(1, len(lin_ir_label)-1):
   label_line += '\t {:15} '.format(lin_ir_label[ind].replace(' ',''))
label_line += '\t {:15}\n'.format(lin_ir_label[-1].replace(' ','')) 
grid_to_write.write(label_line)
#Values:
for ind_val in index_OH_NO_U_sorted:
   val_line = '{:7.7f} '.format(lin_ir_col[0][ind_val])
   for ind2 in range(1, len(lin_ir_label)-1):
      val_line += '\t {:7.7f} '.format(lin_ir_col[ind2][ind_val])
   val_line += '\t {:7.7f}\n'.format(lin_ir_col[-1][ind_val])
   grid_to_write.write(val_line)        
grid_to_write.close()

#Opening sorted grid of models
n_comments = 0
with open('Libraries_ir/'+file_lib, 'r') as file12:
   for line in file12:
      if line[0] == '#':
         n_comments += 1  
grid_aux = np.genfromtxt('Libraries_ir/'+file_lib, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

################################################
###### CONSTRAINTS FOR THE GRID OF MODELS ######
################################################
    
#Reading constraints and creating library with constraints
print ('')
print ('Select a file with the constraints to be used to limit the grid of models when the measurement of a quantity is impossible without any relation.')
print ('')

#Displayig options for the user
print ('')
while question7:
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
   if const == '1' or const == '2' or const == '3' or const == '4' or const == '5' or const == '6': question7 = False 
print ('')
const = int(const)

#Particular file introduced by the user
if const == 6:
   while question8:
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
         with open('Constraints/'+new_const, 'r') as file:
            for line in file:
               if line[0] == '#':
                  n_comments += 1
         const_user = np.genfromtxt('Constraints/'+new_const, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         print(' ')
         print('Loading constraint file '+new_const+'. Checking correct format of the file.')
         question8 = False
      except:
         print(' ')
         print('File was not found in folder "Constraints" or file does not exist.')
   question9 = True
   while question9:
      try:
         #Counting comments:
         n_comments = 0
         with open('Constraints/'+new_const, 'r') as file:
            for line in file:
               if line[0] == '#':
                  n_comments += 1
         const_user = np.genfromtxt('Constraints/'+new_const, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
         #Checking correct format:
         #Counting comments:
         n_comments = 0
         with open('Constraints/template_OH.dat', 'r') as file:
            for line in file:
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
            question9 = False
         else:
            print('File does not present the correct format. The following columns are missing:')
            for need_label in missing_labels:
               print('- '+need_label)
            print('More details on the correct format for the constraint file are found in readme file.')
            if interactive == True:
               print(' ')
               print('Reintroduce name of the file with fixed format.')
               print(' ')
               if int(sys.version[0]) < 3:
                  new_const = raw_input('Name of file: ')
               else:
                  new_const = input('Name of file: ')
            elif interactive == False:
               print (' ')
               print ('Re-run the program. Please correct the format of the constraint file before running the code again')
               sys.exit()
      except:
         if interactive == True:
            print('Something went wrong while reading file. Please, reintroduce name of the file.')
            print(' ')
            if int(sys.version[0]) < 3:
               new_const = raw_input('Name of file: ')
            else:
               new_const = input('Name of file: ')
         elif interactive == False:
            print ('Something went wrong while reading file. Please, check again the file introduced, before running the code again.')
            sys.exit()

#Generation of grids with constraints laws:
if const == 1 or const == 2 or const == 3 or const == 6:
   #First grid does not change
   grid1 = grid_aux
   file_lib_2 = file_lib
elif const == 4 or const == 5:
   lin_ir_agn = []
   #The initial grid need to be constrained in the ionization parameter
   if const == 4:
      U_max = 0.0
      U_min = -2.5
      tag = 'high'
   if const == 5:
      U_max = -2.5
      U_min = -4.0
      tag = 'low'
   #Retrieving each column of the grid
   for label in ir_lin:
      aux_col = grid_aux[label].tolist()
      lin_ir_agn.append(aux_col)
   #Creation of the grid
   file_lib_2 = '.'.join(file_lib.split('.')[0:-1])+'_'+tag+'.'+file_lib.split('.')[-1]
   file_open = open('Libraries_ir/'+file_lib_2, 'w')
   file_open.write('#Library constrained for '+tag+' ionization AGNs\n')
   #Header line
   label_line = '{:15} '.format(lin_ir_label[0].replace(' ',''))
   for ind in range(1, len(lin_ir_label)-1):
      label_line += '\t {:15} '.format(lin_ir_label[ind].replace(' ',''))
   label_line += '\t {:15}\n'.format(lin_ir_label[-1].replace(' ','')) 
   file_open.write(label_line)
   #Values:
   for ind_val in range(0, len(lin_ir_agn[0])):
      if lin_ir_agn[2][ind_val] <= U_max and lin_ir_agn[2][ind_val] >= U_min:
         val_line = '{:7.7f} '.format(lin_ir_agn[0][ind_val])
         for ind2 in range(1, len(lin_ir_label)-1):
            val_line += '\t {:7.7f} '.format(lin_ir_agn[ind2][ind_val])
         val_line += '\t {:7.7f}\n'.format(lin_ir_agn[-1][ind_val])
         file_open.write(val_line)        
   file_open.close()
   #Counting comments:
   n_comments = 0
   with open('Libraries_ir/'+file_lib_2, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid1 = np.genfromtxt('Libraries_ir/'+file_lib_2,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
    
#Generating libraries for the constraints in the files
if const == 1: #Star-Forming Galaxies
   const_file = 'template_OH.dat'
   name_const = 'Constraints/template_OH.dat'
   n_comments = 0
   with open(name_const, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == 2:
   const_file = 'template_OH_eelg.dat'
   name_const = 'Constraints/template_OH_eelg.dat'
   n_comments = 0
   with open(name_const, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == 3:
   name_const = 'Constraints/template_OH_agn.dat'
   const_file = 'template_OH_agn.dat'
   n_comments = 0
   with open(name_const, 'r') as file18:
      for line in file18:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == 4:
   name_const = 'Constraints/template_OH_agn_high.dat'
   const_file = 'template_OH_agn_high.dat'
   n_comments = 0
   with open(name_const, 'r') as file19:
      for line in file19:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == 5:
   name_const = 'Constraints/template_OH_agn_low.dat'
   const_file = 'template_OH_agn_low.dat'
   n_comments = 0
   with open(name_const, 'r') as file20:
      for line in file20:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == 6:
   const_file = new_const
   name_const = 'Constraints/'+new_const
   n_comments = 0
   with open(name_const, 'r') as file21:
      for line in file21:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

#Limiting the grids:
lin_ir_val = []
#The initial grid need to be constrained in the ionization parameter
#Retrieving each column of the grid
for label in ir_lin:
   aux_col = grid1[label].tolist()
   lin_ir_val.append(aux_col)
   #Creation of the grids
   name_OH_U = '.'.join(file_lib_2.split('.')[0:-1])+'_OH_U_constrained.'+file_lib.split('.')[-1]
   name_OH_U_NO = '.'.join(file_lib_2.split('.')[0:-1])+'_OH_U_NO_constrained.'+file_lib.split('.')[-1]
   file_open = open('Libraries_ir/'+ name_OH_U, 'w') #OH and U relation
   file_open_2 = open('Libraries_ir/'+name_OH_U_NO, 'w') #OH, NO and U relation
   file_open.write('#Constrained by relation between 12+log(O/H) and log(U)\n')
   file_open_2.write('#Constrained by relation between 12+log(O/H), log(U) and log(N/O)\n')
   #Header line
   label_line = '{:15} '.format(lin_ir_label[0].replace(' ',''))
   for ind in range(1, len(lin_ir_label)-1):
      label_line += '\t {:15} '.format(lin_ir_label[ind].replace(' ',''))
   label_line += '\t {:15}\n'.format(lin_ir_label[-1].replace(' ','')) 
   file_open.write(label_line)
   file_open_2.write(label_line)
#Values:
for ind_val in range(0, len(lin_ir_val[0])):
   index_desired = np.where(const_data['12logOH'] == lin_ir_val[0][ind_val])[0][0] #Searching for constrain in given value of O/H
   if lin_ir_val[2][ind_val] <= const_data['logU_max'][index_desired] and lin_ir_val[2][ind_val] >= const_data['logU_min'][index_desired]:
      val_line = '{:7.7f} '.format(lin_ir_val[0][ind_val])
      for ind2 in range(1, len(lin_ir_label)-1):
         val_line += '\t {:7.7f} '.format(lin_ir_val[ind2][ind_val])
      val_line += '\t {:7.7f}\n'.format(lin_ir_val[-1][ind_val])
      file_open.write(val_line)
   if lin_ir_val[2][ind_val] <= const_data['logU_max'][index_desired] and lin_ir_val[2][ind_val] >= const_data['logU_min'][index_desired] and lin_ir_val[1][ind_val] <= const_data['logNO_max'][index_desired] and lin_ir_val[1][ind_val] >= const_data['logNO_min'][index_desired]:
      val_line = '{:7.7f} '.format(lin_ir_val[0][ind_val])
      for ind2 in range(1, len(lin_ir_label)-1):
         val_line += '\t {:7.7f} '.format(lin_ir_val[ind2][ind_val])
      val_line += '\t {:7.7f}\n'.format(lin_ir_val[-1][ind_val])
      file_open_2.write(val_line)
file_open.close()
file_open_2.close()
#Counting comments:
n_comments = 0
with open('Libraries_ir/'+name_OH_U, 'r') as file:
   for line in file:
      if line[0] == '#':
         n_comments += 1
grid2 = np.genfromtxt('Libraries_ir/'+name_OH_U,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

n_comments = 0
with open('Libraries_ir/'+name_OH_U_NO, 'r') as file:
   for line in file:
      if line[0] == '#':
         n_comments += 1
grid3 = np.genfromtxt('Libraries_ir/'+name_OH_U_NO,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

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
print (' ')

#################################################
###### CREATING ARRAY TO STORE ESTIMATIONS ######
#################################################

grids = []
OHffs = []
eOHffs = []
SHffs = []
eSHffs = []
NOffs = []
eNOffs = []
logUffs = []
elogUffs = []

#Labels to check information provided in the input file
Label_ID = False
Label_HI_4m = False
Label_eHI_4m = False
Label_Hbeta = False
Label_eHbeta = False
Label_ArII = False
Label_eArII = False
Label_HI_7m = False
Label_eHI_7m = False
Label_ArV_8m = False
Label_eArV_8m = False
Label_ArIII = False
Label_eArIII = False
Label_SIV = False
Label_eSIV = False
Label_HI_12m = False
Label_eHI_12m = False
Label_NeII = False
Label_eNeII = False
Label_ArV_13m = False
Label_eArV_13m = False
Label_NeV_14m = False
Label_eNeV_14m = False
Label_NeIII = False
Label_eNeIII = False
Label_SIII_18m = False
Label_eSIII_18m = False
Label_NeV_24m = False
Label_eNeV_24m = False
Label_SIII_33m = False
Label_eSIII_33m = False
Label_OIV = False
Label_eOIV = False
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

#Checking input information
for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == 'Hbeta':
      Label_Hbeta = True
   if input1.dtype.names[col] == 'eHbeta':
      Label_eHbeta = True
   if input1.dtype.names[col] == 'HI_4m':
      Label_HI_4m = True
   if input1.dtype.names[col] == 'eHI_4m':
      Label_eHI_4m = True
   if input1.dtype.names[col] == 'ArII_7m':
      Label_ArII = True
   if input1.dtype.names[col] == 'eArII_7m':
      Label_eArII = True
   if input1.dtype.names[col] == 'HI_7m':
      Label_HI_7m = True
   if input1.dtype.names[col] == 'eHI_7m':
      Label_eHI_7m = True
   if input1.dtype.names[col] == 'ArV_8m':
      Label_ArV_8m = True
   if input1.dtype.names[col] == 'eArV_8m':
      Label_eArV_8m = True
   if input1.dtype.names[col] == 'ArIII_9m':
      Label_ArIII = True
   if input1.dtype.names[col] == 'eArIII_9m':
      Label_eArIII = True
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
   if input1.dtype.names[col] == 'ArV_13m':
      Label_ArV_13m = True
   if input1.dtype.names[col] == 'eArV_13m':
      Label_eArV_13m = True
   if input1.dtype.names[col] == 'NeV_14m':
      Label_NeV_14m = True
   if input1.dtype.names[col] == 'eNeV_14m':
      Label_eNeV_14m = True
   if input1.dtype.names[col] == 'NeIII_15m':
      Label_NeIII = True
   if input1.dtype.names[col] == 'eNeIII_15m':
      Label_eNeIII = True
   if input1.dtype.names[col] == 'SIII_18m':
      Label_SIII_18m = True
   if input1.dtype.names[col] == 'eSIII_18m':
      Label_eSIII_18m = True
   if input1.dtype.names[col] == 'NeV_24m':
      Label_NeV_24m = True
   if input1.dtype.names[col] == 'eNeV_24m':
      Label_eNeV_24m = True 
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

#Adapting final output with information from given input
if Label_ID == False:
   Names = np.arange(1,input1.size+1,1)
else:
   Names = input1['ID']
if Label_Hbeta == False:
   Hbeta = np.zeros(input1.size)
else:
   Hbeta = input1['Hbeta']
if Label_eHbeta == False:
   eHbeta = np.zeros(input1.size)
else:
   eHbeta = input1['eHbeta']
if Label_HI_4m == False:
   HI_4m = np.zeros(input1.size)
else:
   HI_4m = input1['HI_4m']
if Label_eHI_4m == False:
   eHI_4m = np.zeros(input1.size)
else:
   eHI_4m = input1['eHI_4m']
if Label_ArII == False:
   ArII_7m = np.zeros(input1.size)
else:
   ArII_7m = input1['ArII_7m']
if Label_eArII == False:
   eArII_7m = np.zeros(input1.size)
else:
   eArII_7m = input1['eArII_7m']
if Label_HI_7m == False:
   HI_7m = np.zeros(input1.size)
else:
   HI_7m = input1['HI_7m']
if Label_eHI_7m == False:
   eHI_7m = np.zeros(input1.size)
else:
   eHI_7m = input1['eHI_7m']
if Label_ArV_8m == False:
   ArV_8m = np.zeros(input1.size)
else:
   ArV_8m = input1['ArV_8m']
if Label_eArV_8m == False:
   eArV_8m = np.zeros(input1.size)
else:
   eArV_8m = input1['eArV_8m']
if Label_ArIII == False:
   ArIII_9m = np.zeros(input1.size)
else:
   ArIII_9m = input1['ArIII_9m']
if Label_eArIII == False:
   eArIII_9m = np.zeros(input1.size)
else:
   eArIII_9m = input1['eArIII_9m']
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
if Label_ArV_13m == False:
   ArV_13m = np.zeros(input1.size)
else:
   ArV_13m = input1['ArV_13m']
if Label_eArV_13m == False:
   eArV_13m = np.zeros(input1.size)
else:
   eArV_13m = input1['eArV_13m']
if Label_NeV_14m == False:
    NeV_14m = np.zeros(input1.size)
else:
    NeV_14m = input1['NeV_14m']
if Label_eNeV_14m == False:
    eNeV_14m = np.zeros(input1.size)
else:
    eNeV_14m = input1['eNeV_14m'] 
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
if Label_NeV_24m == False:
    NeV_24m = np.zeros(input1.size)
else:
    NeV_24m = input1['NeV_24m']
if Label_eNeV_24m == False:
    eNeV_24m = np.zeros(input1.size)
else:
    eNeV_24m = input1['eNeV_24m'] 
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


################################################################
###### OUTPUT FORMAT AND INFORMATION: ONLY EMISSION LINES ######
################################################################

#Creation of output only with information from inputs
aux_list = []
aux_list.append(('ID','U12'))
if Label_Hbeta == True:
   aux_list.append(('Hbeta', float))
if Label_eHbeta == True:
   aux_list.append(('eHbeta', float))
if Label_HI_4m == True:
   aux_list.append(('HI_4m', float))
if Label_eHI_4m == True:
   aux_list.append(('eHI_4m', float))
if Label_ArII == True:
   aux_list.append(('ArII_7m', float))
if Label_eArII == True:
   aux_list.append(('eArII_7m', float))
if Label_HI_7m == True:
   aux_list.append(('HI_7m', float))
if Label_eHI_7m == True:
   aux_list.append(('eHI_7m', float))
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
if Label_HI_12m == True:
   aux_list.append(('HI_12m', float))
if Label_eHI_12m == True:
   aux_list.append(('eHI_12m', float))
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

aux_list.append(('grid', int))
aux_list.append(('OH', float))
aux_list.append(('eOH', float))
aux_list.append(('SH', float))
aux_list.append(('eSH', float))
aux_list.append(('NO', float))
aux_list.append(('eNO', float))
aux_list.append(('logU', float))
aux_list.append(('elogU', float))
output = np.zeros(input1.size, dtype=aux_list)

output['ID'] = Names
if Label_Hbeta == True:
   output['Hbeta'] = Hbeta
if Label_eHbeta == True:
   output['eHbeta'] = eHbeta
if Label_HI_4m == True:
   output['HI_4m'] = HI_4m
if Label_eHI_4m == True:
   output['eHI_4m'] = eHI_4m
if Label_ArII == True:
   output['ArII_7m'] = ArII_7m
if Label_eArII == True:
   output['eArII_7m'] = eArII_7m
if Label_HI_7m == True:
   output['HI_7m'] = HI_7m
if Label_eHI_7m == True:
   output['eHI_7m'] = eHI_7m
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
if Label_HI_12m == True:
   output['HI_12m'] = HI_12m
if Label_eHI_12m == True:
   output['eHI_12m'] = eHI_12m
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
###### ESTIMATIONS OF CHEMICAL ABUNDANCES ######
################################################

#Display for the user
print ('Calculating...')
print ('')
print ('')
print ('----------------------------------------------------------------')
print ('(%)   ID    Grid  12+log(O/H)  12+log(S/H)  log(N/O)    log(U)')
print ('----------------------------------------------------------------')


#Beginning of loop of calculation
count = 0
for tab in range(0,len(input1),1):
   try:
      count = count + 1
      OH_mc = []
      NO_mc = []
      logU_mc = []
      SH_mc = []
      OHe_mc = []
      NOe_mc = []
      logUe_mc = []
      SHe_mc = []  

      #Selection of grid   
      if NIII_57m[tab] > 0 and (OIII_52m[tab] > 0 or OIII_88m[tab] > 0 or SIII_18m[tab] > 0 or SIII_33m[tab] > 0 or SIV_10m[tab] > 0):
         grid = grid2
         grid_type = 2
         grids.append(2)
      else:
         grid = grid3
         grid_type = 3
         grids.append(3)      

      ######################
      # Calculation of N/O #
      ######################

      if NIII_57m[tab] <= 0 or (OIII_52m[tab] <= 0 and OIII_88m[tab] <= 0 and SIII_18m[tab] <= 0 and SIII_33m[tab] <= 0 and SIV_10m[tab] <= 0):
         NOff = -10
         eNOff = 0
      else:
         for monte in range(0,n,1):
            NO_p = 0
            den_NO = 0
            NO_e = 0
            den_NO_e = 0
            tol_max = 1e2

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
            SIV_10m_obs = 0
            if SIV_10m[tab] > 0:
               while SIV_10m_obs <= 0:
                  SIV_10m_obs = np.random.normal(SIV_10m[tab], eSIV_10m[tab]+1e-5)
            SIII_18m_obs = 0
            if SIII_18m[tab] > 0:
               while SIII_18m_obs <= 0:
                  SIII_18m_obs = np.random.normal(SIII_18m[tab], eSIII_18m[tab]+1e-5)
            SIII_33m_obs = 0
            if SIII_33m[tab] > 0:
               while SIII_33m_obs <= 0:
                  SIII_33m_obs = np.random.normal(SIII_33m[tab], eSIII_33m[tab]+1e-5)

            N3O3a_obs = -10
            N3O3b_obs = -10
            N3S34a_obs = -10
            N3S34b_obs = -10

            if OIII_52m_obs > 0:
               N3O3a_obs = np.log10(NIII_57m_obs / OIII_52m_obs)
            if (OIII_88m_obs > 0 and OIII_52m_obs == 0) or (sed !=3 and OIII_88m_obs > 0):
               N3O3b_obs = np.log10(NIII_57m_obs / OIII_88m_obs)
            if (SIII_18m_obs > 0 or SIV_10m_obs > 0):
               N3S34a_obs = np.log10(NIII_57m_obs / (SIII_18m_obs + SIV_10m_obs) )
            if ( (SIII_33m_obs > 0 or SIV_10m_obs > 0) and SIII_18m_obs == 0) or (sed !=3 and (SIII_33m_obs > 0 or SIV_10m_obs > 0) ):
               N3S34b_obs = np.log10(NIII_57m_obs / (SIII_33m_obs + SIV_10m_obs) )

            CHI_N3O3a = 0
            CHI_N3O3b = 0
            CHI_N3S34a = 0
            CHI_N3S34b = 0
            CHI_NO = 0

            for index in grid:
               if N3O3a_obs == -10: 
                  CHI_N3O3a = 0
               elif index['NIII_57m'] == 0 or index['OIII_52m'] == 0:
                  CHI_N3O3a = tol_max
               else:   
                  CHI_N3O3a = (np.log10(index['NIII_57m']/index['OIII_52m'])- N3O3a_obs)**2/np.log10(index['NIII_57m']/index['OIII_52m'])          
               if N3O3b_obs == -10: 
                  CHI_N3O3b = 0
               elif index['OIII_88m'] == 0 or index['NIII_57m'] == 0:
                  CHI_N3O3b = tol_max
               else:   
                  CHI_N3O3b = (np.log10(index['NIII_57m']/index['OIII_88m'])- N3O3b_obs)**2/np.log10(index['NIII_57m']/index['OIII_88m'])
               if N3S34a_obs == -10: 
                  CHI_N3S34a = 0
               elif (index['SIII_19m'] == 0 and index['SIV_10m'] == 0) or index['NIII_57m'] == 0:
                  CHI_N3S34a = tol_max
               else:   
                  CHI_N3S34a = (np.log10(index['NIII_57m']/(index['SIII_19m'] + index['SIV_10m']) )- N3S34a_obs)**2/np.log10(index['NIII_57m']/(index['SIII_19m'] + index['SIV_10m']) )
               
               if N3S34b_obs == -10: 
                  CHI_N3S34b = 0
               elif (index['SIII_33m'] == 0 and index['SIV_10m'] == 0) or index['NIII_57m'] == 0:
                  CHI_N3S34b = tol_max
               else:   
                  CHI_N3S34b = (np.log10(index['NIII_57m']/(index['SIII_33m'] + index['SIV_10m']) )- N3S34b_obs)**2/np.log10(index['NIII_57m']/(index['SIII_33m'] + index['SIV_10m']) )

               CHI_NO = (CHI_N3O3a**2 + CHI_N3O3b**2 + CHI_N3S34a**2 + CHI_N3S34b**2 )**0.5

               if CHI_NO == 0:
                  NO_p = NO_p
                  den_NO = den_NO
               else:
                  NO_p = index['logNO'] / (CHI_NO) + NO_p
                  den_NO = 1 / (CHI_NO) + den_NO

            NO = NO_p / den_NO 

            #Calculation of N/O error
            CHI_N3O3a = 0
            CHI_N3O3b = 0
            CHI_N3S34a = 0
            CHI_N3S34b = 0
            CHI_NO = 0
            for index in grid:
               if N3O3a_obs == -10: 
                  CHI_N3O3a = 0
               elif index['OIII_52m'] == 0 or index['NIII_57m'] == 0:
                  CHI_N3O3a = tol_max
               else:   
                  CHI_N3O3a = (np.log10(index['NIII_57m']/index['OIII_52m'])- N3O3a_obs)**2/np.log10(index['NIII_57m']/index['OIII_52m'])            
               if N3O3b_obs == -10: 
                  CHI_N3O3b = 0
               elif index['OIII_88m'] == 0 or index['NIII_57m'] == 0:
                  CHI_N3O3b = tol_max
               else:   
                  CHI_N3O3b = (np.log10(index['NIII_57m']/index['OIII_88m'])- N3O3b_obs)**2/np.log10(index['NIII_57m']/index['OIII_88m'])          
               if N3S34a_obs == -10: 
                  CHI_N3S34a = 0
               elif (index['SIII_19m'] == 0 and index['SIV_10m'] == 0) or index['NIII_57m'] == 0:
                  CHI_N3S34a = tol_max
               else:   
                  CHI_N3S34a = (np.log10(index['NIII_57m']/(index['SIII_19m'] + index['SIV_10m']))- N3S34a_obs)**2/np.log10(index['NIII_57m']/(index['SIII_19m'] + index['SIV_10m']))
               if N3S34b_obs == -10: 
                  CHI_N3S34b = 0
               elif index['NIII_57m'] == 0 or (index['SIII_33m'] == 0 and index['SIV_10m'] == 0):
                  CHI_N3S34b = tol_max
               else:   
                  CHI_N3S34b = (np.log10(index['NIII_57m']/(index['SIII_33m'] + index['SIV_10m']))- N3S34b_obs)**2/np.log10(index['NIII_57m']/(index['SIII_33m'] + index['SIV_10m']))
               
               CHI_NO = (CHI_N3O3a**2 + CHI_N3O3b**2 + CHI_N3S34a**2 + CHI_N3S34b**2)**0.5

               if CHI_NO == 0:
                  NO_e = NO_e
                  den_NO_e = den_NO_e
               else:
                  NO_e = (index['logNO'] - NO)**2 / (CHI_NO) + NO_e
                  den_NO_e = 1 / (CHI_NO) + den_NO_e

            eNO = NO_e / den_NO_e

            #Iterations for the interpolation mode
            if inter == 0 or NO == -10:
               NOf = NO
            elif inter == 1:
               igrid = grid[np.lexsort((grid['12logOH'],grid['logU']))]
               igrid = interpolate(igrid,1,NO-eNO-0.125,NO+eNO+0.125,10)

               CHI_N3O3a = 0
               CHI_N3O3b = 0
               CHI_N3S34a = 0
               CHI_N3S34b = 0
               CHI_NO = 0
               NO_p = 0
               den_NO = 0

               for index in igrid:
                  if N3O3a_obs == -10: 
                     CHI_N3O3a = 0
                  elif index['OIII_52m'] == 0 or index['NIII_57m'] == 0:
                     CHI_N3O3a = tol_max
                  else:   
                     CHI_N3O3a = (np.log10(index['NIII_57m']/index['OIII_52m'])- N3O3a_obs)**2/np.log10(index['NIII_57m']/index['OIII_52m'])
                  if N3O3b_obs == -10: 
                     CHI_N3O3b = 0
                  elif index['NIII_57m'] == 0 or index['OIII_88m'] == 0:
                     CHI_N3O3b = tol_max
                  else:   
                     CHI_N3O3b = (np.log10(index['NIII_57m']/index['OIII_88m'])- N3O3b_obs)**2/np.log10(index['NIII_57m']/index['OIII_88m']+1e-5)
                  if N3S34a_obs == -10: 
                     CHI_N3S34a = 0
                  elif (index['SIII_19m'] == 0 and index['SIV_10m']) or index['NIII_57m'] == 0:
                     CHI_N3S34a = tol_max
                  else:   
                     CHI_N3S34a = (np.log10(index['NIII_57m']/(index['SIII_19m'] + index['SIV_10m']))- N3S34a_obs)**2/np.log10(index['NIII_57m']/(index['SIII_19m'] + index['SIV_10m']))
                  if N3S34b_obs == -10: 
                     CHI_N3S34b = 0
                  elif (index['SIII_33m'] == 0 and index['SIV_10m'] == 0) or index['NIII_57m'] == 0:
                     CHI_N3S34b = tol_max
                  else:   
                     CHI_N3S34b = (np.log10(index['NIII_57m']/(index['SIII_33m'] + index['SIV_10m']))- N3S34b_obs)**2/np.log10(index['NIII_57m']/(index['SIII_33m'] + index['SIV_10m'])+1e-5)

                  CHI_NO = (CHI_N3O3a**2 + CHI_N3O3b**2  + CHI_N3S34a**2 + CHI_N3S34b**2 )**0.5

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

      #Creation of a constrained grid on N/O
      if NOff == -10:
         grid_c = grid
      else:
         grid_mac = []
         for index in grid:
            if np.abs(index['logNO'] - NOff) > np.abs(eNOff+res_NO):
               continue
            else:
               grid_mac.append(index[[name_col for name_col in grid_aux.dtype.names]])
               format_array = []
               for name_col in grid_aux.dtype.names:
                  format_array.append((name_col, grid_aux[0][name_col].dtype))
            
               #grid_c_1 = np.reshape(grid_mac,(int(len(grid_mac)/len(grid_aux.dtype.names[:-3])),len(grid_aux.dtype.names[:-3]))) 
               #Old version
               grid_c = np.array(grid_mac, dtype=format_array)
   
      ####################################
      # Calculation of O/H, S/H and logU #
      ####################################

      for monte in range(0,n,1):

         OH_p = 0
         logU_p = 0
         SH_p = 0
         den_OH = 0
         den_SH = 0
         OH_e = 0
         SH_e = 0
         logU_e = 0
         den_OH_e = 0
         den_SH_e = 0
         tol_max = 1e2
         
         Hbeta_obs = 0
         if Hbeta[tab] > 0:
            while Hbeta_obs <= 0:
               Hbeta_obs = np.random-normal(Hbeta[tab], eHbeta[tab]+1e-5)
         HI_4m_obs = 0
         if HI_4m[tab] > 0:
            while HI_4m_obs <= 0:
               HI_4m_obs = np.random.normal(HI_4m[tab],eHI_4m[tab]+1e-5)
         ArII_7m_obs = 0
         if ArII_7m[tab] > 0:
            while ArII_7m_obs <= 0:
               ArII_7m_obs = np.random.normal(ArII_7m[tab],eArII_7m[tab]+1e-5)
         HI_7m_obs = 0
         if HI_7m[tab] > 0:
            while HI_7m_obs <= 0:
               HI_7m_obs = np.random.normal(HI_7m[tab],eHI_7m[tab]+1e-5)
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
         HI_12m_obs = 0
         if HI_12m[tab] > 0:
            while HI_12m_obs <= 0:
               HI_12m_obs = np.random.normal(HI_12m[tab],eHI_12m[tab]+1e-5)
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
                  NeV_14m_obs = np.random.normal(NeV_14m[tab], eNeV_14m[tab]+1e-3)
         NeIII_15m_obs = 0
         if NeIII_15m[tab] > 0:
            while NeIII_15m_obs <= 0:
               NeIII_15m_obs = np.random.normal(NeIII_15m[tab],eNeIII_15m[tab]+1e-3)
         SIII_18m_obs = 0
         if SIII_18m[tab] > 0:
            while SIII_18m_obs <= 0:
               SIII_18m_obs = np.random.normal(SIII_18m[tab],eSIII_18m[tab]+1e-3)
         NeV_24m_obs = 0
         if NeV_24m[tab] > 0:
               while NeV_24m_obs <= 0:
                  NeV_24m_obs = np.random.normal(NeV_24m[tab], eNeV_24m[tab]+1e-3)
         OIV_26m_obs = 0
         if OIV_26m[tab] > 0:
               while OIV_26m_obs <= 0:
                  OIV_26m_obs = np.random.normal(OIV_26m[tab], eOIV_26m[tab]+1e-3)
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

         #Estimators based on Ne lines      
         if HI_4m_obs == 0 or (NeII_12m_obs == 0 and NeIII_15m_obs == 0 and NeV_14m_obs == 0):
            Ne235a_obs = -10
         else:
            Ne235a_obs = np.log10((NeII_12m_obs + NeIII_15m_obs + NeV_14m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (NeII_12m_obs == 0 and NeIII_15m_obs == 0 and NeV_14m_obs == 0):
            Ne235b_obs = -10
         else:
            Ne235b_obs = np.log10((NeII_12m_obs + NeIII_15m_obs + NeV_14m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (NeII_12m_obs == 0 and NeIII_15m_obs == 0 and NeV_14m_obs == 0):
            Ne235c_obs = -10
         else:
            Ne235c_obs = np.log10((NeII_12m_obs + NeIII_15m_obs + NeV_14m_obs)/HI_12m_obs)
         if HI_4m_obs == 0 or (NeII_12m_obs == 0 and NeIII_15m_obs == 0 and NeV_24m_obs == 0):
            Ne235d_obs = -10
         else:
            Ne235d_obs = np.log10((NeII_12m_obs + NeIII_15m_obs + NeV_24m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (NeII_12m_obs == 0 and NeIII_15m_obs == 0 and NeV_24m_obs == 0):
            Ne235e_obs = -10
         else:
            Ne235e_obs = np.log10((NeII_12m_obs + NeIII_15m_obs + NeV_24m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (NeII_12m_obs == 0 and NeIII_15m_obs == 0 and NeV_24m_obs == 0):
            Ne235f_obs = -10
         else:
            Ne235f_obs = np.log10((NeII_12m_obs + NeIII_15m_obs + NeV_24m_obs)/HI_12m_obs)       
         if NeII_12m_obs != 0 and NeIII_15m_obs != 0 and NeV_14m_obs == 0 and NeV_24m_obs == 0 and sed!=3:
            Ne2Ne3_obs = np.log10(NeII_12m_obs / NeIII_15m_obs)
         else:
            Ne2Ne3_obs = -10      
         if ((NeIII_15m_obs == 0 and NeII_12m_obs == 0) or NeV_14m_obs == 0 ):
            Ne23Ne5a_obs = -10
         else:
            Ne23Ne5a_obs = np.log10(( (NeIII_15m_obs + NeII_12m_obs) / NeV_14m_obs))
         if ((NeIII_15m_obs == 0 and NeII_12m_obs == 0) or NeV_24m_obs == 0  ):
            Ne23Ne5b_obs = -10
         else:
            Ne23Ne5b_obs = np.log10(( (NeIII_15m_obs + NeII_12m_obs) / NeV_24m_obs))
         if Hbeta_obs == 0 or (NeII_12m_obs and NeIII_15m_obs == 0 and NeV_14m_obs == 0):
            Ne235g_obs = -10
         else:
            Ne235g_obs = np.log10(( (NeII_12m_obs + NeIII_15m_obs + NeV_24m_obs)/Hbeta_obs ))
         if Hbeta_obs == 0 or (NeII_12m_obs and NeIII_15m_obs == 0 and NeV_24m_obs == 0):
            Ne235h_obs = -10
         else:
            Ne235h_obs = np.log10(( (NeII_12m_obs + NeIII_15m_obs + NeV_24m_obs)/Hbeta_obs ))
         #Estimators based on Argon lines
         if HI_4m_obs == 0 or (ArII_7m_obs == 0 and ArIII_9m_obs == 0 and ArV_8m_obs == 0):
            Ar235a_obs = -10
         else:
            Ar235a_obs = np.log10((ArII_7m_obs + ArIII_9m_obs + ArV_8m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (ArII_7m_obs == 0 and ArIII_9m_obs == 0 and ArV_8m_obs == 0):
            Ar235b_obs = -10
         else:
            Ar235b_obs = np.log10((ArII_7m_obs + ArIII_9m_obs + ArV_8m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (ArII_7m_obs == 0 and ArIII_9m_obs == 0 and ArV_8m_obs == 0):
            Ar235c_obs = -10
         else:
            Ar235c_obs = np.log10((ArII_7m_obs + ArIII_9m_obs + ArV_8m_obs)/HI_12m_obs)
         if HI_4m_obs == 0 or (ArII_7m_obs == 0 and ArIII_9m_obs == 0 and ArV_13m_obs == 0):
            Ar235d_obs = -10
         else:
            Ar235d_obs = np.log10((ArII_7m_obs + ArIII_9m_obs + ArV_13m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (ArII_7m_obs == 0 and ArIII_9m_obs == 0 and ArV_13m_obs == 0):
            Ar235e_obs = -10
         else:
            Ar235e_obs = np.log10((ArII_7m_obs + ArIII_9m_obs + ArV_13m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (ArII_7m_obs == 0 and ArIII_9m_obs == 0 and ArV_13m_obs == 0):
            Ar235f_obs = -10
         else:
            Ar235f_obs = np.log10((ArII_7m_obs + ArIII_9m_obs + ArV_13m_obs)/HI_12m_obs)       
         if ArII_7m_obs != 0 and ArIII_9m_obs != 0 and ArV_8m_obs == 0 and ArV_8m_obs == 0 and sed!=3:
            Ar2Ar3_obs = np.log10(ArII_7m_obs / ArIII_9m_obs)
         else:
            Ar2Ar3_obs = -10      
         if ((ArIII_9m_obs == 0 and ArII_7m_obs == 0) or ArV_8m_obs == 0 ):
            Ar23Ar5a_obs = -10
         else:
            Ar23Ar5a_obs = np.log10(( (ArIII_9m_obs + ArII_7m_obs) / ArV_8m_obs))
         if ((ArIII_9m_obs == 0 and ArII_7m_obs == 0) or ArV_13m_obs == 0  ):
            Ar23Ar5b_obs = -10
         else:
            Ar23Ar5b_obs = np.log10(( (ArIII_9m_obs + ArII_7m_obs) / ArV_13m_obs))
         if Hbeta_obs == 0 or (ArII_7m_obs and ArIII_9m_obs == 0 and ArV_8m_obs == 0):
            Ar235g_obs = -10
         else:
            Ar235g_obs = np.log10(( (ArII_7m_obs + ArIII_9m_obs + ArV_8m_obs)/Hbeta_obs ))
         if Hbeta_obs == 0 or (ArII_7m_obs and ArIII_9m_obs == 0 and ArV_13m_obs == 0):
            Ar235h_obs = -10
         else:
            Ar235h_obs = np.log10(( (ArII_7m_obs + ArIII_9m_obs + ArV_13m_obs)/Hbeta_obs ))
         #Estimators based on Oxygen lines  
         if HI_4m_obs == 0 or (OIV_26m_obs == 0 and OIII_52m_obs == 0):
            O34a_obs = -10
         else:
            O34a_obs = np.log10((OIV_26m_obs + OIII_52m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (OIV_26m_obs == 0 and OIII_52m_obs == 0):
            O34b_obs = -10
         else:
            O34b_obs = np.log10((OIV_26m_obs + OIII_52m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (OIV_26m_obs == 0 and OIII_52m_obs == 0):
            O34c_obs = -10
         else:
            O34c_obs = np.log10((OIV_26m_obs + OIII_52m_obs)/HI_12m_obs)    
         if HI_4m_obs == 0 or (OIV_26m_obs == 0 and OIII_88m_obs == 0) or (OIII_52m_obs > 0 and sed==3 and HI_4m_obs > 0):
            O34d_obs = -10
         else:
            O34d_obs = np.log10((OIV_26m_obs + OIII_88m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (OIV_26m_obs == 0 and OIII_88m_obs == 0) or (OIII_52m_obs > 0 and sed==3 and HI_7m_obs > 0):
            O34e_obs = -10
         else:
            O34e_obs = np.log10((OIV_26m_obs + OIII_88m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (OIV_26m_obs == 0 and OIII_88m_obs == 0) or (OIII_52m_obs > 0 and sed==3 and HI_12m_obs > 0):
            O34f_obs = -10
         else:
            O34f_obs = np.log10((OIV_26m_obs + OIII_88m_obs)/HI_12m_obs)
         if OIV_26m_obs == 0  or OIII_52m_obs == 0:
            O3O4a_obs = -10
         else:
            O3O4a_obs = np.log10((OIII_52m_obs / OIV_26m_obs))
         if OIV_26m_obs == 0  or OIII_88m_obs == 0 or (OIII_52m_obs > 0 and OIV_26m_obs > 0 and sed==3):
            O3O4b_obs = -10
         else:
            O3O4b_obs = np.log10((OIII_88m_obs / OIV_26m_obs))
         if Hbeta_obs == 0 or (OIV_26m_obs == 0 and OIII_52m_obs == 0):
            O34g_obs = -10
         else:
            O34g_obs = np.log10((OIV_26m_obs + OIII_52m_obs)/Hbeta_obs)    
         if Hbeta_obs == 0 or (OIV_26m_obs == 0 and OIII_88m_obs == 0) or (OIII_52m_obs > 0 and sed==3 and Hbeta_obs > 0):
            O34h_obs = -10
         else:
            O34h_obs = np.log10((OIV_26m_obs + OIII_88m_obs)/Hbeta_obs)
         #Estimators based on Sulphur lines
         if HI_4m_obs == 0 or (SIV_10m_obs == 0 and SIII_18m_obs == 0):
            S34a_obs = -10
         else:
            S34a_obs = np.log10((SIV_10m_obs + SIII_18m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (SIV_10m_obs == 0 and SIII_18m_obs == 0):
            S34b_obs = -10
         else:
            S34b_obs = np.log10((SIV_10m_obs + SIII_18m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (SIV_10m_obs == 0 and SIII_18m_obs == 0):
            S34c_obs = -10
         else:
            S34c_obs = np.log10((SIV_10m_obs + SIII_18m_obs)/HI_12m_obs)
         if HI_4m_obs == 0 or (SIV_10m_obs == 0 and SIII_33m_obs == 0) or (SIII_18m_obs > 0 and sed==3 and HI_4m_obs > 0) :
            S34d_obs = -10
         else:
            S34d_obs = np.log10((SIV_10m_obs + SIII_33m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (SIV_10m_obs == 0 and SIII_33m_obs == 0) or (SIII_18m_obs > 0 and sed==3 and HI_7m_obs > 0) :
            S34e_obs = -10
         else:
            S34e_obs = np.log10((SIV_10m_obs + SIII_33m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (SIV_10m_obs == 0 and SIII_33m_obs == 0) or (SIII_18m_obs > 0 and sed==3 and HI_12m_obs > 0) :
            S34f_obs = -10
         else:
            S34f_obs = np.log10((SIV_10m_obs + SIII_33m_obs)/HI_12m_obs)
         if SIV_10m_obs == 0  or SIII_18m_obs == 0:
            S3S4a_obs = -10
         else:
            S3S4a_obs = np.log10((SIII_18m_obs / SIV_10m_obs))
         if SIV_10m_obs == 0  or SIII_33m_obs == 0 or (SIII_18m_obs > 0 and SIV_10m_obs > 0 and sed==3):
            S3S4b_obs = -10
         else:
            S3S4b_obs = np.log10(SIII_33m_obs / SIV_10m_obs)
         if Hbeta_obs == 0 or (SIV_10m_obs == 0 and SIII_18m_obs == 0):
            S34g_obs = -10
         else:
            S34g_obs = np.log10((SIV_10m_obs + SIII_18m_obs)/Hbeta_obs)
         if Hbeta_obs == 0 or (SIV_10m_obs == 0 and SIII_33m_obs == 0) or (SIII_18m_obs > 0 and sed==3 and HI_4m_obs > 0) :
            S34h_obs = -10
         else:
            S34h_obs = np.log10((SIV_10m_obs + SIII_33m_obs)/Hbeta_obs)
         #Estimators based on Nitrogen lines
         if HI_4m_obs == 0 or (NIII_57m_obs == 0 and NII_122m_obs == 0) or sed==3:
            N23a_obs = -10
         else:
            N23a_obs = np.log10((NII_122m_obs + NIII_57m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (NIII_57m_obs == 0 and NII_122m_obs == 0) or sed==3:
            N23b_obs = -10
         else:
            N23b_obs = np.log10((NII_122m_obs + NIII_57m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (NIII_57m_obs == 0 and NII_122m_obs == 0) or sed==3:
            N23c_obs = -10
         else:
            N23c_obs = np.log10((NII_122m_obs + NIII_57m_obs)/HI_12m_obs)
         if Hbeta_obs == 0 or (NIII_57m_obs == 0 and NII_122m_obs == 0) or sed==3:
            N23g_obs = -10
         else:
            N23g_obs = np.log10((NII_122m_obs + NIII_57m_obs)/Hbeta_obs)             
         if (HI_4m_obs == 0 or NIII_57m_obs == 0) or sed!=3:
            N3a_obs = -10
         else:
            N3a_obs = np.log10((NIII_57m_obs)/HI_4m_obs)
         if (HI_7m_obs == 0 or NIII_57m_obs == 0) or sed!=3:
            N3b_obs = -10
         else:
            N3b_obs = np.log10((NIII_57m_obs)/HI_7m_obs)
         if (HI_12m_obs == 0 or NIII_57m_obs == 0) or sed!=3:
            N3c_obs = -10
         else:
            N3c_obs = np.log10((NIII_57m_obs)/HI_12m_obs)
         if (Hbeta_obs == 0 or NIII_57m_obs == 0) or sed!=3:
            N3d_obs = -10
         else:
            N3d_obs = np.log10((NIII_57m_obs)/Hbeta_obs) 
         if NII_122m_obs == 0 or NIII_57m_obs == 0 or sed==3:
            N2N3a_obs = -10
         else:
            N2N3a_obs = np.log10((NII_122m_obs / NIII_57m_obs))
         if HI_4m_obs == 0 or (NIII_57m_obs == 0 and NII_205m_obs == 0) or sed==3:
            N23d_obs = -10
         else:
            N23d_obs = np.log10((NII_205m_obs + NIII_57m_obs)/HI_4m_obs)
         if HI_7m_obs == 0 or (NIII_57m_obs == 0 and NII_205m_obs == 0) or sed==3:
            N23e_obs = -10
         else:
            N23e_obs = np.log10((NII_205m_obs + NIII_57m_obs)/HI_7m_obs)
         if HI_12m_obs == 0 or (NIII_57m_obs == 0 and NII_205m_obs == 0) or sed==3:
            N23f_obs = -10
         else:
            N23f_obs = np.log10((NII_205m_obs + NIII_57m_obs)/HI_12m_obs)
         if Hbeta_obs == 0 or (NIII_57m_obs == 0 and NII_205m_obs == 0) or sed==3:
            N23h_obs = -10
         else:
            N23h_obs = np.log10((NII_205m_obs + NIII_57m_obs)/Hbeta_obs)
         if NII_205m_obs == 0 or NIII_57m_obs == 0 or sed==3:
            N2N3b_obs = -10
         else:
            N2N3b_obs = np.log10((NII_205m_obs / NIII_57m_obs))
         #Estimator O3N2
         if NII_122m_obs == 0 or OIII_52m_obs == 0 or sed == 3:
            O3N2a_obs = -10
         else:
            O3N2a_obs = np.log10((OIII_52m_obs/NII_122m_obs))
         if NII_122m_obs == 0 or OIII_88m_obs == 0 or sed == 3:
            O3N2b_obs = -10
         else:
            O3N2b_obs = np.log10((OIII_88m_obs/NII_122m_obs))
         if NII_205m_obs == 0 or OIII_52m_obs == 0 or sed == 3:
            O3N2c_obs = -10
         else:
            O3N2c_obs = np.log10((OIII_52m_obs/NII_205m_obs))
         if NII_205m_obs == 0 or OIII_88m_obs == 0 or sed == 3:
            O3N2d_obs = -10
         else:
            O3N2d_obs = np.log10((OIII_88m_obs/NII_205m_obs))
               
         if S34a_obs == -10 and S34b_obs == -10 and S34c_obs == -10 and S34d_obs == -10 and S34e_obs == -10 and S34f_obs == -10 and Ne2Ne3_obs == -10 and S3S4a_obs == -10 and S3S4b_obs == -10 and N23a_obs == -10 and N23b_obs == -10 and N23c_obs == -10 and N2N3a_obs == -10 and O3N2a_obs == -10 and O3N2b_obs == -10 and N23d_obs == -10 and N23e_obs == -10 and N23f_obs == -10 and N2N3b_obs == -10 and O3N2c_obs == -10 and O3N2d_obs == -10 and Ne235a_obs == -10 and Ne235b_obs == -10 and Ne235c_obs == -10 and Ne235d_obs == -10 and Ne235e_obs == -10 and Ne235f_obs == -10 and Ne23Ne5a_obs == -10 and Ne23Ne5b_obs == -10 and O34a_obs == -10 and O34b_obs == -10 and O34c_obs == -10 and O34d_obs == -10 and O34e_obs == -10 and O34f_obs == -10 and O3O4a_obs == -10 and O3O4b_obs == -10 and N3a_obs == -10 and N3b_obs == -10 and N3c_obs == -10 and Ar235a_obs == -10 and Ar235b_obs == -10 and Ar235c_obs == -10 and Ar235d_obs == -10 and Ar235e_obs == -10 and Ar235f_obs == -10 and Ar2Ar3_obs == -10 and Ar23Ar5a_obs == -10 and Ar23Ar5b_obs == -10 and Ne235g_obs == -10 and Ne235h_obs == -10 and Ar235g_obs == -10 and Ar235h_obs == -10 and S34g_obs == -10 and S34h_obs == -10 and O34g_obs == -10 and O34h_obs == -10 and N23g_obs == -10 and N23h_obs == -10 and N3d_obs == -10:
            OH = 0
            logU = 0
            SH = 0
         else:
            CHI_Ne2Ne3 = 0   
            CHI_Ne235a = 0
            CHI_Ne235b = 0
            CHI_Ne235c = 0
            CHI_Ne235d = 0
            CHI_Ne235e = 0
            CHI_Ne235f = 0
            CHI_Ne235g = 0
            CHI_Ne235h = 0
            CHI_Ne23Ne5a = 0
            CHI_Ne23Ne5b = 0

            CHI_Ar2Ar3 = 0   
            CHI_Ar235a = 0
            CHI_Ar235b = 0
            CHI_Ar235c = 0
            CHI_Ar235d = 0
            CHI_Ar235e = 0
            CHI_Ar235f = 0
            CHI_Ar235g = 0
            CHI_Ar235h = 0
            CHI_Ar23Ar5a = 0
            CHI_Ar23Ar5b = 0
         
            CHI_O34a = 0
            CHI_O34b = 0
            CHI_O34c = 0
            CHI_O34d = 0
            CHI_O34e = 0
            CHI_O34f = 0
            CHI_O34g = 0
            CHI_O34h = 0
            CHI_O3O4a = 0
            CHI_O3O4b = 0
         
            CHI_S34a = 0
            CHI_S34b = 0
            CHI_S34c = 0
            CHI_S34d = 0
            CHI_S34e = 0
            CHI_S34f = 0
            CHI_S34g = 0
            CHI_S34h = 0
            CHI_S3S4a = 0
            CHI_S3S4b = 0
               
            CHI_N23a = 0
            CHI_N23b = 0
            CHI_N23c = 0      
            CHI_N23d = 0
            CHI_N23e = 0
            CHI_N23f = 0
            CHI_N23g = 0
            CHI_N23h = 0
            CHI_N3a = 0
            CHI_N3b = 0
            CHI_N3c = 0
            CHI_N3d = 0
            CHI_N2N3a = 0
            CHI_N2N3b = 0

            CHI_O3N2a = 0
            CHI_O3N2b = 0
            CHI_O3N2c = 0
            CHI_O3N2d = 0

            CHI_OH = 0
            CHI_SH = 0

            for index in grid_c:
               #Constraints in N/O and high-ionization emission lines
               if NOff > -10 and np.abs(index['logNO'] - NOff) > np.abs(eNOff+res_NO):
                  continue
               elif SIV_10m_obs > 0 and index['SIV_10m'] == 0:
                  continue
               elif OIV_26m_obs > 0 and index['OIV_26m'] == 0:
                  continue
               elif ArV_8m_obs > 0 and index['ArV_8m'] == 0:
                  continue
               elif ArV_13m_obs > 0 and index['ArV_13m'] == 0:
                  continue
               elif NeV_14m_obs > 0 and index['NeV_14m'] == 0:
                  continue
               elif NeV_24m_obs > 0 and index['NeV_24m'] == 0:
                  continue
               else:            
                  #Neon     
                  if Ne2Ne3_obs == -10: 
                     CHI_Ne2Ne3 = 0
                  elif index['NeII_12m'] == 0 or index['NeIII_15m'] == 0:
                     CHI_Ne2Ne3 = tol_max
                  else:   
                     CHI_Ne2Ne3 = (np.log10((index['NeII_12m']/index['NeIII_15m']))- Ne2Ne3_obs)**2/np.log10((index['NeII_12m']/index['NeIII_15m']))
                  if Ne235a_obs == -10: 
                     CHI_Ne235a = 0
                  elif index['HI_4m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235a = tol_max
                  else:   
                     CHI_Ne235a = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_4m'])- Ne235a_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_4m'])
                  if Ne235b_obs == -10: 
                     CHI_Ne235b = 0
                  elif index['HI_7m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235b = tol_max
                  else:   
                     CHI_Ne235b = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_7m'])- Ne235b_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_7m'])
                  if Ne235c_obs == -10: 
                     CHI_Ne235c = 0
                  elif index['HI_12m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235c = tol_max
                  else:   
                     CHI_Ne235c = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_12m'])- Ne235c_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_12m'])
                  if Ne235d_obs == -10: 
                     CHI_Ne235d = 0
                  elif index['HI_4m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235d = tol_max
                  else:   
                     CHI_Ne235d = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_4m'])- Ne235d_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_4m'])
                  if Ne235e_obs == -10: 
                     CHI_Ne235e = 0
                  elif index['HI_7m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235e = tol_max
                  else:   
                     CHI_Ne235e = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_7m'])- Ne235e_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_7m'])
                  if Ne235f_obs == -10: 
                     CHI_Ne235f = 0
                  elif index['HI_12m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235f = tol_max
                  else:   
                     CHI_Ne235f = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_12m'])- Ne235f_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_12m'])
                  if Ne235g_obs == -10:
                     CHI_Ne235g = 0
                  elif (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235g = tol_max
                  else:
                     CHI_Ne235g = (np.log10(index['NeII_12m']+index['NeIII_15m'] + index['NeV_14m']) - Ne235g_obs)**2 / np.log10(index['NeII_12m'] + index['NeIII_15m'] + index['NeV_14m'])
                  if Ne235h_obs == -10:
                     CHI_Ne235h = 0
                  elif (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235h = tol_max
                  else:
                     CHI_Ne235h = (np.log10(index['NeII_12m'] + index['NeIII_15m'] + index['NeV_24m']) - Ne235h_obs)**2 / np.log10(index['NeII_12m'] + index['NeIII_15m'] + index['NeV_24m'])
                  if Ne23Ne5a_obs == -10: 
                     CHI_Ne23Ne5a = 0
                  elif (index['NeIII_15m'] == 0 and index['NeII_12m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne23Ne5a = tol_max
                  else:   
                     CHI_Ne23Ne5a = (np.log10(( ( index['NeIII_15m'] + index['NeII_12m'])/index['NeV_14m']))- Ne23Ne5a_obs)**2/np.log10(( (index['NeIII_15m']+index['NeII_12m'] )/index['NeV_14m']))
                  if Ne23Ne5b_obs == -10: 
                     CHI_Ne23Ne5b = 0
                  elif (index['NeIII_15m'] == 0 and index['NeII_12m'] == 0) or index['NeV_24m'] == 0:
                     CHI_Ne23Ne5b = tol_max
                  else:   
                     CHI_Ne23Ne5b = (np.log10(( ( index['NeIII_15m'] + index['NeII_12m'])/index['NeV_24m']))- Ne23Ne5b_obs)**2/np.log10(( (index['NeIII_15m']+index['NeII_12m'] )/index['NeV_24m']))
                  #Argon 
                  if Ar2Ar3_obs == -10: 
                     CHI_Ar2Ar3 = 0
                  elif index['ArII_7m'] == 0 or index['ArIII_9m'] == 0:
                     CHI_Ar2Ar3 = tol_max
                  else:   
                     CHI_Ar2Ar3 = (np.log10((index['ArII_7m']/index['ArIII_9m']))- Ar2Ar3_obs)**2/np.log10((index['ArII_7m']/index['ArIII_9m']))
                  if Ar235a_obs == -10: 
                     CHI_Ar235a = 0
                  elif index['HI_4m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235a = tol_max
                  else:   
                     CHI_Ar235a = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_4m'])- Ar235a_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_4m'])
                  if Ar235b_obs == -10: 
                     CHI_Ar235b = 0
                  elif index['HI_7m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235b = tol_max
                  else:   
                     CHI_Ar235b = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_7m'])- Ar235b_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_7m'])
                  if Ar235c_obs == -10: 
                     CHI_Ar235c = 0
                  elif index['HI_12m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235c = tol_max
                  else:   
                     CHI_Ar235c = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_12m'])- Ar235c_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_12m'])
                  if Ne235d_obs == -10: 
                     CHI_Ar235d = 0
                  elif index['HI_4m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235d = tol_max
                  else:   
                     CHI_Ar235d = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_4m'])- Ar235d_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_4m'])
                  if Ar235e_obs == -10: 
                     CHI_Ar235e = 0
                  elif index['HI_7m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235e = tol_max
                  else:   
                     CHI_Ar235e = (np.log10((index['ArII_7m']+index['NeIII_15m']+index['ArV_13m'])/index['HI_7m'])- Ne235e_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_7m'])
                  if Ar235f_obs == -10: 
                     CHI_Ar235f = 0
                  elif index['HI_12m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235f = tol_max
                  else:   
                     CHI_Ar235f = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_12m'])- Ar235f_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_12m'])
                  if Ar235g_obs == -10:
                     CHI_Ar235g = 0
                  elif (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235g = tol_max
                  else:
                     CHI_Ar235g = (np.log10(index['ArII_7m'] + index['ArIII_9m'] + index['ArV_8m']) - Ar235g_obs)**2/np.log10(index['ArII_7m'] + index['ArIII_9m'] + index['ArV_8m'])
                  if Ar235h_obs == -10:
                     CHI_Ar235h = 0
                  elif (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235h = tol_max
                  else:
                     CHI_Ar235h = (np.log10(index['ArII_7m'] + index['ArIII_9m'] + index['ArV_13m']) - Ar235h_obs)**2/np.log10(index['ArII_7m'] + index['ArIII_9m'] + index['ArV_13m'])
                  if Ar23Ar5a_obs == -10: 
                     CHI_Ar23Ar5a = 0
                  elif (index['ArIII_9m'] == 0 and index['ArII_7m'] == 0) or index['ArV_8m'] == 0:
                     CHI_Ar23Ar5a = tol_max
                  else:   
                     CHI_Ar23Ar5a = (np.log10(( ( index['ArIII_9m'] + index['ArII_7m'])/index['ArV_8m']))- Ar23Ar5a_obs)**2/np.log10(( (index['ArIII_9m']+index['ArII_7m'] )/index['ArV_8m']))
                  if Ar23Ar5b_obs == -10: 
                     CHI_Ar23Ar5b = 0
                  elif (index['ArIII_9m'] == 0 and index['ArII_7m'] == 0) or index['ArV_13m'] == 0:
                     CHI_Ar23Ar5b = tol_max
                  else:   
                     CHI_Ar23Ar5b = (np.log10(( ( index['ArIII_9m'] + index['ArII_7m'])/index['ArV_13m']))- Ar23Ar5b_obs)**2/np.log10(( (index['ArIII_9m']+index['ArII_7m'] )/index['ArV_13m']))
                  #Sulphur
                  if S34a_obs == -10: 
                     CHI_S34a = 0
                  elif index['HI_4m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0):
                     CHI_S34a = tol_max
                  else:   
                     CHI_S34a = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_4m'])- S34a_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_4m'])
                  if S34b_obs == -10: 
                     CHI_S34b = 0
                  elif index['HI_7m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0):
                     CHI_S34b = tol_max
                  else:   
                     CHI_S34b = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_7m'])- S34b_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_7m'])
                  if S34c_obs == -10: 
                     CHI_S34c = 0
                  elif index['HI_12m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0): 
                     CHI_S34c = tol_max
                  else:   
                     CHI_S34c = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_12m'])- S34c_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_12m'])
                  if S34d_obs == -10: 
                     CHI_S34d = 0
                  elif index['HI_4m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34d = tol_max
                  else:   
                     CHI_S34d = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_4m'])- S34d_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_4m'])
                  if S34e_obs == -10: 
                     CHI_S34e = 0
                  elif index['HI_7m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34e = tol_max
                  else:   
                     CHI_S34e = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_7m'])- S34e_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_7m'])
                  if S34f_obs == -10: 
                     CHI_S34f = 0
                  elif index['HI_12m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34f = tol_max
                  else:   
                     CHI_S34f = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_12m'])- S34f_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_12m'])
                  if S34g_obs == -10:
                     CHI_S34g = 0
                  elif (index['SIV_10m'] == 0 and index['SIII_19m'] == 0):
                     CHI_S34g = tol_max
                  else:
                     CHI_S34g = (np.log10(index['SIV_10m'] + index['SIII_19m']) - S34g_obs)**2/np.log10(index['SIV_10m'] + index['SIII_19m'])
                  if S34h_obs == -10:
                     CHI_S34h = 0
                  elif (index['SIV_10m'] == 0 and index['SIII_33m'] == 0):
                     CHI_S34h = tol_max
                  else:
                     CHI_S34h = (np.log10(index['SIV_10m'] + index['SIII_33m']) - S34h_obs)**2/np.log10(index['SIV_10m'] + index['SIII_33m']) 
                  if S3S4a_obs == -10: 
                     CHI_S3S4a = 0
                  elif index['SIV_10m'] == 0 or index['SIII_19m'] == 0 : 
                     CHI_S3S4a = tol_max
                  else:   
                     CHI_S3S4a = (np.log10(index['SIII_19m']/index['SIV_10m'])- S3S4a_obs)**2/np.log10((index['SIII_19m']/index['SIV_10m']))
                  if S3S4b_obs == -10: 
                     CHI_S3S4b = 0
                  elif index['SIV_10m'] == 0 or index['SIII_33m'] == 0 : 
                     CHI_S3S4b = tol_max
                  else:   
                     CHI_S3S4b = (np.log10(index['SIII_33m']/index['SIV_10m'])- S3S4b_obs)**2/np.log10((index['SIII_33m']/index['SIV_10m']))
                  #Oxygen
                  if O34a_obs == -10: 
                     CHI_O34a = 0
                  elif index['HI_4m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0):
                     CHI_O34a = tol_max
                  else:   
                     CHI_O34a = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_4m'])- O34a_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_4m'])
                  if O34b_obs == -10: 
                     CHI_O34b = 0
                  elif index['HI_7m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0):
                     CHI_O34b = tol_max
                  else:   
                     CHI_O34b = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_7m'])- O34b_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_7m'])
                  if O34c_obs == -10: 
                     CHI_O34c = 0
                  elif index['HI_12m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0): 
                     CHI_O34c = tol_max
                  else:   
                     CHI_O34c = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_12m'])- O34c_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_12m'])
                  if O34d_obs == -10: 
                     CHI_O34d = 0
                  elif index['HI_4m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34d = tol_max
                  else:   
                     CHI_O34d = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_4m'])- O34d_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_4m'])
                  if O34e_obs == -10: 
                     CHI_O34e = 0
                  elif index['HI_7m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34e = tol_max
                  else:   
                     CHI_O34e = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_7m'])- O34e_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_7m'])
                  if O34f_obs == -10: 
                     CHI_O34f = 0
                  elif index['HI_12m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34f = tol_max
                  else:   
                     CHI_O34f = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_12m'])- O34f_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_12m'])
                  if O34g_obs == -10:
                     CHI_O34g = 0
                  elif (index['OIV_26m'] == 0 and index['OIII_52m'] == 0):
                     CHI_O34g = tol_max
                  else:
                     CHI_O34g = (np.log10(index['OIV_26m'] + index['OIII_52m']) - O34g_obs)**2/np.log10(index['OIV_26m'] + index['OIII_52m'])
                  if O34h_obs == -10:
                     CHI_O34h = 0
                  elif (index['OIV_26m'] == 0 and index['OIII_88m'] == 0):
                     CHI_O34h = tol_max
                  else:
                     CHI_O34h = (np.log10(index['OIV_26m'] + index['OIII_88m']) - O34h_obs)**2/np.log10(index['OIV_26m'] + index['OIII_88m'])
                  if O3O4a_obs == -10: 
                     CHI_O3O4a = 0
                  elif index['OIV_26m'] == 0 or index['OIII_52m'] == 0 : 
                     CHI_O3O4a = tol_max
                  else:   
                     CHI_O3O4a = (np.log10(index['OIII_52m']/index['OIV_26m'])- O3O4a_obs)**2/np.log10((index['OIII_52m']/index['OIV_26m']))
                  if O3O4b_obs == -10: 
                     CHI_O3O4b = 0
                  elif index['OIV_26m'] == 0 or index['OIII_88m'] == 0 : 
                     CHI_O3O4b = tol_max
                  else:   
                     CHI_O3O4b = (np.log10(index['OIII_88m']/index['OIV_26m'])- O3O4b_obs)**2/np.log10((index['OIII_88m']/index['OIV_26m']))
                  #Nitrogen   
                  if N23a_obs == -10: 
                     CHI_N23a = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23a = tol_max
                  else:   
                     CHI_N23a = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_4m'])- N23a_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_4m'])
                  if N23b_obs == -10: 
                     CHI_N23b = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23b = tol_max
                  else:   
                     CHI_N23b = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_7m'])- N23b_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_7m'])
                  if N23c_obs == -10: 
                     CHI_N23c = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23c = tol_max
                  else:   
                     CHI_N23c = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_12m'])- N23c_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_12m'])
                  if N23g_obs == -10:
                     CHI_N23g = 0
                  elif (index['NIII_57m'] == 0 and index['NII_122m'] == 0):
                     CHI_N23g = tol_max
                  else:
                     CHI_N23g = (np.log10(index['NIII_57m'] + index['NII_122m']) - N23g_obs)**2/np.log10(index['NIII_57m'] + index['NII_122m'])         
                  if N3a_obs == -10: 
                     CHI_N3a = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3a = tol_max
                  else:   
                     CHI_N3a = (np.log10((index['NIII_57m'])/index['HI_4m'])- N3a_obs)**2/np.log10((index['NIII_57m'])/index['HI_4m'])
                  if N3b_obs == -10: 
                     CHI_N3b = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3b = tol_max
                  else:   
                     CHI_N3b = (np.log10((index['NIII_57m'])/index['HI_7m'])- N3b_obs)**2/np.log10((index['NIII_57m'])/index['HI_7m'])
                  if N3c_obs == -10: 
                     CHI_N3c = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3c = tol_max
                  else:   
                     CHI_N3c = (np.log10((index['NIII_57m'])/index['HI_12m'])- N3c_obs)**2/np.log10((index['NIII_57m'])/index['HI_12m'])
                  if N3d_obs == -10:
                     CHI_N3d = 0
                  elif index['NIII_57m'] == 0:
                     CHI_N3d = tol_max
                  else:
                     CHI_N3d = (np.log10(index['NIII_57m']) - N3d_obs)**2/np.log10(index['NIII_57m'])          
                  if N2N3a_obs == -10: 
                     CHI_N2N3a = 0
                  elif index['NIII_57m'] == 0 or index['NII_122m'] == 0: 
                     CHI_N2N3a = tol_max
                  else:   
                     CHI_N2N3a = (np.log10(index['NII_122m']/index['NIII_57m'])- N2N3a_obs)**2/np.log10((index['NII_122m']/index['NIII_57m']))
                  if N23d_obs == -10: 
                     CHI_N23d = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23d = tol_max
                  else:   
                     CHI_N23d = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_4m'])- N23d_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_4m'])
                  if N23e_obs == -10: 
                     CHI_N23e = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23e = tol_max
                  else:   
                     CHI_N23e = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_7m'])- N23e_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_7m'])
                  if N23f_obs == -10: 
                     CHI_N23f = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23f = tol_max
                  else:   
                     CHI_N23f = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_12m'])- N23f_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_12m'])
                  if N23h_obs == -10:
                     CHI_N23h = 0
                  elif (index['NIII_57m'] == 0 and index['NII_205m'] == 0):
                     CHI_N23h = tol_max
                  else:
                     CHI_N23h = (np.log10(index['NIII_57m']+index['NII_205m']) - N23h_obs)**2/np.log10(index['NIII_57m']+index['NII_205m']) 
                  if N2N3b_obs == -10: 
                     CHI_N2N3b = 0
                  elif index['NIII_57m'] == 0 or index['NII_205m'] == 0: 
                     CHI_N2N3b = tol_max
                  else:   
                     CHI_N2N3b = (np.log10(index['NII_205m']/index['NIII_57m'])- N2N3b_obs)**2/np.log10((index['NII_205m']/index['NIII_57m']))
                  #O3N2
                  if O3N2a_obs == -10: 
                     CHI_O3N2a = 0
                  elif index['OIII_52m'] == 0 or index['NII_122m'] == 0: 
                     CHI_O3N2a = tol_max
                  else:   
                     CHI_O3N2a = (np.log10(index['OIII_52m']/index['NII_122m'])- O3N2a_obs)**2/np.log10((index['OIII_52m']/index['NII_122m']))
                  if O3N2b_obs == -10: 
                     CHI_O3N2b = 0
                  elif index['OIII_88m'] == 0 or index['NII_122m'] == 0: 
                     CHI_O3N2b = tol_max
                  else:   
                     CHI_O3N2b = (np.log10(index['OIII_88m']/index['NII_122m'])- O3N2b_obs)**2/np.log10((index['OIII_88m']/index['NII_122m']))
                  if O3N2c_obs == -10: 
                     CHI_O3N2c = 0
                  elif index['OIII_52m'] == 0 or index['NII_205m'] == 0: 
                     CHI_O3N2c = tol_max
                  else:   
                     CHI_O3N2c = (np.log10(index['OIII_52m']/index['NII_205m'])- O3N2c_obs)**2/np.log10((index['OIII_52m']/index['NII_205m']))
                  if O3N2d_obs == -10: 
                     CHI_O3N2d = 0
                  elif index['OIII_88m'] == 0 or index['NII_205m'] == 0: 
                     CHI_O3N2d = tol_max
                  else:   
                     CHI_O3N2d = (np.log10(index['OIII_88m']/index['NII_205m'])- O3N2d_obs)**2/np.log10((index['OIII_88m']/index['NII_205m']))

                  CHI_OH = (CHI_Ne2Ne3**2 + CHI_S3S4a**2 + CHI_S3S4b**2 +CHI_N23a**2+ CHI_N23b**2 +CHI_N23c**2 + CHI_N2N3a**2 + CHI_O3N2a**2 + CHI_O3N2b**2 + CHI_N23d**2 + CHI_N23e**2 + CHI_N23f**2 + CHI_N2N3b**2 + CHI_O3N2c**2 + CHI_O3N2d**2 + CHI_O34a**2 + CHI_O34b**2 + CHI_O34c**2 + CHI_O34d**2 + CHI_O34e**2 + CHI_O34f**2 + CHI_O3O4a**2 + CHI_O3O4b**2 + CHI_Ne235a**2 + CHI_Ne235b**2 + CHI_Ne235c**2 + CHI_Ne235d**2 + CHI_Ne235e**2 + CHI_Ne235f**2 + CHI_Ne23Ne5a**2 + CHI_Ne23Ne5b**2 + CHI_N3a**2 + CHI_N3b**2 + CHI_N3c**2 + CHI_Ar2Ar3**2 + CHI_Ar235a**2 + CHI_Ar235b**2 + CHI_Ar235c**2 + CHI_Ar235d**2 + CHI_Ar235e**2 + CHI_Ar235f**2 + CHI_Ar23Ar5a**2 +CHI_Ar23Ar5b**2 + CHI_Ne235g**2 + CHI_Ne235h**2 + CHI_Ar235g**2 + CHI_Ar235h**2 +  CHI_O34g**2 + CHI_O34h**2 + CHI_N23g**2 + CHI_N23h**2 + CHI_N3d**2)**0.5
                  CHI_SH = (CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S34d**2 + CHI_S34e**2 + CHI_S34f**2 +CHI_S3S4a**2 + CHI_S3S4b**2 + CHI_S34g**2 + CHI_S34h**2)**0.5

                  if CHI_OH == 0:
                     OH_p = OH_p
                     logU_p = logU_p
                     den_OH = den_OH
                  else:
                     OH_p = index['12logOH'] / (CHI_OH) + OH_p
                     logU_p = index['logU'] / (CHI_OH) + logU_p
                     den_OH = 1 / (CHI_OH) + den_OH

                  if CHI_SH == 0:
                     SH_p = SH_p
                     den_SH = den_SH
                  else:
                     SH_p = index['12logOH'] / (CHI_SH) + SH_p
                     den_SH = 1/(CHI_SH) + den_SH
   

            if sed >= 3 and HI_4m[tab] == 0 and HI_7m[tab] == 0 and HI_12m[tab] == 0 and Hbeta[tab] == 0:
               OH = 0
            elif OH_p == 0:
               OH = 0
            else:
               OH = OH_p / den_OH
            if sed >= 3 and HI_4m[tab] == 0 and HI_7m[tab] == 0 and HI_12m[tab] == 0 and Hbeta[tab] == 0:
               SH = 0
            elif SH_p == 0:
               SH = 0
            else:
               SH = SH_p / den_SH - 1.57
            if logU_p == 0:
               logU = 0
            else:
               logU = logU_p / den_OH


         #Calculation of error of O/H, S/H and logU
         if S34a_obs == -10 and S34b_obs == -10 and S34c_obs == -10 and S34d_obs == -10 and S34e_obs == -10 and S34f_obs == -10 and Ne2Ne3_obs == -10 and S3S4a_obs == -10 and S3S4b_obs == -10 and N23a_obs == -10 and N23b_obs == -10 and N23c_obs == -10 and N2N3a_obs == -10 and O3N2a_obs == -10 and O3N2b_obs == -10 and N23d_obs == -10 and N23e_obs == -10 and N23f_obs == -10 and N2N3b_obs == -10 and O3N2c_obs == -10 and O3N2d_obs == -10 and Ne235a_obs == -10 and Ne235b_obs == -10 and Ne235c_obs == -10 and Ne235d_obs == -10 and Ne235e_obs == -10 and Ne235f_obs == -10 and Ne23Ne5a_obs == -10 and Ne23Ne5b_obs == -10 and O34a_obs == -10 and O34b_obs == -10 and O34c_obs == -10 and O34d_obs == -10 and O34e_obs == -10 and O34f_obs == -10 and O3O4a_obs == -10 and O3O4b_obs == -10 and N3a_obs == -10 and N3b_obs == -10 and N3c_obs == -10 and Ar2Ar3_obs == -10 and Ar235a_obs == -10 and Ar235b_obs == -10 and Ar235c_obs == -10 and Ar235d_obs == -10 and Ar235e_obs == -10 and Ar235f_obs == -10 and Ar23Ar5a_obs == -10 and Ar23Ar5b_obs == -10 and Ne235g_obs == -10 and Ne235h_obs == -10 and Ar235g_obs == -10 and Ar235h_obs == -10 and S34g_obs == -10 and S34h_obs == -10 and O34g_obs == -10 and O34h_obs == -10 and N23g_obs == -10 and N23h_obs == -10 and N3d_obs == -10: 
            eOH = 0
            elogU = 0
            eSH = 0
         else:
            CHI_Ne2Ne3 = 0
            CHI_Ne235a = 0
            CHI_Ne235b = 0
            CHI_Ne235c = 0
            CHI_Ne235d = 0
            CHI_Ne235e = 0
            CHI_Ne235f = 0
            CHI_Ne235g = 0
            CHI_Ne235h = 0
            CHI_Ne23Ne5a = 0
            CHI_Ne23Ne5b = 0
               
            CHI_Ar2Ar3 = 0
            CHI_Ar235a = 0
            CHI_Ar235b = 0
            CHI_Ar235c = 0
            CHI_Ar235d = 0
            CHI_Ar235e = 0
            CHI_Ar235f = 0
            CHI_Ar235g = 0
            CHI_Ar235h = 0
            CHI_Ar23Ar5a = 0
            CHI_Ar23Ar5b = 0

            CHI_O34a = 0
            CHI_O34b = 0
            CHI_O34c = 0
            CHI_O34d = 0
            CHI_O34e = 0
            CHI_O34f = 0
            CHI_O34g = 0
            CHI_O34h = 0
            CHI_O3O4a = 0
            CHI_O3O4b = 0

            CHI_S34a = 0
            CHI_S34b = 0
            CHI_S34c = 0
            CHI_S34d = 0
            CHI_S34e = 0
            CHI_S34f = 0
            CHI_S34g = 0
            CHI_S34h = 0
            CHI_S3S4a = 0
            CHI_S3S4b = 0

            CHI_N23a = 0
            CHI_N23b = 0
            CHI_N23c = 0
            CHI_N23d = 0
            CHI_N23e = 0
            CHI_N23f = 0
            CHI_N23g = 0
            CHI_N23h = 0
            CHI_N3a = 0
            CHI_N3b = 0
            CHI_N3c = 0
            CHI_N3d = 0
            CHI_N2N3a = 0
            CHI_N2N3b = 0

            CHI_O3N2a = 0
            CHI_O3N2b = 0
            CHI_O3N2c = 0
            CHI_O3N2d = 0

            CHI_OH = 0
            CHI_SH = 0
            #Constraints in N/O and high-ionization emission lines
            for index in grid_c:
               if NOff > -10 and np.abs(index['logNO'] - NOff) > np.abs(eNOff+res_NO):
                  continue
               elif SIV_10m_obs > 0 and index['SIV_10m'] == 0:
                  continue
               elif OIV_26m_obs > 0 and index['OIV_26m'] == 0:
                  continue
               elif ArV_8m_obs > 0 and index['ArV_8m'] == 0:
                  continue
               elif ArV_13m_obs > 0 and index['ArV_13m'] == 0:
                  continue
               elif NeV_14m_obs > 0 and index['NeV_14m'] == 0:
                  continue
               elif NeV_24m_obs > 0 and index['NeV_24m'] == 0:
                  continue
               else:            
                  #Neon     
                  if Ne2Ne3_obs == -10: 
                     CHI_Ne2Ne3 = 0
                  elif index['NeII_12m'] == 0 or index['NeIII_15m'] == 0:
                     CHI_Ne2Ne3 = tol_max
                  else:   
                     CHI_Ne2Ne3 = (np.log10((index['NeII_12m']/index['NeIII_15m']))- Ne2Ne3_obs)**2/np.log10((index['NeII_12m']/index['NeIII_15m']))
                  if Ne235a_obs == -10: 
                     CHI_Ne235a = 0
                  elif index['HI_4m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235a = tol_max
                  else:   
                     CHI_Ne235a = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_4m'])- Ne235a_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_4m'])
                  if Ne235b_obs == -10: 
                     CHI_Ne235b = 0
                  elif index['HI_7m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235b = tol_max
                  else:   
                     CHI_Ne235b = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_7m'])- Ne235b_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_7m'])
                  if Ne235c_obs == -10: 
                     CHI_Ne235c = 0
                  elif index['HI_12m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235c = tol_max
                  else:   
                     CHI_Ne235c = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_12m'])- Ne235c_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_12m'])
                  if Ne235d_obs == -10: 
                     CHI_Ne235d = 0
                  elif index['HI_4m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235d = tol_max
                  else:   
                     CHI_Ne235d = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_4m'])- Ne235d_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_4m'])
                  if Ne235e_obs == -10: 
                     CHI_Ne235e = 0
                  elif index['HI_7m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235e = tol_max
                  else:   
                     CHI_Ne235e = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_7m'])- Ne235e_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_7m'])
                  if Ne235f_obs == -10: 
                     CHI_Ne235f = 0
                  elif index['HI_12m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235f = tol_max
                  else:   
                     CHI_Ne235f = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_12m'])- Ne235f_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_12m']) 
                  if Ne23Ne5a_obs == -10: 
                     CHI_Ne23Ne5a = 0
                  elif (index['NeIII_15m'] == 0 and index['NeII_12m'] == 0) or index['NeV_14m'] == 0:
                     CHI_Ne23Ne5a = tol_max
                  else:   
                     CHI_Ne23Ne5a = (np.log10(( ( index['NeIII_15m'] + index['NeII_12m'])/index['NeV_14m']))- Ne23Ne5a_obs)**2/np.log10(( (index['NeIII_15m']+index['NeII_12m'] )/index['NeV_14m']))
                  if Ne23Ne5b_obs == -10: 
                     CHI_Ne23Ne5b = 0
                  elif (index['NeIII_15m'] == 0 and index['NeII_12m'] == 0) or index['NeV_24m'] == 0:
                     CHI_Ne23Ne5b = tol_max
                  else:   
                     CHI_Ne23Ne5b = (np.log10(( ( index['NeIII_15m'] + index['NeII_12m'])/index['NeV_24m']))- Ne23Ne5b_obs)**2/np.log10(( (index['NeIII_15m']+index['NeII_12m'] )/index['NeV_24m']))
                  #Argon 
                  if Ar2Ar3_obs == -10: 
                     CHI_Ar2Ar3 = 0
                  elif index['ArII_7m'] == 0 or index['ArIII_9m'] == 0:
                     CHI_Ar2Ar3 = tol_max
                  else:   
                     CHI_Ar2Ar3 = (np.log10((index['ArII_7m']/index['ArIII_9m']))- Ar2Ar3_obs)**2/np.log10((index['ArII_7m']/index['ArIII_9m']))
                  if Ar235a_obs == -10: 
                     CHI_Ar235a = 0
                  elif index['HI_4m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235a = tol_max
                  else:   
                     CHI_Ar235a = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_4m'])- Ar235a_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_4m'])
                  if Ar235b_obs == -10: 
                     CHI_Ar235b = 0
                  elif index['HI_7m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235b = tol_max
                  else:   
                     CHI_Ar235b = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_7m'])- Ar235b_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_7m'])
                  if Ar235c_obs == -10: 
                     CHI_Ar235c = 0
                  elif index['HI_12m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235c = tol_max
                  else:   
                     CHI_Ar235c = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_12m'])- Ar235c_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_12m'])
                  if Ne235d_obs == -10: 
                     CHI_Ar235d = 0
                  elif index['HI_4m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235d = tol_max
                  else:   
                     CHI_Ar235d = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_4m'])- Ar235d_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_4m'])
                  if Ar235e_obs == -10: 
                     CHI_Ar235e = 0
                  elif index['HI_7m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235e = tol_max
                  else:   
                     CHI_Ar235e = (np.log10((index['ArII_7m']+index['NeIII_15m']+index['ArV_13m'])/index['HI_7m'])- Ne235e_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_7m'])
                  if Ar235f_obs == -10: 
                     CHI_Ar235f = 0
                  elif index['HI_12m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235f = tol_max
                  else:   
                     CHI_Ar235f = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_12m'])- Ar235f_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_12m']) 
                  if Ar23Ar5a_obs == -10: 
                     CHI_Ar23Ar5a = 0
                  elif (index['ArIII_9m'] == 0 and index['ArII_7m'] == 0) or index['ArV_8m'] == 0:
                     CHI_Ar23Ar5a = tol_max
                  else:   
                     CHI_Ar23Ar5a = (np.log10(( ( index['ArIII_9m'] + index['ArII_7m'])/index['ArV_8m']))- Ar23Ar5a_obs)**2/np.log10(( (index['ArIII_9m']+index['ArII_7m'] )/index['ArV_8m']))
                  if Ar23Ar5b_obs == -10: 
                     CHI_Ar23Ar5b = 0
                  elif (index['ArIII_9m'] == 0 and index['ArII_7m'] == 0) or index['ArV_13m'] == 0:
                     CHI_Ar23Ar5b = tol_max
                  else:   
                     CHI_Ar23Ar5b = (np.log10(( ( index['ArIII_9m'] + index['ArII_7m'])/index['ArV_13m']))- Ar23Ar5b_obs)**2/np.log10(( (index['ArIII_9m']+index['ArII_7m'] )/index['ArV_13m']))
                  #Sulphur
                  if S34a_obs == -10: 
                     CHI_S34a = 0
                  elif index['HI_4m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0):
                     CHI_S34a = tol_max
                  else:   
                     CHI_S34a = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_4m'])- S34a_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_4m'])
                  if S34b_obs == -10: 
                     CHI_S34b = 0
                  elif index['HI_7m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0):
                     CHI_S34b = tol_max
                  else:   
                     CHI_S34b = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_7m'])- S34b_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_7m'])
                  if S34c_obs == -10: 
                     CHI_S34c = 0
                  elif index['HI_12m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0): 
                     CHI_S34c = tol_max
                  else:   
                     CHI_S34c = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_12m'])- S34c_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_12m'])
                  if S34d_obs == -10: 
                     CHI_S34d = 0
                  elif index['HI_4m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34d = tol_max
                  else:   
                     CHI_S34d = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_4m'])- S34d_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_4m'])
                  if S34e_obs == -10: 
                     CHI_S34e = 0
                  elif index['HI_7m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34e = tol_max
                  else:   
                     CHI_S34e = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_7m'])- S34e_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_7m'])
                  if S34f_obs == -10: 
                     CHI_S34f = 0
                  elif index['HI_12m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34f = tol_max
                  else:   
                     CHI_S34f = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_12m'])- S34f_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_12m'])
                  if S3S4a_obs == -10: 
                     CHI_S3S4a = 0
                  elif index['SIV_10m'] == 0 or index['SIII_19m'] == 0 : 
                     CHI_S3S4a = tol_max
                  else:   
                     CHI_S3S4a = (np.log10(index['SIII_19m']/index['SIV_10m'])- S3S4a_obs)**2/np.log10((index['SIII_19m']/index['SIV_10m']))
                  if S3S4b_obs == -10: 
                     CHI_S3S4b = 0
                  elif index['SIV_10m'] == 0 or index['SIII_33m'] == 0 : 
                     CHI_S3S4b = tol_max
                  else:   
                     CHI_S3S4b = (np.log10(index['SIII_33m']/index['SIV_10m'])- S3S4b_obs)**2/np.log10((index['SIII_33m']/index['SIV_10m']))
                  #Oxygen
                  if O34a_obs == -10: 
                     CHI_O34a = 0
                  elif index['HI_4m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0):
                     CHI_O34a = tol_max
                  else:   
                     CHI_O34a = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_4m'])- O34a_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_4m'])
                  if O34b_obs == -10: 
                     CHI_O34b = 0
                  elif index['HI_7m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0):
                     CHI_O34b = tol_max
                  else:   
                     CHI_O34b = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_7m'])- O34b_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_7m'])
                  if O34c_obs == -10: 
                     CHI_O34c = 0
                  elif index['HI_12m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0): 
                     CHI_O34c = tol_max
                  else:   
                     CHI_O34c = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_12m'])- O34c_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_12m'])
                  if O34d_obs == -10: 
                     CHI_O34d = 0
                  elif index['HI_4m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34d = tol_max
                  else:   
                     CHI_O34d = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_4m'])- O34d_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_4m'])
                  if O34e_obs == -10: 
                     CHI_O34e = 0
                  elif index['HI_7m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34e = tol_max
                  else:   
                     CHI_O34e = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_7m'])- O34e_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_7m'])
                  if O34f_obs == -10: 
                     CHI_O34f = 0
                  elif index['HI_12m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34f = tol_max
                  else:   
                     CHI_O34f = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_12m'])- O34f_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_12m'])
                  if O3O4a_obs == -10: 
                     CHI_O3O4a = 0
                  elif index['OIV_26m'] == 0 or index['OIII_52m'] == 0 : 
                     CHI_O3O4a = tol_max
                  else:   
                     CHI_O3O4a = (np.log10(index['OIII_52m']/index['OIV_26m'])- O3O4a_obs)**2/np.log10((index['OIII_52m']/index['OIV_26m']))
                  if O3O4b_obs == -10: 
                     CHI_O3O4b = 0
                  elif index['OIV_26m'] == 0 or index['OIII_88m'] == 0 : 
                     CHI_O3O4b = tol_max
                  else:   
                     CHI_O3O4b = (np.log10(index['OIII_88m']/index['OIV_26m'])- O3O4b_obs)**2/np.log10((index['OIII_88m']/index['OIV_26m']))
                  #Nitrogen   
                  if N23a_obs == -10: 
                     CHI_N23a = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23a = tol_max
                  else:   
                     CHI_N23a = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_4m'])- N23a_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_4m'])
                  if N23b_obs == -10: 
                     CHI_N23b = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23b = tol_max
                  else:   
                     CHI_N23b = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_7m'])- N23b_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_7m'])
                  if N23c_obs == -10: 
                     CHI_N23c = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23c = tol_max
                  else:   
                     CHI_N23c = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_12m'])- N23c_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_12m'])          
                  if N3a_obs == -10: 
                     CHI_N3a = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3a = tol_max
                  else:   
                     CHI_N3a = (np.log10((index['NIII_57m'])/index['HI_4m'])- N3a_obs)**2/np.log10((index['NIII_57m'])/index['HI_4m'])
                  if N3b_obs == -10: 
                     CHI_N3b = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3b = tol_max
                  else:   
                     CHI_N3b = (np.log10((index['NIII_57m'])/index['HI_7m'])- N3b_obs)**2/np.log10((index['NIII_57m'])/index['HI_7m'])
                  if N3c_obs == -10: 
                     CHI_N3c = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3c = tol_max
                  else:   
                     CHI_N3c = (np.log10((index['NIII_57m'])/index['HI_12m'])- N3c_obs)**2/np.log10((index['NIII_57m'])/index['HI_12m'])           
                  if N2N3a_obs == -10: 
                     CHI_N2N3a = 0
                  elif index['NIII_57m'] == 0 or index['NII_122m'] == 0: 
                     CHI_N2N3a = tol_max
                  else:   
                     CHI_N2N3a = (np.log10(index['NII_122m']/index['NIII_57m'])- N2N3a_obs)**2/np.log10((index['NII_122m']/index['NIII_57m']))
                  if N23d_obs == -10: 
                     CHI_N23d = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23d = tol_max
                  else:   
                     CHI_N23d = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_4m'])- N23d_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_4m'])
                  if N23e_obs == -10: 
                     CHI_N23e = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23e = tol_max
                  else:   
                     CHI_N23e = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_7m'])- N23e_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_7m'])
                  if N23f_obs == -10: 
                     CHI_N23f = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23f = tol_max
                  else:   
                     CHI_N23f = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_12m'])- N23f_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_12m'])
                  if N2N3b_obs == -10: 
                     CHI_N2N3b = 0
                  elif index['NIII_57m'] == 0 or index['NII_205m'] == 0: 
                     CHI_N2N3b = tol_max
                  else:   
                     CHI_N2N3b = (np.log10(index['NII_205m']/index['NIII_57m'])- N2N3b_obs)**2/np.log10((index['NII_205m']/index['NIII_57m']))
                  #O3N2
                  if O3N2a_obs == -10: 
                     CHI_O3N2a = 0
                  elif index['OIII_52m'] == 0 or index['NII_122m'] == 0: 
                     CHI_O3N2a = tol_max
                  else:   
                     CHI_O3N2a = (np.log10(index['OIII_52m']/index['NII_122m'])- O3N2a_obs)**2/np.log10((index['OIII_52m']/index['NII_122m']))
                  if O3N2b_obs == -10: 
                     CHI_O3N2b = 0
                  elif index['OIII_88m'] == 0 or index['NII_122m'] == 0: 
                     CHI_O3N2b = tol_max
                  else:   
                     CHI_O3N2b = (np.log10(index['OIII_88m']/index['NII_122m'])- O3N2b_obs)**2/np.log10((index['OIII_88m']/index['NII_122m']))
                  if O3N2c_obs == -10: 
                     CHI_O3N2c = 0
                  elif index['OIII_52m'] == 0 or index['NII_205m'] == 0: 
                     CHI_O3N2c = tol_max
                  else:   
                     CHI_O3N2c = (np.log10(index['OIII_52m']/index['NII_205m'])- O3N2c_obs)**2/np.log10((index['OIII_52m']/index['NII_205m']))
                  if O3N2d_obs == -10: 
                     CHI_O3N2d = 0
                  elif index['OIII_88m'] == 0 or index['NII_205m'] == 0: 
                     CHI_O3N2d = tol_max
                  else:   
                     CHI_O3N2d = (np.log10(index['OIII_88m']/index['NII_205m'])- O3N2d_obs)**2/np.log10((index['OIII_88m']/index['NII_205m']))

                  CHI_OH = (CHI_Ne2Ne3**2 + CHI_S3S4a**2 + CHI_S3S4b**2 +CHI_N23a**2+ CHI_N23b**2 +CHI_N23c**2 + CHI_N2N3a**2 + CHI_O3N2a**2 + CHI_O3N2b**2 + CHI_N23d**2 + CHI_N23e**2 + CHI_N23f**2 + CHI_N2N3b**2 + CHI_O3N2c**2 + CHI_O3N2d**2 + CHI_O34a**2 + CHI_O34b**2 + CHI_O34c**2 + CHI_O34d**2 + CHI_O34e**2 + CHI_O34f**2 + CHI_O3O4a**2 + CHI_O3O4b**2 + CHI_Ne235a**2 + CHI_Ne235b**2 + CHI_Ne235c**2 + CHI_Ne235d**2 + CHI_Ne235e**2 + CHI_Ne235f**2 + CHI_Ne23Ne5a**2 + CHI_Ne23Ne5b**2 + CHI_N3a**2 + CHI_N3b**2 + CHI_N3c**2 + CHI_Ar2Ar3**2 + CHI_Ar235a**2 + CHI_Ar235b**2 + CHI_Ar235c**2 + CHI_Ar235d**2 + CHI_Ar235e**2 + CHI_Ar235f**2 + CHI_Ar23Ar5a**2 +CHI_Ar23Ar5b**2 + CHI_Ne235g**2 + CHI_Ne235h**2 + CHI_Ar235g**2 + CHI_Ar235h**2 + CHI_O34g**2 + CHI_O34h**2 + CHI_N3d**2 + CHI_N23g**2 + CHI_N23h**2)**0.5

                  CHI_SH = (CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S34d**2 + CHI_S34e**2 + CHI_S34f**2 +CHI_S3S4a**2 + CHI_S3S4b**2 + CHI_S34g**2 + CHI_S34h**2)**0.5

                  if CHI_OH == 0:
                     OH_e = OH_e
                     logU_e = logU_e
                     den_OH_e = den_OH_e
                  else:
                     OH_e = (index['12logOH'] - OH)**2 / (CHI_OH) + OH_e
                     logU_e = (index['logU'] - logU)**2 / (CHI_OH) + logU_e
                     den_OH_e = 1 / (CHI_OH) + den_OH_e
                  if CHI_SH == 0:
                     SH_e = SH_e
                     den_SH_e = den_SH_e
                  else:
                     SH_e = (index['12logOH'] - 1.57 - SH)**2 / (CHI_SH) + SH_e
                     den_SH_e = 1 / (CHI_SH) + den_SH_e

            if sed >= 3 and HI_4m[tab] == 0 and HI_7m[tab] == 0 and HI_12m[tab] == 0 and Hbeta[tab] == 0:
               eOH = 0
            elif OH_e == 0:
               eOH = 0
            else:
               eOH = OH_e / den_OH_e
            if OH_e == 0:
               elogU = 0
            else:
               elogU = logU_e / den_OH_e
            if sed >= 3 and HI_4m[tab] == 0 and HI_7m[tab] == 0 and HI_12m[tab] == 0 and Hbeta[tab] == 0:
               eSH = 0
            elif SH_e == 0:
               eSH = 0
            else:
               eSH = SH_e / den_SH_e

         #Iterations for interpolated models
         if inter == 0 or OH == 0:
            OHf = OH
            logUf = logU
            SHf = SH
         elif inter == 1:
            igrid = interpolate(grid_c,2,logU-elogU-0.25,logU+elogU+0.25,10)
            igrid = igrid[np.lexsort((igrid['12logOH'],igrid['logU']))]
            igrid = interpolate(igrid,0,OH-eOH-0.1,OH+eOH+0.1,10)
            igrid = igrid[np.lexsort((igrid['12logOH'],igrid['logU']))]

            CHI_Ne235a = 0
            CHI_Ne235b = 0
            CHI_Ne235c = 0
            CHI_Ne235d = 0
            CHI_Ne235e = 0
            CHI_Ne235f = 0 
            CHI_Ne2Ne3 = 0
            CHI_Ne23Ne5a = 0
            CHI_Ne23Ne5b = 0

            CHI_Ar235a = 0
            CHI_Ar235b = 0
            CHI_Ar235c = 0
            CHI_Ar235d = 0
            CHI_Ar235e = 0
            CHI_Ar235f = 0 
            CHI_Ar2Ar3 = 0
            CHI_Ar23Ar5a = 0
            CHI_Ar23Ar5b = 0
         
            CHI_S34a = 0
            CHI_S34b = 0
            CHI_S34x = 0
            CHI_S34d = 0
            CHI_S34e = 0
            CHI_S34f = 0
            CHI_S3S4a = 0
            CHI_S3S4b = 0

            CHI_O34a = 0
            CHI_O34b = 0
            CHI_O34c = 0
            CHI_O34d = 0
            CHI_O34e = 0
            CHI_O34f = 0
            CHI_O3O4a = 0
            CHI_O3O4b = 0
      
            CHI_N23a = 0
            CHI_N23b = 0
            CHI_N23c = 0
            CHI_N3a = 0
            CHI_N3b = 0
            CHI_N3c = 0   
            CHI_N2N3a = 0
            CHI_N23d = 0
            CHI_N23e = 0
            CHI_N23f = 0
            CHI_N2N3b = 0

            CHI_O3N2a = 0
            CHI_O3N2b = 0
            CHI_O3N2c = 0
            CHI_O3N2d = 0

            CHI_OH = 0
            CHI_SH = 0

            OH_p = 0
            logU_p = 0
            SH_p = 0
            den_OH = 0
            den_SH = 0
         
            for index in igrid:
               if NOff > -10 and np.abs(index['logNO'] - NOff) > np.abs(eNOff+res_NO):
                  continue
               elif SIV_10m_obs > 0 and index['SIV_10m'] == 0:
                  continue
               elif OIV_26m_obs > 0 and index['OIV_26m'] == 0:
                  continue
               elif ArV_8m_obs > 0 and index['ArV_8m'] == 0:
                  continue
               elif ArV_13m_obs > 0 and index['ArV_13m'] == 0:
                  continue
               elif NeV_14m_obs > 0 and index['NeV_14m'] == 0:
                  continue
               elif NeV_24m_obs > 0 and index['NeV_24m'] == 0:
                  continue
               else:            
                  #Neon     
                  if Ne2Ne3_obs == -10: 
                     CHI_Ne2Ne3 = 0
                  elif index['NeII_12m'] == 0 or index['NeIII_15m'] == 0:
                     CHI_Ne2Ne3 = tol_max
                  else:   
                     CHI_Ne2Ne3 = (np.log10((index['NeII_12m']/index['NeIII_15m']))- Ne2Ne3_obs)**2/np.log10((index['NeII_12m']/index['NeIII_15m']))
                  if Ne235a_obs == -10: 
                     CHI_Ne235a = 0
                  elif index['HI_4m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235a = tol_max
                  else:   
                     CHI_Ne235a = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_4m'])- Ne235a_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_4m'])
                  if Ne235b_obs == -10: 
                     CHI_Ne235b = 0
                  elif index['HI_7m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235b = tol_max
                  else:   
                     CHI_Ne235b = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_7m'])- Ne235b_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_7m'])
                  if Ne235c_obs == -10: 
                     CHI_Ne235c = 0
                  elif index['HI_12m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235c = tol_max
                  else:   
                     CHI_Ne235c = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_12m'])- Ne235c_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_14m'])/index['HI_12m'])
                  if Ne235d_obs == -10: 
                     CHI_Ne235d = 0
                  elif index['HI_4m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235d = tol_max
                  else:   
                     CHI_Ne235d = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_4m'])- Ne235d_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_4m'])
                  if Ne235e_obs == -10: 
                     CHI_Ne235e = 0
                  elif index['HI_7m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235e = tol_max
                  else:   
                     CHI_Ne235e = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_7m'])- Ne235e_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_7m'])
                  if Ne235f_obs == -10: 
                     CHI_Ne235f = 0
                  elif index['HI_12m'] == 0 or (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235f = tol_max
                  else:   
                     CHI_Ne235f = (np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_12m'])- Ne235f_obs)**2/np.log10((index['NeII_12m']+index['NeIII_15m']+index['NeV_24m'])/index['HI_12m'])
                  if Ne235g_obs == -10:
                     CHI_Ne235g = 0
                  elif (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_14m'] == 0):
                     CHI_Ne235g = tol_max
                  else:
                     CHI_Ne235g = (np.log10(index['NeII_12m'] + index['NeIII_15m'] + index['NeV_14m']) - Ne235g_obs)**2/np.log10(index['NeII_12m'] + index['NeIII_15m'] + index['NeV_14m'])
                  if Ne235h_obs == -10:
                     CHI_Ne235h = 0
                  elif (index['NeII_12m'] == 0 and index['NeIII_15m'] == 0 and index['NeV_24m'] == 0):
                     CHI_Ne235h = tol_max
                  else:
                     CHI_Ne235h = (np.log10(index['NII_12m'] + index['NeIII_15m'] + index['NeV_24m']) - Ne235h_obs)**2/np.log10(index['NeII_12m'] + index['NeIII_15m'] + index['NeV_24m'])
                  if Ne23Ne5a_obs == -10: 
                     CHI_Ne23Ne5a = 0
                  elif (index['NeIII_15m'] == 0 and index['NeII_12m'] == 0) or index['NeV_14m'] == 0:
                     CHI_Ne23Ne5a = tol_max
                  else:   
                     CHI_Ne23Ne5a = (np.log10(( ( index['NeIII_15m'] + index['NeII_12m'])/index['NeV_14m']))- Ne23Ne5a_obs)**2/np.log10(( (index['NeIII_15m']+index['NeII_12m'] )/index['NeV_14m']))
                  if Ne23Ne5b_obs == -10: 
                     CHI_Ne23Ne5b = 0
                  elif (index['NeIII_15m'] == 0 and index['NeII_12m'] == 0) or index['NeV_24m'] == 0:
                     CHI_Ne23Ne5b = tol_max
                  else:   
                     CHI_Ne23Ne5b = (np.log10(( ( index['NeIII_15m'] + index['NeII_12m'])/index['NeV_24m']))- Ne23Ne5b_obs)**2/np.log10(( (index['NeIII_15m']+index['NeII_12m'] )/index['NeV_24m']))
                  #Argon 
                  if Ar2Ar3_obs == -10: 
                     CHI_Ar2Ar3 = 0
                  elif index['ArII_7m'] == 0 or index['ArIII_9m'] == 0:
                     CHI_Ar2Ar3 = tol_max
                  else:   
                     CHI_Ar2Ar3 = (np.log10((index['ArII_7m']/index['ArIII_9m']))- Ar2Ar3_obs)**2/np.log10((index['ArII_7m']/index['ArIII_9m']))
                  if Ar235a_obs == -10: 
                     CHI_Ar235a = 0
                  elif index['HI_4m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235a = tol_max
                  else:   
                     CHI_Ar235a = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_4m'])- Ar235a_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_4m'])
                  if Ar235b_obs == -10: 
                     CHI_Ar235b = 0
                  elif index['HI_7m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235b = tol_max
                  else:   
                     CHI_Ar235b = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_7m'])- Ar235b_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_7m'])
                  if Ar235c_obs == -10: 
                     CHI_Ar235c = 0
                  elif index['HI_12m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235c = tol_max
                  else:   
                     CHI_Ar235c = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_12m'])- Ar235c_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_8m'])/index['HI_12m'])
                  if Ne235d_obs == -10: 
                     CHI_Ar235d = 0
                  elif index['HI_4m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235d = tol_max
                  else:   
                     CHI_Ar235d = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_4m'])- Ar235d_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_4m'])
                  if Ar235e_obs == -10: 
                     CHI_Ar235e = 0
                  elif index['HI_7m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235e = tol_max
                  else:   
                     CHI_Ar235e = (np.log10((index['ArII_7m']+index['NeIII_15m']+index['ArV_13m'])/index['HI_7m'])- Ne235e_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_7m'])
                  if Ar235f_obs == -10: 
                     CHI_Ar235f = 0
                  elif index['HI_12m'] == 0 or (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235f = tol_max
                  else:   
                     CHI_Ar235f = (np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_12m'])- Ar235f_obs)**2/np.log10((index['ArII_7m']+index['ArIII_9m']+index['ArV_13m'])/index['HI_12m'])
                  if Ar235g_obs == -10:
                     CHI_Ar235g = 0
                  elif (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_8m'] == 0):
                     CHI_Ar235g = tol_max
                  else:
                     CHI_Ar235g = (np.log10(index['ArII_7m'] + index['ArIII_9m'] + index['ArV_8m']) - Ar235g_obs)**2/np.log10(index['ArII_7m'] + index['ArIII_9m'] + index['ArV_8m'])
                  if Ar235h_obs == -10:
                     CHI_Ar235h = 0
                  elif (index['ArII_7m'] == 0 and index['ArIII_9m'] == 0 and index['ArV_13m'] == 0):
                     CHI_Ar235h = tol_max
                  else:
                     CHI_Ar235h = (np.log10(index['ArII_7m'] + index['ArIII_9m'] + index['ArV_13m']) - Ar235h_obs)**2/np.log10(index['ArII_7m'] + index['ArIII_9m'] + index['ArV_13m'])
                  if Ar23Ar5a_obs == -10: 
                     CHI_Ar23Ar5a = 0
                  elif (index['ArIII_9m'] == 0 and index['ArII_7m'] == 0) or index['ArV_8m'] == 0:
                     CHI_Ar23Ar5a = tol_max
                  else:   
                     CHI_Ar23Ar5a = (np.log10(( ( index['ArIII_9m'] + index['ArII_7m'])/index['ArV_8m']))- Ar23Ar5a_obs)**2/np.log10(( (index['ArIII_9m']+index['ArII_7m'] )/index['ArV_8m']))
                  if Ar23Ar5b_obs == -10: 
                     CHI_Ar23Ar5b = 0
                  elif (index['ArIII_9m'] == 0 and index['ArII_7m'] == 0) or index['ArV_13m'] == 0:
                     CHI_Ar23Ar5b = tol_max
                  else:   
                     CHI_Ar23Ar5b = (np.log10(( ( index['ArIII_9m'] + index['ArII_7m'])/index['ArV_13m']))- Ar23Ar5b_obs)**2/np.log10(( (index['ArIII_9m']+index['ArII_7m'] )/index['ArV_13m']))
                  #Sulphur
                  if S34a_obs == -10: 
                     CHI_S34a = 0
                  elif index['HI_4m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0):
                     CHI_S34a = tol_max
                  else:   
                     CHI_S34a = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_4m'])- S34a_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_4m'])
                  if S34b_obs == -10: 
                     CHI_S34b = 0
                  elif index['HI_7m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0):
                     CHI_S34b = tol_max
                  else:   
                     CHI_S34b = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_7m'])- S34b_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_7m'])
                  if S34c_obs == -10: 
                     CHI_S34c = 0
                  elif index['HI_12m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_19m'] == 0): 
                     CHI_S34c = tol_max
                  else:   
                     CHI_S34c = (np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_12m'])- S34c_obs)**2/np.log10((index['SIV_10m']+index['SIII_19m'])/index['HI_12m'])
                  if S34d_obs == -10: 
                     CHI_S34d = 0
                  elif index['HI_4m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34d = tol_max
                  else:   
                     CHI_S34d = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_4m'])- S34d_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_4m'])
                  if S34e_obs == -10: 
                     CHI_S34e = 0
                  elif index['HI_7m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34e = tol_max
                  else:   
                     CHI_S34e = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_7m'])- S34e_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_7m'])
                  if S34f_obs == -10: 
                     CHI_S34f = 0
                  elif index['HI_12m'] == 0 or (index['SIV_10m'] == 0 and index['SIII_33m'] == 0): 
                     CHI_S34f = tol_max
                  else:   
                     CHI_S34f = (np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_12m'])- S34f_obs)**2/np.log10((index['SIV_10m']+index['SIII_33m'])/index['HI_12m'])
                  if S34g_obs == -10:
                     CHI_S34g = 0
                  elif (index['SIV_10m'] == 0 and index['SIII_19m'] == 0):
                     CHI_S34g = tol_max
                  else:
                     CHI_S34g = (np.log10(index['SIV_10m'] + index['SIII_19m']) - S34g_obs)**2/np.log10(index['SIV_10m'] + index['SIII_19m'])
                  if S34h_obs == -10:
                     CHI_S34h = 0
                  elif (index['SIV_10m'] == 0 and index['SIII_33m'] == 0):
                     CHI_S34h = 0
                  else:
                     CHI_S34h = (np.log10(index['SIV_10m'] + index['SIII_33m']) - S34h_obs)**2/np.log10(index['SIV_10m'] + index['SIII_33m'])
                  if S3S4a_obs == -10: 
                     CHI_S3S4a = 0
                  elif index['SIV_10m'] == 0 or index['SIII_19m'] == 0 : 
                     CHI_S3S4a = tol_max
                  else:   
                     CHI_S3S4a = (np.log10(index['SIII_19m']/index['SIV_10m'])- S3S4a_obs)**2/np.log10((index['SIII_19m']/index['SIV_10m']))
                  if S3S4b_obs == -10: 
                     CHI_S3S4b = 0
                  elif index['SIV_10m'] == 0 or index['SIII_33m'] == 0 : 
                     CHI_S3S4b = tol_max
                  else:   
                     CHI_S3S4b = (np.log10(index['SIII_33m']/index['SIV_10m'])- S3S4b_obs)**2/np.log10((index['SIII_33m']/index['SIV_10m']))
                  #Oxygen
                  if O34a_obs == -10: 
                     CHI_O34a = 0
                  elif index['HI_4m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0):
                     CHI_O34a = tol_max
                  else:   
                     CHI_O34a = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_4m'])- O34a_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_4m'])
                  if O34b_obs == -10: 
                     CHI_O34b = 0
                  elif index['HI_7m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0):
                     CHI_O34b = tol_max
                  else:   
                     CHI_O34b = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_7m'])- O34b_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_7m'])
                  if O34c_obs == -10: 
                     CHI_O34c = 0
                  elif index['HI_12m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_52m'] == 0): 
                     CHI_O34c = tol_max
                  else:   
                     CHI_O34c = (np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_12m'])- O34c_obs)**2/np.log10((index['OIV_26m']+index['OIII_52m'])/index['HI_12m'])
                  if O34d_obs == -10: 
                     CHI_O34d = 0
                  elif index['HI_4m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34d = tol_max
                  else:   
                     CHI_O34d = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_4m'])- O34d_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_4m'])
                  if O34e_obs == -10: 
                     CHI_O34e = 0
                  elif index['HI_7m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34e = tol_max
                  else:   
                     CHI_O34e = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_7m'])- O34e_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_7m'])
                  if O34f_obs == -10: 
                     CHI_O34f = 0
                  elif index['HI_12m'] == 0 or (index['OIV_26m'] == 0 and index['OIII_88m'] == 0): 
                     CHI_O34f = tol_max
                  else:   
                     CHI_O34f = (np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_12m'])- O34f_obs)**2/np.log10((index['OIV_26m']+index['OIII_88m'])/index['HI_12m'])
                  if O34g_obs == -10:
                     CHI_O34g = 0
                  elif (index['OIV_26m'] == 0 and index['OIII_52m'] == 0):
                     CHI_O34g = tol_max
                  else:
                     CHI_O34g = (np.log10(index['OIV_26m'] + index['OIII_52m']) - O34g_obs)**2/np.log10(index['OIV_26m'] + index['OIII_52m'])
                  if O34h_obs == -10:
                     CHI_O34h = 0
                  elif (index['OIV_26m'] == 0 and index['OIII_88m'] == 0):
                     CHI_O34h = tol_max
                  else:
                     CHI_O34h = (np.log10(index['OIV_26m'] + index['OIII_88m']) - O34h_obs)**2/np.log10(index['OIV_26m'] + index['OIII_88m'])
                  if O3O4a_obs == -10: 
                     CHI_O3O4a = 0
                  elif index['OIV_26m'] == 0 or index['OIII_52m'] == 0 : 
                     CHI_O3O4a = tol_max
                  else:   
                     CHI_O3O4a = (np.log10(index['OIII_52m']/index['OIV_26m'])- O3O4a_obs)**2/np.log10((index['OIII_52m']/index['OIV_26m']))
                  if O3O4b_obs == -10: 
                     CHI_O3O4b = 0
                  elif index['OIV_26m'] == 0 or index['OIII_88m'] == 0 : 
                     CHI_O3O4b = tol_max
                  else:   
                     CHI_O3O4b = (np.log10(index['OIII_88m']/index['OIV_26m'])- O3O4b_obs)**2/np.log10((index['OIII_88m']/index['OIV_26m']))
                  #Nitrogen   
                  if N23a_obs == -10: 
                     CHI_N23a = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23a = tol_max
                  else:   
                     CHI_N23a = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_4m'])- N23a_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_4m'])
                  if N23b_obs == -10: 
                     CHI_N23b = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23b = tol_max
                  else:   
                     CHI_N23b = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_7m'])- N23b_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_7m'])
                  if N23c_obs == -10: 
                     CHI_N23c = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0 and index['NII_122m'] == 0): 
                     CHI_N23c = tol_max
                  else:   
                     CHI_N23c = (np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_12m'])- N23c_obs)**2/np.log10((index['NIII_57m']+index['NII_122m'])/index['HI_12m'])          
                  if N3a_obs == -10: 
                     CHI_N3a = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3a = tol_max
                  else:   
                     CHI_N3a = (np.log10((index['NIII_57m'])/index['HI_4m'])- N3a_obs)**2/np.log10((index['NIII_57m'])/index['HI_4m'])
                  if N3b_obs == -10: 
                     CHI_N3b = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3b = tol_max
                  else:   
                     CHI_N3b = (np.log10((index['NIII_57m'])/index['HI_7m'])- N3b_obs)**2/np.log10((index['NIII_57m'])/index['HI_7m'])
                  if N3c_obs == -10: 
                     CHI_N3c = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0): 
                     CHI_N3c = tol_max
                  else:   
                     CHI_N3c = (np.log10((index['NIII_57m'])/index['HI_12m'])- N3c_obs)**2/np.log10((index['NIII_57m'])/index['HI_12m'])
                  if N3d_obs == -10:
                     CHI_N3d = 0
                  elif index['NIII_57m'] == 0:
                     CHI_N3d = tol_max
                  else:
                     CHI_N3d = (np.log10(index['NIII_57m']) - N3d_obs)**2/np.log10(index['NIII_57m'])           
                  if N2N3a_obs == -10: 
                     CHI_N2N3a = 0
                  elif index['NIII_57m'] == 0 or index['NII_122m'] == 0: 
                     CHI_N2N3a = tol_max
                  else:   
                     CHI_N2N3a = (np.log10(index['NII_122m']/index['NIII_57m'])- N2N3a_obs)**2/np.log10((index['NII_122m']/index['NIII_57m']))
                  if N23d_obs == -10: 
                     CHI_N23d = 0
                  elif index['HI_4m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23d = tol_max
                  else:   
                     CHI_N23d = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_4m'])- N23d_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_4m'])
                  if N23e_obs == -10: 
                     CHI_N23e = 0
                  elif index['HI_7m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23e = tol_max
                  else:   
                     CHI_N23e = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_7m'])- N23e_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_7m'])
                  if N23f_obs == -10: 
                     CHI_N23f = 0
                  elif index['HI_12m'] == 0 or (index['NIII_57m'] == 0 and index['NII_205m'] == 0): 
                     CHI_N23f = tol_max
                  else:   
                     CHI_N23f = (np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_12m'])- N23f_obs)**2/np.log10((index['NIII_57m']+index['NII_205m'])/index['HI_12m'])
                  if N23g_obs == -10:
                     CHI_N23g = 0
                  elif (index['NIII_57m'] == 0 and index['NII_122m'] == 0):
                     CHI_N23g = tol_max
                  else:
                     CHI_N23g = (np.log10(index['NIII_57m']+index['NII_122m']) - N23g_obs)**2/np.log10(index['NIII_57m']+index['NII_122m'])
                  if N23h_obs == -10:
                     CHI_N23h = 0
                  elif (index['NIII_57m'] == 0 and index['NII_205m'] == 0):
                     CHI_N23h = tol_max
                  else:
                     CHI_N23h = (np.log10(index['NIII_57m']+index['NII_205m']) - N23h_obs)**2/np.log10(index['NIII_57m']+index['NII_205m'])
                  if N2N3b_obs == -10: 
                     CHI_N2N3b = 0
                  elif index['NIII_57m'] == 0 or index['NII_205m'] == 0: 
                     CHI_N2N3b = tol_max
                  else:   
                     CHI_N2N3b = (np.log10(index['NII_205m']/index['NIII_57m'])- N2N3b_obs)**2/np.log10((index['NII_205m']/index['NIII_57m']))
                  #O3N2
                  if O3N2a_obs == -10: 
                     CHI_O3N2a = 0
                  elif index['OIII_52m'] == 0 or index['NII_122m'] == 0: 
                     CHI_O3N2a = tol_max
                  else:   
                     CHI_O3N2a = (np.log10(index['OIII_52m']/index['NII_122m'])- O3N2a_obs)**2/np.log10((index['OIII_52m']/index['NII_122m']))
                  if O3N2b_obs == -10: 
                     CHI_O3N2b = 0
                  elif index['OIII_88m'] == 0 or index['NII_122m'] == 0: 
                     CHI_O3N2b = tol_max
                  else:   
                     CHI_O3N2b = (np.log10(index['OIII_88m']/index['NII_122m'])- O3N2b_obs)**2/np.log10((index['OIII_88m']/index['NII_122m']))
                  if O3N2c_obs == -10: 
                     CHI_O3N2c = 0
                  elif index['OIII_52m'] == 0 or index['NII_205m'] == 0: 
                     CHI_O3N2c = tol_max
                  else:   
                     CHI_O3N2c = (np.log10(index['OIII_52m']/index['NII_205m'])- O3N2c_obs)**2/np.log10((index['OIII_52m']/index['NII_205m']))
                  if O3N2d_obs == -10: 
                     CHI_O3N2d = 0
                  elif index['OIII_88m'] == 0 or index['NII_205m'] == 0: 
                     CHI_O3N2d = tol_max
                  else:   
                     CHI_O3N2d = (np.log10(index['OIII_88m']/index['NII_205m'])- O3N2d_obs)**2/np.log10((index['OIII_88m']/index['NII_205m']))

                  CHI_OH = (CHI_Ne2Ne3**2 + CHI_N23a**2+ CHI_N23b**2 +CHI_N23c**2 + CHI_N2N3a**2 + CHI_O3N2a**2 + CHI_O3N2b**2 + CHI_N23d**2 + CHI_N23e**2 + CHI_N23f**2 + CHI_N2N3b**2 + CHI_O3N2c**2 + CHI_O3N2d**2 + CHI_O34a**2 + CHI_O34b**2 + CHI_O34c**2 + CHI_O34d**2 + CHI_O34e**2 + CHI_O34f**2 + CHI_O3O4a**2 + CHI_O3O4b**2 + CHI_Ne235a**2 + CHI_Ne235b**2 + CHI_Ne235c**2 + CHI_Ne235d**2 + CHI_Ne235e**2 + CHI_Ne235f**2 + CHI_Ne23Ne5a**2 + CHI_Ne23Ne5b**2 + CHI_N3a**2 + CHI_N3b**2 + CHI_N3c**2 + CHI_Ar2Ar3**2 + CHI_Ar235a**2 + CHI_Ar235b**2 + CHI_Ar235c**2 + CHI_Ar235d**2 + CHI_Ar235e**2 + CHI_Ar235f**2 + CHI_Ar23Ar5a**2 +CHI_Ar23Ar5b**2 + CHI_Ne235g**2 + CHI_Ne235h**2 + CHI_Ar235g**2 + CHI_Ar235h**2 + CHI_O34g**2 + CHI_O34h**2 + CHI_N3d**2 + CHI_N23g**2 + CHI_N23h**2)**0.5

                  CHI_SH = (CHI_S34a**2 + CHI_S34b**2 + CHI_S34c**2 + CHI_S34d**2 + CHI_S34e**2 + CHI_S34f**2 +CHI_S3S4a**2 + CHI_S3S4b**2 + CHI_S34g**2 + CHI_S34h**2 )**0.5           

                  OH_p = index['12logOH'] / CHI_OH**2 + OH_p
                  logU_p = index['logU'] / CHI_OH**2 + logU_p
                  den_OH = 1 / CHI_OH**2 + den_OH
                  SH_p = index['12logOH'] / CHI_SH**2 + SH_p
                  den_SH = 1 / CHI_SH**2 + den_SH

            if OH == 0:
               OHf = OH
               logUf = logU
            else:         
               OHf = OH_p / den_OH
               logUf = logU_p / den_OH
            if SH == 0:
               SHf = SH
            else:
               SHf = SH_p / den_SH - 1.57

         OH_mc.append(OHf)
         logU_mc.append(logUf)
         SH_mc.append(SHf)
         OHe_mc.append(eOH)
         logUe_mc.append(elogU)
         SHe_mc.append(eSH)
      
      OHff = np.mean(OH_mc)
      eOHff = (np.std(OH_mc)**2+np.mean(OHe_mc)**2)**0.5
      logUff = np.mean(logU_mc)
      elogUff = (np.std(logU_mc)**2+np.mean(logUe_mc)**2)**0.5
      SHff = np.mean(SH_mc)
      eSHff = (np.std(SH_mc)**2 + np.mean(SHe_mc)**2)**0.5

   except:
      OHff = 9999
      eOHff = 9999
      NOff = 9999
      eNOff = 9999
      logUff = 9999
      elogUff = 9999
      SHff = 9999
      eSHff = 9999
      
   OHffs.append(OHff)
   eOHffs.append(eOHff)
   NOffs.append(NOff)
   eNOffs.append(eNOff)
   logUffs.append(logUff)
   elogUffs.append(elogUff)
   SHffs.append(SHff)
   eSHffs.append(eSHff)
         
   ##################################
   # Displaying results in terminal #
   ##################################

   if input0.size == 1 and tab==0: continue
   print (round(100*(count)/float(len(input1)),1),'%',Names[tab],grid_type,'', round(OHff,2), round(eOHff,2),'',round(SHff,2), round(eSHff,2),'',round(NOff,2), round(eNOff,2), '',round(logUff,2), round(elogUff,2))

####################################################
###### OUTPUT FORMAT AND INFORMATION: RESULTS ######
####################################################

#Grid used and results from the free parameters
output['grid'] = grids
output['OH'] = OHffs
output['eOH'] = eOHffs
output['SH'] = SHffs
output['eSH'] = eSHffs
output['NO'] = NOffs
output['eNO'] = eNOffs
output['logU'] = logUffs
output['elogU'] = elogUffs

#Header comments for the file
lineas_header = [' HII-CHI-mistry-IR v.3.2 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'']
#Labels for columns (emission lines)
line_label = '{:30}  '.format(output.dtype.names[0])
for ind2 in range(1, len(output.dtype.names)):
   line_label += '{:30}  '.format(output.dtype.names[ind2])

#Labels for columns 
lineas_header.append(line_label)
header = '\n'.join(lineas_header)

#Results
np.savetxt('.'.join(input00.split('.')[:-1])+'_hcm-ir-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*(len(output.dtype.names)-10)+['%i']+['%.2f']*8),header=header)

lines_stor = []
with open('.'.join(input00.split('.')[:-1])+'_hcm-ir-output.dat', 'r+') as output_file:
   for line in output_file:
      lines_stor.append(line)

#Reformating output for better reading of the table
file_overwrite = open('.'.join(input00.split('.')[:-1])+'_hcm-ir-output.dat', 'r+')
file_overwrite.seek(0)
for line_n in lines_stor:  
   if line_n[0] == '#' and line_n[2:4] == 'ID':
      file_overwrite.write(line_n[2:])
   else:
      file_overwrite.write(line_n)
file_overwrite.truncate()   
file_overwrite.close()

print ('-------------------------------------------------')
print ('Results are stored in ' + '.'.join(input00.split('.')[:-1])+'_hcm-ir-output.dat')
print ('-------------------------------------------------')

#############################################
###### INFORMATION AND CONTACT DETAILS ######
#############################################

# Enrique Perez-Montero, epm@iaa.es
# Borja Perez-Diaz, bperez@iaa.es


#################
###### END ######
#################
