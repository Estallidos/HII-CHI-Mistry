# Filename: HCm_UV_v5.3.py

#####################
###### IMPORTS ######
#####################

import string
import numpy as np
import sys
sys.stderr = open('errorlog.txt', 'w')
import warnings
warnings.filterwarnings("ignore")

#######################
###### FUNCTIONS ######
#######################

#Function for interpolation of grids

def interpolate(grid,z,zmin,zmax,n):
   #Columns of the library
   n_comments = 0
   with open('Libraries_uv/C17_POPSTAR_1myr_uv.dat', 'r') as file1:
      for line in file1:
         if line[0] == '#':
            n_comments += 1
   auxiliar_labels = np.genfromtxt('Libraries_uv/C17_POPSTAR_1myr_uv.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
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
print ('This is HII-CHI-mistry for UV version 5.3')
print ('See Perez-Montero, & Amorin (2017) for details')
print ('Insert the name of your input text file with some or all of the following columns:')
print (' Lya 1216')
print (' NV] 1239')
print (' CIV] 1549')
print (' HeII 1640')
print (' OIII] 1665')
print (' CIII] 1909')
print (' [OII] 3727+3729')
print (' [NeIII] 3869')
print (' [OIII] 4363')
print (' [OIII] 4959')
print (' Hb 4861')
print (' OIII 5007')
print ('in arbitrary units and reddening corrected. Each column must be given with labels for the lines and their corresponding flux errors.')
print ('-------------------------------------------------')


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

#############################################
###### NON INTERACTIVE OR INTERACTIVE CODE ######
#############################################

interactive = True #Change this value to False to run the code in non-interactive mode

#Questions (inputs fromt terminal)
question1 = interactive #Question to select the grids of models
question2 = interactive #Question for the value of alphaOX (only required if "sed" is set to 3)
question2a = interactive #Question for the value of effective temperature (only required if "sed" is set to 4)
question3 = interactive #Question for the value of efrac (only required if "sed" is set to 3), (4), or (5)
question3a = interactive #Question for the consideration of dust within photoionization models for sed = (3) and alpha = (1) or (3)
question4 = interactive #Question required to introduce a particular file by the user (only required if "sed" is set 4)
question6 = interactive #Question to use or not interpolation for the grids
question7 = interactive #Question regarding the constrain law to be assumed by the code
question8 = interactive #Question required to tntroduced a particular file by the user (only if "const" is set to 6)

#Set values to operate
if question1 == False:
   sed = 1 #Choose grid of model: (1) for POPSTAR models; (2) for BPASS models; (3) for AGN models; (4) pAGB models assuming different effective temperatures; (5) advected-dominated models; (6) models introduced by the user (will require file). Replace value for the given option
if question2 == False and sed == 3:
   alpha = 2 #Choose value of alpha_OX: (1) for alpha_OX = -0.8; (2) for alpha_OX = -1.0; (3) for alpha_OX = -1.2; (4) for alpha_OX = -1.4; (5) for alpha_OX = -1.6; (6) for alpha_OX = -1.8; (7) for alpha_OX = -2.0. Replace value with any of these options.
if question2a == False and sed == 4:
   teff_pagb = 1 #Choose value of the effective temperature: (1) T_eff = 50 kK; (2) T_eff = 100 kK; (3) T_eff = 150 kK
if question3 == False and (sed == 3 or sed == 4 or sed == 5):
   efrac = 1 #Choose value for the stopping criteria, i.e., the fraction of free electrons: (1) for efrac = 0.02; (2) for efrac = 0.98
if question3a == False and sed == 3 and (alpha == 1 or alpha == 3):
   grains = 1 #Choose value for dust: (1) grains are considered in the photoionization models; (2) photoionization models without dust
if question4 == False and sed == 6:
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



#Interface with the user
print ('')
while question1:
   print ('-------------------------------------------------')
   print ('Default SEDs')
   print ('------------')
   print ('(1) POPSTAR with Chabrier IMF, age = 1 Myr')
   print ('(2) BPASS v.2.1 a_IMF = 1.35, Mup = 300, age = 1Myr')
   print ('(3) AGN, double component,  a(UV) = -1.0')
   print ('(4) pAGB NLTE models (Rauch 2003)')
   print ('(5) Advected-dominated Accretion Flow for AGN (ADAF, Nemmen et al. 2014)')
   print (' ')
   print ('Other SED')
   print ('---------')
   print ('(6) Different library')
   print ('-------------------------------------------------')
   if int(sys.version[0]) < 3:
      sed = raw_input('Choose SED of the models: ')
   else:
      sed = input('Choose SED of the models: ')
   if sed == '1' or sed == '2' or sed == '3' or sed == '4' or sed == '5' or sed == '6': question1 = False 
print ('')

# Further questions on the AGN models
if sed == '3':
   #SLOPE ALPHA
   while question2:
      if int(sys.version[0]) < 3:
         alpha = raw_input('Choose value for alpha(OX) in the AGN models: [1] -0.8; [2] -1.0; [3] -1.2; [4] -1.4; [5] -1.6; [6] -1.8; [7] -2.0: ')
      else:
         alpha = input('Choose value for alpha(OX) in the AGN models: [1] -0.8; [2] -1.0; [3] -1.2; [4] -1.4; [5] -1.6; [6] -1.8; [7] -2.0:  ')
      if alpha == '1' or alpha == '2' or alpha == '3' or alpha == '4' or alpha == '5' or alpha == '6' or alpha == '7': question2 = False
      alpha = int(alpha)
   print ('')
   #Fraction of free electrons (stopping criteria in the models)
   while question3:
      if int(sys.version[0]) < 3:
         efrac = raw_input('Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons: ')
      else:
         efrac = input('Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons:  ')
      if efrac == '1' or efrac == '2': question3 = False
      efrac = int(efrac)
   print ('')
   #Presence or absence of dust in the models
   if alpha == 1 or alpha == 3:
      while question3a:
         if int(sys.version[0]) < 3:
            grains = raw_input('Choose AGN models with [1] or without [2] dust grains: ')
         else:
            grains = input('Choose AGN models with [1] or without [2] dust grains:  ')
         if grains == '1' or grains == '2': question3a = False
   else:
      grains = 1
   print ('')


# Further questions on the pAGB models
if sed == '4':
   #Effective temperature in pAGB models
   while question2a:
      if int(sys.version[0]) < 3:
         teff_pagb = raw_input('Choose value for T_effective in pAGB models: [1] 50 kK; [2] 100 kK; [3] 150 kK: ')
      else:
         teff_pagb = input('Choose value for T_effective in pAGB models: [1] 50 kK; [2] 100 kK; [3] 150 kK:  ')
      if teff_pagb == '1' or teff_pagb == '2' or teff_pagb == '3': question2a = False
      teff_pagb = int(teff_pagb)
   print ('')
   #Fraction of free electrons (stopping criteria in the models)
   while question3:
      if int(sys.version[0]) < 3:
         efrac = raw_input('Choose stop criterion in the pAGB models: [1] 2% free electrons [2] 98% free electrons: ')
      else:
         efrac = input('Choose stop criterion in the pAGB models: [1] 2% free electrons [2] 98% free electrons:  ')
      if efrac == '1' or efrac == '2': question3 = False
      efrac = int(efrac)
   print ('')

# Further questions on the ADAF model
if sed == '5':
   #Fraction of free electrons (stopping criteria in the models)
   while question3:
      if int(sys.version[0]) < 3:
         efrac = raw_input('Choose stop criterion in the pAGB models: [1] 2% free electrons [2] 98% free electrons: ')
      else:
         efrac = input('Choose stop criterion in the pAGB models: [1] 2% free electrons [2] 98% free electrons:  ')
      if efrac == '1' or efrac == '2': question3 = False
      efrac = int(efrac)
   print ('')
 

#Particular file introduced by the user
if sed == '6':
   while question4:
      print ('Introduce name of the file containing the models. It must be located in the folder "Libraries_uv".')
      print (' ')
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
         print (' ')
         print ('Loading library '+new_library+'. Checking correct format of the file.')
         question4 = False
      except:
         print (' ')
         print ('Library was not found in folder "Libraries_uv" or file does not exist.')
   question5 = True
   while question5:
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
         with open('Libraries_uv/C17_POPSTAR_1myr_uv.dat', 'r') as file5:
            for line in file5:
               if line[0] == '#':
                  n_comments += 1
         auxiliar_labels = np.genfromtxt('Libraries_uv/C17_POPSTAR_1myr_uv.dat', dtype=None, names=True, encoding = 'ascii', skip_header=n_comments).dtype.names
         missing_labels = []
         for label in auxiliar_labels:
            if label in library_user.dtype.names:
               continue
            else:
               missing_labels.append(label)
         #Displaying message for the user:
         print('Succesfully reading of the file')
         if len(missing_labels) == 0:
            print ('File presents the correct format')
            question5 = False
         else:
            print ('File does not present the correct format. The following columns are missing:')
            for need_label in missing_labels:
               print('- '+need_label)
            print ('More details on the correct format for the library are found in readme file.')
            print (' ')
            print ('Reintroduce name of the file with fixed format:')
            print (' ')
            if int(sys.version[0]) < 3:
               new_library = raw_input('Name of file: ')
            else:
               new_library = input('Name of file: ')
      except:
         print ('Something went wrong while reading file. Please, reintroduce name of the file:')
         print ('')
         if int(sys.version[0]) < 3:
            new_library = raw_input('Name of file: ')
         else:
            new_library = input('Name of file: ')

#Interpolation in the grid of models
print ('')
while question6:
   if int(sys.version[0]) < 3:
      inter = raw_input('Choose models [0] No interpolated [1] Interpolated: ')
   else:
      inter = input('Choose models [0] No interpolated [1] Interpolated: ')
   if inter == '0' or inter == '1': question6 = False
print ('')

sed = int(sed)
inter = int(inter)
if sed == 3:
   alpha = int(alpha)
   efrac = int(efrac)
   grains = int(grains)


#POPSTAR MODEL
if sed==1:
   file_lib = 'C17_POPSTAR_1myr_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file6:
      for line in file6:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation.'
      print ('No interpolation for the POPSTAR models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      print ('')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. Interpolation.'
      print ('Interpolation for the POPSTAR models is going to be used.')
      print ('The grid has a resolution of 0.01dex for O/H and 0.0125dex for C/O.')
      print ('')
      res_CO = 0.125
        
#BPASS MODEL
elif sed==2:
   file_lib = 'C17_BPASS_IMF135_mup300_1myr_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file7:
      for line in file7:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)  
   if inter == 0:
      sed_type = 'BPASS a_IMF = 1.35, M_up = 300, age = 1Myr, with binaries. No interpolation.'
      print ('No interpolation for the BPASS models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      print ('')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'BPASS v.2.1, a_IMF = 1.35, M_up = 300, age = 1Myr. Interpolation.'
      print ('Interpolation for the BPASS  models is going to be used.')
      print ('The grid has a resolution of 0.01dex for O/H and 0.0125dex for C/O.')
      print ('')
      res_CO = 0.125

#AGN MODEL FOR alpha_OX = -0.8, efrac = 2%, with dust grains
elif sed==3 and alpha ==1 and efrac == 1 and grains == 1:
   file_lib = 'C17_AGN_alpha08_efrac02_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 2% with dust grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -0.8 with 2% free electrons and dust grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -0.8, free electron fraction = 2% and with dust grains. Interpolation.'
      print ('Interpolation for the AGN a(ox) = -0.8, 2% free electrons and with dust models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O.')
      res_CO = 0.125

#AGN MODEL FOR alpha_OX = -0.8, efrac = 2%, without dust grains
elif sed==3 and alpha ==1 and efrac == 1 and grains == 2:
   file_lib = 'C17_AGN_alpha08_efrac02_CNfix_nograins_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 2% without dust grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -0.8 with 2% free electrons models without grains is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -0.8, free electron fraction = 2% and without dust grains. Interpolation.'
      print ('Interpolation for the AGN a(ox) = -0.8, 2% free electrons and without dust models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O.')
      res_CO = 0.125
        
#AGN MODEL FOR alpha_OX = -0.8, efrac = 98%, with dust grains
elif sed==3 and alpha ==1 and efrac == 2 and grains == 1:
   file_lib = 'C17_AGN_alpha08_efrac98_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 98% with dust grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -0.8 with 98% free electrons and dust grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -0.8, free electron fraction = 98% and with dust grains. Interpolation.'
      print ('Interpolation for the AGN a(ox) = -0.8, 98% free electrons and with dust models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O.')
      res_CO = 0.125

#AGN MODEL FOR alpha_OX = -0.8, efrac = 98%, without dust grains
elif sed==3 and alpha ==1 and efrac == 2 and grains == 2:
   file_lib = 'C17_AGN_alpha08_efrac98_CNfix_nograins_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 98% without dust grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -0.8 with 98% free electrons models without grains is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -0.8, free electron fraction = 98% and without dust grains. Interpolation.'
      print ('Interpolation for the AGN a(ox) = -0.8, 98% free electrons and without dust models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O.')
      res_CO = 0.125

#AGN MODEL FOR alpha_OX = -1.0, efrac = 2% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==2 and efrac == 1:
   file_lib = 'C17_AGN_alpha10_efrac02_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.0 and free electron fraction = 2%. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.0 with 2% free electrons models with grains is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.0 and free electron fraction = 2% with grains interpolated.'
      print ('Interpolation for the AGN a(ox) = -1.0 with 2% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#AGN MODEL FOR alpha_OX = -1.0, efrac = 98% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==2 and efrac == 2:
   file_lib = 'C17_AGN_alpha10_efrac98_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.0 and free electron fraction = 98% with grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.0 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.0 and free electron fraction = 98% interpolated.'
      print ('Interpolation for the AGN a(ox) = -1.0 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125



#AGN MODEL FOR alpha_OX = -1.2, efrac = 2%, with dust grains
elif sed==3 and alpha ==2 and efrac == 1 and grains == 1:
   file_lib = 'C17_AGN_alpha12_efrac02_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 2% with dust grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.2 with 2% free electrons and dust grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.2, free electron fraction = 2% and with dust grains. Interpolation.'
      print ('Interpolation for the AGN a(ox) = -1.2, 2% free electrons and with dust models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O.')
      res_CO = 0.125

#AGN MODEL FOR alpha_OX = -1.2, efrac = 2%, without dust grains
elif sed==3 and alpha ==2 and efrac == 1 and grains == 2:
   file_lib = 'C17_AGN_alpha12_efrac02_CNfix_nograins_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 2% without dust grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.2 with 2% free electrons models without grains is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.2, free electron fraction = 2% and without dust grains. Interpolation.'
      print ('Interpolation for the AGN a(ox) = -1.2, 2% free electrons and without dust models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O.')
      res_CO = 0.125
        
#AGN MODEL FOR alpha_OX = -1.2, efrac = 98%, with dust grains
elif sed==3 and alpha ==2 and efrac == 2 and grains == 1:
   file_lib = 'C17_AGN_alpha12_efrac98_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 98% with dust grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.2 with 98% free electrons and dust grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.2, free electron fraction = 98% and with dust grains. Interpolation.'
      print ('Interpolation for the AGN a(ox) = -1.2, 98% free electrons and with dust models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O.')
      res_CO = 0.125

#AGN MODEL FOR alpha_OX = -1.2, efrac = 98%, without dust grains
elif sed==3 and alpha ==2 and efrac == 2 and grains == 2:
   file_lib = 'C17_AGN_alpha12_efrac98_CNfix_nograins_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file8:
      for line in file8:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 98% without dust grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.2 with 98% free electrons models without grains is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O.')
      res_NO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.2, free electron fraction = 98% and without dust grains.  Interpolation.'
      print ('Interpolation for the AGN a(ox) = -1.2, 98% free electrons and without dust models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O.')
      res_NO = 0.125


#AGN MODEL FOR alpha_OX = -1.4, efrac = 2% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==4 and efrac == 1:
   file_lib = 'C17_AGN_alpha14_efrac02_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.4 and free electron fraction = 2% with grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.4 with 2% free electrons and grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.4 and free electron fraction = 2% with grains interpolated.'
      print ('Interpolation for the AGN a(ox) = -1.4 with 2% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#AGN MODEL FOR alpha_OX = -1.4, efrac = 98% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==4 and efrac == 2:
   file_lib = 'C17_AGN_alpha14_efrac98_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.4 and free electron fraction = 98% with grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.4 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.4 and free electron fraction = 98% with grains interpolated.'
      print ('Interpolation for the AGN a(ox) = -1.4 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')

#AGN MODEL FOR alpha_OX = -1.6, efrac = 2% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==5 and efrac == 1:
   file_lib = 'C17_AGN_alpha16_efrac02_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.6 and free electron fraction = 2% with grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.6 with 2% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.6 and free electron fraction = 2% with grains interpolated.'
      print ('Interpolation for the AGN a(ox) = -1.6 with 2% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#AGN MODEL FOR alpha_OX = -1.6, efrac = 98% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==5 and efrac == 2:
   file_lib = 'C17_AGN_alpha16_efrac98_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.6 and free electron fraction = 98% with grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.6 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.6 and free electron fraction = 98% with grains interpolated.'
      print ('Interpolation for the AGN a(ox) = -1.6 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')

#AGN MODEL FOR alpha_OX = -1.8, efrac = 2% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==6 and efrac == 1:
   file_lib = 'C17_AGN_alpha18_efrac02_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.8 and free electron fraction = 2%. No with grains interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.8 with 2% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.8 and free electron fraction = 2% interpolated with grains.'
      print ('Interpolation for the AGN a(ox) = -1.8 with 2% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#AGN MODEL FOR alpha_OX = -1.8, efrac = 98% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==6 and efrac == 2:
   file_lib = 'C17_AGN_alpha18_efrac98_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -1.8 and free electron fraction = 98% with grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -1.8 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -1.8 and free electron fraction = 98% with grains interpolated.'
      print ('Interpolation for the AGN a(ox) = -1.8 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')

#AGN MODEL FOR alpha_OX = -2.0, efrac = 2% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==7 and efrac == 1:
   file_lib = 'C17_AGN_alpha20_efrac02_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -2.0 and free electron fraction = 2%. No with grains interpolation.'
      print ('No interpolation for the AGN a(ox) = -2.0 with 2% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -2.0 and free electron fraction = 2% with grains interpolated.'
      print ('Interpolation for the AGN a(ox) = -2.0 with 2% free electrons models with grains is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#AGN MODEL FOR alpha_OX = -2.0, efrac = 98% and logU in [-4.0, -0.5]
elif sed==3 and alpha ==7 and efrac == 2:
   file_lib = 'C17_AGN_alpha20_efrac98_CNfix_grains_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'Double composite AGN, a(OX) = -2.0 and free electron fraction = 98% with grains. No interpolation.'
      print ('No interpolation for the AGN a(ox) = -2.0 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'Double composite AGN, a(OX) = -2.0 and free electron fraction = 98% with grains interpolated.'
      print ('Interpolation for the AGN a(ox) = -2.0 with 98% free electrons with grains models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')

#pAGB MODEL FOR T_eff = 50 kK, efrac = 2% and logU in [-4.0, -0.5]
elif sed==4 and teff_pagb==1 and efrac == 1:
   file_lib = 'C17_PAGB_T50kk_efrac02_CNfix_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'pAGB model T_eff = 50 kK and free electron fraction = 2%. No interpolation.'
      print ('No interpolation for the pAGB model T_eff = 50 kK with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'pAGB model T_eff = 50 kK and free electron fraction = 2% interpolated.'
      print ('Interpolation for the pAGB model T_eff = 50 kK with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#pAGB MODEL FOR T_eff = 50 kK, efrac = 98% and logU in [-4.0, -0.5]
elif sed==4 and teff_pagb==1 and efrac == 2:
   file_lib = 'C17_PAGB_T50kk_efrac98_CNfix_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'pAGB model T_eff = 50 kK and free electron fraction = 98%. No interpolation.'
      print ('No interpolation for the pAGB model T_eff = 50 kK with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'pAGB model T_eff = 50 kK and free electron fraction = 98% interpolated.'
      print ('Interpolation for the pAGB model T_eff = 50 kK with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')

#pAGB MODEL FOR T_eff = 100 kK, efrac = 2% and logU in [-4.0, -0.5]
elif sed==4 and teff_pagb==2 and efrac == 1:
   file_lib = 'C17_PAGB_T100kk_efrac02_CNfix_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'pAGB model T_eff = 100 kK and free electron fraction = 2%. No interpolation.'
      print ('No interpolation for the pAGB model T_eff = 100 kK with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'pAGB model T_eff = 100 kK and free electron fraction = 2% interpolated.'
      print ('Interpolation for the pAGB model T_eff = 100 kK with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#pAGB MODEL FOR T_eff = 100 kK, efrac = 98% and logU in [-4.0, -0.5]
elif sed==4 and teff_pagb==2 and efrac == 2:
   file_lib = 'C17_PAGB_T100kk_efrac98_CNfix_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'pAGB model T_eff = 100 kK and free electron fraction = 98%. No interpolation.'
      print ('No interpolation for the pAGB model T_eff = 100 kK with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'pAGB model T_eff = 100 kK and free electron fraction = 98% interpolated.'
      print ('Interpolation for the pAGB model T_eff = 100 kK with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')

#pAGB MODEL FOR T_eff = 150 kK, efrac = 2% and logU in [-4.0, -0.5]
elif sed==4 and teff_pagb==3 and efrac == 1:
   file_lib = 'C17_PAGB_T150kk_efrac02_CNfix_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'pAGB model T_eff = 150 kK and free electron fraction = 2%. No interpolation.'
      print ('No interpolation for the pAGB model T_eff = 150 kK with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'pAGB model T_eff = 150 kK and free electron fraction = 2% interpolated.'
      print ('Interpolation for the pAGB model T_eff = 150 kK with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#pAGB MODEL FOR T_eff = 150 kK, efrac = 98% and logU in [-4.0, -0.5]
elif sed==4 and teff_pagb==3 and efrac == 2:
   file_lib = 'C17_PAGB_T150kk_efrac98_CNfix_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'pAGB model T_eff = 150 kK and free electron fraction = 98%. No interpolation.'
      print ('No interpolation for the pAGB model T_eff = 150 kK with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'pAGB model T_eff = 150 kK and free electron fraction = 98% interpolated.'
      print ('Interpolation for the pAGB model T_eff = 150 kK with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')

#ADAF model, efrac = 2% and logU in [-4.0, -0.5]
elif sed==5 and efrac == 1:
   file_lib = 'C17_ADAF_efrac02_CNfix_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file10:
      for line in file10:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'ADAF model and free electron fraction = 2%. No interpolation.'
      print ('No interpolation for the ADAF model with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'ADAF model and free electron fraction = 2% interpolated.'
      print ('Interpolation for the ADAF model with 2% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')
      res_CO = 0.125
 
#ADAF model, efrac = 98% and logU in [-4.0, -0.5]
elif sed==5 and efrac == 2:
   file_lib = 'C17_ADAF_efrac98_CNfix_uv.dat'
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib, 'r') as file11:
      for line in file11:
         if line[0] == '#':
            n_comments += 1
   grid_aux = np.genfromtxt('Libraries_uv/'+file_lib,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
   if inter == 0:
      sed_type = 'ADAF model and free electron fraction = 98%. No interpolation.'
      print ('No interpolation for the ADAF model with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.1dex for O/H and 0.125dex for C/O')
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'ADAF model and free electron fraction = 98% interpolated.'
      print ('Interpolation for the ADAF model with 98% free electrons models is going to be used.')
      print ('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for C/O')


#Different library
elif sed==6:
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
      res_CO = 0.125
   elif inter == 1:
      sed_type = 'User file ' + new_library + ' used as library for the models interpolated'
      print ('Interpolation for the library '+new_library)
      res_CO = 0.125

#Valuable columns of the files
uv_lin = ['12logOH', 'logCO', 'logU', 'Lya_1216', 'NV_1239', 'CIV_1549', 'HeII_1640', 'OIII_1665', 'CIII_1909', 'OII_3727', 'NeIII_3869', 'OIII_4363', 'OIII_5007']
lin_uv_label = ['12+log(O/H)', 'log(C/O)', 'log(U)', 'Lya_1216', 'NV_1239', 'CIV_1549', 'HeII_1640', 'OIII_1665', 'CIII_1909', 'OII_3727', 'NeIII_3869', 'OIII_4363', 'OIII_5007']

########################################
###### SORTING THE GRID OF MODELS ######
########################################

print (' ')
print ('Sorting the grid of models')
print (' ')

index_OH_CO_U_sorted = [] #storing the correct order of the indexes

#Sorting abundances 12+log(O/H)
OH_values = grid_aux['12logOH'] #Oxygen abundances
if len(OH_values) != 1:
   sorted_list_OH = sorted(range(len(OH_values)),key=OH_values.__getitem__)
if len(OH_values) == 1:
   sorted_list_OH = [0]

#Sorting abundance ratios log(C/O)
OH_values_diff = list(set(OH_values[sorted_list_OH]))
OH_values_diff.sort() #It is necessary to sort again the list of different elements
for OH_num in OH_values_diff:
   index_OH_fix = np.where(OH_values == OH_num)[0] #Index(es) for a particular abundance 12+log(O/H)
   CO_values = grid_aux['logCO'][index_OH_fix]
   if len(CO_values) != 1:
      sorted_list_CO = sorted(range(len(CO_values)), key=CO_values.__getitem__)
   if len(CO_values) == 1:
      sorted_list_CO = [0]
   CO_values_diff = list(set(CO_values[sorted_list_CO]))
   CO_values_diff.sort() #It s necessary to sort again the list of different elements
   for CO_num in CO_values_diff:
      index_OH_CO_fix = np.where(CO_values == CO_num)[0] #Index(es) for particular abundances 12+log(O/H) and log(C/O)
      #Sorting ionization parameters
      U_values = grid_aux['logU'][index_OH_fix[index_OH_CO_fix]]
      if len(U_values) != 1:
         sorted_list_U = sorted(range(len(U_values)), key=U_values.__getitem__)
      if len(U_values) == 1:
         sorted_list_U = [0]
      index_OH_CO_U = index_OH_fix[index_OH_CO_fix[sorted_list_U]] #Sorted index(es) for U at fixed O/H and C/O
      for index_sort in index_OH_CO_U:
         index_OH_CO_U_sorted.append(index_sort) #Adding index in the correct order

#Generating new library file
list_comments = [] #Storing comments in the file:
with open('Libraries_uv/'+file_lib, 'r') as file_aux:
   for line in file_aux:
      if line[0] == '#':
         list_comments.append(line)

#Storing columns:
lin_uv_col = []
#Retrieving each column of the grid
for label in uv_lin:
   aux_col = grid_aux[label].tolist()
   lin_uv_col.append(aux_col)

#Comments
grid_to_write = open('Libraries_uv/'+file_lib, 'w')
for line_com in list_comments:
   grid_to_write.write(line_com)
#Header line
label_line = '{:15} '.format(lin_uv_label[0].replace(' ',''))
for ind in range(1, len(lin_uv_label)-1):
   label_line += '\t {:15} '.format(lin_uv_label[ind].replace(' ',''))
label_line += '\t {:15}\n'.format(lin_uv_label[-1].replace(' ','')) 
grid_to_write.write(label_line)
#Values:
for ind_val in index_OH_CO_U_sorted:
   val_line = '{:7.7f} '.format(lin_uv_col[0][ind_val])
   for ind2 in range(1, len(lin_uv_label)-1):
      val_line += '\t {:7.7f} '.format(lin_uv_col[ind2][ind_val])
   val_line += '\t {:7.7f}\n'.format(lin_uv_col[-1][ind_val])
   grid_to_write.write(val_line)        
grid_to_write.close()

#Opening sorted grid of models
n_comments = 0
with open('Libraries_uv/'+file_lib, 'r') as file12:
   for line in file12:
      if line[0] == '#':
         n_comments += 1  
grid_aux = np.genfromtxt('Libraries_uv/'+file_lib, dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

################################################
###### CONSTRAINTS FOR THE GRID OF MODELS ######
################################################

#Reading constraints and creating library with constraints
print (' ')
print ('Select a file with the constraint laws to be used to limit the grid of models when the measurement of a quantity is impossible without any relation.')
print (' ')

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

#Particular file introduced by the user
if const == '6':
   while question8:
      print ('Introduce name of the file containing the constraints for the grids. It must be located in the folder "Constraints".')
      print (' ')
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
         print (' ')
         print ('Loading constraint file '+new_const+'. Checking correct format of the file.')
         question8 = False
      except:
         print (' ')
         print ('File was not found in folder "Constraints" or file does not exist.')
   question9 = True
   while question9:
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
         print ('Succesfully reading of the file')
         if len(missing_labels) == 0:
            print ('File presents the correct format')
            question9 = False
         else:
            print ('File does not present the correct format. The following columns are missing:')
            for need_label in missing_labels:
               print('- '+need_label)
            print ('More details on the correct format for the library are found in readme file.')
            print (' ')
            print ('Reintroduce name of the file with fixed format:')
            print (' ')
            if int(sys.version[0]) < 3:
               new_const = raw_input('Name of file: ')
            else:
               new_const = input('Name of file: ')
      except:
         print ('Something went wrong while reading file. Please, reintroduce name of the file:')
         print (' ')
         if int(sys.version[0]) < 3:
            new_const = raw_input('Name of file: ')
         else:
            new_const = input('Name of file: ')

const = int(const)

#Generation of grids with constraints laws:
if const == 1 or const == 2 or const == 3 or const == 6:
   #First grid does not change
   grid1 = grid_aux
   file_lib_2 = file_lib
elif const == 4 or const == 5:
   lin_uv_agn = []
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
   for label in uv_lin:
      aux_col = grid_aux[label].tolist()
      lin_uv_agn.append(aux_col)
   #Creation of the grid
   file_lib_2 = '.'.join(file_lib.split('.')[0:-1])+'_'+tag+'.'+file_lib.split('.')[-1]
   file_open = open('Libraries_uv/'+file_lib_2, 'w')
   file_open.write('#Library constrained for '+tag+' ionization AGNs\n')
   #Header line
   label_line = '{:15} '.format(lin_uv_label[0].replace(' ',''))
   for ind in range(1, len(lin_uv_label)-1):
      label_line += '\t {:15} '.format(lin_uv_label[ind].replace(' ',''))
   label_line += '\t {:15}\n'.format(lin_uv_label[-1].replace(' ','')) 
   file_open.write(label_line)
   #Values:
   for ind_val in range(0, len(lin_uv_agn[0])):
      if lin_uv_agn[2][ind_val] <= U_max and lin_uv_agn[2][ind_val] >= U_min:
         val_line = '{:7.7f} '.format(lin_uv_agn[0][ind_val])
         for ind2 in range(1, len(lin_uv_label)-1):
            val_line += '\t {:7.7f} '.format(lin_uv_agn[ind2][ind_val])
         val_line += '\t {:7.7f}\n'.format(lin_uv_agn[-1][ind_val])
         file_open.write(val_line)        
   file_open.close()
   #Counting comments:
   n_comments = 0
   with open('Libraries_uv/'+file_lib_2, 'r') as file:
      for line in file:
         if line[0] == '#':
            n_comments += 1
   grid1 = np.genfromtxt('Libraries_uv/'+file_lib_2,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)


#Generating libraries for the constraints in the files
if const == 1: #Star-Forming Galaxies
   const_file = 'template_OH.dat'
   name_const = 'Constraints/template_OH.dat'
   n_comments = 0
   with open(name_const, 'r') as file12:
      for line in file12:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)
if const == 2:
   const_file = 'template_OH_eelg.dat'
   name_const = 'Constraints/template_OH_eelg.dat'
   n_comments = 0
   with open(name_const, 'r') as file13:
      for line in file13:
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
   with open(name_const, 'r') as file14:
      for line in file14:
         if line[0] == '#':
            n_comments += 1
   const_data = np.genfromtxt(name_const,dtype=None,names=True, encoding = 'ascii', skip_header=n_comments)

#Limiting the grids:
lin_uv_val = []
#The initial grid need to be constrained in the ionization parameter

#Retrieving each column of the grid
for label in uv_lin:
   aux_col = grid1[label].tolist()
   lin_uv_val.append(aux_col)
   #Creation of the grids
   name_OH_U = '.'.join(file_lib_2.split('.')[0:-1])+'_OH_U_constrained.'+file_lib.split('.')[-1]
   name_OH_U_CO = '.'.join(file_lib_2.split('.')[0:-1])+'_OH_U_CO_constrained.'+file_lib.split('.')[-1]
   file_open = open('Libraries_uv/'+ name_OH_U, 'w') #OH and U relation
   file_open_2 = open('Libraries_uv/'+name_OH_U_CO, 'w') #OH, CO and U relation
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

###########################################
###### SUMMARY OF THE GRID OF MODELS ######
###########################################
  
print ('-------------------------------------------------')
print ('Summary of the models')
print ('---------------------')
print ('Libraries generated with the constraints. The following grids are going to be used:')
print ('- Full library (Grid#1): '+file_lib_2)
print ('       Total number of models: ' + str(len(grid1)))
print ('- Library constrained by 12+log(O/H) - log(U) relation (Grid#2): '+name_OH_U)
print ('       Total number of models: ' + str(len(grid2)))
print ('- Library constrained by 12+log(O/H) - log(U) - log(C/O) relation (Grid#3): '+name_OH_U_CO)
print ('       Total number of models: ' + str(len(grid3)))
print ('-------------------------------------------------')
print (' ')

#################################################
###### CREATING ARRAY TO STORE ESTIMATIONS ######
#################################################

grids = []
OHffs = []
eOHffs = []
COffs = []
eCOffs = []
logUffs = []
elogUffs = []

#Labels to check information provided in the input file

Label_ID = False
Label_Lya = False
Label_eLya = False
Label_NV = False
Label_eNV = False
Label_CIV = False
Label_eCIV = False
Label_HeII = False
Label_eHeII = False
Label_OIII_1665 = False
Label_eOIII_1665 = False
Label_CIII = False
Label_eCIII = False
Label_OII_3727 = False
Label_eOII_3727 = False
Label_NeIII_3869 = False
Label_eNeIII_3869 = False
Label_OIII_4363 = False
Label_eOIII_4363 = False
Label_OIII_5007 = False
Label_eOIII_4959 = False
Label_OIII_4959 = False
Label_eOIII_5007 = False
Label_Hbeta = False
Label_eHbeta = False

#Checking input information
for col in range(0,len(input1.dtype.names),1):
   if input1.dtype.names[col] == 'ID':
      Label_ID = True
   if input1.dtype.names[col] == 'Lya_1216':
      Label_Lya = True
   if input1.dtype.names[col] == 'eLya_1216':
      Label_eLya = True
   if input1.dtype.names[col] == 'NV_1239':
      Label_NV = True
   if input1.dtype.names[col] == 'eNV_1239':
      Label_eNV = True
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
   if input1.dtype.names[col] == 'OII_3727':
      Label_OII_3727 = True
   if input1.dtype.names[col] == 'eOII_3727':
      Label_eOII_3727 = True
   if input1.dtype.names[col] == 'NeIII_3869':
      Label_NeIII_3869 = True
   if input1.dtype.names[col] == 'eNeIII_3869':
      Label_eNeIII_3869 = True
   if input1.dtype.names[col] == 'OIII_4363':
      Label_OIII_4363 = True
   if input1.dtype.names[col] == 'eOIII_4363':
      Label_eOIII_4363 = True
   if input1.dtype.names[col] == 'OIII_4959':
      Label_OIII_4959 = True
   if input1.dtype.names[col] == 'eOIII_4959':
      Label_eOIII_4959 = True
   if input1.dtype.names[col] == 'Hb_4861':
      Label_Hbeta = True
   if input1.dtype.names[col] == 'eHb_4861':
      Label_eHbeta = True
   if input1.dtype.names[col] == 'OIII_5007':
      Label_OIII_5007 = True
   if input1.dtype.names[col] == 'eOIII_5007':
      Label_eOIII_5007 = True

#Adapting final output with information from given input
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
if Label_NV == False:
   NV_1239 = np.zeros(input1.size)
else:
   NV_1239 = input1['NV_1239']
if Label_eNV == False:
   eNV_1239 = np.zeros(input1.size)
else:
   eNV_1239 = input1['eNV_1239']
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
if Label_OII_3727 == False:
   OII_3727 = np.zeros(input1.size)
else:
   OII_3727 = input1['OII_3727']
if Label_eOII_3727 == False:
   eOII_3727 = np.zeros(input1.size)
else:
   eOII_3727 = input1['eOII_3727']
if Label_NeIII_3869 == False:
   NeIII_3869 = np.zeros(input1.size)
else:
   NeIII_3869 = input1['NeIII_3869']
if Label_eNeIII_3869 == False:
   eNeIII_3869 = np.zeros(input1.size)
else:
   eNeIII_3869 = input1['eNeIII_3869']
if Label_OIII_4363 == False:
   OIII_4363 = np.zeros(input1.size)
else:
   OIII_4363 = input1['OIII_4363']
if Label_eOIII_4363 == False:
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

################################################################
###### OUTPUT FORMAT AND INFORMATION: ONLY EMISSION LINES ######
################################################################

#Creation of output only with information from inputs
aux_list = []
aux_list.append(('ID','U12'))
if Label_Lya == True:
   aux_list.append(('Lya_1216', float))
if Label_eLya == True:
   aux_list.append(('eLya_1216', float))
if Label_NV == True:
   aux_list.append(('NV_1239', float))
if Label_eNV == True:
   aux_list.append(('eNV_1239', float))
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
if Label_OII_3727 == True:
   aux_list.append(('OII_3727', float))
if Label_eOII_3727 == True:
   aux_list.append(('eOII_3727', float))
if Label_NeIII_3869 == True:
   aux_list.append(('NeIII_3869', float))
if Label_eNeIII_3869 == True:
   aux_list.append(('eNeIII_3869', float))
if Label_OIII_4363 == True:
   aux_list.append(('OIII_4363', float))
if Label_eOIII_4363 == True:
   aux_list.append(('eOIII_4363', float))
if Label_OIII_4959 == True:
   aux_list.append(('OIII_4959', float))
if Label_eOIII_4959 == True:
   aux_list.append(('eOIII_4959', float))
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

output['ID'] = Names
if Label_Lya == True:
   output['Lya_1216'] = Lya_1216
if Label_eLya == True:
   output['eLya_1216'] = eLya_1216
if Label_NV == True:
   output['NV_1239'] = NV_1239
if Label_eNV == True:
   output['eNV_1239'] = eNV_1239
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
   output['eCIII_1909'] = eCIII_1909
if Label_OII_3727 == True:
   output['OII_3727'] = OII_3727
if Label_eOII_3727 == True:
   output['eOII_3727'] = eOII_3727
if Label_NeIII_3869 == True:
   output['NeIII_3869'] = NeIII_3869
if Label_eNeIII_3869 == True:
   output['eNeIII_3869'] = eNeIII_3869
if Label_OIII_4363 == True:
   output['OIII_4363'] = OIII_4363
if Label_eOIII_4363 == True:
   output['eOIII_4363'] = eOIII_4363
if Label_OIII_4959 == True:
   output['OIII_4959'] = OIII_4959
if Label_eOIII_4959 == True:
   output['eOIII_4959'] = eOIII_4959
if Label_Hbeta == True:
   output['Hb_4861'] = Hb_4861
if Label_eHbeta == True:
   output['eHb_4861'] = eHb_4861
if Label_OIII_5007 == True:
   output['OIII_5007'] = OIII_5007
if Label_eOIII_5007 == True:
   output['eOIII_5007'] = eOIII_5007

################################################
###### ESTIMATIONS OF CHEMICAL ABUNDANCES ######
################################################

#Display for the user
print ('Calculating....')
print ('')
print ('')
print ('----------------------------------------------------------------')
print ('(%)   ID    Grid  12+log(O/H)  log(C/O)    log(U)')
print ('----------------------------------------------------------------')

# Beginning of loop of calculation
count = 0
for tab in range(0,len(input1),1):
   try:
      count = count + 1
      OH_mc = []
      CO_mc = []
      logU_mc = []
      OHe_mc = []
      COe_mc = []
      logUe_mc = []  


      #Starting Montecarlo
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
         tol_max = 1e3
         
         #Generating observable values for emission lines
         Lya_1216_obs = 0
         if Lya_1216[tab] <= 0:
            Lya_1216_obs = 0
         else:
            while Lya_1216_obs <= 0:
               Lya_1216_obs = np.random.normal(Lya_1216[tab],eLya_1216[tab]+1e-5)
         NV_1239_obs = 0
         if NV_1239[tab]<= 0:
            NV_1239_obs = 0
         else:
            while NV_1239_obs <= 0:
               NV_1239_obs = np.random.normal(NV_1239[tab],eNV_1239[tab]+1e-5)
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
         OII_3727_obs = 0
         if OII_3727[tab] <= 0:
            OII_3727_obs = 0
         else:
            while OII_3727_obs <= 0:
               OII_3727_obs = np.random.normal(OII_3727[tab],eOII_3727[tab]+1e-5)
         NeIII_3869_obs = 0
         if NeIII_3869[tab] <= 0:
            NeIII_3869_obs = 0
         else:
            while NeIII_3869_obs <= 0:
               NeIII_3869_obs = np.random.normal(NeIII_3869[tab],eNeIII_3869[tab]+1e-5)
         OIII_4363_obs = 0
         if OIII_4363[tab] <= 0:
            OIII_4363_obs = 0
         else:
            while OIII_4363_obs <= 0:
               OIII_4363_obs = np.random.normal(OIII_4363[tab],eOIII_4363[tab]+1e-5)
         Hb_4861_obs = 0
         if Hb_4861[tab] <= 0:
            Hb_4861_obs = 0
         else:
            while Hb_4861_obs <= 0:
               Hb_4861_obs = np.random.normal(Hb_4861[tab],eHb_4861[tab]+1e-5)
         OIII_5007_obs = 0
         if OIII_4959[tab] <= 0 and OIII_5007[tab] <= 0:
            OIII_5007_obs = 0
         elif OIII_4959[tab] <= 0 and OIII_5007[tab] > 0:
            while OIII_5007_obs <= 0:
               OIII_5007_obs = np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-5)
         elif OIII_4959[tab] > 0 and OIII_5007[tab] <= 0:
            while OIII_5007_obs <= 0:
               OIII_5007_obs = 3*np.random.normal(OIII_4959[tab],eOIII_4959[tab]+1e-5)
         else:
            while OIII_5007_obs <= 0:
               OIII_5007_obs = (np.random.normal(OIII_5007[tab],eOIII_5007[tab]+1e-5)+np.random.normal(OIII_4959[tab],eOIII_4959[tab]+1e-5))/1.33

         #Observables
         if OIII_1665_obs <= 0 or OIII_5007_obs <= 0:
            ROIII_obs = 0
         else:
            ROIII_obs = OIII_5007_obs/OIII_1665_obs
         if OIII_1665_obs <= 0 or OIII_4363_obs <= 0:
            R2OIII_obs = 0
         else:
            R2OIII_obs = OIII_4363_obs/OIII_1665_obs
         if OIII_4363_obs <= 0 or OIII_5007_obs <= 0:
            R3OIII_obs = 0
         else:
            R3OIII_obs = OIII_5007_obs/OIII_4363_obs
         if Lya_1216_obs == 0 or NV_1239_obs == 0:
            N5_obs = 0
         else:
            N5_obs = (NV_1239_obs ) / (Lya_1216_obs)
         if HeII_1640_obs == 0 or NV_1239_obs == 0:
            N5He2_obs = 0
         else:
            N5He2_obs = (NV_1239_obs) / (HeII_1640_obs)
         if Lya_1216_obs <= 0 or CIII_1909_obs <= 0 or CIV_1549_obs <= 0:
            C34_obs = 0
         else:
            C34_obs = (CIII_1909_obs + CIV_1549_obs) / (Lya_1216_obs)
         if HeII_1640_obs <= 0 or CIII_1909_obs <= 0  or CIV_1549_obs <= 0:
            C34He2_obs = 0
         else:
            C34He2_obs = (CIII_1909_obs + CIV_1549_obs) / (HeII_1640_obs)
         if CIII_1909_obs <= 0 or OIII_1665_obs <= 0:
            C3O3_obs = -10
         else:   
            C3O3_obs = np.log10((CIII_1909_obs) / (OIII_1665_obs))
         if CIII_1909_obs <= 0 or OIII_5007_obs <= 0:
            C3O5007_obs = -10
         else:   
            C3O5007_obs = np.log10((CIII_1909_obs) / (OIII_5007_obs))
         if CIII_1909_obs <= 0 or NeIII_3869_obs <= 0 :
            C3Ne3_obs = -10
         else:   
            C3Ne3_obs = np.log10((CIII_1909_obs) / (NeIII_3869_obs))
         if CIII_1909_obs <= 0 or CIV_1549_obs <= 0:
            C3C4_obs = 0
         else:
            C3C4_obs = (CIII_1909_obs/CIV_1549_obs)
         if CIII_1909_obs <= 0 or Hb_4861_obs <= 0:
            C34Hb_obs = 0
         else:
            C34Hb_obs = (CIII_1909_obs + CIV_1549_obs) / Hb_4861_obs
         if OII_3727_obs <= 0 or OIII_5007_obs <= 0 or Hb_4861_obs <= 0:
            R23_obs = -10
         else:
            R23_obs = np.log10((OII_3727_obs + OIII_5007_obs) / Hb_4861_obs)
         if OII_3727_obs <= 0 or OIII_5007_obs <= 0:
            O2O3_obs = 0
         else:
            O2O3_obs = (OII_3727_obs / OIII_5007_obs) 
         if OII_3727_obs <= 0 or NeIII_3869_obs <= 0 or Hb_4861_obs <= 0:
            R2Ne3_obs = -10
         else:
            R2Ne3_obs = np.log10((OII_3727_obs + NeIII_3869_obs) / Hb_4861_obs)
         if OII_3727_obs <= 0 or NeIII_3869_obs <= 0:
            O2Ne3_obs = 0
         else:
            O2Ne3_obs = (OII_3727_obs / NeIII_3869_obs)




               
         # Selection of grid
         if (OIII_1665[tab] > 0 and OIII_5007[tab] > 0) or (OIII_1665[tab] > 0 and OIII_4363[tab] > 0) or (OIII_5007[tab] > 0 and OIII_4363[tab] > 0):
            grid = grid1
            if monte == n-1: grids.append(1)
            grid_type = 1
         elif (OIII_1665[tab] > 0 and CIII_1909[tab] > 0) or (CIII_1909[tab] > 0 and NeIII_3869[tab] > 0) or (CIII_1909[tab] > 0 and OIII_5007[tab] > 0):
            grid = grid2
            if monte == n-1: grids.append(2)
            grid_type = 2
         else:
            grid = grid3
            if monte == n-1: grids.append(3)
            grid_type = 3

         ######################
         # Calculation of C/O #
         ######################
         if C3O3_obs == -10 and C3O5007_obs == -10 and C3Ne3_obs == -10:
            CO = -10
         else:
            CHI_ROIII = 0
            CHI_R2OIII = 0
            CHI_R3OIII = 0
            CHI_C3O3 = 0   
            CHI_C3O5007 = 0   
            CHI_C3Ne3 = 0   
            CHI_CO = 0
            for index in grid:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index['OIII_1665'] == 0 or index['OIII_5007'] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index['OIII_5007']/index['OIII_1665'] - ROIII_obs)**2/(index['OIII_5007']/index['OIII_1665'])
               if R2OIII_obs == 0:
                  CHI_R2OIII = 0
               elif index['OIII_1665'] == 0 or index['OIII_4363'] == 0:
                  CHI_R2OIII = tol_max
               else:
                  CHI_R2OIII = (index['OIII_4363']/index['OIII_1665'] - R2OIII_obs)**2/(index['OIII_4363']/index['OIII_1665'])
               if R3OIII_obs == 0:
                  CHI_R3OIII = 0
               elif index['OIII_4363'] == 0 or index['OIII_4363'] == 0:
                  CHI_R3OIII = tol_max
               else:
                  CHI_R3OIII = (index['OIII_5007']/index['OIII_4363'] - R3OIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
               if C3O3_obs == -10:
                  CHI_C3O3 = 0          
               elif index['CIII_1909'] == 0 or index['OIII_1665'] == 0:
                  CHI_C3O3 = tol_max
               else:
                  CHI_C3O3 =(np.log10((index['CIII_1909'])/index['OIII_1665']) - C3O3_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_1665']+1e-5))
               if C3O5007_obs == -10:
                  CHI_C3O5007 = 0          
               elif index['CIII_1909'] == 0 or index['OIII_5007'] == 0:
                  CHI_C3O5007 = tol_max
               else:
                  CHI_C3O5007 =(np.log10((index['CIII_1909'])/index['OIII_5007']) - C3O5007_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_5007']+1e-5))
               if C3Ne3_obs == -10:
                  CHI_C3Ne3 = 0          
               elif index['CIII_1909'] == 0 or index['NeIII_3869'] == 0:
                  CHI_C3Ne3 = tol_max
               else:
                  CHI_C3Ne3 =(np.log10((index['CIII_1909'])/index['NeIII_3869']) - C3Ne3_obs)**2/np.log10((index['CIII_1909'])/(index['NeIII_3869']+1e-5))
               CHI_CO = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 +CHI_C3O3**2 + CHI_C3O5007**2 + CHI_C3Ne3**2)**0.5

               if CHI_CO == 0:
                  CO_p = CO_p
                  den_CO = den_CO
               else:
                  CO_p = index['logCO'] /(CHI_CO)**2 + CO_p
                  den_CO = 1 / (CHI_CO)**2 + den_CO

            CO = CO_p / den_CO 

         # Calculation of C/O error
         if C3O3_obs == -10 and C3O5007_obs == -10 and C3Ne3_obs == -10:
            eCO = 0
         else:
            CHI_ROIII = 0
            CHI_R2OIII = 0
            CHI_R3OIII = 0
            CHI_C3O3 = 0   
            CHI_C3O5007 = 0   
            CHI_C3Ne3 = 0   
            CHI_CO = 0
            for index in grid:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index['OIII_1665'] == 0 or index['OIII_5007'] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index['OIII_5007']/index['OIII_1665'] - ROIII_obs)**2/(index['OIII_5007']/index['OIII_1665'])
               if R2OIII_obs == 0:
                  CHI_R2OIII = 0
               elif index['OIII_1665'] == 0 or index['OIII_4363'] == 0:
                  CHI_R2OIII = tol_max
               else:
                  CHI_R2OIII = (index['OIII_4363']/index['OIII_1665'] - R2OIII_obs)**2/(index['OIII_4363']/index['OIII_1665'])
               if R3OIII_obs == 0:
                  CHI_R3OIII = 0
               elif index['OIII_4363'] == 0 or index['OIII_4363'] == 0:
                  CHI_R3OIII = tol_max
               else:
                  CHI_R3OIII = (index['OIII_5007']/index['OIII_4363'] - R3OIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
               if C3O3_obs == -10:
                  CHI_C3O3 = 0          
               elif index['CIII_1909'] == 0 or index['OIII_1665'] == 0:
                  CHI_C3O3 = tol_max
               else:
                  CHI_C3O3 =(np.log10((index['CIII_1909'])/index['OIII_1665']) - C3O3_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_1665']+1e-5))
               if C3O5007_obs == -10:
                  CHI_C3O5007 = 0          
               elif index['CIII_1909'] == 0 or index['OIII_5007'] == 0:
                  CHI_C3O5007 = tol_max
               else:
                  CHI_C3O5007 =(np.log10((index['CIII_1909'])/index['OIII_5007']) - C3O5007_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_5007']+1e-5))
               if C3Ne3_obs == -10:
                  CHI_C3Ne3 = 0          
               elif index['CIII_1909'] == 0 or index['NeIII_3869'] == 0:
                  CHI_C3Ne3 = tol_max
               else:
                  CHI_C3Ne3 =(np.log10((index['CIII_1909'])/index['NeIII_3869']) - C3Ne3_obs)**2/np.log10((index['CIII_1909'])/(index['NeIII_3869']+1e-5))
               CHI_CO = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 +CHI_C3O3**2 + CHI_C3O5007**2 + CHI_C3Ne3**2)**0.5


               if CHI_CO == 0:
                  CO_e = CO_e
                  den_CO_e = den_CO_e  
               else:
                  CO_e = (index['logCO'] - CO)**2 / (CHI_CO)**2 + CO_e
                  den_CO_e = 1 /(CHI_CO)**2 + den_CO_e  


            eCO = CO_e / den_CO_e 

         
         ###############################
         # Calculation of O/H and logU #
         ###############################
         if C34_obs == 0 and ROIII_obs == 0 and R2OIII_obs == 0 and R3OIII_obs == 0 and C34Hb_obs == 0 and C34He2_obs == 0 and N5_obs == 0 and N5He2_obs == 0 and R23_obs == -10 and O2O3_obs == 0 and O2Ne3_obs == 0:
            OH = 0
            logU = 0
         else:
            CHI_ROIII = 0 
            CHI_R2OIII = 0 
            CHI_R3OIII = 0 
            CHI_C3C4 = 0
            CHI_C34He2 = 0
            CHI_C34 = 0
            CHI_C34Hb = 0
            CHI_N5 = 0
            CHI_N5He2 = 0
            CHI_R23 = 0
            CHI_O2O3 = 0
            CHI_R2Ne3 = 0
            CHI_O2Ne3 = 0
            CHI_OH = 0
            for index in grid:
               if CO > -10 and np.abs(index['logCO'] - CO) > np.abs(eCO+0.125):
                  continue
               if NV_1239_obs > 0 and index['NV_1239'] == 0:
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
                  if R2OIII_obs == 0:
                     CHI_R2OIII = 0
                  elif index['OIII_1665'] == 0 or index['OIII_4363'] == 0:
                     CHI_R2OIII = tol_max
                  else:
                     CHI_R2OIII = (index['OIII_4363']/index['OIII_1665'] - R2OIII_obs)**2/(index['OIII_4363']/index['OIII_1665'])
                  if R3OIII_obs == 0:
                     CHI_R3OIII = 0
                  elif index['OIII_4363'] == 0 or index['OIII_5007'] == 0:
                     CHI_R3OIII = tol_max
                  else:
                     CHI_R3OIII = (index['OIII_5007']/index['OIII_4363'] - R3OIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
                  if N5_obs == 0:
                     CHI_N5 = 0
                  elif index['Lya_1216'] == 0 or index['NV_1239'] == 0:
                     CHI_N5 = tol_max
                  else:
                     CHI_N5 = ((index['NV_1239'])/index['Lya_1216'] - N5_obs)**2/((index['NV_1239'])/index['Lya_1216'])
                  if N5He2_obs == 0:
                     CHI_N5He2 = 0
                  elif index['HeII_1640'] == 0 or index['NV_1239'] == 0:
                     CHI_N5He2 = tol_max
                  else:
                     CHI_N5He2 = ((index['NV_1239'])/index['HeII_1640'] - N5He2_obs)**2/((index['NV_1239'])/index['HeII_1640'])
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
                  if R23_obs == -10:
                     CHI_R23 = 0
                  elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                     CHI_R23 = tol_max
                  else:
                     CHI_R23 = (np.log10(index['OII_3727']+index['OIII_5007']) - R23_obs)**2/np.log10(index['OII_3727']+index['OIII_5007'])
                  if O2O3_obs == 0:
                     CHI_O2O3 = 0
                  elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                     CHI_O2O3 = tol_max
                  else:
                     CHI_O2O3 = ((index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/(index['OII_3727']/index['OIII_5007'])
                  if R2Ne3_obs == -10:
                     CHI_R2Ne3 = 0
                  elif index['OII_3727'] == 0 or index['NeIII_3869'] == 0:
                     CHI_R2Ne3 = tol_max
                  else:
                     CHI_R2Ne3 = (np.log10(index['OII_3727']+index['NeIII_3869']) - R2Ne3_obs)**2/np.log10(index['OII_3727']+index['NeIII_3869'])
                  if O2Ne3_obs == 0:
                     CHI_O2Ne3 = 0
                  elif index['OII_3727'] == 0 or index['NeIII_3869'] == 0:
                     CHI_O2Ne3 = tol_max
                  else:
                     CHI_O2Ne3 = ((index['OII_3727']/index['NeIII_3869']) - O2Ne3_obs)**2/(index['OII_3727']/index['NeIII_3869'])
                  if R23_obs > -10:
                     CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34Hb**2  + CHI_C3C4**2 + CHI_R23**2 + CHI_O2O3**2)**0.5
                  elif R23_obs == -10 and R2Ne3_obs > -10:
                     CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34Hb**2  + CHI_C3C4**2 +  CHI_R2Ne3**2 + CHI_O2Ne3**2)**0.5
                  else:	
                     CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_N5**2 + CHI_N5He2**2 + CHI_C3C4**2 + CHI_O2Ne3**2)**0.5

                  if CHI_OH == 0:
                     OH_p = OH_p
                     logU_p = logU_p
                     den_OH = den_OH
                  else:
                     OH_p = index['12logOH'] / (CHI_OH)**2 + OH_p
                     logU_p = index['logU'] / (CHI_OH)**2 + logU_p
                     den_OH = 1 /(CHI_OH)**2 + den_OH

            if OH_p == 0:
               OH = 0
            else:
               OH = OH_p / den_OH
            if logU_p == 0:
               logU = 0
            else:
               logU = logU_p / den_OH

      #Impossibility for AGN in the estimation
            if sed == 3 and Lya_1216[tab] == 0 and HeII_1640[tab] == 0 and Hb_4861[tab] == 0:
               OH = 0

         # Calculation of error of O/H and logU
         if C34_obs == 0 and ROIII_obs == 0 and R2OIII_obs == 0 and R3OIII_obs == 0 and  C34Hb_obs == 0  and C34He2_obs == 0 and N5_obs == 0 and N5He2_obs == 0 and O2O3_obs == 0 and O2Ne3_obs == 0:
            eOH = 0
            elogU = 0
         else:
            CHI_ROIII = 0 
            CHI_R2OIII = 0
            CHI_R3OIII = 0
            CHI_N5 = 0
            CHI_N5He2 = 0
            CHI_C3C4 = 0
            CHI_C34 = 0
            CHI_C34He2 = 0
            CHI_C34Hb = 0
            CHI_R23 = 0
            CHI_O2O3 = 0
            CHI_R2Ne3 = 0
            CHI_O2Ne3 = 0
            CHI_OH = 0
            for index in grid:
               if CO > -10 and np.abs(index['logCO'] - CO) > np.abs(eCO+res_CO):
                  continue
               if NV_1239_obs > 0 and index['NV_1239'] == 0:
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
                  if R2OIII_obs == 0:
                     CHI_R2OIII = 0
                  elif index['OIII_1665'] == 0 or index['OIII_4363'] == 0:
                     CHI_R2OIII = tol_max
                  else:
                     CHI_R2OIII = (index['OIII_4363']/index['OIII_1665'] - R2OIII_obs)**2/(index['OIII_4363']/index['OIII_1665'])
                  if R3OIII_obs == 0:
                     CHI_R3OIII = 0
                  elif index['OIII_4363'] == 0 or index['OIII_5007'] == 0:
                     CHI_R3OIII = tol_max
                  else:
                     CHI_R3OIII = (index['OIII_5007']/index['OIII_4363'] - R3OIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
                  if N5_obs == 0:
                     CHI_N5 = 0
                  elif index['Lya_1216'] == 0 or index['NV_1239'] == 0:
                     CHI_N5 = tol_max
                  else:
                     CHI_N5 = ((index['NV_1239'])/index['Lya_1216'] - N5_obs)**2/((index['NV_1239'])/index['Lya_1216'])
                  if N5He2_obs == 0:
                     CHI_N5He2 = 0
                  elif index['HeII_1640'] == 0 or index['NV_1239'] == 0:
                     CHI_N5He2 = tol_max
                  else:
                     CHI_N5He2 = ((index['NV_1239'])/index['HeII_1640'] - N5He2_obs)**2/((index['NV_1239'])/index['HeII_1640'])
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
                  if R23_obs == -10:
                     CHI_R23 = 0
                  elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                     CHI_R23 = tol_max
                  else:
                     CHI_R23 = (np.log10(index['OII_3727']+index['OIII_5007']) - R23_obs)**2/np.log10(index['OII_3727']+index['OIII_5007'])
                  if O2O3_obs == 0:
                     CHI_O2O3 = 0
                  elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                     CHI_O2O3 = tol_max
                  else:
                     CHI_O2O3 = ((index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/(index['OII_3727']/index['OIII_5007'])
                  if R2Ne3_obs == -10:
                     CHI_R2Ne3 = 0
                  elif index['OII_3727'] == 0 or index['NeIII_3869'] == 0:
                     CHI_R2Ne3 = tol_max
                  else:
                     CHI_R2Ne3 = (np.log10(index['OII_3727']+index['NeIII_3869']) - R2Ne3_obs)**2/np.log10(index['OII_3727']+index['NeIII_3869'])
                  if O2Ne3_obs == 0:
                     CHI_O2Ne3 = 0
                  elif index['OII_3727'] == 0 or index['NeIII_3869'] == 0:
                     CHI_O2Ne3 = tol_max
                  else:
                     CHI_O2Ne3 = ((index['OII_3727']/index['NeIII_3869']) - O2Ne3_obs)**2/(index['OII_3727']/index['NeIII_3869'])

                  if R23_obs > -10:
                     CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34Hb**2  + CHI_C3C4**2 + CHI_R23**2 + CHI_O2O3**2)**0.5
                  elif R23_obs == -10 and R2Ne3_obs > -10:
                     CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34Hb**2  + CHI_C3C4**2 +  CHI_R2Ne3**2 + CHI_O2Ne3**2)**0.5
                  else:	
                     CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_N5**2 + CHI_N5He2**2 + CHI_C3C4**2 + CHI_O2Ne3**2)**0.5


               if CHI_OH == 0:
                  OH_e = OH_e
                  logU_e = logU_e
                  den_OH_e = den_OH_e
               else:
                  OH_e = (index['12logOH'] - OH)**2 /(CHI_OH)**2 + OH_e
                  logU_e = (index['logU'] - logU)**2 /(CHI_OH)**2 + logU_e
                  den_OH_e = 1 /(CHI_OH)**2 + den_OH_e 

            if OH_e == 0:
               eOH = 0
            else:
               eOH = OH_e / den_OH_e
            if logU_e == 0:
               elogU = 0
            else:
               elogU = logU_e / den_OH_e

            #Impossiiblity in AGNs to determine O/H without recombination lines
            if sed == 3 and Lya_1216[tab] == 0 and HeII_1640[tab] == 0 and Hb_4861[tab] == 0:
               eOH = 0

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
            CHI_R2OIII = 0
            CHI_R3OIII = 0
            CHI_C3O3 = 0   
            CHI_C3O5007 = 0   
            CHI_C3Ne3 = 0   
            CHI_C3C4 = 0
            CHI_N5 = 0
            CHI_N5He2 = 0
            CHI_C34He2 = 0
            CHI_C34 = 0
            CHI_C34Hb = 0
            CHI_R23 = 0
            CHI_O2O3 = 0
            CHI_R2Ne3 = 0
            CHI_O2Ne3 = 0
            CHI_OH = 0
            CHI_CO = 0


            for index in igrid:
               if ROIII_obs == 0:
                  CHI_ROIII = 0
               elif index['OIII_1665'] == 0 or index['OIII_5007'] == 0:
                  CHI_ROIII = tol_max
               else:
                  CHI_ROIII = (index['OIII_5007']/index['OIII_1665'] - ROIII_obs)**2/(index['OIII_5007']/index['OIII_1665'])
               if R2OIII_obs == 0:
                  CHI_R2OIII = 0
               elif index['OIII_1665'] == 0 or index['OIII_4363'] == 0:
                  CHI_R2OIII = tol_max
               else:
                  CHI_R2OIII = (index['OIII_4363']/index['OIII_1665'] - R2OIII_obs)**2/(index['OIII_4363']/index['OIII_1665'])
               if R3OIII_obs == 0:
                  CHI_R3OIII = 0
               elif index['OIII_4363'] == 0 or index['OIII_5007'] == 0:
                  CHI_R3OIII = tol_max
               else:
                  CHI_R3OIII = (index['OIII_5007']/index['OIII_4363'] - R3OIII_obs)**2/(index['OIII_5007']/index['OIII_4363'])
               if N5_obs == 0:
                  CHI_N5 = 0
               elif index['Lya_1216'] == 0 or index['NV_1239'] == 0:
                  CHI_N5 = tol_max
               else:
                  CHI_N5 = ((index['NV_1239'])/index['Lya_1216'] - N5_obs)**2/((index['NV_1239'])/index['Lya_1216'])
               if N5He2_obs == 0:
                  CHI_N5He2 = 0
               elif index['HeII_1640'] == 0 or index['NV_1239'] == 0:
                  CHI_N5He2 = tol_max
               else:
                  CHI_N5He2 = ((index['NV_1239'])/index['HeII_1640'] - N5He2_obs)**2/((index['NV_1239'])/index['HeII_1640'])
               if C3O3_obs == -10:
                  CHI_C3O3 = 0          
               elif index['CIII_1909'] == 0 or index['OIII_1665'] == 0:
                  CHI_C3O3 = tol_max
               else:
                  CHI_C3O3 =(np.log10((index['CIII_1909'])/index['OIII_1665']) - C3O3_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_1665']+1e-5))
               if C3O5007_obs == -10:
                  CHI_C3O5007 = 0          
               elif index['CIII_1909'] == 0 or index['OIII_5007'] == 0:
                  CHI_C3O5007 = tol_max
               else:
                  CHI_C3O5007 =(np.log10((index['CIII_1909'])/index['OIII_5007']) - C3O5007_obs)**2/np.log10((index['CIII_1909'])/(index['OIII_5007']+1e-5))
               if C3Ne3_obs == -10:
                  CHI_C3Ne3 = 0          
               elif index['CIII_1909'] == 0 or index['NeIII_3869'] == 0:
                  CHI_C3Ne3 = tol_max
               else:
                  CHI_C3Ne3 =(np.log10((index['CIII_1909'])/index['NeIII_3869']) - C3Ne3_obs)**2/np.log10((index['CIII_1909'])/(index['NeIII_3869']+1e-5))
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
               if R23_obs == -10:
                  CHI_R23 = 0
               elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                  CHI_R23 = tol_max
               else:
                  CHI_R23 = (np.log10(index['OII_3727']+index['OIII_5007']) - R23_obs)**2/np.log10(index['OII_3727']+index['OIII_5007'])
               if O2O3_obs == 0:
                  CHI_O2O3 = 0
               elif index['OII_3727'] == 0 or index['OIII_5007'] == 0:
                  CHI_O2O3 = tol_max
               else:
                  CHI_O2O3 = ((index['OII_3727']/index['OIII_5007']) - O2O3_obs)**2/(index['OII_3727']/index['OIII_5007'])
               if R2Ne3_obs == -10:
                  CHI_R2Ne3 = 0
               elif index['OII_3727'] == 0 or index['NeIII_3869'] == 0:
                  CHI_R2Ne3 = tol_max
               else:
                  CHI_R2Ne3 = (np.log10(index['OII_3727']+index['NeIII_3869']) - R2Ne3_obs)**2/np.log10(index['OII_3727']+index['NeIII_3869'])
               if O2Ne3_obs == 0:
                  CHI_O2Ne3 = 0
               elif index['OII_3727'] == 0 or index['NeIII_3869'] == 0:
                  CHI_O2Ne3 = tol_max
               else:
                  CHI_O2Ne3 = ((index['OII_3727']/index['NeIII_3869']) - O2Ne3_obs)**2/(index['OII_3727']/index['NeIII_3869'])

               if R23_obs > -10:
                  CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34Hb**2  + CHI_C3C4**2 + CHI_R23**2 + CHI_O2O3**2)**0.5
               elif R23_obs == -10 and R2Ne3_obs > -10:
                  CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34Hb**2  + CHI_C3C4**2 +  CHI_R2Ne3**2 + CHI_O2Ne3**2)**0.5
               else:	
                  CHI_OH = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C34**2 + CHI_C34He2**2 + CHI_N5**2 + CHI_N5He2**2 + CHI_C3C4**2 + CHI_O2Ne3**2)**0.5


               if CHI_OH == 0:
                  OH_p = OH_p
                  logU_p = logU_p
                  den_OH = den_OH
               else:
                  OH_p = index['12logOH'] /(CHI_OH)**2 + OH_p
                  logU_p = index['logU'] /(CHI_OH)**2 + logU_p
                  den_OH = 1 /(CHI_OH)**2 + den_OH

               CHI_CO = (CHI_ROIII**2 + CHI_R2OIII**2 + CHI_R3OIII**2 + CHI_C3O3**2 + CHI_C3O5007**2 + CHI_C3Ne3**2)**0.5

               if CHI_CO == 0:
                  CO_p = CO_p
                  den_CO = den_CO
               else:
                  CO_p = index['logCO'] /(CHI_CO)**2**2 + CO_p
                  den_CO = 1 /(CHI_CO)**2**2 + den_CO

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


         if OHf > 0: OH_mc.append(OHf)
         if COf > -10: CO_mc.append(COf)
         if logUf < 0: logU_mc.append(logUf)
         if OHf > 0: OHe_mc.append(eOH)
         if COf > -10: COe_mc.append(eCO)
         if logUf < 0: logUe_mc.append(elogU)


         if len(OH_mc) > 0:
            OHff = np.mean(OH_mc)
            eOHff = (np.std(OH_mc)**2+np.mean(OHe_mc)**2)**0.5
         else:
            OHff = 0
            eOHff = 0
         if len(logU_mc) > 0:
            logUff = np.mean(logU_mc)
            elogUff = (np.std(logU_mc)**2+np.mean(logUe_mc)**2)**0.5
         else:
            elogUff = 0
            logUff = 0
         if len(CO_mc) > 0:
            COff = np.mean(CO_mc)
            eCOff = (np.std(CO_mc)**2+np.mean(COe_mc)**2)**0.5
         else:
            COff = -10
            eCOff = 0
            
   except:
      grids = 9999
      grid_type = 0
      OHff = 9999
      eOHff = 9999
      COff = 9999
      eCOff = 9999
      logUff = 9999
      elogUff = 9999

      
        
   OHffs.append(OHff)
   eOHffs.append(eOHff)
   COffs.append(COff)
   eCOffs.append(eCOff)
   logUffs.append(logUff)
   elogUffs.append(elogUff)
         
   ##################################
   # Displaying results in terminal #
   ##################################

   if input0.size == 1 and tab==0: continue

   print (round(100*(count)/float(input1.size),1),'%',Names[tab],grid_type,'', round(OHff,2), round(eOHff,2),'',round(COff,2), round(eCOff,2), '',round(logUff,2), round(elogUff,2))

####################################################
###### OUTPUT FORMAT AND INFORMATION: RESULTS ######
####################################################

#Grid used and results from the free parameters
output['grid'] = grids
output['OH'] = OHffs
output['eOH'] = eOHffs
output['CO'] = COffs
output['eCO'] = eCOffs
output['logU'] = logUffs
output['elogU'] = elogUffs

if input0.size == 1:  output = np.delete(output,obj=1,axis=0)

#Header comments for the file
lineas_header = [' HII-CHI-mistry_UV v.5.3 output file', 'Input file:'+input00,'Iterations for MonteCarlo: '+str(n),'Used models: '+sed_type,'Library file used : '+file_lib_2, 'Template used to constraint grid of models: '+const_file,'']

#Labels for columns (emission lines)
line_label = '{:3}  '.format(output.dtype.names[0])
for ind2 in range(1, len(output.dtype.names)):
   line_label += '{:3}  '.format(output.dtype.names[ind2])

#Labels for columns
lineas_header.append(line_label)
header = '\n'.join(lineas_header)

#Results
np.savetxt('.'.join(input00.split('.')[:-1])+'_hcm-uv-output.dat',output,fmt=' '.join(['%s']*1+['%.3f']*(len(output.dtype.names)-8)+['%i']+['%.2f']*6), header=header)

lines_stor = []
with open('.'.join(input00.split('.')[:-1])+'_hcm-uv-output.dat', 'r+') as output_file:
   for line in output_file:
      lines_stor.append(line)

#Reformating output for better reading of the table
file_overwrite = open('.'.join(input00.split('.')[:-1])+'_hcm-uv-output.dat', 'r+')
file_overwrite.seek(0)
for line_n in lines_stor:  
   if line_n[0] == '#' and line_n[2:4] == 'ID':
      file_overwrite.write(line_n[2:])
   else:
      file_overwrite.write(line_n)
file_overwrite.truncate()   
file_overwrite.close()
print ('-------------------------------------------------')
print ('Results are stored in ' + '.'.join(input00.split('.')[:-1]) + '_hcm-uv-output.dat')
print ('-------------------------------------------------')

#############################################
###### INFORMATION AND CONTACT DETAILS ######
#############################################

# Enrique Perez-Montero, epm@iaa.es
# Borja Perez-Diaz, bperez@iaa.es


#################
###### END ######
#################
