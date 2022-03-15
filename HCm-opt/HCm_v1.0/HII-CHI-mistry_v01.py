#!/usr/bin/python
# Filename: clasificacion.py

print ' ---------------------------------------------------------------------'
print ' This is HII-CHI-mistry version 0.1'
print ' See Perez-Montero, E. (2014) for details'
print  ' Insert the name of your input text file with following columns:'
print ' 3727 [OII], 4363 [OIII], 5007 [OIII], 6584 [NII], 6725 [SII]'
print '---------------------------------------------------------------------'


import string
import asciidata
from math import log10

# Reading of models grids. These can be changed

grid1 = asciidata.open('C13_cha_1Myr.dat')
grid2 = asciidata.open('C13_cha_1Myr_logU_adapted_emp.dat')
grid3 = asciidata.open('C13_cha_1Myr_logU-NO_adapted_emp.dat')

# Input file reading

input = asciidata.open(raw_input('Insert input file name:'))

print 'Reading grids ....'
print ''
print ''
print '----------------------------------------------------------------'
print '(%)   Grid  12+log(O/H)  log(N/O)    log(U)'
print '-----------------------------------------------------------------'





# Initial conditions. These can be edited. See readme file for an explanation

nmod = 8  #number of models with the lowest tolerance to be considered in the mean
ntol = 0.5  #discrete step of tolerance to be added to the minimum tolerance if the number of models is not reached..
tol_max = 20 #maximum tolerance taken into account in the process
mod_min = 3 #minimum number of models of calculate a mean if the maxium tolerance is reached.


# Beginning of loop of calculation

for tab in range(input.nrows):

   OH_p = 0
   logU_p = 0
   NO_p = 0
   den_OH = 0
   den_NO = 0
   OH_e = 0
   NO_e = 0
   logU_e = 0
   den_OH_e = 0
   den_NO_e = 0

# definicion de variables iniciales de los bucles

   tol_OH = 0
   tol_NO = 0
   tol_OH_e = 0
   tol_NO_e = 0
   modelos_OH = 0
   modelos_NO = 0
   modelos_OH_e = 0
   modelos_NO_e = 0

   

   OII_3727_obs = input[0][tab]
   OIII_4363_obs = input[1][tab]
   if OIII_4363_obs == 0:
      ROIII_obs = 0
   else:
      ROIII_obs = input[2][tab]/input[1][tab]
   OIII_5007_obs = input[2][tab]
   if input[3][tab] == 0:
      NII_6584_obs = 0
   else:
      NII_6584_obs = input[3][tab]
   SII_6725_obs = input[4][tab]
   if NII_6584_obs == 0 or OII_3727_obs == 0:
      N2O2_obs = -10
   else:
      N2O2_obs = log10(input[3][tab] / input[0][tab])
   if NII_6584_obs == 0 or SII_6725_obs == 0:
      N2S2_obs = -10
   else:   
      N2S2_obs = log10(input[3][tab] / input[4][tab])
   if OII_3727_obs == 0 or OIII_5007_obs == 0:
      O2O3_obs = 0
      R23_obs = 0
   else:
      O2O3_obs = (input[0][tab] / input[2][tab])
      R23_obs = log10(input[0][tab] + input[2][tab])
   if OIII_5007_obs == 0 or NII_6584_obs == 0:
      O3N2_obs = 0
   else:
      O3N2_obs = (input[2][tab] / input[3][tab])
   if OIII_5007_obs == 0 or SII_6725_obs == 0:
      O3S2_obs = 0
   else:
      O3S2_obs = (input[2][tab] / input[4][tab])


   if ROIII_obs > 0:
      grid = grid1
      input['grid'][tab] = 1
   elif N2O2_obs > -10 or N2S2_obs > -10:
      grid = grid2
      input['grid'][tab] = 2
   else:
      grid = grid3
      input['grid'][tab] = 3

# definicion de variables iniciales de los bucles


   while modelos_NO < nmod and tol_NO <= tol_max:
      tol_NO = tol_NO + ntol
      modelos_NO = 0
      NO_p = 0
      den_NO = 0



      for index in range(grid.nrows):
         if ROIII_obs == 0:
            CHI_ROIII = 0
         else:
           CHI_ROIII = (grid[5][index]/grid[4][index] - ROIII_obs)**2/(grid[5][index]/grid[4][index])
         if NII_6584_obs == 0 or OII_3727_obs == 0:
            CHI_N2O2 = 0
         else:
            CHI_N2O2 =(log10(grid[6][index]/grid[3][index]+0.001) - N2O2_obs)**2/(log10(grid[6][index]/grid[3][index])+0.001)
         if NII_6584_obs == 0 or SII_6725_obs == 0:
            CHI_N2S2 = 0
         elif grid[6][index] < 0.001 and grid[7][index] < 0.001:
            CHI_N2S2 = tol_max
         else:
            CHI_N2S2 =(log10(grid[6][index]/grid[7][index]) +0.001 - N2S2_obs)**2/log10(grid[6][index]/grid[7][index] + 0.001) 
             

         CHI_NO = (CHI_ROIII**2 + CHI_N2O2**2 + CHI_N2S2**2)**0.5
         
         if CHI_NO == 0:
            NO_p = NO_p
            den_NO = den_NO
            modelos_NO = modelos_NO
         elif CHI_NO < tol_NO:
            modelos_NO = modelos_NO + 1
            NO_p = grid[1][index] / CHI_NO + NO_p
            den_NO = 1/CHI_NO + den_NO
#            print grid[1][index], CHI_NO, NO_p/den_NO
         
      if tol_NO >= tol_max and modelos_NO < mod_min:
         NO_p = -10
         den_NO = 1
         break

   
    

   
        
   while modelos_OH < nmod and tol_OH <= tol_max:
      tol_OH = tol_OH + ntol
      OH_p = 0
      logU_p = 0
      den_OH = 0
      modelos_OH = 0
      

      for index in range(grid.nrows):
         if ROIII_obs == 0:
            CHI_ROIII = 0
         else:
            CHI_ROIII = (grid[5][index]/grid[4][index] - ROIII_obs)**2/(grid[5][index]/grid[4][index])
         if NII_6584_obs == 0:
            CHI_NII = 0
         else:
            CHI_NII = ((grid[6][index]) - (NII_6584_obs))**2/grid[6][index]
         if SII_6725_obs == 0:
            CHI_SII = 0
         else:
             CHI_SII = ((grid[7][index]) - (SII_6725_obs))**2/(grid[7][index])
         if OII_3727_obs == 0:
            CHI_OII = 0
         else:
             CHI_OII=((grid[3][index])-(OII_3727_obs))**2/(grid[3][index])
         if OIII_5007_obs == 0:
            CHI_OIII = 0
         else:
             CHI_OIII=((grid[5][index])-(OIII_5007_obs))**2/(grid[5][index])
             
         if OII_3727_obs == 0 or OIII_5007_obs == 0:
            CHI_O2O3 = 0
            CHI_R23 = 0
         else:
             CHI_O2O3 =((grid[3][index]/grid[5][index]) - O2O3_obs)**2/log10(grid[3][index]/grid[5][index]) 
             CHI_R23 =(log10(grid[3][index]+grid[5][index]) - R23_obs)**2/log10(grid[3][index]+grid[5][index]) 
         if OIII_5007_obs == 0 or NII_6584_obs == 0:
            CHI_O3N2 = 0
         elif grid[5][index] < 0.001 and grid[6][index] < 0.001:
            CHI_O3N2 = tol_max
         else:
            CHI_O3N2=((grid[5][index]/grid[6][index]) - O3N2_obs)**2/(grid[5][index]/grid[6][index])
         if OIII_5007_obs == 0 or SII_6725_obs == 0:
            CHI_O3S2 = 0
         elif grid[5][index] < 0.001 and grid[7][index] < 0.001:
            CHI_O3S2 = tol_max
         else:
            CHI_O3S2=((grid[5][index]/grid[7][index]) - O3S2_obs)**2/(grid[5][index]/grid[7][index])

                   
         if OIII_4363_obs > 0:
            CHI_OH = (CHI_ROIII**2 + CHI_OII**2 + CHI_OIII**2 + CHI_NII**2 + CHI_SII**2)**0.5
         else:
            if OII_3727_obs == 0:
               if NII_6584_obs == 0:
                  OH_p = 0
                  logU_p = 0
                  den_OH = 1
                  break
               else:
                  CHI_OH = (CHI_NII**2 + CHI_OIII**2 + CHI_O3N2**2)**0.5
            else:
               CHI_OH = (CHI_NII**2 + CHI_R23**2 + CHI_O2O3**2 )**0.5
         

         if CHI_OH < tol_OH and CHI_OH > 0 and abs(NO_p/den_NO-grid[1][index])<0.13:
            modelos_OH = modelos_OH + 1
            OH_p = grid[0][index] / CHI_OH + OH_p
            logU_p = grid[2][index] / CHI_OH + logU_p
            den_OH = 1/CHI_OH + den_OH
         elif CHI_OH < tol_OH and CHI_OH > 0 and NO_p == -10:
            modelos_OH = modelos_OH + 1
            OH_p = grid[0][index] / CHI_OH + OH_p
            logU_p = grid[2][index]/CHI_OH + logU_p
            den_OH = 1/CHI_OH + den_OH
         

      if tol_OH >= tol_max and modelos_OH < mod_min:
         OH_p = 0
         logU_p = 0
         den_OH = 1
         break


  


   while modelos_OH_e < nmod and tol_OH_e <= tol_max:
      tol_OH_e = tol_OH_e + ntol
      OH_e = 0.05**2/0.1
      logU_e = 0
      den_OH_e == 0



      for index in range(grid.nrows):
         if ROIII_obs == 0:
            CHI_ROIII = 0
         else:
            CHI_ROIII = (grid[5][index]/grid[4][index] - ROIII_obs)**2/(grid[5][index]/grid[4][index])
         if NII_6584_obs == 0:
            CHI_NII = 0
         else:
            CHI_NII = ((grid[6][index]) - (NII_6584_obs))**2/grid[6][index]
         if SII_6725_obs == 0:
            CHI_SII = 0
         else:
             CHI_SII = ((grid[7][index]) - (SII_6725_obs))**2/(grid[7][index])
         if OII_3727_obs == 0:
            CHI_OII = 0
         else:
             CHI_OII=((grid[3][index])-(OII_3727_obs))**2/(grid[3][index])
         if OIII_5007_obs == 0:
            CHI_OIII = 0
            CHI_R23 = 0
         else:
             CHI_OIII=((grid[5][index])-(OIII_5007_obs))**2/(grid[5][index])
             
         if OII_3727_obs == 0:
            CHI_O2O3 = 0
         else:
             CHI_O2O3 =((grid[3][index]/grid[5][index]) - O2O3_obs)**2/(grid[3][index]/grid[6][index]) 
             CHI_R23 =(log10(grid[3][index]+grid[5][index]) - R23_obs)**2/log10(grid[3][index]+grid[6][index]) 
         if OIII_5007_obs == 0 or NII_6584_obs == 0:
            CHI_O3N2 = 0
         elif grid[6][index] < 0.001 and grid[7][index] < 0.001:
            CHI_O3N2 = tol_max
         else:
            CHI_O3N2=((grid[5][index]/grid[6][index]) - O3N2_obs)**2/(grid[5][index]/grid[6][index])

         if OIII_5007_obs == 0 or SII_6725_obs == 0:
            CHI_O3S2 = 0
         elif grid[6][index] < 0.001 and grid[7][index] < 0.001:
            CHI_O3S2 = tol_max
         else:
            CHI_O3S2=((grid[5][index]/grid[7][index]) - O3S2_obs)**2/(grid[5][index]/grid[7][index])


         if OIII_4363_obs > 0:
            CHI_OHe = (CHI_ROIII**2 + CHI_OIII**2 + CHI_OII**2 + CHI_NII**2 + CHI_SII**2)**0.5
         else:
            if OII_3727_obs == 0:
               if NII_6584_obs == 0:
                  OH_e = 0
                  logU_e = 0
                  den_OH_e = 1
                  break
               else:
                  CHI_OHe = (CHI_NII**2 + CHI_OIII**2 + CHI_O3N2**2)**0.5
            else:
               CHI_OHe = (CHI_NII**2 + CHI_R23**2 + CHI_O2O3**2 )**0.5

           

         if CHI_OHe < tol_OH_e and CHI_OHe > 0 and abs(NO_p/den_NO-grid[1][index])<0.13:
            OH_e = (OH_p/den_OH - grid[0][index])**2/CHI_OHe + OH_e
            logU_e = (logU_p/den_OH - grid[2][index])**2/CHI_OHe + logU_e
            den_OH_e = 1/CHI_OHe + den_OH_e
            modelos_OH_e = modelos_OH_e + 1 
         elif CHI_OHe < tol_OH_e and CHI_OHe > 0 and NO_p == -10:
            OH_e = (OH_p/den_OH - grid[0][index])**2/CHI_OHe + OH_e
            logU_e = (logU_p/den_OH - grid[2][index])**2/CHI_OHe + logU_e
            den_OH_e = 1/CHI_OHe + den_OH_e
            modelos_OH_e = modelos_OH_e + 1 



      if tol_OH_e >= tol_max and modelos_OH_e < mod_min:
         OH_e = 0
         logU_e = 0
         den_OH_e = 1
         break





              
 

   while modelos_NO_e < nmod and tol_NO_e <= tol_max:
      tol_NO_e = tol_NO_e + ntol
      NO_e = 0
      den_NO_e = 0
      modelos_NO_e = 0


      for index in range(grid.nrows):
         if ROIII_obs == 0:
            CHI_ROIII = 0
         else:
           CHI_ROIII = (grid[5][index]/grid[4][index] - ROIII_obs)**2/(grid[5][index]/grid[4][index])
         if NII_6584_obs == 0 or OII_3727_obs == 0:
            CHI_N2O2 = 0
         else:
             CHI_N2O2 =(log10(grid[6][index]/grid[3][index]+0.001) - (N2O2_obs))**2/(log10(grid[6][index]/grid[3][index])+0.001 )
         if NII_6584_obs == 0 or SII_6725_obs == 0:
            CHI_N2S2 = 0
         elif grid[6][index] < 0.001 and grid[7][index] < 0.001:
            CHI_N2S2 = tol_max
         else:
             CHI_N2S2 =(log10(grid[6][index]/grid[7][index]+0.001) - (N2S2_obs))**2/(log10(grid[6][index]/grid[7][index])+0.001 )
             

         CHI_NOe = (CHI_N2O2**2 + CHI_N2S2**2)**0.5 + 1e-3



         
#         if CHI_NOe == 0:
#            NO_e = 1
#            den_NO_e = 1
#            modelos_NO_e = modelos_NO_e + 1
         if CHI_NOe < tol_NO_e:
            modelos_NO_e = modelos_NO_e + 1
            NO_e = (NO_p/den_NO - grid[1][index])**2/ CHI_NOe + NO_e
            den_NO_e = 1/CHI_NOe + den_NO_e



      if tol_NO_e >= tol_max and modelos_NO_e < mod_min or (N2O2_obs == -10 and N2S2_obs == -10):
         NO_e = 0
         den_NO_e = 1
         break



   input['OH'][tab] = OH_p/den_OH
   input['OH_error'][tab] = (OH_e/den_OH_e)**0.5
   input['NO'][tab] = NO_p/den_NO  
   input['NO_error'][tab] = (NO_e/den_NO_e)**0.5
   input['logU'][tab] = logU_p/den_OH   
   input['logU_error'][tab] = (logU_e/den_OH_e)**0.5

   print round(100*tab/float(input.nrows),1),'%', input['grid'][tab],'   ', round(OH_p/den_OH,3), round((OH_e/den_OH_e)**0.5,3), round(NO_p/den_NO,2), round((NO_e/den_NO_e)**0.5,3), round(logU_p/den_OH,2), round((logU_e/den_OH_e)**0.5,2)
   

#print modelos, modelos_2, tol, tol_2
#print '------------------------------'
#print '12+log(O/H) =',round(OH,2),'+/-',round(OH_error,2)
#print 'log(N/O) =',round(NO,2),'+/-',round(NO_error,2)
#print 'log U =',round(logU,2),'+/-',round(logU_error,2)
#print'-------------------------------'

grid.flush()
input.flush()
