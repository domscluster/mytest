# "=================================================================================="
#"                                                                                   "
#"                   PROGRAM TO FIT ANY KIND OF DATA WITH LM ALGORITHM               "
#                    WRITTEN :Milan Hazra(5th SEPT 2016)                             "
#                    some functions are built in there aditional model functions     "
#                        you can incorporate in the specific section                 "
#"==================================================================================="
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit import  Model
import math
import csv
import os
import glob

max_col_ss = 2
max_col_se = 3
max_col_ee = 5
lenset     = 76
nset_ee       =343
nset_se       = 532
n_set_ss      = 1

mypath='/home/milan/FOUR.DIMENSION.POTENTIAL/Ellpot_PHIPSI_Scheragalist_assym_fin_new'
os.listdir(mypath)
dirs = os.listdir(mypath)
#print dirs[1]

def zomb(p, a, b ,c ):
    return 4*a * (((b/p)**12) -((c/p)**6))








path = '/home/milan/FOUR.DIMENSION.POTENTIAL/Ellpot_PHIPSI_Scheragalist_assym_fin_new/*.dat'
files=glob.glob(path)
c_file = 0
for file in files:
    mm = 0    
    list_of_lists = []
    with open(file) as filein:
         c_file = c_file+1
         for line in filein:
             line = line.strip()
             if not  (line.startswith("S")):
                if not  (line.startswith("@")):
                   columns = line.split()
                   new_list = [float(i) for i in columns]
                   len_col  = len(new_list)
                   list_of_lists.append((new_list))
    filein.close()
    print ' <==> '*10
    print  'I ate it'
    print ' <==> '*10
    theta1=[]
    theta2=[]
    phi   =[]
    x=[]
    y=[]
    var=np.zeros((lenset , 4))
    print len(list_of_lists)
    
    nset = len(list_of_lists)/lenset
    print (nset)

    for ii in range(1 , nset):
        for jj in range ( lenset):
            
            mm = (ii)*lenset +jj
            if(len(new_list) == max_col_ee) :
              theta1_num=list_of_lists[mm][0]
              theta2_num=list_of_lists[mm][1]
              phi_num   =list_of_lists[mm][2]
              rij = list_of_lists[mm][3]
              Eij = list_of_lists[mm][4]
              theta1.append(theta1_num)
              theta2.append(theta2_num)
              phi.append(phi_num)
              x.append(rij)
              y.append(Eij)

            elif(len(new_list) == max_col_se):
              theta1_num=list_of_lists[mm][0]
              rij = list_of_lists[mm][1]
              Eij = list_of_lists[mm][2]
              theta1.append(theta1_num)
              x.append(rij)
              y.append(Eij)
            elif(len(new_list) == max_col_ss):
              rij = list_of_lists[mm][0]
              Eij = list_of_lists[mm][1]
              x.append(rij)
              y.append(Eij)

    if(len(new_list) == max_col_ee) :
        for i in range(0, len(y), lenset):
          y_plot =np.array(y[i+25:i + lenset])
          x_plot =np.array(x[i+25:i + lenset])
          theta1_plot = theta1[i]
          theta2_plot = theta2[i]
          phi_plot    = phi[i]
          #fit curve here
          index_min = np.argmin(y_plot)
          minval    = min(y_plot)
          min_x     =x_plot[index_min]
          popt, pcov = curve_fit(zomb, x_plot, y_plot  , p0=( minval , min_x , min_x))
          plt.plot(x_plot,y_plot, '.')
          plt.plot(x_plot,zomb(x_plot , popt[0] , popt[1] , popt[2]), 'r-')
          plt.show()
    elif(len(new_list) == max_col_se) :
        for i in range(0, len(y), lenset):
          y_plot =np.array(y[i+25:i + lenset])
          x_plot =np.array(x[i+25:i + lenset])
          theta1_plot = theta1[i]
          #fit curve here
          index_min = np.argmin(y_plot)
          minval    = min(y_plot)
          min_x     =x_plot[index_min]
          popt, pcov = curve_fit(zomb, x_plot, y_plot  , p0=( minval , min_x , min_x))
          plt.plot(x_plot,y_plot, '.')
          plt.plot(x_plot,zomb(x_plot , popt[0] , popt[1] , popt[2]), 'r-')
          plt.show()

    elif(len(new_list) == max_col_ss) :
        for i in range(0, len(y), lenset):
          y_plot =np.array(y[i+25:i + lenset])
          x_plot =np.array(x[i+25:i + lenset])
          #fit curve here
          index_min = np.argmin(y_plot)
          minval    = min(y_plot)
          min_x     =x_plot[index_min]
          popt, pcov = curve_fit(zomb, x_plot, y_plot  , p0=( minval , min_x , min_x))
          plt.plot(x_plot,y_plot, '.')
          plt.plot(x_plot,zomb(x_plot , popt[0] , popt[1] , popt[2]), 'r-')
          plt.show()
       
    
              
    del list_of_lists[:]
    del theta1[:]
    del theta2[:]
    del phi[:]
    del x[:]
    del y[:]
#iprint c_file
