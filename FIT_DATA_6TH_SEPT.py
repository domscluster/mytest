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

path = '/home/milan/FOUR.DIMENSION.POTENTIAL/Ellpot_PHIPSI_Scheragalist_assym_fin_new/*.dat'
files=glob.glob(path)

#define the functions here
def gb(var1,k1,k2,miu,niu,sigma_o,epsilon_o):
    lenset =76
    x=np.zeros(lenset)
    for ii in range(lenset):
        x[ii] = var1[ii][0]
    theta1 = var1[0][1]
    theta2 = var1[0][2]
    dih    = var1[0][3]
    #print x 
    theta1= theta1*0.017444444
    theta2= theta2*0.017444444
    dih   = dih *0.017444444
    chi = ((k1**2)-1)*(((k1**2)+1))**(-1)
    #print chi
    term = k2**(0.5)
    chi_1=(term-1)*(term+1)**(-1)
    chi_by2     =chi/2
    chi_1_by2   = chi_1/2
    e_ref_1 = [1 , 0 , 0]
    e_ref_2 =[1 , 0 , 0]
    e1_theta=[e_ref_1[0]*np.cos(theta1)-e_ref_1[2]*np.sin(theta1) , e_ref_1[1] , e_ref_1[2]*np.cos(theta1) + e_ref_1[0]*np.sin(theta1)]
    e2_theta=[e_ref_2[0]*np.cos(theta2)-e_ref_2[2]*np.sin(theta2) , e_ref_2[1] , e_ref_2[2]*np.cos(theta2) + e_ref_2[0]*np.sin(theta2)]
    e2_phi=[e2_theta[0] , e2_theta[1]*np.cos(dih) + e2_theta[2]*np.sin(dih) , e2_theta[2]*np.cos(dih) - e2_theta[1]*np.sin(dih)]
    #print e1_theta, e2_theta,e2_phi
    rij=[1 , 0 , 0]
    epsilon_1 =1-((chi**2)*(np.dot(e1_theta , e2_phi)))**2
    epsilon_1=epsilon_1**(-0.5)
    term1     = (np.dot(e1_theta , rij) + np.dot(e2_theta , rij) )**2
    term2     = 1+chi*np.dot(e1_theta , e2_phi)
    term3     = (np.dot(e1_theta , rij) - np.dot(e2_theta , rij) )**2
    term4     = 1-chi*np.dot(e1_theta , e2_phi)
    term5     = 1+chi_1*np.dot(e1_theta , e2_phi)
    term6     = 1-chi_1*np.dot(e1_theta , e2_phi)
    epsilon_2 = 1-((chi_1_by2)*(-((term1)/term5)+((term3)/term6)))
    sigma    = (1-((chi_by2)*(((term1)/term2)+((term3)/term4))))
    sigma    = sigma**(-0.5)
    epsil     = epsilon_o * (epsilon_1**niu) *(epsilon_2**miu)
    #print term1 , term2 , term3 , term4
    r2i       = (1/(x-sigma +sigma_o))**2
    r6i       = r2i**3
    return 4.0*epsil*r6i*(r6i - 1)
    

def lj(x_dum, a, b ,c , d , e ):
    test1 = b*(x_dum+c)**(-12)
    test2 = d*(x_dum+e)**(-6)
    return a*(test1-test2)

def xifi(x, a, b ,c ):
    return (a * (((b/(x))**12) -((c/(x))**6)))
def zomb(p, a, b ,c ):
    return 4*a * (((b/p)**12) -((c/p)**6))



def exp_f(x , a0 , a1):
    return a0 * np.exp(-a1*x)


def exp_lj_f(x , a0 , a1 , a2 , a3):
    return a0 * np.exp(-a1*x) - 4*a2 * np.power((a3/(x+0.002)) , 6)


def exp_db_exp(x , a0 , a1 , a2 , a3):
    return a0 * np.exp(-a1*x) + a2 * np.exp(-a3*x)

#=========================================================================================#
#    if you need you can incorporate function here                                        #
#=========================================================================================#



print "==================================================================================="
print "      Hey man you successfully incorporated the model potential congrats           "
print "                   GET A SCOTCH                                                    "
print "==================================================================================="






#Now reading a file with maximum 5 coloumn dimensions a trick 
#pure python way



#filein = open('sample.dat' , 'r')    #opening the file through this standard file open statement
#data = filein.read()                 # Reads the file as a complete string  #optional
#for line in filein:
#     print repr(line)
#     line = line.strip()
#     print line
     
#filein.close()
#print data


#optional for others in case required  :)

#with open ('sample.dat' , 'r') as filein:
#     reader = csv.reader(filein)
#     for line in reader:
#         print line
        
max_col_ss = 2 
max_col_se = 3
max_col_ee = 5 
lenset     = 76
nset       =343
 
#list_of_lists_ss = []
#list_of_lists_se = []
list_of_lists = []


for file in files:
 with open(file) as filein:
    for line in filein:
        line = line.strip()
        if not  (line.startswith("S")):
           if not  (line.startswith("@")):
              columns = line.split()
              new_list = [float(i) for i in columns]
#              if len(line.split()) == max_col_ee:
              list_of_lists.append((new_list))
#              elif len(line.split()) == max_col_se:
#                 list_of_lists_se.append((new_list))
#              else:
#                 list_of_lists_ss.append((new_list))
 filein.close()
#print list_of_lists_se       


 print ' <==> '*10
 print  'I ate it'
 print ' <==> '*10

 print "number of lines ", len(list_of_lists)




 print "==================================================================================="
 print "      Hey man I ate it now I will chew the meal and break it into pieces           "
 print "==================================================================================="

 theta1=[]
 theta2=[]
 phi   =[]
 x=[]
 y=[]
 var=np.zeros((lenset , 4)) 
 for ii in range(nset):
    for jj in range (lenset):
        mm = (ii)*lenset +jj
        theta1_num=list_of_lists[mm][0]
#        print theta1_num        
        theta2_num=list_of_lists[mm][1]
        phi_num   =list_of_lists[mm][2]
        rij = list_of_lists[mm][3]
        Eij = list_of_lists[mm][4]
        theta1.append(theta1_num)
        theta2.append(theta2_num)
        phi.append(phi_num)
        x.append(rij)
        y.append(Eij)
#for ii in range(5) :       
#    print theta1[ii]    
#print theta1
 count_n = 0 
 for i in range(0, len(y), lenset):                  #change lenset to len(y)
        count_n = count_n +1 
        y_plot =np.array(y[i:i + lenset])
        x_plot =np.array(x[i:i + lenset])            #made the data into chunks
        
        theta1_plot = theta1[i]
        theta2_plot = theta2[i]
        phi_plot    = phi[i]
        for jj in range( lenset):
            var[jj][0] = x_plot[jj]            
            var[jj][1] = theta1_plot
            var[jj][2] = theta2_plot    
            var[jj][3] = phi_plot   


        #sigma_plot =   1.0
        #epsilon_plot = 1.0
        #k1_plot      = 3
        #k2_plot      = 5
        #miu_plot     = 1
        #niu_plot     =2
#       plt.ylim((-20 , 20))
        index_min = np.argmin(y_plot)
        minval    = min(y_plot)
        min_x     =x_plot[index_min]
#        print (min_x , minval)
#        plt.plot(x_plot,y_plot, '.')
#        xfine = np.linspace(0.13, 10, num=100)
#        plt.ylim((-10 , 100))
#        yfine = exp_f(x_plot, 20 , 0.5 )
#        print yfine
#        popt, pcov = curve_fit(exp_lj_f, x_plot, y_plot  , p0=(3000 , 0.5 , minval , min_x))
        popt, pcov = curve_fit(exp_f, x_plot, y_plot  )                                       #fitting the curve 
#        gmod = Model(zomb)
        plt.plot(x_plot, y_plot, '.')
        plt.plot(x_plot, exp_f(x_plot, popt[0], popt[1] ), 'r-')
        plt.show()
        print count_n  , popt
        
 print count_n





#delete all the lists used 
 del x[:]
 del y[:]
 del list_of_lists[:]
 del theta1[:]
 del theta2[:]
 del phi[:]

