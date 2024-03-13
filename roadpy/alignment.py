#1. IMPORT PACKAGES
import math
import pandas as pd
import haversine as hs
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

#2. DEFINING VARIABLES, FUNCTIONS
# variables:
# a = Design Speed, Vd(Kph) 
# b = emax  
# c = fmax 
# d = R design, Rd (m)     
# e = T, Transition curve travel time, tmax=3,00”
# Dd =
# ed = e_design
# f = delta_i, The rate of change in slope: 0.035 (Vr<80) atau 0,025 (Vr>80), (ref SE Dirjen Bina Marga Nomor 20/SE/Db/2021, page:100)
# g = Angle PI, derajat
# Ls = Length of spiral design
# Angle Spiral
# Angle Circle
# P =
# K =
# Lc = Length of circle, Ls = Length of Spiral
# u, v = the angle of vectors u, v


def hitungFmax(a):
    return (0.192-(0.00065*(a)))

def hitungRmin(a,b,c):
    return (a**2)/(127*(b+c))

def hitungDmax(a,b,c):
    return ((181913.53*(b+c))/(a**2))

def hitungDd(d):
    return 1432.4/(d)

def e_Hitung_def(b, Dd, Dmax):
    return (((-b*Dd**2)/(Dmax**2))+((2*b*Dd)/(Dmax)))

def hitung_Lst(a,e):
    return (a/3.6)*e

def hitung_Lsms(a,d):   # modifikasi shortt: return (((0.022*a**3)/(d*0.4))-(2.727*(a*ed)/0.4))
    return (((0.0214*a**3)/(d*1.2)))  # C=1,2, ref: SE 20/2021 hal 125

def hitung_Lsi(a,f):       # Diantara 2 Ls, diambil yang terpanjang
    return (0.1-0.02)*a/(3.6*f)

def sudutSpiral_def(d, Ls):
    return (90*Ls/((math.pi)*d))

def sudutCircle(g, sudutSpiral_def):
    return (g-(2*(sudutSpiral_def)))

def hitung_Lc(sudutCircle, d):
    return (sudutCircle*math.pi*d)/180

def hitung_Xs(Ls, d):
    return Ls*(1-((Ls**2)/(40*d**2)))

def hitung_Ys(Ls, d):
    return (Ls**2)/(6*d)

def hitung_P(Ls,d):
    return Ls**2/(24*d)

def hitung_K(Ls, d, sudutSpiral):
    return Ls-(Ls**2/(40*(d**2)))-d*math.sin(sudutSpiral*math.pi/180)

def hitung_Tt(d, P, g, K):
    return (d+P)*math.tan(0.5*g*math.pi/180)+K

def hitung_Et(d, hitung_P, g):
    return ((d+hitung_P)/(math.cos(0.5*g*math.pi/180)))-d

def hitung_Lt(hitung_Lc, Ls):
    return hitung_Lc+2*Ls

def angle(u, v):
    """return the angle between two vectors in any dimension space, indegrees"""
    return np.degrees(
        math.acos(np.dot(u,v)/ (np.linalg.norm(u)*np.linalg.norm(v))))
    #Ref:angle between three point, iTechNote, https://itecnote.com/tecnote/python-code-to-calculate-angle-between-three-points-lat-long-coordinates/

#3.CALCULATING HORIZONTAL ALIGNMENT ELEMENTS
def alignmentHorizontal():
    #3.1. INPUT DATA KOORDINAT
    print('(In this example the file name is provided in the folder: coordinate_input_data.csv)')
    filepath = input('Enter the coordinate data point file (latitude, longitude) in CSV format: ')
    data = pd.read_csv(filepath)
    print()
    print('1. Coordinate Data: ')
    print(data)
    #3.2. CALCULATE THE DISTANCE OF 2 POINTS
    print()
    print('2. Distance Calculation Results:')
    outputj=[]
    for i in range(len(data)-1):
        loc1=(data.iloc[i,1], data.iloc[i,2])
        loc2=(data.iloc[i+1,1], data.iloc[i+1,2])
        a=hs.haversine(loc1,loc2)
        b='{:.2f}'.format(a)
        outputj.append(b)
        df=pd.DataFrame(outputj, columns=(['distance']))
        print(i, 'Distance ' +data.iloc[i,0]+ '-' +data.iloc[i+1,0]+ ' = {:.3f}'.format(a),  ' km')
        
    data1=data.assign(jarak=df) #using assign() method to add a new column
    print()
    #print(data1)
    
    #3.3.CALCULATE ANGLES,
    print('3. Angle Calculation Results:')
    outputs=[]
    for i in range(len(data)-2):    
        loc1=(data.iloc[i,1], data.iloc[i,2])
        loc2=(data.iloc[i+1,1], data.iloc[i+1,2])
        loc3=(data.iloc[i+2,1], data.iloc[i+2,2])
        a=np.radians(np.array(loc1))
        b=np.radians(np.array(loc2))
        c=np.radians(np.array(loc3))
        #vector inlat/lon space
        avec=a-b
        cvec=c-b
        #Adjust vectors for changed longitude scale at given lat into 2D
        lat = b[0]
        avec[1] *= math.cos(lat)
        cvec[1] *= math.cos(lat)
        #angle between the vectors in 2 D space
        PI_angle = angle(avec, cvec)
        c='{:.2f}'.format(PI_angle)
        outputs.append(c)
        df1=pd.DataFrame(outputs, columns=(['PI_angle']))
        print(i, 'Angle '+data.iloc[i,0] + '-' +data.iloc[i+1,0] + '-' +data.iloc[i+2,0] + ' = {:.3f}'.format(PI_angle), ' degree')
    data2=data1.assign(PI_angle=df1) #join using assign() method to add a new column
    data3=(data2.iloc[0:3,4]) 
    df=pd.DataFrame(data3)

    print()
    print('(In this example the file name is provided in the folder: geometric_input_data.csv)')
    filepath = input('Enter the planning data file, CSV format: ')
    data_geometric = pd.read_csv(filepath)
    df=df.join(data_geometric, lsuffix="left") 
    print()
    print('4. Calculation Results as Geometric Planning Input Data:')
    print(df)
    df.to_csv("design_input_data.csv") 

    #3.4. CALCULATING GEOMETRIC ELEMENTS 
    output=[]
    for i in range(len(df)):
         Fmax=hitungFmax(df.iloc[i,1])                                # (Vd) 
         Rmin=hitungRmin(df.iloc[i,1], df.iloc[i,2], df.iloc[i,3])    # (Vd, emax, fmax )
         Dmax=hitungDmax(df.iloc[i,1], df.iloc[i,2], df.iloc[i,3])    
         Dd=hitungDd(df.iloc[i,4])                                     
         e_hitung=e_Hitung_def(df.iloc[i,2], Dd, Dmax)                   
         Lst=hitung_Lst(df.iloc[i,1],df.iloc[i,6])
         Lsms=hitung_Lsms(df.iloc[i,1],df.iloc[i,4])
         Lsi=hitung_Lsi(df.iloc[i,1],df.iloc[i,7])
         L=(Lst,Lsms,Lsi)
         Lsort=sorted(L)
         Ls_hitung=Lsort[2]
         Ls_tabel=0.556*df.iloc[i,1]-0.029 #Vr vs Ls, tabel AASHTO, 2011, SE 20/2021 hal 125
         sudutSpiral=sudutSpiral_def(df.iloc[i,4],df.iloc[i,8])
         sudutCircle=float(df.iloc[i,0])-2*sudutSpiral
         Lc=hitung_Lc(sudutCircle, df.iloc[i,4])
         Xs=hitung_Xs(df.iloc[i,8], df.iloc[i,4])
         Ys=hitung_Ys(df.iloc[i,8], df.iloc[i,4])
         P=hitung_P(df.iloc[i,8],df.iloc[i,4]) #SE 20/2021 p.271 if p<0.25 m, full circle, Ls runoff, hal 104 rumus 14
         K=hitung_K(df.iloc[i,8], df.iloc[i,4], sudutSpiral)
         Tt=hitung_Tt(float(df.iloc[i,4]), P, float(df.iloc[i,0]), K)
         Et=hitung_Et(float(df.iloc[i,4]), P, float(df.iloc[i,0]))
         Lt=hitung_Lt(Lc, df.iloc[i,8])
         a=('{:,.3f}'.format(Rmin))
         b=('{:,.3f}'.format(df.iloc[i,4])) # Rd 
         c=('{:,.3f}'.format(e_hitung))
         d=('{:,.3f}'.format(df.iloc[i,5])) #ed
         e=('{:,.3f}'.format(Ls_tabel))
         f=('{:,.3f}'.format(Ls_hitung))
         g=('{:,.3f}'.format(df.iloc[i,8])) #Ls design
         h=('{:,.3f}'.format(Lc))           #L circle
         i=('{:,.3f}'.format(Lt))           #L total
         j=('{:,.3f}'.format(sudutSpiral))
         k=('{:,.3f}'.format(sudutCircle))
         l=('{:,.3f}'.format(Xs))
         m=('{:,.3f}'.format(Ys))
         n=('{:,.3f}'.format(P))
         o=('{:,.3f}'.format(K))
         p=('{:,.3f}'.format(Tt))
         q=('{:,.3f}'.format(Et))
     
         output.append((a,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q))
         # Create the pandas DataFrame
         df1 = (pd.DataFrame(output, columns=(['Rmin', 'e-hit', 'e-design','Ls-tabel','Ls-hit','Ls-renc','Lc','Lt',
                                                'spiral angle','circle angle','Xs','Ys','P','K','Tt','Et'])))  
           
    #3.5 TABEL GABUNGAN OUTPUT
    df1=df.join(df1, lsuffix="_left")

    output=df1.T
 
    print('5. Final Results of Road Geometric Calculations:')
    print(output)
    print('Continued to check the condition of the field and further provisions regarding Horizontal alignment design')
    print()    
    # 3.6. SAVE OUTPUT TO CSV FILE
    output.to_csv("geometric_element_calculation_results.csv")
    print('6. The calculation results are saved in a file: geometric_element_calculation_results.csv')
    print()
    print('7. Image of Road Horizontal Alignment Elements: ')
    img1 = mpimg.imread('Figure Detail of Alignment Horizontal.png') 
    plt.figure(figsize=(20,10))
    plt.axis("off")
    plt.imshow(img1)
    plt.show()
    
    


