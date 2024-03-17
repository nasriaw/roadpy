''' Roadpy, Geometric Design Horizontal Alignment of Road/Highway'''
#1. IMPORT PACKAGES
import math
import pandas as pd
import haversine as hs
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

#2. DEFINING VARIABLES, FUNCTIONS
def hitung_fmax(_a):
    """a = Design Speed, Vd(Kph)"""
    return (0.192-(0.00065*(_a)))

def hitung_rmin(_a,_b,_c):
    """b = emax, c = fmax"""
    return (_a**2)/(127*(_b+_c))

def hitung_d_max(_a,_b,_c):
    """a = Design Speed, Vd(Kph),b = emax, c = fmax"""
    return ((181913.53*(_b+_c))/(_a**2))

def hitung_dd(_d):
    """d = R design, Rd (m)"""
    return 1432.4/(_d)

def e_hitung_def(_b, _dd, _d_max):
    """b = emax, Dd, d_max"""
    return (((-_b*_dd**2)/(_d_max**2))+((2*_b*_dd)/(_d_max)))

def hitung_lst(_a,_e):
    """e = T, Transition curve travel time, tmax=3,00” """
    return (_a/3.6)*_e

def hitung_lsms(_a,_d):
    '''modifikasi shortt: return (((0.022*a**3)/(d*0.4))-(2.727*(a*ed)/0.4))
    C=1,2, ref: SE 20/2021 hal 125'''
    return (((0.0214*_a**3)/(_d*1.2)))

def hitung_lsi(_a,_f):
    """a = Design Speed, Vd(Kph)
    f = delta_i, : 0.035 (Vr<80) atau 0,025 (Vr>80) """
    return (0.1-0.02)*_a/(3.6*_f)

def angle_spiral_def(_d, _ls):
    """d = R design, Rd (m), Ls"""
    return 90*_ls/((math.pi)*_d)

def angle_circle_def(_g, _h):
    """g = Angle PI, derajat"""
    _h = angle_spiral_def
    return (_g-(2*(_h)))

_i = angle_circle_def
def hitung_lc(_i, _d):
    """d = R design, Rd (m)"""
    return (_i)*(math.pi*(_d)/180)

def hitung_xs(_ls, _d):
    """d = R design, Rd (m)"""
    return _ls*(1-((_ls**2)/(40*_d**2)))

def hitung_ys(_ls, _d):
    """d = R design, Rd (m)"""
    return (_ls**2)/(6*_d)

def hitung_p(_ls,_d):
    """d = R design, Rd (m)"""
    return _ls**2/(24*_d)

def hitung_k(_ls, _d, _angle_spiral):
    """d = R design, Rd (m)"""
    return _ls-(_ls**2/(40*(_d**2)))-_d*math.sin(_angle_spiral*math.pi/180)

def hitung_tt(_d, _p, _g, _k):
    """d = R design, Rd (m)"""
    return (_d+_p)*math.tan(0.5*_g*math.pi/180)+_k

_j = hitung_p
def hitung_et(_d, _j, _g):
    """d = R design, Rd (m)"""
    return ((_d+_j)/(math.cos(0.5*_g*math.pi/180)))-_d

_k = hitung_lc
def hitung_lt(_k, _ls):
    """ls = length of spiral)"""
    return _k+2*_ls

def angle(_u, _v):
    """return the angle between two vectors in any dimension space, indegrees
    Ref:angle between three point, iTechNote,"""
    return np.degrees(
        math.acos(np.dot(_u,_v)/ (np.linalg.norm(_u)*np.linalg.norm(_v))))

# pylint: disable=C0301
# pylint: disable=C0103
#3.CALCULATING HORIZONTAL ALIGNMENT ELEMENTS
def Horizontal():
    """ to calculate elements oh horizontal alignment"""
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
    data1=data.assign(jarak=df)
    print()
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
        avec=a-b
        cvec=c-b
        lat = b[0]
        avec[1] *= math.cos(lat)
        cvec[1] *= math.cos(lat)
        #angle between the vectors in 2 D space
        PI_angle = angle(avec, cvec)
        c=f"{(PI_angle):0.3f}"
        outputs.append(c)
        df1=pd.DataFrame(outputs, columns=(['PI_angle']))
        print(i,'Angle '+data.iloc[i,0]+'-'
              +data.iloc[i+1,0]+'-'+data.iloc[i+2,0]
              + ' = {:.3f}'.format(PI_angle), ' degree')
    data2=data1.assign(PI_angle=df1)
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
        Rmin=hitung_rmin(df.iloc[i,1],df.iloc[i,2],df.iloc[i,3]) # (Vd,emax,fmax )
        d_max=hitung_d_max(df.iloc[i,1],df.iloc[i,2],df.iloc[i,3])
        Dd=hitung_dd(df.iloc[i,4])
        e_hitung=e_hitung_def(df.iloc[i,2], Dd, d_max)
        lst=hitung_lst(df.iloc[i,1],df.iloc[i,6])
        lsms=hitung_lsms(df.iloc[i,1],df.iloc[i,4])
        lsi=hitung_lsi(df.iloc[i,1],df.iloc[i,7])
        L=(lst,lsms,lsi)
        lsort=sorted(L)
        Ls_hitung=lsort[2]
        Ls_tabel=0.556*df.iloc[i,1]-0.029
        angle_spiral=angle_spiral_def(df.iloc[i,4],df.iloc[i,8])
        j=pd.to_numeric(df.iloc[i,0])
        angle_circle=(j)-2*angle_spiral
        Lc=hitung_lc(angle_circle, df.iloc[i,4])
        Xs=hitung_xs(df.iloc[i,8], df.iloc[i,4])
        Ys=hitung_ys(df.iloc[i,8], df.iloc[i,4])
        P=hitung_p(df.iloc[i,8],df.iloc[i,4])
        K=hitung_k(df.iloc[i,8], df.iloc[i,4], angle_spiral)
        Tt=hitung_tt(float(df.iloc[i,4]), P, float(df.iloc[i,0]), K)
        Et=hitung_et(float(df.iloc[i,4]), P, float(df.iloc[i,0]))
        Lt=hitung_lt(Lc, df.iloc[i,8])
        a=('{:,.3f}'.format(Rmin))
        b=('{:,.3f}'.format(df.iloc[i,4])) # R rencana
        c=('{:,.3f}'.format(e_hitung))
        d=('{:,.3f}'.format(df.iloc[i,5])) #e-design
        e=('{:,.3f}'.format(Ls_tabel))
        f=('{:,.3f}'.format(Ls_hitung))
        g=('{:,.3f}'.format(df.iloc[i,8])) #Ls rencana
        h=('{:,.3f}'.format(Lc))           #L circle
        i=('{:,.3f}'.format(Lt))           #L total
        j=('{:,.3f}'.format(angle_spiral))
        k=('{:,.3f}'.format(angle_circle))
        l=('{:,.3f}'.format(Xs))
        m=('{:,.3f}'.format(Ys))
        n=('{:,.3f}'.format(P))
        o=('{:,.3f}'.format(K))
        p=('{:,.3f}'.format(Tt))
        q=('{:,.3f}'.format(Et))
        output.append((a,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q))
         # Create the pandas DataFrame
        df1 = (pd.DataFrame(output, columns=(['Rmin', 'e-hit',
                                              'e-design','ls-tabel','ls-hit',
                                              'ls-renc','Lc','Lt',
                                              'spiral angle','circle angle',
                                              'Xs','Ys','P','K','Tt','Et'])))
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
# End-of-file (EOF)
