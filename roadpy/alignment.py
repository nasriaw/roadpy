''' Roadpy, Geometric Design Horizontal Alignment of Road/Highway
    Author: m nasri aw, nasriaw@gmail.com, Januari 2024
    '''
import math
import pandas as pd
import haversine as hs
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

def calculate_fmax(_a):
    """a = Design Speed, Vd(Kph)"""
    return (0.192-(0.00065*(_a)))

def calculate_rmin(_a,_b,_c):
    """b = emax, c = fmax"""
    return (_a**2)/(127*(_b+_c))

def calculate_d_max(_a,_b,_c):
    """a = Design Speed, Vd(Kph),b = emax, c = fmax"""
    return ((181913.53*(_b+_c))/(_a**2))

def calculate_dd(_d):
    """d = R design, Rd (m)"""
    return 1432.4/(_d)

def e_calculate_def(_b, _dd, _d_max):
    """b = emax, Dd, d_max"""
    return (((-_b*_dd**2)/(_d_max**2))+((2*_b*_dd)/(_d_max)))

def calculate_lst(_a,_e):
    """e = T, Transition curve travel time, tmax=3,00” """
    return (_a/3.6)*_e

def calculate_lsms(_a,_d):
    '''modifikasi shortt: return (((0.022*a**3)/(d*0.4))-(2.727*(a*ed)/0.4))
    C=1,2, ref: SE 20/2021 hal 125'''
    return (((0.0214*_a**3)/(_d*1.2)))

def calculate_lsi(_a,_f):
    """a = Design Speed, Vd(Kph)
    f = delta_i, 0.035 (Vr<80) atau 0,025 (Vr>80) """
    return (0.1-0.02)*_a/(3.6*_f)

def angle_spiral_def(_d, _ls):
    """d = R design, Rd (m), Ls"""
    return 90*_ls/((math.pi)*_d)

def angle_circle_def(_g, _h):
    """g = Angle PI, derajat"""
    _h = angle_spiral_def
    return (_g-(2*(_h)))

_i = angle_circle_def
def calculate_lc(_i, _d):
    """d = R design, Rd (m)"""
    return (_i)*(math.pi*(_d)/180)

def calculate_xs(_ls, _d):
    """d = R design, Rd (m)"""
    return _ls*(1-((_ls**2)/(40*_d**2)))

def calculate_ys(_ls, _d):
    """d = R design, Rd (m)"""
    return (_ls**2)/(6*_d)

def calculate_p(_ls,_d):
    """d = R design, Rd (m)"""
    return _ls**2/(24*_d)

def calculate_k(_ls, _d, _angle_spiral):
    """d = R design, Rd (m)"""
    return _ls-(_ls**2/(40*(_d**2)))-_d*math.sin(_angle_spiral*math.pi/180)

def calculate_tt(_d, _p, _g, _k):
    """d = R design, Rd (m)"""
    return (_d+_p)*math.tan(0.5*_g*math.pi/180)+_k

_j = calculate_p
def calculate_et(_d, _j, _g):
    """d = R design, Rd (m)"""
    return ((_d+_j)/(math.cos(0.5*_g*math.pi/180)))-_d

_k = calculate_lc
def calculate_lt(_k, _ls):
    """ls = length of spiral)"""
    return _k+2*_ls

def angle(_u, _v):
    """return the angle between two vectors in any dimension space, indegrees
    Ref:angle between three point, iTechNote,"""
    return np.degrees(
        math.acos(np.dot(_u,_v)/(np.linalg.norm(_u)*np.linalg.norm(_v))))

def distance(coordinate):
    """Calculate distances and angles between 3 or more points
        input file longitude latitude data coordinate_input_data.csv as coordinate
    """
    output=[]
    for i in range(len(coordinate)-1):
        loc1=(coordinate.iloc[i,1], coordinate.iloc[i,2])
        loc2=(coordinate.iloc[i+1,1], coordinate.iloc[i+1,2])
        a=hs.haversine(loc1,loc2)
        #b='{:.3f}'.format(a)
        #output.append(b)
        df=pd.DataFrame(output, columns=(['distance, km']))
        print(i, 'Distance ' +coordinate.iloc[i,0]+ '-' +coordinate.iloc[i+1,0]+ ' = {:.3f}'.format(a),  ' km')
    #data1=coordinate.assign(distance=df)
    return 
    
def angles(coordinate):
    """Calculate Angles (Point of Intersection, PI)"""
    output=[]
    for i in range(len(coordinate)-2):
        loc1=(coordinate.iloc[i,1], coordinate.iloc[i,2])
        loc2=(coordinate.iloc[i+1,1], coordinate.iloc[i+1,2])
        loc3=(coordinate.iloc[i+2,1], coordinate.iloc[i+2,2])
        a=np.radians(np.array(loc1))
        b=np.radians(np.array(loc2))
        c=np.radians(np.array(loc3))
        avec=a-b
        cvec=c-b
        lat = b[0]
        avec[1] *= math.cos(lat)
        cvec[1] *= math.cos(lat)
        PI_angle = angle(avec, cvec)
        c=f"{(PI_angle):0.3f}"
        output.append(c)
        df1=pd.DataFrame(output, columns=(['PI_angle(degree)']))
        print(i,'Angle '+coordinate.iloc[i,0]+'-'
               +coordinate.iloc[i+1,0]+'-'+coordinate.iloc[i+2,0]
               + ' = {:.3f}'.format(PI_angle), ' degree')
    data2=coordinate.assign(PI_angle=df1)
    #data3=(data2.iloc[0:3,4])
    df1=pd.DataFrame(data2)
    return df1
    
def Horizontal(df):
    """Calculating the geometric elements (alignment) of the highway
        using above functions, file input geometric_input_data.csv join angles, as df
    """
    output=[]
    for i in range(len(df)):
        Rmin=calculate_rmin(df.iloc[i,1],df.iloc[i,2],df.iloc[i,3]) # (Vd,emax,fmax )
        d_max=calculate_d_max(df.iloc[i,1],df.iloc[i,2],df.iloc[i,3])
        Dd=calculate_dd(df.iloc[i,4])
        e_calculate=e_calculate_def(df.iloc[i,2], Dd, d_max)
        lst=calculate_lst(df.iloc[i,1],df.iloc[i,6])
        lsms=calculate_lsms(df.iloc[i,1],df.iloc[i,4])
        lsi=calculate_lsi(df.iloc[i,1],df.iloc[i,7])
        L=(lst,lsms,lsi)
        lsort=sorted(L)
        Ls_calculate=lsort[2]
        Ls_tabel=0.556*df.iloc[i,1]-0.029
        angle_spiral=angle_spiral_def(df.iloc[i,4],df.iloc[i,8])
        j=pd.to_numeric(df.iloc[i,0])
        angle_circle=(j)-2*angle_spiral
        Lc=calculate_lc(angle_circle, df.iloc[i,4])
        Xs=calculate_xs(df.iloc[i,8], df.iloc[i,4])
        Ys=calculate_ys(df.iloc[i,8], df.iloc[i,4])
        P=calculate_p(df.iloc[i,8],df.iloc[i,4])
        K=calculate_k(df.iloc[i,8], df.iloc[i,4], angle_spiral)
        Tt=calculate_tt(float(df.iloc[i,4]), P, float(df.iloc[i,0]), K)
        Et=calculate_et(float(df.iloc[i,4]), P, float(df.iloc[i,0]))
        Lt=calculate_lt(Lc, df.iloc[i,8])
        a=('{:,.3f}'.format(Rmin))
        b=('{:,.3f}'.format(df.iloc[i,4])) 
        c=('{:,.3f}'.format(e_calculate))
        d=('{:,.3f}'.format(df.iloc[i,5])) 
        e=('{:,.3f}'.format(Ls_tabel))
        f=('{:,.3f}'.format(Ls_calculate))
        g=('{:,.3f}'.format(df.iloc[i,8])) 
        h=('{:,.3f}'.format(Lc))           
        i=('{:,.3f}'.format(Lt))           
        j=('{:,.3f}'.format(angle_spiral))
        k=('{:,.3f}'.format(angle_circle))
        l=('{:,.3f}'.format(Xs))
        m=('{:,.3f}'.format(Ys))
        n=('{:,.3f}'.format(P))
        o=('{:,.3f}'.format(K))
        p=('{:,.3f}'.format(Tt))
        q=('{:,.3f}'.format(Et))
        output.append((a,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q))
        df1 = (pd.DataFrame(output, columns=(['Rmin, m', 'e-hit',
                                              'e-design','ls-tabel, m','ls-hit, m',
                                              'ls-renc, m','Lc, m','Lt, m',
                                              'spiral angle','circle angle',
                                              'Xs, m','Ys, m','P, m','K, m','Tt, m','Et, m'])))
    df1=df.join(df1, lsuffix=" ")
    return df1

""" RUN THE PROGRAM """
'''Inputing data files '''
print('(In this example the file name is provided in the folder: coordinate_input_data.csv)')
coordinate=input('Enter the coordinate data point file (latitude, longitude) in CSV format: ')
coordinate = pd.read_csv(coordinate)
print()
print("1. Coordinate input data : ")
print(coordinate)
print()
print("2. Distance between points : ")
a=distance(coordinate)
print(a)
print()
print("3. Point of Intersections, PI : ")
b=angles(coordinate)
""" save data PI"""
b.to_csv("PI_data.csv")
print()
print('(In this example the file name is provided in the folder: geometric_input_data.csv)')
geometric = input('Enter the planning data file, CSV format: ')
df_geometric = pd.read_csv(geometric)
print()
print("4. Horizontal Element Calculation Input Data : ")
c=pd.DataFrame(b.iloc[0:3,3])
df=c.join(df_geometric,  rsuffix=" ")
print(df)
print()
print('5. Results of Geometric Calculations of Highway Horizontal Alignment:')
d=Horizontal(df)
print(d.T)
"""output save a file"""
d.to_csv("geometric_element_calculation_results.csv")
print('Continued to check the condition of the field and further provisions regarding Horizontal alignment design')
print()
print('6. The calculation results are saved in a file: geometric_element_calculation_results.csv')
print()
print('7. Image of Road Horizontal Alignment Elements: ')
img1 = mpimg.imread('Figure Detail of Alignment Horizontal.png')
plt.figure(figsize=(20,10))
plt.axis("off")
plt.imshow(img1)
plt.show()
# End-of-file (EOF)
