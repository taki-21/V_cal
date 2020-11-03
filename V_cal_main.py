import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d
import time


#データ読み込みの関数(現在txtデータのみ適応)
def dateraeder_pandas(pointfile1,pointfile2,delimiter):
    filename1=pointfile1.split('\\')[len(pointfile1.split('\\'))-1]
    filename2=pointfile2.split('\\')[len(pointfile2.split('\\'))-1]
    ext1=filename1.split('.')[len(filename1.split('.'))-1]
    ext2=filename2.split('.')[len(filename2.split('.'))-1]
    if ext1=='txt' and ext2=='txt':
        time1=pd.read_table(pointfile1,sep=delimiter).values
        time2=pd.read_table(pointfile2,sep=delimiter).values
    else:
        print("This extention is not cover yet.")
    return time1,time2

def triarea_calculation(ptz):
    return np.abs( (ptz[:,0,0]-ptz[:,2,0]) * (ptz[:,1,0]-ptz[:,2,0]) - (ptz[:,0,1]-ptz[:,2,1]) * (ptz[:,1,1]-ptz[:,2,1]))/2

def Heron_formula(ptz):
    a=np.sqrt(np.square(ptz[:,2,0]-ptz[:,0,0])+np.square(ptz[:,2,1]-ptz[:,0,1]))
    b=np.sqrt(np.square(ptz[:,1,0]-ptz[:,0,0])+np.square(ptz[:,1,1]-ptz[:,0,1]))
    c=np.sqrt(np.square(ptz[:,2,0]-ptz[:,1,0])+np.square(ptz[:,2,1]-ptz[:,1,1]))
    s=(a+b+c)/2
    return np.sqrt(s*(s-a)*(s-b)*(s-c))

def height_calculation(ptz):
    return np.min(ptz[:,:,2],axis=1)

def Triangular_pyramid(ptz,i):
    if i==2:
        vol=Heron_formula(ptz)*(np.max(ptz[:,:,2],axis=1)-np.min(ptz[:,:,2],axis=1))/3
    else:
        vol=0
    return vol

def Volume_calculation(times,i):
    points=times[:,:3]
    #print(points)
    x=points[:,0]
    y=points[:,1]
    z=points[:,2]
    #　ドロネー三角分割（三次元）
    print("Delaunay split starts ..")
    start = time.time()
    tri = Delaunay(np.array([x,y]).T)
    #　面を構成する点の座標
    ptz = points[tri.simplices]
    #　面の数
    len = points[tri.simplices].shape[0]
    elapsed_time = time.time() - start
    print ("elapsed_time of Delaunay split:{0}".format(elapsed_time) + "[sec]")
    #体積計算
    vol1 = np.sum(triarea_calculation(ptz)*height_calculation(ptz))
    vol2 = np.sum(Heron_formula(ptz)*height_calculation(ptz))
    vol3 = vol2 + np.sum(Triangular_pyramid(ptz,i))
    # fig = delaunay_plot_2d(tri)
    # fig.savefig('scipy_matplotlib_delaunay'+str(i)+'.png')
    return vol1,vol2,vol3




def main():
    #1時期目のファイル名
    pointfile1="n2_06.txt"
    #2時期目のファイル名
    pointfile2="a2_05_21.txt"
    #各時期の区切り文字
    delimiter=" "

    #データの読み込み
    print("Date Read Start (pandas) ..")
    start = time.time()
    time1,time2=dateraeder_pandas(pointfile1,pointfile2,delimiter)
    elapsed_time = time.time() - start
    print ("elapsed_time of data reder:{0}".format(elapsed_time) + "[sec]")

    #体積計算
    print("Volume calculation starts ..")
    start = time.time()
    vol11,vol12,vol13=Volume_calculation(time1,1)
    vol21,vol22,vol23=Volume_calculation(time2,2)
    elapsed_time = time.time() - start
    print ("elapsed_time of Volume calculation:{0}".format(elapsed_time) + "[sec]")


    print("体積1(sakai):"+str(format(vol11,"f")))
    print("体積1(Heron):"+str(format(vol12,"f")))
    print("体積1(Heron+pyramid):"+str(format(vol13,"f")))
    print("体積2(sakai):"+str(format(vol21,"f")))
    print("体積2(Heron):"+str(format(vol22,"f")))
    print("体積2(Heron+pyramid):"+str(format(vol23,"f")))
    print("差分値(sakai):"+str(format((vol21-vol11),"f")))
    print("差分値(Heron):"+str(format((vol22-vol12),"f")))
    print("差分値(Heron+pyramid):"+str(format((vol23-vol11),"f")))



if __name__ == "__main__":
    main()
