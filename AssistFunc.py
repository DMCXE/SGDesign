import numpy as np
class Newton:
    def __init__(self,arr1):
        self.arr1 = arr1
        self.arr1_x = arr1[:,0]
        self.arr1_y = arr1[:,1]
        self.lenth = len(arr1)
        self.fr= self.f()[0]

    def f(self):
        list = [self.arr1_y] #list可以包容不同长度的向量，以区分不同阶
        fx = np.array([self.arr1_y[0]])
        for j in range(0,self.lenth-1):
            list2 = []
            long = len(list[j])
            for i in range(0,long-1):
                l2 = (list[j][i]-list[j][i+1])/(self.arr1_x[i]-self.arr1_x[j+i+1])
                list2.append(l2)
            list.append(list2)
            fx = np.append(fx,list2[0])
        return fx,list

    def num(self,x):
        num = self.arr1_y[0]
        for i in range(1,self.lenth):
            prod = 1
            for j in range(0,i):
                eq = x-self.arr1_x[j]
                prod = prod*eq
            num = num + self.fr[i]*prod
        return num

class CubicSplineFree:
    def __init__(self,arr1):
        self.arr1 = arr1
        self.arr1_x = arr1[:,0]
        self.arr1_y = arr1[:,1]
        self.lenth = len(arr1)
    #hn为x之间的间隔
    def hn(self):
        hnn = np.array([])
        for i in range(0,self.lenth-1):
            hnn =np.append(hnn,self.arr1_x[i+1]-self.arr1_x[i])
        return hnn

    def mu(self):
        mu = np.zeros(1)
        hn = self.hn()
        for i in range(1,len(hn)):
            mu = np.append(mu,hn[i-1]/(hn[i-1]+hn[i]))
        return mu

    def lam(self):
        lam = np.zeros(1)
        hn = self.hn()
        for i in range(1,len(hn)):
            lam = np.append(lam,hn[i]/(hn[i-1]+hn[i]))
        return lam
    #fm为余项，定义与牛顿插值相同
    def fm(self,i):
        return (self.arr1_y[i]-self.arr1_y[i+1])/(self.arr1_x[i]-self.arr1_x[i+1])\
               -(self.arr1_y[i]-self.arr1_y[i-1])/(self.arr1_x[i]-self.arr1_x[i-1])

    def dn(self):
        dn = np.zeros(1)
        hn = self.hn()
        for i in range(1,len(hn)):
            dn = np.append(dn,6*self.fm(i)/(hn[i-1]+hn[i]))
        return dn

    def TDMA(self,a, b, c, d):
        try:
            n = len(d)  #确定长度以生成矩阵
            # 通过输入的三对角向量a,b,c以生成矩阵A
            A = np.array([[0] * n] * n, dtype='float64')
            for i in range(n):
                A[i, i] = b[i]
                if i > 0:
                    A[i, i - 1] = a[i]
                if i < n - 1:
                    A[i, i + 1] = c[i]
            # 初始化代计算矩阵
            c_1 = np.array([0] * n)
            d_1 = np.array([0] * n)
            for i in range(n):
                if not i:
                    c_1[i] = c[i] / b[i]
                    d_1[i] = d[i] / b[i]
                else:
                    c_1[i] = c[i] / (b[i] - c_1[i - 1] * a[i])
                    d_1[i] = (d[i] - d_1[i - 1] * a[i]) / (b[i] - c_1[i - 1] * a[i])
            # x: Ax=d的解
            x = np.array([0] * n)
            for i in range(n - 1, -1, -1):
                if i == n - 1:
                    x[i] = d_1[i]
                else:
                    x[i] = d_1[i] - c_1[i] * x[i + 1]
            #x = np.array([round(_, 4) for _ in x])
            return x
        except Exception as e:
            return e

    def Mn(self):
        a = np.append(self.mu(),0)
        c = np.append(self.lam(),0)
        b = 2*np.ones(self.lenth)
        d = np.append(self.dn(),0)
        Mn = self.TDMA(a,b,c,d)
        return Mn

    def zone(self,x):
        if x < np.min(self.arr1_x): zone = 0
        if x > np.max(self.arr1_x): zone = self.lenth-2
        for i in range(0,self.lenth-1):
            if x-self.arr1_x[i]>=0 and x-self.arr1_x[i+1]<=0:
                zone = i
        return zone

    def num(self,x):
        j = self.zone(x) #zone函数的作用为确定输入量x处于的区间
        M = self.Mn()
        h = self.hn()
        S = M[j]*((self.arr1_x[j+1]-x)**3)/(6*h[j]) \
            + M[j+1]*((x-self.arr1_x[j])**3)/(6*h[j]) \
            + (self.arr1_y[j]-(M[j]*(h[j]**2))/6)*(self.arr1_x[j+1]-x)/h[j] \
            + (self.arr1_y[j+1]-M[j+1]*h[j]**2/6)*(x-self.arr1_x[j])/h[j]
        return S


def Narrow(x):
    arrC = np.array([[0.01,0.618],[0.1,0.632],[0.2,0.644],[0.3,0.659]
                    ,[0.4,0.676],[0.5,0.696],[0.6,0.717],[0.7,0.744]
                    ,[0.8,0.784],[0.9,0.890],[1,1.0]])
    arre = np.array([[0.01, 0.5], [0.1, 0.469], [0.2, 0.431], [0.3, 0.387]
                        , [0.4, 0.343], [0.5, 0.298], [0.6, 0.257], [0.7, 0.212]
                        , [0.8, 0.161], [0.9, 0.079], [1, 0]])
    C = Newton(arrC)
    e = Newton(arre)
    return C.num(x),e.num(x)

def support(x):
    arr = np.array([[0.08561151079136689, 219.22330097087385],
    [0.09136690647482013, 199.4174757281554],
    [0.09712230215827336, 180.7766990291263],
    [0.1028776978417266, 167.96116504854376],
    [0.11151079136690648, 152.23300970873794],
    [0.12302158273381295, 128.34951456310685],
    [0.14028776978417265, 105.04854368932044],
    [0.16618705035971224, 76.50485436893209] ,
    [0.18633093525179859, 60.77669902912626],
    [0.21510791366906476, 44.466019417475735],
    [0.243884892086331, 33.39805825242726],
    [0.28705035971223025, 23.495145631068],
    [0.36187050359712236, 13.009708737864145],
    [0.40791366906474835, 10.097087378640822],
    [0.47410071942446064, 7.184466019417528],
    [0.5201438848920865, 5.436893203883528],
    [0.5633093525179859, 4.854368932038881],
    [0.6007194244604318, 4.271844660194205],
    [0.6323741007194246, 3.1067961165048814],
    [0.6640287769784174, 3.1067961165048814],
    [0.6956834532374102, 2.524271844660234],
    [0.7388489208633096, 1.9417475728155864],
    [0.9978417266187053, 1.359223300970882], [1, 0]])
    num = CubicSplineFree(arr)
    return num.num(x)

