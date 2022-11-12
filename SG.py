from iapws import iapws97,IAPWS97
import math
from PipeArrange import Pip_arrangement
class matr():
    '''设计书给出的参数'''
    '''一二回路压力参数'''
    P_loop_1_0 = 15.0               #一回路额定工作压力，MPa
    P_loop_1 = P_loop_1_0*1.25      #一回路设计压力，MPa
    P_loop_2_0 = 5.0                #二回路额定工作压力，MPa
    P_loop_2_s = P_loop_2_0*1.25    #二回路设计压力，MPa
    '''一二回路温度参数'''
    T_loop_1_i = 310  # 冷却剂进口温度，摄氏度
    T_loop_1_o = 290  # 冷却剂出口温度，摄氏度
    Delta_T_loop_1 = T_loop_1_i-T_loop_1_o  # 冷却剂进出口温差，摄氏度
    T_2_s = iapws97._TSat_P(P_loop_2_s)-273.15 #二回路饱和温度，摄氏度
    T_2_in = 220       #二回路给水温度，摄氏度
    '''传热管参数'''
    d_i = 17 * (1e-3)  # 传热管内径
    d_o = 19 * (1e-3)  # 传热管外径
    lamb_w = 17.4  # 管壁导热系数，W/mC
    R_f = 2.60*1e-5     #传热管管壁污垢热阻,不锈钢：(5.2-6.9)*1e-5;镍基合金：2.60*1e-5;
                        #污垢厚度0.05mm
    t = 1.40*d_o        #传热管最小节距（采用正方形排列）
    s3 = 32 * (1e-3)    #两侧管间距
    L = 9               #换热管长度

    '''蒸汽参数'''
    D = 126             #蒸汽产量，kg/s
    x = 0.99            #蒸汽干度
    D_d = 0.01*D        #排污量

    SG_yita = 0.99      #蒸汽发生器热效率
    C = 1.09            #传热面积余量系数
    '''假设的参数'''

    u_f = 1         #传热管内水的流动速度，m/s

'''蒸汽发生器热力计算'''

DTln = matr.Delta_T_loop_1/(math.log((matr.T_loop_1_i-matr.T_2_s)/(matr.T_loop_1_o-matr.T_2_s)))
'''第一步：计算单相对流换热系数'''
AT_loop_1 = 0.5 * (matr.T_loop_1_i+matr.T_loop_1_o) #一回路冷却剂平均温度

water_loop_1_in = IAPWS97(T=matr.T_loop_1_i+273.15,P = matr.P_loop_1)
water_loop_1_out = IAPWS97(T=matr.T_loop_1_o+273.15,P = matr.P_loop_1)
i1_i = water_loop_1_in.h    #一回路入口焓
i1_o = water_loop_1_out.h   #一回路出口焓

water_loop_2 = IAPWS97(P = matr.P_loop_2_s , x=0.00001)
r = water_loop_2.Hvap * 1e3 #二回路水的汽化潜热
i_s = water_loop_2.Liquid.h * 1e3 #二回路饱和水比焓
i_f = IAPWS97(P = matr.P_loop_2_s , T = matr.T_2_in+273.15).Liquid.h * 1e3
Q = matr.D * r + (matr.D + matr.D_d)*(i_s-i_f)

G_loop1 = Q/(matr.SG_yita*(i1_i-i1_o)) #一回路质量流量


#实例化IAPWS97
water_loop_1 = IAPWS97(T = AT_loop_1+273.15,P = matr.P_loop_1)
lamb_f = water_loop_1.Liquid.k  #传热系数k
mu_f = water_loop_1.Liquid.mu   #mu
Pr_f = water_loop_1.Pr          #普朗克数
v_f = water_loop_1.Liquid.v
c = 0.023 * (lamb_f*Pr_f**(0.4))/((v_f*mu_f)**(0.8))
h_i = c * matr.d_i**(-0.2) * matr.u_f**(0.8)


'''第二步：计算管壁的导热热阻'''
R_w = (matr.d_o/(2*matr.lamb_w))*(math.log(matr.d_o/matr.d_i))

'''第三步：取污垢导热热阻'''
R_f = matr.R_f

'''第四步：设定初步k'''
k = 1000
es = 10
while(es>0.01):
    F = Q/(k*DTln)
    F_w = F * matr.C
    ess = 10
    N = (math.sqrt(matr.L ** 2 + 2 * F_w) - matr.L) / (2 * math.pi * matr.d_o)
    while(ess>10):
        pipe = Pip_arrangement(matr.t, matr.t, matr.s3, matr.s3, matr.d_o, N, 'Squar')
        pipe.arrangement()
        Nc = pipe.N
        R_part = pipe.R_part
        F1 = Nc * math.pi * matr.d_o * matr.L + R_part*math.pi * matr.d_o
        N = Nc + 1
        ess = abs(F_w-F1)
        #print(ess)
    '''第五步：计算传热热流密度'''
    print("ttttee",N)
    q = k * DTln

    '''第六步：计算沸腾换热系数'''
    h_o = 0.557 * (matr.P_loop_2_s*1e6)**(0.15) * q**(0.7)

    '''第七步：计算k_c'''
    k_c = ((matr.d_o/matr.d_i)*(1/h_i)+R_w+R_f+(1/h_o))**(-1)

    '''第八步：比较、迭代'''
    es = abs((k-k_c)/k)
    k = k_c
    print(k)

'''第九步：确定换热面积'''
F = Q/(k*DTln)
F_w = F * matr.C
print("k",k,"Q",Q,"F",F,"F_w",F_w)
print(Q/F,q)

'''换热直径'''
'''
L = 5
N = (math.sqrt(L**2+2*F)-L)/(2*math.pi*matr.d_o)
'''
print(N)
print(i_s,i_f)
print(DTln)
print("v_f",v_f,"mu_f",mu_f,"vmu",(v_f*mu_f)**(-1)*matr.d_i)
Re = (matr.u_f*matr.d_i)/(v_f*mu_f)
print("雷诺数",Re)
print(matr.T_2_s)

