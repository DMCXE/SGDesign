from iapws import iapws97,IAPWS97
import math
import numpy as np
from PipeArrange import Pip_arrangement
import AssistFunc as AFunc
'''
计算过程中换热管根数N指的是截面的全部圆截面个数，
实际的传热管数量应当为N/2
'''
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
    L = 6               #换热管长度

    '''蒸汽参数'''
    D = 126             #蒸汽产量，kg/s
    x = 0.99            #蒸汽干度
    D_d = 0.01*D        #排污量

    SG_yita = 0.99      #蒸汽发生器热效率
    C = 1.09            #传热面积余量系数

    steam = IAPWS97(P = P_loop_2_s,x = x)
    nu_l = steam.Liquid.nu
    nu_g = steam.Vapor.nu

    '''主管道参数'''
    u_1_0_c = 10          #主管道计算流速
    u_2_c = 38            #蒸汽管计算流速
    u_3_c = 4.5           #给水计算流速

    '''下架空间'''
    epsilon_in = 1
    epsilon_out = 1
    epsilon_poi = 1
    delta_down = 0.15 * (1e-3) #mm,下降空间绝对粗燥度

    '''上升空间'''
    '''支撑板开孔采用四叶梅花孔'''
    N_support = 6
    P_c = 221.15*1e5 #临界压力

    '''流量分配孔板'''
    A_u_distribution = 533 * (1e-6)
    a_u_distribution = 216  * (1e-6)

    '''机械强度'''
    sigma1 = 18        #kg/mm2 传热管壁许用应力
    sigma2 = 18        #kg/mm2 下筒体许用应力
    sigma3 = 18        #kg/mm2 上筒体许用应力
    sigma4 = 14.5      #kg/mm2 球形下封头 许用应力
    sigma5 = 1800  # kg/mm2 传热管壁许用应力




class pip():
    #s1 = 32 * 1e-3                      #管心距
    #s2 = s1 * math.sqrt(3)/2 * 1e-3     #形成三角形的高
    s3 = 4 * matr.t                     #两侧间距,最小节圆直径
    l = 6                               #单程直管段长度
    a = 0.25 * math.pi * matr.d_i ** 2

'''在假设管子数量中计算管子的总长度'''
def PipL_tot(N_setting):
    N = N_setting
    pipe = Pip_arrangement(matr.t, matr.t, matr.s3, matr.s3, 0.5*matr.d_o, N, 'Squar')
    pipe.arrangement()
    R_part = pipe.R_part
    L_tot = matr.L*N + R_part
    N_arr = pipe.PipeNum
    R_max = pipe.R_max
    R_min = pipe.R_min
    D_tb = pipe.R
    return L_tot,N_arr,R_part,R_max,R_min,D_tb

def flowstate(Rel,Reg):
    a = 't' if Rel>1000 else 'l'
    b = 't' if Reg > 1000 else 'l'
    return a+b

def Chisholm(X,StateState):
    if StateState == 'tt':   C = 20
    elif StateState == 'lt': C = 15
    elif StateState == 'tl': C = 10
    elif StateState == 'll': C = 5
    Phi_l = 1+C/X+1/(X**2)
    Phi_g = 1+C*X+X**2
    return Phi_l,Phi_g


DTln = matr.Delta_T_loop_1/(math.log((matr.T_loop_1_i-matr.T_2_s)/(matr.T_loop_1_o-matr.T_2_s)))
'''第一步：计算单相对流换热系数'''
AT_loop_1 = 0.5*(matr.T_loop_1_i+matr.T_loop_1_o)
water_loop_1 = IAPWS97(T = AT_loop_1+273.15,P = matr.P_loop_1)
lamb_f = water_loop_1.Liquid.k
mu_f = water_loop_1.Liquid.mu
Pr_f = water_loop_1.Pr
v_f = water_loop_1.Liquid.v
#Re_f = (matr.u_f*matr.d_i)/(v_f*mu_f)

water_loop_1_in = IAPWS97(T=matr.T_loop_1_i+273.15,P = matr.P_loop_1)
water_loop_1_out = IAPWS97(T=matr.T_loop_1_o+273.15,P = matr.P_loop_1)
i1_i = water_loop_1_in.h    #一回路入口焓
i1_o = water_loop_1_out.h   #一回路出口焓

water_loop_2 = IAPWS97(P = matr.P_loop_2_s , x=0.00001)
r = water_loop_2.Hvap           #二回路水的汽化潜热
i_s = water_loop_2.Liquid.h     #二回路饱和水比焓
i_f = IAPWS97(P = matr.P_loop_2_s , T = matr.T_2_in+273.15).Liquid.h    #二回路给水比焓
v_f_2 = IAPWS97(P = matr.P_loop_2_s , T = matr.T_2_in+273.15).Liquid.v    #二回路给水比容
rho_f_2 = 1/v_f_2
v_f_s_2 = IAPWS97(P = matr.P_loop_2_s , x=0.99999).Vapor.v
rho_f_s_2 = 1/v_f_s_2
Q = matr.D * r + (matr.D + matr.D_d)*(i_s-i_f)
G_loop1 = Q / (matr.SG_yita * (i1_i - i1_o))  # 一回路质量流量
#print(G_loop1,i_f,Q)


'''
#D-B公式
c = 0.023 * (lamb_f*Pr_f**(0.3))/((v_f*mu_f)**(0.8))
h_i = c * matr.d_i**(-0.2) * matr.u_f**(0.8)
'''
N = 5000
ess = 1000
while ess > 200:
    A = 0.5* N * pip.a
    u_f = G_loop1 * v_f / A
    Piptot = PipL_tot(N)
    L_tot = Piptot[0]
    N = Piptot[1]
    R_part = Piptot[2]
    #print("N=",N)
    #print("Ltot=",L_tot)
    l_aver = L_tot / N
    '''Gnielinski公式'''
    Re_f = (u_f * matr.d_i) / (v_f * mu_f)
    c = 0.023 * (lamb_f * Pr_f ** (0.3)) / ((v_f * mu_f) ** (0.8))
    h_i = c * matr.d_i ** (-0.2) * u_f ** (0.8)
    '''
    f = ( 1.8 * np.log(Re_f) - 1.5 )**(-2)
    Nu_f = (((f/8)*(Re_f-1000)*Pr_f)/(1+12.7*np.sqrt(f/8)*(Pr_f**(2/3)-1)))*(1+(matr.d_i/l_aver)**(2/3))
    h_i = Nu_f*lamb_f/matr.d_i
    '''
    #print(AT_loop_1)
    '''第二步：计算管壁的导热热阻'''
    R_w = (matr.d_o/(2*matr.lamb_w))*(math.log(matr.d_o/matr.d_i))

    '''第三步：取污垢导热热阻'''
    R_f = matr.R_f

    '''第四步：设定初步k'''
    k = 1000
    es = 10
    while(es>0.01):
        '''第五步：计算传热热流密度'''
        q = k * DTln

        '''第六步：计算沸腾换热系数'''
        h_o = 0.557 * (matr.P_loop_2_s*1e6)**(0.15) * q**(0.7)

        '''第七步：计算k_c'''
        k_c = ((matr.d_o/matr.d_i)*(1/h_i)+R_w+R_f+(1/h_o))**(-1)

        '''第八步：比较、迭代'''
        es = abs((k-k_c)/k)
        k = k_c

    '''第九步：确定换热面积'''
    #print("k=",k)
    F = (Q*10**3)/(k*DTln)
    F_w = F * matr.C
    #print("Fw=", F_w)
    '''迭代确定管长和管数'''
    L_tot_c = F_w / (np.pi * matr.d_o)
    #print("L_tot_c",L_tot_c)
    ess = abs((L_tot - L_tot_c))
    N = N + 2

    print(ess,N,u_f)
'''管束的排布'''
#目前尚却管板高度
R_max = Piptot[3]
R_min = Piptot[4]
D_tb = 2*(Piptot[5] - pip.s3)   #仅为最外侧管子对应直径
L_max = 2*matr.L+np.pi*R_max
L_min = 2*matr.L+np.pi*R_min
L_str = N*matr.L
H_pips = matr.L+R_max


'''主要管道参数'''
d_1_i_c = np.sqrt((4*G_loop1*v_f)/(np.pi*matr.u_1_0_c))
d_1_i = d_1_i_c                         #选定的设计内径
u_1_0 = 4*G_loop1*v_f/(np.pi*d_1_i**2)  #设计流速

G_loop2_o = matr.D              #二回路蒸汽质量流量（产量）
G_loop2_i = matr.D + matr.D_d   #二回路给水质量流量

steam = IAPWS97(P = matr.P_loop_2_s , x=matr.x)
v_2 = steam.Liquid.v*(1-matr.x)+steam.Vapor.v*matr.x

d_2_i_c = np.sqrt(4*G_loop2_o*v_2/(np.pi*matr.u_2_c))
d_2_i = d_2_i_c             #选定
u_2 = 4*G_loop2_o*v_2/(np.pi*d_2_i**2)  #设计流速

water_loop_2_in = IAPWS97(T=matr.T_2_in+273.14,P = matr.P_loop_2_s)
v_3 = water_loop_2_in.v

d_3_i_c = np.sqrt(4*G_loop2_i*v_3/(np.pi*matr.u_3_c))
d_3_i = d_3_i_c             #选定
u_3 = 4*G_loop2_i*v_3/(np.pi*d_3_i**2)  #设计流速

'''U型管摩擦阻力'''
u_f_pip = 1.05*u_f
lamb = 0.3164*Re_f**(-0.25) #一回路水阻力系数

mu_f_wall = IAPWS97(T = 0.5*(AT_loop_1+matr.T_2_s)+273.15,P = matr.P_loop_1).Liquid.mu
phi = (mu_f/mu_f_wall)**(0.14)

#U型管摩擦阻力
P_f = lamb * (matr.L/matr.d_i)*(u_f_pip**2)/(phi*v_f)


D_down = 2*Piptot[5]        #下封头内径 == 管束考虑2e后度直径
F_c = np.pi*(D_down**2)/8   #水室截面积
A_1_i = 0.25*np.pi*d_1_i**2 #进口管截面积
ApF = A_1_i/F_c
epsilon1 = (1-ApF)**2       #突扩系数
v_1_i = IAPWS97(P=matr.P_loop_1,T = matr.T_loop_1_i).v  #一回路入口处比容
rho_1_i = 1/v_1_i                                       #一回路入口处密度

P_1 = 0.5*epsilon1*rho_1_i*u_1_0**2

epsilon2 = (45/90)*(0.131+0.163*2**3.5) #45度弯管局部阻力系数，弯管近似
P_2 = 0.5*epsilon2*rho_1_i*u_1_0**2

A_pips = A/1.05
ApippF = A_pips/F_c
Nar = AFunc.Narrow(ApippF)
#epsilon3 = Nar[1]/Nar[0]**2 + (1/Nar[0]-1)**2      #认为突缩插值
epsilon3 = Nar[1]
P_3 = 0.5*epsilon3*rho_1_i*u_f_pip**2

epsilon4 = 0.5  #180度，来自核动力设备
P_4 = 0.5*epsilon4*(1/v_f)*u_f_pip**2

epsilon5 = (1-ApippF)**2
v_1_o = IAPWS97(P=matr.P_loop_1,T = matr.T_loop_1_o).v  #一回路出口处比容
rho_1_o = 1/v_1_o
P_5 = 0.5*epsilon5*rho_1_o*u_f_pip**2

u_1_o = G_loop1*v_1_o/A_1_i

epsilon6 = epsilon2
P_6 = 0.5*epsilon6*rho_1_o*u_1_o**2

Narr = AFunc.Narrow(ApF)
epsilon7 = Nar[1]
P_7 = 0.5*epsilon7*rho_1_o*u_1_o**2

P_tot_c = P_f + P_1 + P_2 + P_3 + P_4 + P_5 + P_6 + P_7
P_tot = 1.1*P_tot_c

'''二回路水循环阻力计算'''
'''下降空间'''
CR = 3
H_down = H_pips + 1
De = 0.18                 #下降空间当量直径
D_w_o = 2*Piptot[5] + 2*(matr.d_o-matr.d_i) + 2*12*(1e-3)     #管束+外推
D_s_i = D_w_o + De               #下筒体内径
lamb_d = (1.74+2*np.log(De/(2*matr.delta_down)))**(-2)
F_d = 0.25 * np.pi * (D_s_i**2-D_w_o**2)
u_d = CR * matr.D * v_f_2 /F_d
P_d = (lamb_d*H_down/De+matr.epsilon_poi+matr.epsilon_out+matr.epsilon_in)*0.5*rho_f_2*u_d**2

'''上升空间'''
D_w_i = 2*Piptot[5]
n_l = 0 #拉杆数量
F_u = 0.25*np.pi*(D_w_i**2-(2*0.5*N+n_l)*matr.d_o**2)
de = 4*F_u/(np.pi*(D_w_i+(2*0.5*N+n_l)*matr.d_o))
u_o = CR*matr.D*v_f_2/F_u
u_o_l = (CR-1)*matr.D*v_f_2/F_u #出口水相折算速度
u_o_l_ave = 0.5*(u_o+u_o_l)     #水相平均折算速度
u_o_g = matr.D*v_f_s_2/F_u      #出口汽相折算速度
u_o_g_ave = 0.5*u_o_g           #汽相平均折算速度
Re_lo = u_o_l_ave*de/matr.nu_l
Re_go = u_o_g_ave*de/matr.nu_g
flow_state = flowstate(Re_lo,Re_go) #ll,lt,tl,tt
lamb_lo = 0.3164*Re_lo**(-0.25)
lamb_go = 0.3164*Re_go**(-0.25)
P_f_lo = lamb_lo*(pip.l/de)*0.5*(u_o_l_ave**2*rho_f_2)
P_f_go = (1/3)*lamb_go*(pip.l/de)*0.5*(u_o_g_ave**2*rho_f_s_2)
X = np.sqrt(P_f_lo/P_f_go)
Phi_l = Chisholm(X,flow_state)[0]
Phi_g = Chisholm(X,flow_state)[1]
P_f_l = Phi_l*P_f_lo
P_f_g = Phi_g*P_f_go
P_f_up = 0.5*(P_f_l+P_f_g)

'''局部阻力'''

A_u = F_u/(2*0.5*N+n_l)
a_u = (np.sqrt(2)*0.5-0.25*np.pi+0.5)*matr.d_o**2 #四瓣梅花，开孔面积：传热管外径内接正方形延伸
                                                #四个长0.5R,高1.414R矩形面积-弓形面积
auA = a_u/A_u           #相当于流道缩小
epsilon_support = AFunc.support(auA)
P_l_lo_support = matr.N_support*epsilon_support*0.5*rho_f_2*u_o_l_ave**2
P_l_go_support = (1/3)*matr.N_support*epsilon_support*0.5*rho_f_s_2*u_o_g_ave**2
X_support = np.sqrt(P_l_lo_support/P_l_go_support)
Z_R = (0.19+0.92*matr.P_loop_2_s/matr.P_c)**(-1)
K = Z_R + 1/Z_R
Phi_l_support = 1 + K/X_support + 1/(X_support**2)
Phi_g_support = 1 + K*X_support + (X_support**2)
P_l_l_support = Phi_l_support*P_l_lo_support
P_l_g_support = Phi_g_support*P_l_go_support
P_l_support = 0.5*(P_l_l_support+P_l_g_support)

'''弯管区阻力'''
d_b = D_tb - matr.d_o
y_s = 0.2122*d_b
N_rush = math.ceil(y_s/matr.t - 1)
x1x2 = matr.t/matr.d_o
n_rush = 0.43+1.13/x1x2
Re_lo_rush = u_o_l*de/matr.nu_l
Re_go_rush = u_o_g*de/matr.nu_g
lamb_lo_rush = 4*(0.044+0.08*x1x2/((x1x2-1)**n_rush))*Re_lo_rush**(-0.15)
lamb_go_rush = 4*(0.044+0.08*x1x2/((x1x2-1)**n_rush))*Re_go_rush**(-0.15)
P_b_lo = N_rush*lamb_lo_rush*0.5*rho_f_2*u_o_l
P_b_go = (1/3)*N_rush*lamb_go_rush*0.5*rho_f_s_2*u_o_g
X_rush = np.sqrt(P_b_lo/P_b_go)
flow_state_rush = flowstate(Re_lo_rush,Re_go_rush) #ll,lt,tl,tt
Phi_lo_rush = Chisholm(X,flow_state_rush)[0]
Phi_go_rush = Chisholm(X,flow_state_rush)[1]
P_b_l = Phi_lo_rush*P_b_lo
P_b_g = Phi_go_rush*P_b_go
P_b = 0.5*(P_b_l+P_b_g)

'''加速阻力'''
x2 = 1/CR
beta2 = (x2/rho_f_s_2)/(x2/rho_f_s_2+(1-x2)/rho_f_2)
C_a = 0.833 + 0.05*np.log(matr.P_loop_2_s)
phi_2_a = C_a*beta2
G_a = u_o*rho_f_2
P_a = (G_a**2)*(((1-x2)**2)/(rho_f_2*(1-phi_2_a))+(x2**2)/(rho_f_s_2*phi_2_a)-1/rho_f_2)

'''流量分配孔阻力'''
audAu = matr.a_u_distribution/matr.A_u_distribution
epsilon_h = AFunc.support(audAu)
P_h = epsilon_h*0.5*rho_f_2*u_o**2

'''上升空间阻力'''
P_r = P_f_up + P_l_support + P_b + P_a + P_h

'''汽水分离器阻力'''
if CR == 3:
    P_s = 12600
elif CR == 4:
    P_s = 14900
elif CR == 5:
    P_s = 17090

'''循环总阻力'''
P_loop2_tot = P_d + P_r + P_s

'''运动压头计算'''
'''预热段高度计算'''
p_low = matr.P_loop_2_s + (9.8*rho_f_2*H_down)/1e6      #Mpa
H_low = IAPWS97(P=p_low,T = iapws97._TSat_P(p_low)).h
diddp = ((H_low-i_s)*1e3)/((p_low-matr.P_loop_2_s)*1e6)
G_loop_water = CR*matr.D
H_p = (((i_s-i_f)*1e3)/CR + diddp*(9.8*rho_f_2*H_down-P_d))\
      /((2*np.pi*matr.d_o*0.5*N*q)/G_loop_water+diddp*9.8*rho_f_2)

'''运动压头计算'''
H_r1 = H_pips - H_p
H_r2 = H_down - H_pips
x_1_a = x2 * 0.5
beta_1_a = (x_1_a/rho_f_s_2)/((x_1_a/rho_f_s_2)+(1-x_1_a)/rho_f_2)
Phi_1_a = C_a*beta_1_a
Phi_2_a =phi_2_a
P_m1 = (rho_f_2-rho_f_s_2)*9.8*Phi_1_a*H_r1
P_m2 = (rho_f_2-rho_f_s_2)*9.8*Phi_2_a*H_r2
P_m = P_m1+P_m2


'''蒸汽发生器强度计算'''
'''传热管'''
P_pip = matr.P_loop_1 * (10**6) * 9.8          #Kg/m2
S_pip_thi = (P_pip*matr.d_o)/(200*(matr.sigma1*1e6)+0.8*P_pip)   #
phi_gongcha = 1.102
phi_R = 1 + matr.d_o/(4*(2*matr.t))
S_cacu = S_pip_thi*phi_gongcha*phi_R

'''下筒体'''
P_down_sp = matr.P_loop_2_s * (10**6) * 9.8
D_i_down = D_w_o + 2 * 88 * (1e-3)
S_down_sp_c = P_down_sp*D_i_down/(200*matr.sigma2*1e6-1.2*P_down_sp)
S_down_sp = S_down_sp_c
D_o_down = D_i_down + 2*S_down_sp
S_I_down = 50 *1e-3

'''上筒体'''
D_i_up = 3.2    #上筒体内径，3.2m
S_up_sp_c = P_down_sp*D_i_up/(200*matr.sigma3*1e6-1.2*P_down_sp)
S_up_sp = S_up_sp_c

'''球形下封头'''
D_o_sp = D_i_down
S_downsphead_c = P_pip*D_o_sp/(400*matr.sigma4*1e6 + 1.6*P_pip)
S_downsphead = S_downsphead_c

'''管板'''
SdD = S_down_sp/D_i_down
F_TEMA = 1.04
S_plate_c = 0.5*F_TEMA*D_i_down*np.sqrt(P_pip/(matr.sigma5*1e6))
S_plate = S_plate_c
S_duihan = 8*(1e-3)


'''输出'''
print("k",k,"Q",Q,"F",F,"F_w",F_w)
print(Q/F,q)
print(N/2)
print(i_s,i_f)
print(DTln)
print("v_f",v_f,"mu_f",mu_f,"vmu",(v_f*mu_f)**(-1)*matr.d_i)
#print("雷诺数",Re_f)
print(matr.T_2_s)
print(N,L_tot_c)
print("流速",u_f)
'''
pipe = Pip_arrangement(matr.t, matr.t, matr.s3, matr.s3, 0.5*matr.d_o, N, 'Squar')
pipe.arrangement()
pipe.visualize()
'''