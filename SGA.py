from iapws import iapws97,IAPWS97
import math
import numpy as np
from PipeArrange import Pip_arrangement
import AssistFunc as AFunc



def print1():
    print("附录一 蒸汽发生器热力计算")
    print("一、热平衡")
    print("1.一回路放热量\t Q\t\t", Q)
    print("2.一回路工作压力\t p1\t\t", matr.P_loop_1)
    print("3.一回路水进口温度\t t1'\t", matr.T_loop_1_i)
    print("4.一回路水出口温度\t t1''\t", matr.T_loop_1_o)
    print("5.一回路水平均温度\t ta\t\t", AT_loop_1)
    print("6.一回路水进口焓\t i1' \t", i1_i)
    print("7.一回路水出口焓\t i1'' \t", i1_o)
    print("8.蒸发器热效率\t η\t\t", matr.SG_yita)
    print("9.一回路水流量\t G1 \t", G_loop1)
    print("10.二回路工作压力\t ps \t", matr.P_loop_2_s)
    print("11.二回路饱和温度\t ts \t", matr.T_2_s)
    print("12.二回路水饱和焓\t is \t", i_s)
    print("13.二回路给水温度\t tf \t", matr.T_2_in)
    print("14.二回路给水焓\t if   \t", i_f)
    print("15.汽化潜热    \t r   \t", r)
    print("16.蒸汽干度    \t x   \t", matr.x)
    print("17.排放系数    \t Cs  \t", "未给定")
    print("18.二回路蒸汽产量\t D  \t", matr.D)
    print("19.二回路排污量\t Dd   \t", matr.D_d)

    print("二、传热计算")
    print("20.传热管外径       \t do     \t", matr.d_o)
    print("21.传热管内径       \t di     \t", matr.d_i)
    print("22.单管流通面积     \t a      \t", pip.a)
    print("23.U型管数目        \t n      \t", 0.5 * N)
    print("24.一回路流通面积    \t A      \t", A)
    print("25.一回路水平均比容  \t nu1     \t", v_f)
    print("26.一回路水流速     \t u1      \t", u_f)
    print("27.一回路水导热系数  \t lambda1 \t", lamb_f)
    print("28.一回路水动力粘度  \t eta1    \t", mu_f)
    print("29.一回路水普朗特数  \t Prf     \t", Pr_f)
    print("30.一回路水雷诺数    \t Ref     \t", Re_f)
    print("31.一回路侧放热系数  \t alpha1   \t", h_i)
    print("32.传热管导热系数    \t lambdaw  \t", matr.lamb_w)
    print("33.传热管壁热阻      \t Rw       \t", R_w)
    print("34.污垢热阻         \t Rf       \t", R_f)
    print("35.二回路侧放热系数  \t alpha2   \t", h_o)
    print("36.传热系数         \t k        \t", k)
    print("37.大端温差         \t A        \t", matr.T_loop_1_i - matr.T_2_s)
    print("38.小端温差         \t A        \t", matr.T_loop_1_o - matr.T_2_s)
    print("39.对数平均温差      \t Dtln     \t", DTln)
    print("40.热负荷           \t q        \t", q)
    print("41.计算传热面积      \t F        \t", F)
    print("42.传热裕度系数      \t         \t", matr.C)
    print("43.设计传热面积      \t F设      \t", F_w)

    print("三、管束结构")
    print("1.传热管总长         \t L总   \t", L_tot)
    print("2.排列方式          \t     \t", "正方形")
    print("3.节距              \t t   \t", matr.t)
    print("4.最小U型管节圆直径   \t D节   \t", 4 * matr.t)
    print("5.实际布管数         \t n   \t", 0.5 * N)
    print("6.管束直径          \t Dtb   \t", D_tb)
    print("7.弯管总长          \t L弯   \t", R_part)
    print("8.直管总长          \t L直   \t", L_str)
    print("9.管束直段高        \t H直   \t", matr.L)
    print("10.管束弯段高       \t H弯   \t", R_max)
    print("11.管束总高         \t Htb   \t", H_pips)
    print("12.传热管实际平均长度 \t l    \t", l_aver)
    print("13.最长管子长        \t Lmax \t", L_max)
    print("14.最短管子长        \t lmin \t", L_min)

    print("四、主要管道内径")
    print("1.冷却剂平均比容      \t nu     \t", v_f)
    print("2.主管道计算流速      \t u'10   \t", matr.u_1_0_c)
    print("3.主管道计算内径      \t d'1i   \t", d_1_i_c)
    print("4.主管道设计内径      \t d1i    \t", d_1_i)
    print("5.主管道设计流速      \t u10    \t", u_1_0)
    print("6.新蒸汽比容          \t v2     \t", v_2)
    print("7.蒸汽管计算流速      \t u'2    \t", matr.u_2_c)
    print("8.蒸汽管计算内径      \t d'2i   \t", d_2_i_c)
    print("9.蒸汽管设计内径      \t d2i    \t", d_2_i)
    print("10.蒸汽管设计流速     \t u2     \t", u_2)
    print("11.二回路给水比容     \t v3     \t", v_3)
    print("12.给水管计算流速     \t u'3    \t", matr.u_3_c)
    print("13.给水管计算内径     \t d'3i   \t", d_3_i_c)
    print("14.给水管设计内径     \t d3i    \t", d_3_i)
    print("15.给水管设计流速     \t u3     \t", u_3)

    print("附录二 蒸汽发生器水力计算")
    print("I 一回路水阻力计算")
    print("一、 U型管内摩擦阻力计算")
    print("1.传热管实际平均长度      \t l      \t", l_aver)
    print("2.当量直径              \t di     \t", matr.d_i)
    print("3.一回路水流量           \t G1     \t", G_loop1)
    print("4.一回路水平均比容        \t v1     \t", v_f)
    print("5.一回路水流速           \t u1     \t", u_f)
    print("6.考虑堵管后流速         \t u'1     \t", u_f_pip)
    print("7.一回路水雷诺数         \t Re      \t", Re_f)
    print("8.摩擦阻力系数           \t lambda  \t", lamb)
    print("9.平均温度下动力粘度      \t eta1    \t", mu_f)
    print("10.壁温下动力粘度        \t eta'1   \t", mu_f_wall)
    print("11.温度修正系数          \t phi     \t", phi)
    print("12.摩擦阻力             \t ΔPf     \t", P_f)

    print("二、 局部阻力计算")
    print("13.下封头内径         \t D1     \t", D_down)
    print("14.水室截面积         \t Fc     \t", F_c)
    print("15.进口管内径         \t d1i    \t", d_1_i)
    print("16.进口管截面积       \t A1     \t", A_1_i)
    print("17.比值              \t        \t", ApF)
    print("18.突扩阻力系数       \t eps1   \t", epsilon1)
    print("19.一回路入水口处比容  \t v1i    \t", v_1_i)
    print("20.一回路入水口处密度  \t rho1i  \t", rho_1_i)
    print("21.入口管内流速       \t u1i    \t", u_1_0)
    print("22.从入口管到水室阻力  \t ΔP1    \t", P_1)
    print("23.水室转弯45度阻力系数\t eps2   \t", epsilon2)
    print("24.水室转弯45度阻力    \t ΔP2    \t", P_2)
    print("25.传热管流道截面     \t A       \t", A)
    print("26.考虑有堵管后截面    \t A'     \t", A_pips)
    print("27.系数              \t        \t", ApippF)
    print("28.传热管入口阻力系数  \t eps3   \t", epsilon3)
    print("29.传热管入口阻力      \t ΔP3    \t", P_3)
    print("30.U型管转180度阻力系数\t eps4   \t", epsilon4)
    print("31.U型管转180度阻力    \t ΔP4    \t", P_4)
    print("32.传热管出口阻力系数   \t eps5   \t", epsilon5)
    print("33.出口处水比容        \t v2     \t", v_1_o)
    print("34.出口处水密度        \t rho2   \t", rho_1_o)
    print("35.传热管出口阻力      \t ΔP5    \t", P_5)
    print("36.出口管内流速       \t u2      \t", u_1_o)
    print("37.水室转弯阻力系数    \t eps6    \t", epsilon6)
    print("38.水室转弯阻力       \t ΔP6      \t", P_6)
    print("39.出口管突缩阻力系数  \t eps7     \t", epsilon7)
    print("40.出口管突缩阻力     \t ΔP7      \t", P_7)

    print("三、 总阻力")
    print("41.总阻力      \t ΔP      \t", P_tot_c)
    print("42.设计阻力    \t ΔP设     \t", P_tot)

def printtable():
    print("II 二回路水循环阻力计算")
    print("一、 下降阻力计算")
    print("1.循环倍率           \t CR     \t",   CR)
    print("2.给水温度           \t tf     \t",   matr.T_2_in)
    print("3.二回路饱和温度      \t ts     \t",   matr.T_2_s)
    print("4.下降空间水比容      \t vd     \t",   v_f_2)
    print("5.下降空间水密度      \t rhod   \t",   rho_f_2)
    print("6.入口阻力系数        \t epsin  \t",   matr.epsilon_in)
    print("7.出口阻力系数        \t epsout \t",   matr.epsilon_out)
    print("8.定位装置阻力系数     \t epsf    \t",   matr.epsilon_poi)
    print("9.下降空间高度        \t H0     \t",   H_down)
    print("10.套筒外径          \t Dw0    \t",   D_w_o)
    print("11.下筒体内径         \t Dsi    \t",   D_s_i)
    print("12.下降空间当量直径    \t De     \t",   De)
    print("13.绝对粗燥度         \t Δ      \t",   matr.delta_down)
    print("14.摩擦系数          \t lambd  \t",   lamb_d)
    print("15.下降空间截面积     \t Fd     \t",   F_d)
    print("16.下降空间水流速     \t ud     \t",   u_d)
    print("17.下降空间阻力       \t ΔPd    \t",   P_d)

    print("二、 上升空间阻力计算")
    print("1、 摩擦阻力")
    print("1.饱和水比容                    \t v'        \t",   v_f_2)
    print("2.饱和水密度                    \t rho'      \t",   rho_f_2)
    print("3.饱和蒸汽比容                  \t v''       \t",   v_f_s_2)
    print("4.饱和蒸汽密度                  \t rho''     \t",   rho_f_s_2)
    print("5.套筒内径                      \t Dw1       \t",   D_w_i)
    print("6.传热管外径                    \t do        \t",   matr.d_o)
    print("7.支撑板定位拉杆数量              \t n'       \t",   n_l)
    print("8.上升空间流通面积               \t Fu        \t",   F_u)
    print("9.上升空间当量直径               \t de        \t",   de)
    print("10.循环速度                     \t uo        \t",   u_o)
    print("11.出口水相折算速度              \t u'o2      \t",   u_o_l)
    print("12.水相平均折算速度              \t u'oa      \t",   u_o_l_ave)
    print("13.出口汽相折算速度              \t u''o2     \t",   u_o_g)
    print("14.汽相平均折算速度              \t u''oa     \t",   u_o_g_ave)
    print("15.水相运动粘度                 \t vl        \t",   matr.nu_l)
    print("16.汽相运动粘度                 \t vg        \t",   matr.nu_g)
    print("17.水相雷诺数                   \t Relo      \t",   Re_lo)
    print("18.汽相雷诺数                   \t Rego      \t",   Re_go)
    print("19.判别流态                     \t          \t",   flow_state)
    print("20.管束直段高                   \t Hs        \t",   pip.l)
    print("21.水相摩擦系数                 \t lamblo    \t",   lamb_lo)
    print("22.汽相摩擦系数                 \t lambgo    \t",   lamb_go)
    print("23.按折算速度计算的水相摩擦阻力    \t ΔPflo     \t",   P_f_lo)
    print("24.按折算速度计算的汽相摩擦阻力    \t ΔPfgo     \t",   P_f_go)
    print("25.参量X                       \t           \t",   X)
    print("26.Phi_l**2                   \t           \t",   Phi_l)
    print("27.Phi_g**2                   \t           \t",   Phi_g)
    print("28.水相摩擦阻力                 \t ΔPfl      \t",   P_f_l)
    print("29.汽相摩擦阻力                 \t ΔPfg      \t",   P_f_g)
    print("30.摩擦阻力                    \t ΔPf       \t",   P_f_up)

    print("2、 局部阻力")
    print("1.支撑板数目                   \t N        \t",   matr.N_support)
    print("2.上升流道单元面积              \t Au      \t",   A_u)
    print("3.支撑板单元开孔面积            \t au       \t",   a_u)
    print("4.面积比                      \t          \t",   auA)
    print("5.局部阻力系数                 \t epsl     \t",   epsilon_support)
    print("6.按折算速度计算的水相局部阻力    \t ΔPllo    \t",   P_l_lo_support)
    print("7.按折算速度计算的汽相局部阻力    \t ΔPlgo    \t",   P_l_go_support)
    print("8.参量X                       \t           \t",   X_support)
    print("9.二回路工作压力               \t Ps        \t",   matr.P_loop_2_s)
    print("10.临界压力                   \t Pc        \t",   matr.P_c)
    print("11.参数ZR                    \t           \t",   Z_R)
    print("12.参数K                     \t           \t",   K)
    print("13.Phi_l                    \t           \t",   Phi_l_support)
    print("14.Phi_g                    \t           \t",   Phi_g_support)
    print("15.水相局部阻力               \t ΔPll      \t",   P_l_l_support)
    print("16.汽相局部阻力               \t ΔPlg      \t",   P_l_g_support)
    print("17.局部阻力                  \t ΔPl      \t",   P_l_support)

    print("3、 弯管区阻力")
    print("1.管束弯头最大节圆直径        \t db        \t",   d_b)
    print("2.弯管区重心至圆心距离        \t ys        \t",   y_s)
    print("3.节距                     \t t         \t",   matr.t)
    print("4.计算冲刷排数              \t N'        \t",   N_rush)
    print("5.系数                     \t x1=x2     \t",   x1x2)
    print("6.系数                     \t n         \t",   n_rush)
    print("7.水相雷诺数                \t Relo      \t",   Re_lo_rush)
    print("8.汽相雷诺数                \t Rego      \t",   Re_go_rush)
    print("9.水相摩擦阻力系数           \t lamb_lo   \t",   lamb_lo_rush)
    print("10.汽相摩擦阻力系数          \t lamb_go   \t",   lamb_go_rush)
    print("11.全水相阻力               \t ΔPblo     \t",   P_b_lo)
    print("12.全汽相阻力               \t ΔPbgo     \t",   P_b_go)
    print("13.参量X                   \t           \t",   X_rush)
    print("14.Phi_l                  \t           \t",   Phi_lo_rush)
    print("15.Phi_g                  \t           \t",   Phi_go_rush)
    print("16.水相阻力                \t ΔPbl      \t",   P_b_g)
    print("17.汽相阻力                \t ΔPbg      \t",   P_b_l)
    print("18.弯管区阻力              \t ΔPb       \t",   P_b)

    print("4、加速阻力")
    print("1.管束出口质量含汽率        \t x2        \t",   x2)
    print("2.管束出口体积含汽率        \t beta2     \t",   beta2)
    print("3.系数                    \t C         \t",   C_a)
    print("4.管束出口截面含汽率        \t phi2      \t",   phi_2_a)
    print("5.质量流速                 \t G        \t",   G_a)
    print("6.加速阻力                 \t ΔPa      \t",   P_a)

    print("5、流量分配阻力")
    print("1.单元面积               \t Au'       \t",   matr.A_u_distribution)
    print("2.单元开孔面积            \t au'       \t",   matr.a_u_distribution)
    print("3.系数                   \t          \t",   audAu)
    print("4.阻力系数               \t epsh      \t",   epsilon_h)
    print("5.孔板局部阻力            \t Ph        \t",   P_h)

    print("6、上升空间阻力")
    print("1.上升空间阻力            \t ΔPr       \t",   P_r)

    print("三、汽水分离器阻力")
    print("1.汽水分离器阻力            \t ΔPs       \t",   P_s)

    print("四、循环总阻力")
    print("1.循环总阻力            \t ΔP总       \t",   P_loop2_tot)


    print("III、运动压头计算")
    print("一、预热段高度计算")
    print("1.循环倍率              \t CR        \t",   CR)
    print("2.二回路给水焓           \t if        \t",   i_f)
    print("3.二回路饱和水焓         \t is        \t",   i_s)
    print("4.液面高度              \t Ho        \t",   H_down)
    print("5.下降空间水密度         \t rhod      \t",   rho_f_2)
    print("6.下降空间下端压力       \t Plow      \t",   p_low)
    print("7.plow压力下的饱和水焓   \t isl       \t",   H_low)
    print("8.Δi/ΔP               \t           \t",   diddp)
    print("9.热负荷               \t q         \t",   q)
    print("10.循环水量             \t G        \t",   G_loop_water)
    print("11.预热段高度           \t Hp       \t",   H_p)

    print("二、运动压头计算")
    print("1.蒸发段高度                \t Hr1        \t",   H_r1)
    print("2.管束上方区段高度           \t Hr2        \t",   H_r2)
    print("3.蒸发段平均质量含汽率        \t x1a        \t",   x_1_a)
    print("4.蒸发段平均体积含汽率        \t beta1a     \t",   beta_1_a)
    print("5.蒸发段平均截面含汽率        \t phi1a      \t",   Phi_1_a)
    print("6.管束上方区段平均截面含汽率   \t phi2a      \t",   Phi_2_a)
    print("7.蒸发段运动压头             \t Pm1        \t",   P_m1)
    print("8.管束上方区段压头           \t Pm2        \t",   P_m2)
    print("9.运动压头                  \t Pm         \t",   P_m)

    print("附录三、蒸汽发生器强度计算")
    print("一、传热管")
    print("1.设计压力            \t P设1      \t",   P_pip)
    print("2.许用应力            \t [σ1]       \t",   matr.sigma1)
    print("3.管子外径            \t do        \t",   matr.d_o)
    print("4.直管计算壁厚         \t S''       \t",   S_pip_thi)
    print("5.负公差修正系数        \t phi      \t",   phi_gongcha)
    print("6.弯曲减薄系数         \t phiR      \t",   phi_R)
    print("7.计算壁厚             \t S'       \t",   S_cacu)
    print("8.设计壁厚             \t S        \t",   S_cacu)

    print("二、下筒体")
    print("1.设计压力            \t P设2      \t",   P_down_sp)
    print("2.许用应力            \t [σ2]      \t",   matr.sigma2)
    print("3.筒体内径            \t Di下      \t",   D_i_down)
    print("4.计算壁厚            \t S'        \t",   S_down_sp_c)
    print("5.设计壁厚            \t S         \t",   S_down_sp)
    print("6.筒体外径            \t Do下      \t",   D_o_down)
    print("7.与管板连接壁厚       \t Si        \t",   S_I_down)

    print("三、上筒体")
    print("1.设计压力            \t P设2      \t",   P_down_sp)
    print("2.许用应力            \t [σ3]      \t",   matr.sigma3)
    print("3.筒体内径            \t Di上      \t",   D_i_up)
    print("4.计算壁厚            \t S'        \t",   S_up_sp_c)
    print("5.设计壁厚            \t S         \t",   S_up_sp)

    print("四、球形下封头")
    print("1.设计压力            \t P设1      \t",   P_down_sp)
    print("2.许用应力            \t [σ4]      \t",   matr.sigma4)
    print("3.球型封头外径         \t D0        \t",   D_o_sp)
    print("4.计算壁厚            \t S'        \t",   S_downsphead_c)
    print("5.设计壁厚            \t S         \t",   S_downsphead)

    print("五、管板")
    print("1.设计压力            \t P设1      \t",   P_pip)
    print("2.许用应力            \t [σ1]      \t",   matr.sigma5)
    print("3.承压部分直径         \t do        \t",   D_i_down)
    print("4.筒体根部壁厚与直径比  \t           \t",   SdD)
    print("5.系数                \t F        \t",   F_TEMA)
    print("6.计算壁厚             \t S'       \t",   S_plate_c)
    print("7.设计壁厚             \t S        \t",   S_plate)
    print("8.堆焊层厚度           \t S'        \t",   S_duihan)


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
    T_2_s = iapws97._TSat_P(P_loop_2_0)-273.15 #二回路饱和温度，摄氏度
    T_2_in = 220       #二回路给水温度，摄氏度
    '''传热管参数'''
    #d_i = 19.4 * (1e-3)  # 传热管内径
    d_o = 22 * (1e-3)  # 传热管外径

    lamb_w = 17.4       # 管壁导热系数，W/mC
    R_f = 2.60*1e-5     #传热管管壁污垢热阻,不锈钢：(5.2-6.9)*1e-5;镍基合金：2.60*1e-5;
                        #污垢厚度0.05mm
    t = 1.40*d_o        #传热管最小节距（采用正方形排列）
    s3 = 32 * (1e-3)    #两侧管间距
    L = 4.4               #换热管长度

    '''蒸汽参数'''
    D = 126             #蒸汽产量，kg/s
    x = 0.99            #蒸汽干度
    D_d = 0.01*D        #排污量

    SG_yita = 0.99      #蒸汽发生器热效率
    C = 1.09            #传热面积余量系数

    steam = IAPWS97(P = P_loop_2_0,x = x)
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

    P_pip = P_loop_1 * (10 ** 6) * 9.8  # Kg/m2
    S_pip_thi = (P_pip * d_o) / (200 * (sigma1 * 1e6) + 0.8 * P_pip)  #
    phi_gongcha = 1.102
    phi_R = 1 + d_o / (4 * (2 * t))
    S_cacu = S_pip_thi * phi_gongcha * phi_R

    d_i = d_o - 2 * S_cacu


class pip():
    #s1 = 32 * 1e-3                      #管心距
    #s2 = s1 * math.sqrt(3)/2 * 1e-3     #形成三角形的高
    s3 = 4 * matr.t                     #两侧间距,最小节圆直径
    l = 6                               #单程直管段长度
    a = 0.25 * math.pi * matr.d_i ** 2

'''在假设管子数量中计算管子的总长度'''
def PipL_tot(N_setting):
    N = N_setting
    pipe = Pip_arrangement(matr.t, matr.t, matr.s3, 20*(1e-3), 0.5*matr.d_o, N, 'Squar')
    pipe.arrangement()
    R_part = pipe.R_part
    L_tot = matr.L*N + R_part
    N_arr = pipe.PipeNum
    R_max = pipe.R_max
    R_min = pipe.R_min
    D_tb = pipe.R
    return L_tot,N_arr,R_part,R_max,R_min,D_tb




DTln = matr.Delta_T_loop_1/(math.log((matr.T_loop_1_i-matr.T_2_s)/(matr.T_loop_1_o-matr.T_2_s)))
'''第一步：计算单相对流换热系数'''
AT_loop_1 = 0.5*(matr.T_loop_1_i+matr.T_loop_1_o)
water_loop_1 = IAPWS97(T = AT_loop_1+273.15,P = matr.P_loop_1_0)
lamb_f = water_loop_1.Liquid.k
mu_f = water_loop_1.Liquid.mu
Pr_f = water_loop_1.Pr
v_f = water_loop_1.Liquid.v
#Re_f = (matr.u_f*matr.d_i)/(v_f*mu_f)

water_loop_1_in = IAPWS97(T=matr.T_loop_1_i+273.15,P = matr.P_loop_1_0)
water_loop_1_out = IAPWS97(T=matr.T_loop_1_o+273.15,P = matr.P_loop_1_0)
i1_i = water_loop_1_in.h    #一回路入口焓
i1_o = water_loop_1_out.h   #一回路出口焓

water_loop_2 = IAPWS97(P = matr.P_loop_2_0 , x=0.00001)
r = water_loop_2.Hvap           #二回路水的汽化潜热
i_s = water_loop_2.Liquid.h     #二回路饱和水比焓
i_f = IAPWS97(P = matr.P_loop_2_0 , T = matr.T_2_in+273.15).Liquid.h    #二回路给水比焓
v_f_2 = IAPWS97(P = matr.P_loop_2_0 , T = matr.T_2_in+273.15).Liquid.v    #二回路给水比容
rho_f_2 = 1/v_f_2
v_f_s_2 = IAPWS97(P = matr.P_loop_2_0 , x=0.99999).Vapor.v
rho_f_s_2 = 1/v_f_s_2
Q = matr.D * r + (matr.D + matr.D_d)*(i_s-i_f)
G_loop1 = Q / (matr.SG_yita * (i1_i - i1_o))  # 一回路质量流量
#print(G_loop1,i_f,Q)


'''
#D-B公式
c = 0.023 * (lamb_f*Pr_f**(0.3))/((v_f*mu_f)**(0.8))
h_i = c * matr.d_i**(-0.2) * matr.u_f**(0.8)
'''
N = 3300
ess = 1000
while ess > 150:
    A = 0.5 *  N * pip.a
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
        h_o = 0.557 * (matr.P_loop_2_0*1e6)**(0.15) * q**(0.7)

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
D_tb = 2*(Piptot[5] - 20*(1e-3))   #仅为最外侧管子对应直径

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

steam = IAPWS97(P = matr.P_loop_2_0 , x=matr.x)
v_2 = steam.Liquid.v*(1-matr.x)+steam.Vapor.v*matr.x

d_2_i_c = np.sqrt(4*G_loop2_o*v_2/(np.pi*matr.u_2_c))
d_2_i = d_2_i_c             #选定
u_2 = 4*G_loop2_o*v_2/(np.pi*d_2_i**2)  #设计流速

water_loop_2_in = IAPWS97(T=matr.T_2_in+273.14,P = matr.P_loop_2_0)
v_3 = water_loop_2_in.v

d_3_i_c = np.sqrt(4*G_loop2_i*v_3/(np.pi*matr.u_3_c))
d_3_i = d_3_i_c             #选定
u_3 = 4*G_loop2_i*v_3/(np.pi*d_3_i**2)  #设计流速

'''U型管摩擦阻力'''
u_f_pip = 1.05*u_f
lamb = 0.3164*Re_f**(-0.25) #一回路水阻力系数

mu_f_wall = IAPWS97(T = 0.5*(AT_loop_1+matr.T_2_s)+273.15,P = matr.P_loop_1_0).Liquid.mu
phi = (mu_f/mu_f_wall)**(0.14)

#U型管摩擦阻力
P_f = lamb * (matr.L/matr.d_i)*(u_f_pip**2)/(phi*v_f)


D_down = 2*Piptot[5]        #下封头内径 == 管束考虑2e后度直径
F_c = np.pi*(D_down**2)/8   #水室截面积
A_1_i = 0.25*np.pi*d_1_i**2 #进口管截面积
ApF = A_1_i/F_c
epsilon1 = (1-ApF)**2       #突扩系数
v_1_i = IAPWS97(P=matr.P_loop_1_0,T = matr.T_loop_1_i).v  #一回路入口处比容
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
v_1_o = IAPWS97(P=matr.P_loop_1_0,T = matr.T_loop_1_o).v  #一回路出口处比容
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

print1()

'''二回路水循环阻力计算'''
'''下降空间'''
P_m_array = np.zeros(3)
Delta_P_array = np.zeros(3)

for i in range(3,6):
    CR = i
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
    flow_state = AFunc.flowstate(Re_lo,Re_go) #ll,lt,tl,tt
    lamb_lo = 0.3164*Re_lo**(-0.25)
    lamb_go = 0.3164*Re_go**(-0.25)
    P_f_lo = lamb_lo*(pip.l/de)*0.5*(u_o_l_ave**2*rho_f_2)
    P_f_go = (1/3)*lamb_go*(pip.l/de)*0.5*(u_o_g_ave**2*rho_f_s_2)
    X = np.sqrt(P_f_lo/P_f_go)
    Phi_l = AFunc.Chisholm(X,flow_state)[0]
    Phi_g = AFunc.Chisholm(X,flow_state)[1]
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
    Z_R = (0.19+0.92*matr.P_loop_2_0/matr.P_c)**(-1)
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
    flow_state_rush = AFunc.flowstate(Re_lo_rush,Re_go_rush) #ll,lt,tl,tt
    Phi_lo_rush = AFunc.Chisholm(X,flow_state_rush)[0]
    Phi_go_rush = AFunc.Chisholm(X,flow_state_rush)[1]
    P_b_l = Phi_lo_rush*P_b_lo
    P_b_g = Phi_go_rush*P_b_go
    P_b = 0.5*(P_b_l+P_b_g)

    '''加速阻力'''
    x2 = 1/CR
    beta2 = (x2/rho_f_s_2)/(x2/rho_f_s_2+(1-x2)/rho_f_2)
    C_a = 0.833 + 0.05*np.log(matr.P_loop_2_0)
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
    p_low = matr.P_loop_2_0 + (9.8*rho_f_2*H_down)/1e6      #Mpa
    H_low = IAPWS97(P=p_low,T = iapws97._TSat_P(p_low)).h
    diddp = ((H_low-i_s)*1e3)/((p_low-matr.P_loop_2_0)*1e6)
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

    printtable()

    P_m_array[i-3] = P_m
    Delta_P_array[i-3] = P_loop2_tot

print(P_m_array)
print(Delta_P_array)

CRxCacu = AFunc.CReq(P_m_array,Delta_P_array)
CRx = CRxCacu.CRx()
print(CRx)
CRxCacu.visualize()

pipe = Pip_arrangement(matr.t, matr.t, matr.s3, 20*(1e-3), 0.5*matr.d_o, N, 'Squar')
pipe.arrangement()
pipe.visualize()




"""
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
"""