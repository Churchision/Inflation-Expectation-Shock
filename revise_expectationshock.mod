//-------------------------------------------------------------------------
//      Fasani, Mumtaz, Rossi (2022)
//      "Monetary Policy Uncertainty and Firm Dynamics"
//      compatible with Dynare 4.4 onwards
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//    VARIABLE DEFINITIONS
//-------------------------------------------------------------------------
parameters  
rhopi,
sigma_pi,
sigma_a,

   beta,                   %贴现因子
sigma,                     %相对风险厌恶系数
chi,                       %chi劳动的相对效用权重
h,                         %消费习惯h
phil,                      %sigma_L劳动供给弹性
deltak,                    %delta^k资本折旧率
gammai,                    %phi_k投资调整成本
gammau1,                   %gamma1资本利用调整函数的参数
gammau2,                   %gamma2资本利用调整函数的参数
 gammaw,                   %phi_w名义工资调整成本系数
   fe,                     %牌照价格
 big_theta,                %Theta^e 企业进入拥堵外部性的系数
big_thetaxc,               %Theta^x 企业退出拥堵外部性的系数
epsw,                      %theta_w
varsigma,                  %zeta_e 企业进入拥堵外部性的弹性
varsigmaxc,                %zeta_x 企业退出拥堵外部性的弹性
tau,                       %1-退出企业向家庭回扣的固定进入成本份额
    xie,                   %生产率Pareto分布的参数
%   xiess,
 zmin,                     %生产率Pareto分布的参数
theta,                     %theta_p中间产品的替代弹性参数，加成率的参数
gammap,                    %phi_p 价格调整成本系数
alfa,                      %alfa资本份额
phii,                      %phi_R以下四个都是taylor规则参数
   phipai,                 %phi_pi
phiy,                      %phi_y
phidy,                     %phi_dy
 taxss,
 nxxcy, 
neecy, 
nen, 
   rhob,
 rhoi, 
rhow, 
rhoe, 
rhox,
 rhoa, 
rhoii, 
%rhosigma, 
rhod,
   rhop, 
V_normalization, 
alphaEZ, 
ckappa
%rhoxie
%rhosigma_xie
%sigmabar_ii
    ;

var 

paie,
zeps_pi,

rk,                  %r^k资本服务的租金
rir,                    %实际利率r
ii,                     %名义利率i
q,                      %Tobin q
cu,                     %资本利用程度
adj_cu,                 %资本利用率调整成本
pai,                    %pi
pai_w,                  %pi^w
pai_is,                 %【？】
vtilde,                 %v(\tilde{z})现存企业的价值
vhat,                   %v(\bar{z})临界企业的价值
eta,                    %企业退出率
lv,                     %企业清算价值
divtilde,               %j_t(\tilde{z})现存企业利润分红
divhat,                 %j_t(\bar{z])临界企业利润分红
div,                    %现存企业的总利润
ne,                     %新进入企业数量，模型中时间下标是t-1，但对应代码中时间下标t
n,                      %t期企业数量
   nx,                  %退出企业数量
ec,                     %进入拥堵外部性
xc,                     %退出拥堵外部性
zhat,                   %退出的临界生产率
ztilde,                 %现存企业的平均生产率
mctilde,                %mc(\tilde{z})现存企业的边际成本
muptilde,               %mu 中间品厂商的加成率
   rhotilde,            %\tilde{\rho}平均企业的最优价格
rho_I,                  %rho^I
tax,                    %税收=-转移支付
ldata, 
invnedata,
   tfp_prd,
lab_share,
pro_share,
P_long,
R_long,
P_free,
R_free,
T_premium,
T_spread,
R_equity
   E_premium,
vtilde_fake, 
R_equity_fake,
E_premium_fake,
vtilde_cons,
   R_equity_cons,
E_premium_cons,
zeps_b,                  %贴现因子冲击【？】
zeps_i,                  %投资效率冲击【？】
zeps_w,     
   zeps_e,
zeps_x,
zeps_a,
zeps_ii,
%sigma_eps_ii,
zeps_p,
zeps_d,
gydata,
gcdata,
   ginvdata,
ginvnedata,
ginvtotdata,
gldata,
gnxdata,
gnedata,
y,                     %大Y
ymean,                 %平均每家企业的产出
   c,                  %消费
inv,                   %投资
k,                     %资本存量
ks,                    %资本服务
ytilde,                %y(\tilde{z})现存企业的产出
ydata,                 %Taylor规则中的y
cdata, 
   invdata,
tfpdata,
invtotdata,
l,                     %劳动
lambda,                %消费边际效用
nu,                    %原文psi
V,
E_t_V_tp1_1_minus_alpha,
M,                     %理性预期算子都会有这个【？】
   mrs,                %消费劳动边际替代率
w,                     %工资
lprod,
A,                     %可能是全要生产率【？】
%xie,                   %生产率Pareto分布的参数
%sigma_xie,
R_equity_var,
R_equity_fake_var
   sdf
 R_capital ,
K_premium,
 R_equityhat,
 AlogR ,
Alogpai ,
AlogRR,
E_y,var_y,
var_pai

   ;

varexo
    eps_ii, eps_b, eps_i, eps_w, eps_e, eps_x, eps_a, eps_p, eps_d, eps_pi %eps_xie, eta_xie,
   ;

//-------------------------------------------------------------------------
//                         CALIBRATION
//-------------------------------------------------------------------------
%gammai=1.43450435425568;
%rhoii=0.000730575183330449;

%gammap=9.88696172422216;
%gammaw=221.069634447268;
%phii=0.531553062896770;
%phipai=3.56725097546153;
%phiy=0;
%phidy=0.0167559402877902;

/*
// uncomment the following calibration for robustness
%gammai=0.289698779974552;%
%rhoii=1.72917746457380e-11;%

%gammap=29.6624746232485;%
%gammaw=693.491504323785;%
%phii=0.5;%
%phipai=2.5;%
%phiy=0;
%phidy=0.05;%
*/

//*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%Ascari(2023)%%%%%%%%%%%%%%%%%%%%%%%%
rhopi=0.9;
rhoa=0.9;
sigma_a=0.1;
sigma_pi=0.1;


gammai=0.5;%

gammap=40;%
gammaw=40;%
phii=0.75;%
phipai=2.5;%
phiy=0.01;
phidy=0.05;%


phil=5.0100;       %5?
epsw=4.3;
theta=4.3;

varsigma=2.0001;
varsigmaxc=1.001;
h=0.6;

gammau2=(1-0.54)/0.54 ;
beta=0.9967;           
xie=6.51;
sigma=1.5;
zmin=1;
alfa=0.33;         
deltak=.0067;
%rhosigma_xie=0; 
alphaEZ=0;      %没有Epstein-Zin效用

ySS = 1;
cuSS = 1;
paiSS = 1.0017;   
nSS = 1;
ASS = 1;
qSS = 1;                                                                                      
pai_wSS = 1;                                                                               
adj_cuSS = 0;                                                                                

rhoii=0.5;%rhoii=1.72917746457380e-11

rhob=0.5;         
rhoi=0.5;         
rhow=0.5;
rhoe=0.5;
rhox=0.5; 
rhop=0.5;
rhod=0.5;
%rhoxie=0.5; 
  
zeps_bSS=0;         
zeps_iSS=0;         
zeps_wSS=0;
zeps_eSS=0;
zeps_xSS=0; 
zeps_aSS=0;
zeps_iiSS=0;
zeps_pSS=0;
zeps_sigmaSS=0;
zeps_dSS=0;
%zeps_xieSS=0; 
sigma_eps_iiSS=0;
%sigma_xieSS = 1;

nen=0.098901; %（是进入率也是退出率 nx/n=eta/(1-eta) Fasani的季度死亡率为3%，故月度退出概率=0.030/3）
neecy=0.016;    %ec/y稳态值
nxxcy=0.0121;   %xc/y稳态值

//STEADY STATE

rirSS=1/beta-1;
iiSS=((1+rirSS)*paiSS)-1;
pai_isSS=paiSS;

invk=deltak;                                                                                  
rkSS=1/beta-(1-deltak);                                                                      
gammau1 = rkSS;                                                                              
muptildeSS=theta/(theta-1); 
    
etaSS=1-1/(1+nen);                                                                                  
neSS=nen*nSS;
nxSS=etaSS*(nSS+neSS);
nxn=nen;
ecSS=neecy*ySS/neSS; 
xcSS=nxxcy*ySS/nxSS;
 
zhatSS=(1-etaSS)^(-1/xie)*zmin; %as eta=1-(zmin/zhat)^xie;
rho_ISS=(theta-1)/theta;
rhotildeSS=nSS^(1/(theta-1))*rho_ISS;     

ztildeSS=(xie/(xie+1-theta))^(1/(theta-1))*zhatSS;
mctildeSS=rhotildeSS/muptildeSS;                                                        

big_theta=ecSS/(nen^varsigma);    
big_thetaxc=xcSS/(nxSS/nSS)^varsigmaxc;

ymeanSS=ySS/nSS;                                                                                   
lks=(rkSS/(ztildeSS*mctildeSS*alfa))^(1/(1-alfa));                                        
wSS=ztildeSS*mctildeSS*(1-alfa)*(lks)^(-alfa);                                              
ksSS=1/(ztildeSS)*(((1-alfa)/alfa)*rkSS/wSS)^(alfa-1)*ySS;                         
lSS=(ySS/(ztildeSS*ksSS^alfa*nSS^(1/(theta-1))))^(1/(1-alfa));                            
ytildeSS=ztildeSS*lSS^(1-alfa)*ksSS^(alfa)/nSS;                                                       
kSS=ksSS/(nSS*cuSS)*ySS;                                                                                     
invSS=invk*kSS;    
cSS=ySS-(invSS+neSS*ecSS+nxSS*xcSS);
divhatSS=(theta/(theta-1)-1)*mctildeSS*(ztildeSS/zhatSS)^(1-theta)*ytildeSS;        
vhatSS=divhatSS/(1-beta*(1-etaSS)) ;                             
lvSS = vhatSS;                                                                                     
divSS=1-wSS*lSS-rkSS*ksSS ;                                                                                
divtildeSS=divSS/nSS;                                                                                 
vtildeSS = beta*((1-etaSS)*divSS+etaSS*lvSS)/(1-beta*(1-etaSS));                        
fe=vtildeSS-ecSS;                                                                                       
tau=1-(lvSS+xcSS)/fe; % as lv=(1-tau)*fe-xc
lab_shareSS = wSS*lSS/ySS;
pro_shareSS = divSS/ySS;

chi=((epsw-1))/(((((cSS-h*cSS)))*lSS^(phil))*(epsw/wSS));
mrsSS =chi*((cSS-h*cSS))*lSS^(phil);   
lambdaSS= (cSS-h*cSS)^(-sigma)*exp(((sigma-1)*chi*(lSS)^(1+phil))/(1+phil));
nuSS=qSS*lambdaSS;                                                                                         

taxss = -fe*neSS/ySS+(1-tau)*fe*nxSS/ySS;
taxSS=taxss*ySS;

lprodSS=(ySS/lSS);
 
cdataSS=cSS/nSS^(1/(1-theta));
ydataSS=ySS/nSS^(1/(1-theta));
invdataSS=invSS/nSS^(1/(1-theta));
divdataSS=divSS/nSS^(1/(1-theta));
divtildedataSS=divtildeSS/nSS^(1/(1-theta));
invnedataSS=(ecSS+fe)*neSS/nSS^(1/(1-theta));
invtotdataSS=(invnedataSS+invdataSS);
ldataSS=lSS/nSS^(1/(theta-1));

tfp_prdSS=ydataSS/(lSS^(1-alfa)*ksSS^(alfa));
tfpdataSS=log(tfp_prdSS/(tfp_prdSS));

gydataSS=0;
gcdataSS=0;
ginvdataSS=0;
ginvnedataSS=0;
ginvtotdataSS=0;
gldataSS=0;
gnxdataSS=0;
gnedataSS=0;

MSS=1;
VSS=1;
E_t_V_tp1_1_minus_alphaSS=VSS^(1-alphaEZ);
V_normalization=(1-beta)*VSS*(1-sigma)/(((cSS-h*cSS)^(1-sigma))*exp(((sigma-1)*chi*(lSS)^(1+phil))/(1+phil)));

ckappa=(1+iiSS)*((40-1)/40);
P_longSS=1/(1-ckappa*beta*MSS/paiSS);
P_freeSS=1/(1-ckappa/(1+iiSS));
R_longSS=log(ckappa*P_longSS/(P_longSS-1));//R_long=log(1+ii);
R_freeSS=log(ckappa*P_freeSS/(P_freeSS-1));//R_free=R_long;
T_premiumSS=(R_longSS-R_freeSS);
T_spreadSS=4*(R_longSS-log(1+iiSS));

R_equitySS=((1-etaSS)*(vtildeSS+divtildeSS)+etaSS*lvSS)/vtildeSS;
R_equity_varSS=(((1-etaSS)*(vtildeSS+divtildeSS)+etaSS*lvSS)/vtildeSS-R_equitySS)^2;

E_premiumSS=4*(log(R_equitySS)-log(1+rirSS));

vtilde_fakeSS=(beta*MSS*(1-etaSS))/(1-beta*MSS*(1-etaSS))*(divtildeSS+(etaSS/(1-etaSS))*(lvSS));

R_equity_fakeSS=(((vtildeSS+divtildeSS)/vtildeSS)+(etaSS/(1-etaSS)*lvSS/vtildeSS));
R_equity_fake_varSS=((((vtildeSS+divtildeSS)/vtildeSS)+(etaSS/(1-etaSS)*lvSS/vtildeSS))-R_equity_fakeSS)^2;
E_premium_fakeSS=4*(log(R_equity_fakeSS)-log(1+rirSS));

vtilde_consSS=beta*MSS*lambdaSS/lambdaSS*((cSS)^3)/(1-beta*MSS*lambdaSS/lambdaSS);
R_equity_consSS=(vtilde_consSS+(cSS)^3)/vtilde_consSS;
E_premium_consSS=4*(log(R_equity_consSS)-log(1+rirSS));

R_capitalSS=((rkSS*cuSS-adj_cuSS)+qSS*(1-deltak)-qSS*(((gammai/2)*(invSS/kSS-deltak)^2)-((gammai)*(invSS/kSS-deltak))*invSS/kSS))/qSS;
K_premiumSS=4*(log(R_capitalSS)-log(1+rirSS));

R_equityhatSS=((vhatSS)/vhatSS);

log_rirSS=log((1+rirSS)/(1+(rirSS)));
log_iiSS=log((1+iiSS)/(1+(iiSS)));

sdfSS=beta*MSS*lambdaSS/lambdaSS;

AlogRSS=(4*(log(1+iiSS)));
AlogpaiSS=(4*(log(paiSS)));
AlogRRSS=(4*(log(1+rirSS)));



zeps_piSS=0;
E_ySS=ydataSS;
var_ySS=0;

var_paiSS=0;
paieSS=paiSS;

  
//-------------------------------------------------------------------------
//                            MODEL EQUATIONS  
//-------------------------------------------------------------------------

model;



%expectation shock
paie=pai(+1)*exp(zeps_pi);



[name='E.1.']          %长期债券价格
P_long=1+ckappa*(P_long(+1))*beta*M/(paie);

[name='E.2.']          %长期债券收益率
R_long=log(ckappa*P_long/(P_long-1));

[name='E.3.']          %短期债券价格
P_free=1+ckappa*(P_free(+1)/(1+ii));

[name='E.4.']          %短期收益率
R_free=log(ckappa*P_free/(P_free-1));

[name='E.5.']          %期限溢价
T_premium=(R_long-R_free);
 
[name='E.6.']          %期限利差（年化）
T_spread=4*(R_long-log(1+ii));

[name='E.7.']          %平均企业的股权回报率
R_equity=((1-eta(+1))*(vtilde(+1)+divtilde(+1))+eta(+1)*lv(+1))/vtilde;

[name='E.8.']          %平均企业股权回报率的方差(VXO)
R_equity_var=(((1-eta(+1))*(vtilde(+1)+divtilde(+1))+eta(+1)*lv(+1))/vtilde-R_equity)^2;

[name='E.9.']          %股权风险溢价（年化）
E_premium=4*(log(R_equity)-log(1+rir));

[name='E.10.']         %假想企业价值（如果不退出的价值）
vtilde_fake=beta*M*lambda(+1)/lambda*(1-eta(+1))
*(vtilde_fake(+1)+divtilde(+1)
+(eta(+1)/(1-eta(+1)))*(lv(+1)));

[name='E.11.']        %假想企业股权回报率
R_equity_fake=((vtilde(+1)+divtilde(+1)+eta(+1)/(1-eta(+1))*lv(+1))/vtilde);

[name='E.12.']        %假想企业股权回报率的方差
R_equity_fake_var=((((vtilde(+1)+divtilde(+1))/vtilde)+(eta(+1)/(1-eta(+1))*lv(+1)/vtilde))-R_equity_fake)^2;

[name='E.13.']        %假想股权溢价
E_premium_fake=4*(log(R_equity_fake)-log(1+rir));

[name='E.14.']        
vtilde_cons=beta*M*lambda(+1)/lambda*(vtilde_cons(+1)+(c(+1))^3);

[name='E.15.']
R_equity_cons=(vtilde_cons(+1)+(c(+1))^3)/vtilde_cons;

[name='E.16.']
E_premium_cons=4*(log(R_equity_cons)-log(1+rir));

[name='E.17.']        %资本回报率
R_capital=(
            (rk(+1)*cu(+1)-adj_cu(+1))+
            q(+1)*(1-deltak)-
            q(+1)*(((gammai/2)*(inv(+1)/k-deltak)^2)-((gammai)*(inv(+1)/k-deltak))*inv(+1)/k)
            )/q;

[name='E.18.']       %资本风险溢价
K_premium=4*(log(R_capital)-log(1+rir));

[name='E.19.']       %随机贴现因子
sdf=beta*M*lambda(+1)/lambda;

[name='E.20.']       %临界企业的股权回报率
R_equityhat=((vhat(+1))/vhat);

[name='E.21.']      %下一期的终身效用风险调整
E_t_V_tp1_1_minus_alpha=V(+1)^(1-alphaEZ);

[name='E.22.']
M=(V(+1)^(-alphaEZ)) / (E_t_V_tp1_1_minus_alpha^(-alphaEZ/(1-alphaEZ)));

[name='E.23.']      %Epstein-Zin递归效用，V_normalization包含了(1-beta)
      V=V_normalization*(((c-h*c(-1))^(1-sigma))*exp(((sigma-1)*chi*(l)^(1+phil))/(1+phil)))/(1-sigma)
+beta*(E_t_V_tp1_1_minus_alpha)^(1/(1-alphaEZ));

[name='E.24.']      %（1）消费边际效用
lambda= (c-h*c(-1))^(-sigma)
*exp(((sigma-1)*chi*(l)^(1+phil))/(1+phil));

[name='E.25.']       %（2）边际替代率
mrs=chi*((c-h*c(-1)))*(l)^(phil);

[name='E.26.']       %（4）欧拉方程，但贴现因子有冲击zeps_b,有M【？】
lambda=beta*M*lambda(+1)*(1+rir)*exp(zeps_b);

[name='E.27.']        %（3）资本运动方程，但-1变成-deltak【？】，投资inv*exp(zeps_i)【难道是投资效率？】
k=(1-deltak)*k(-1)+inv*exp(zeps_i)
-((gammai/2)*(inv/k(-1)-deltak)^2)*k(-1);

[name='E.28.']        %（8）投资的一阶条件，但1*exp(zeps_i)【？】
1=q*(1*exp(zeps_i)-gammai*(inv/k(-1)-deltak));

[name='E.29.']        % (9)Tobin q
q=nu/lambda;  

[name='E.30.']        %（7）资本的欧拉方程同除了lambda
q=beta*M*((lambda(+1))/(lambda)*(rk(+1)*cu(+1)-adj_cu(+1))+
         ((lambda(+1))/(lambda))*q(+1)*(1-deltak)-
         ((lambda(+1))/(lambda))*q(+1)*
(((gammai/2)*(inv(+1)/k-deltak)^2)-((gammai)*(inv(+1)/k-deltak))*inv(+1)/k));

[name='E.31.']        %（10）资本利用率的一阶条件
rk=gammau1+gammau2*(cu-1);

[name='E.32.']        %（11）资本利用率调整成本
adj_cu=gammau1*(cu-1)+gammau2/2*(cu-1)^2;

[name='E.33.']        %（5）现存企业的欧拉方程
vtilde=beta*M*lambda(+1)/lambda
*((1-eta(+1))*(vtilde(+1)+divtilde(+1))+eta(+1)*lv(+1));

[name='E.34.']         %（13）工资NKPC，同乘(epsw-1)
mrs*(epsw/w)+(1-epsw)-gammaw/(l*w(-1))*(pai_w*pai-STEADY_STATE(pai))*pai*y+
beta*M*lambda(+1)/lambda*   %期望算子又出现了M
(gammaw/(l*w)*(pai_w(+1)*paie-STEADY_STATE(pai))*pai_w(+1)*paie*y(+1))+zeps_w=0;

[name='E.35.']         %（14）工资通胀率，少输了pi，后续使用时要*pi
pai_w=w/w(-1);

[name='E.36.']         %（12）企业数量运动方程，注意eta(+1)【？】
n=(1-eta(+1))*(n(-1)+ne);

[name='E.37.']         %（21）退出企业数量，注意eta和n(-2)和ne(-1)【？】
nx=eta*(n(-2)+ne(-1));
 
[name='E.38.']         %（6）进入企业的欧拉方程
vtilde=(ec+fe)*exp(zeps_e);

[name='E.39.']         %（19）进入的拥堵外部性，ne和n的时间下标有问题【？】
ec=big_theta*(ne/n(-1))^varsigma;
 
[name='E.40.']         %（22）清算价值，多了个exp(zeps_z)【？】
lv=((1-tau)*fe-xc)*exp(zeps_x);

[name='E.41']          %（20）退出拥堵外部性，时间下标【？】
xc=big_thetaxc*(nx/n(-1))^varsigmaxc;

[name='E.42.']         %（23）退出概率
eta=1-(zmin/zhat)^xie;

[name='E.43.']         %（24）平均生产率
ztilde=(xie/(xie+1-(theta)))^(1/((theta)-1))*zhat;

[name='E.44.']         %（25）临界企业的价值
vhat=divhat+
beta*M*lambda(+1)/lambda*(1-eta(+1))*(vhat(+1));

[name='E.45.']         %（26）退出的临界条件
vhat=lv;

[name='E.46.']         %（27）临界企业的利润
divhat=(theta/(theta-1)-1)*mctilde*(ztilde/zhat)^(1-(theta))*ytilde;

[name='E.47.']         %（28）平均企业的最优价格
rhotilde=muptilde*mctilde;
    
[name='E.48,']         %（29）中间品企业的加成率【定义式】
muptilde=(theta)/(((theta)-1));

[name='E.49.']         %（16）价格NKPC，同乘了theta-1。其中STEADY_STATE(pai)对应原文1。但带入随机贴现率少了个(1-eta),还多了个zps_P
(1-(theta))+(theta)*rho_I -(1-(theta))*(gammap/2)*(pai-STEADY_STATE(pai))^2
-gammap*(pai-STEADY_STATE(pai))*pai+beta*M*(lambda(+1)/lambda)*gammap*(paie-STEADY_STATE(pai))*paie*((y(+1))/(y))+zeps_p=0;

[name='E.50.']         %（17）中间品厂商的价格聚合
rhotilde=n(-1)^(1/((theta)-1))*rho_I;

[name='E.51.']         %（31）劳动需求，多了个A，可能是全要素生产率【？】
w=A*ztilde*mctilde*(1-alfa)*(l/ks)^(-alfa);

[name='E.52.']         %（32）资本服务需求，多个A
rk=A*ztilde*mctilde*alfa*(l/ks)^(1-alfa);

[name='E.53.']         %（34）资本服务
ks=cu*k(-1);

[name='E.54.']         %（33）平均企业的产出，多个A，还要除以n(-1)【？】
ytilde=A*ztilde*l^(1-alfa)*ks^alfa/n(-1);

[name='E.55.']         %【】平均每家企业的产出
ymean=y/n(-1);
 
[name='E.56.']         %【】pai_is定义式，类比（14）【？】
pai_is=rhotilde/rhotilde(-1)*pai;
 
[name='E.57.']         %（36）总量资源约束
y=c+inv+adj_cu*k(-1)
+ne*ec*exp(zeps_e)     %*exp(zps_e)
+nx*xc*exp(zeps_x)     %*exp(zps_x)外生冲击
+(gammap/2)*(pai-STEADY_STATE(pai))^2*y   %PAC
+(gammaw/2)*(w*pai/w(-1)-STEADY_STATE(pai))^2*y  %WAC
+zeps_d;               %【？】

[name='E.58.']         %（37）总产出
y=A*n(-1)^(1/((theta)-1))*ztilde*l^(1-alfa)*ks^(alfa);

[name='E.59.']         %（30）现存企业的总利润
div=y-w*l-rk*ks;

[name='E.60.']         %（30）平均企业的利润
divtilde=div/n(-1);

[name='E.61.']        %（38）政府预算约束，同除y，此处tax，原文是转移支付
0=
-(tax/y)
-(fe*ne*exp(zeps_e))/y  %*exp(zeps_e)  
+((1-tau)*fe*nx*exp(zeps_x))/y;     %*exp(zeps_x)

[name='E.62.']          %（40）费雪方程
(1+rir)=(1+ii)/paie;

[name='E.63.']          %（39）Taylor规则（参考Asacri2023原文，与附录不完全相同）
log((1+ii)/(1+STEADY_STATE(ii)))    %i用TEADY_STATE(ii)代替
=(phii)*log(((1+ii(-1))/(1+STEADY_STATE(ii))))+
//(1+ii)/(1+ii(-1))=
(phipai*(1-phii))*log((pai/STEADY_STATE(pai)))+
(phiy*(1-phii))*log((ydata/(STEADY_STATE(ydata))))+    %这个ydata和y不一样
(phidy*(1-phii))*log((ydata/(ydata(-1))))+
(zeps_ii);              %冲击

[name='E.64.'] 
tfp_prd=A*ztilde;

[name='E.65.']         
lprod=(ydata/l);

[name='E.66.']          %劳动收入份额
lab_share = w*l/y;

[name='E.67.']          %利润收入份额
pro_share = div/y;

[name='E.68.']          %全要素生产率
A = exp(zeps_a);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%                          SHOCKS                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[name='E.69.']
zeps_ii=rhoii*zeps_ii(-1)+eps_ii;

[name='E.71.']
zeps_b=rhob*zeps_b(-1)+eps_b;

[name='E.72.']
zeps_i=rhoi*zeps_i(-1)+eps_i;

[name='E.73.']
zeps_w=rhow*zeps_w(-1)+eps_w;

[name='E.74.']
zeps_e=rhoe*zeps_e(-1)+eps_e;

[name='E.75.']
zeps_x=rhox*zeps_x(-1)+eps_x;

[name='E.76.']
zeps_a=rhoa*zeps_a(-1)+sigma_a*eps_a;

[name='E.77.']
zeps_p=rhop*zeps_p(-1)+eps_p; 

[name='E.78.']
zeps_d=rhod*zeps_d(-1)+eps_d;

%expectation shock
zeps_pi=rhopi*zeps_pi(-1)+sigma_pi*eps_pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%                        var w/o LOVE                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[name='E.79.']
ydata=y/n(-1)^(1/(theta-1));

[name='E.80']
cdata=c/n(-1)^(1/(theta-1));

[name='E.81.'] 
invdata=inv/n(-1)^(1/(theta-1));

[name='E.82.'] 
invnedata=((ec+fe)*exp(zeps_e))*ne/n(-1)^(1/(theta-1));

[name='E.83.'] 
invtotdata=(invnedata+invdata);

[name='E.84.'] 
ldata=l/n(-1)^(1/(theta-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%             Growth rates              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[name='E.85.'] 
gydata= log(ydata/ydata(-1));

[name='E.86.'] 
gcdata= log(cdata/cdata(-1));

[name='E.87.'] 
ginvdata= log(invdata/invdata(-1));

[name='E.88.'] 
ginvnedata= log(invnedata/invnedata(-1));

[name='E.89.'] 
ginvtotdata=log(invtotdata/invtotdata(-1));

[name='E.90.'] 
gldata = log(l/l(-1));

[name='E.91.']      
gnedata = log(ne/ne(-1));

[name='E.92.']         
gnxdata = log(nx/nx(-1));

[name='E.93.'] 
tfpdata= log(tfp_prd/tfp_prd(-1));

%[name='E.94.'] 
%log(xie/STEADY_STATE(xie))=log((xie(-1)/STEADY_STATE(xie))^rhoxie)+sigma_xie*eps_xie;

%[name='E.95.']
%log(sigma_xie)=(rhosigma_xie)*log(sigma_xie(-1))+eta_xie;

[name='E.96.']
AlogR=(4*(log(1+ii)));

[name='E.97.']
Alogpai=(4*(log(pai)));

[name='E.98.']
AlogRR=(4*(log(1+rir)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%                          Uncertainty                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E_y=ydata(+1);
var_y=(ydata(+1)-E_y)^2;

var_pai=(pai(+1)-paie)^2;


%VXO即为R_equity_var



end;

//-------------------------------------------------------------------------
//                            INITIAL VALUES  
//-------------------------------------------------------------------------

initval;
rk=rkSS;
rir=rirSS;
ii=iiSS;
q=qSS;
cu=cuSS;
adj_cu=adj_cuSS;
pai=paiSS;
pai_w=pai_wSS;
pai_is=pai_isSS;
vtilde=vtildeSS;
vhat=vhatSS;
eta=etaSS;
lv=lvSS;
divtilde=divtildeSS;
divhat=divhatSS;
div=divSS;
ne=neSS;
n=nSS;
nx=nxSS;
ec=ecSS;
xc=xcSS;
zhat=zhatSS;
ztilde=ztildeSS;
mctilde=mctildeSS;
muptilde=muptildeSS;
rhotilde=rhotildeSS;
rho_I=rho_ISS;
tax=taxSS;
ldata=ldataSS;
invnedata=invnedataSS;
tfp_prd=tfp_prdSS;
lab_share=lab_shareSS;
pro_share=pro_shareSS;
P_long=P_longSS;
R_long=R_longSS;
P_free=P_freeSS;
R_free=R_freeSS;
T_premium=T_premiumSS;
T_spread=T_spreadSS;
R_equity=R_equitySS;
E_premium=E_premiumSS;
vtilde_fake=vtilde_fakeSS; 
R_equity_fake=R_equity_fakeSS;
E_premium_fake=E_premium_fakeSS;
vtilde_cons=vtilde_consSS;
R_equity_cons=R_equity_consSS;
E_premium_cons=E_premium_consSS;
zeps_b=zeps_bSS;

zeps_i=zeps_iSS;
zeps_w=zeps_wSS;    
zeps_e=zeps_eSS;
zeps_x=zeps_xSS;
zeps_a=zeps_aSS;
zeps_ii=zeps_iiSS;

zeps_p=zeps_pSS;
zeps_d=zeps_dSS;
gydata=gydataSS;
gcdata=gcdataSS;
ginvdata=ginvdataSS;
ginvnedata=ginvnedataSS;
ginvtotdata=ginvtotdataSS;
gldata=gldataSS;
gnxdata=gnxdataSS;
gnedata=gnedataSS;
y=ySS;
ymean=ymeanSS;
c=cSS;
inv=invSS;
k=kSS;
ks=ksSS;
ytilde=ytildeSS;
ydata=ydataSS;
cdata=cdataSS;
invdata=invdataSS;
tfpdata=tfpdataSS;
invtotdata=invtotdataSS;
l=lSS;
lambda=lambdaSS;
nu=nuSS;
V=VSS;
E_t_V_tp1_1_minus_alpha=E_t_V_tp1_1_minus_alphaSS;
M=MSS;
mrs=mrsSS;
w=wSS;
lprod=lprodSS;
A=ASS;
%xie=xiess;
%sigma_xie=sigma_xieSS;
R_equity_var=R_equity_varSS;
R_equity_fake_var=R_equity_fake_varSS;
sdf=sdfSS;
R_capital=R_capitalSS; 
K_premium=K_premiumSS;
R_equityhat=R_equityhatSS; 
AlogR=AlogRSS; 
Alogpai=AlogpaiSS; 
AlogRR=AlogRRSS;

zeps_pi=zeps_piSS;
E_y=E_ySS;
var_y=var_ySS;

paie=paieSS;
var_pai=var_paiSS;



end;

//-------------------------------------------------------------------------
//                            SHOCKS STDERR  
//-------------------------------------------------------------------------

shocks;
var eps_pi=1;

var   eps_ii=0;

%var   eps_xie=0;
%var   eta_xie=0;

var   eps_b=0;
var   eps_i=0;
var   eps_w=0;
var   eps_e=0;
var   eps_x=0;
var   eps_a=0;
var   eps_p=0;
var   eps_d=0;

end;

//-------------------------------------------------------------------------
//                          CHECK AND SIMULATION  
//-------------------------------------------------------------------------

resid(1);
check;
steady;

stoch_simul(order=3,periods=0,irf=0,noprint,nograph,nomoments,nofunctions,nocorr,pruning);


IRF_periods=60;
burnin=5000; %periods for convergence

reorderidx(1)=(strmatch('ydata',M_.endo_names,'exact'));
reorderidx(2)=(strmatch('cdata',M_.endo_names,'exact'));
reorderidx(3)=(strmatch('invdata',M_.endo_names,'exact'));
reorderidx(4)=(strmatch('vtilde',M_.endo_names,'exact'));
reorderidx(5)=(strmatch('ne',M_.endo_names,'exact'));
reorderidx(6)=(strmatch('nx',M_.endo_names,'exact'));
reorderidx(7)=(strmatch('ldata',M_.endo_names,'exact'));
reorderidx(8)=(strmatch('rk',M_.endo_names,'exact'));
reorderidx(9)=(strmatch('rir',M_.endo_names,'exact'));
reorderidx(10)=(strmatch('pai',M_.endo_names,'exact'));
reorderidx(11)=(strmatch('pai_w',M_.endo_names,'exact'));
reorderidx(12)=(strmatch('ii',M_.endo_names,'exact'));
reorderidx(13)=(strmatch('tfp_prd',M_.endo_names,'exact'));

reorderidx(14)=(strmatch('T_spread',M_.endo_names,'exact'));
reorderidx(15)=(strmatch('E_premium',M_.endo_names,'exact'));
reorderidx(16)=(strmatch('E_premium_fake',M_.endo_names,'exact'));
reorderidx(17)=(strmatch('E_premium_cons',M_.endo_names,'exact'));
reorderidx(18)=(strmatch('K_premium',M_.endo_names,'exact'));

%reorderidx(19)=(strmatch('sigma_eps_ii',M_.endo_names,'exact'));
reorderidx(19)=(strmatch('paie',M_.endo_names,'exact'));

reorderidx(20)=(strmatch('R_equity',M_.endo_names,'exact'));
reorderidx(21)=(strmatch('R_equity_fake',M_.endo_names,'exact'));
reorderidx(22)=(strmatch('R_equity_var',M_.endo_names,'exact'));
reorderidx(23)=(strmatch('R_equity_fake_var',M_.endo_names,'exact'));
reorderidx(24)=(strmatch('AlogR',M_.endo_names,'exact'));
reorderidx(25)=(strmatch('Alogpai',M_.endo_names,'exact'));
reorderidx(26)=(strmatch('AlogRR',M_.endo_names,'exact'));

//-------------------------------------------------------------------------
//                          VOLATILITY SHOCK  
//-------------------------------------------------------------------------

shock_mat_with_zeros=zeros(burnin+IRF_periods,M_.exo_nbr); %shocks set to 0 to simulate without uncertainty
IRF_no_shock_mat = simult_(oo_.dr.ys,oo_.dr,shock_mat_with_zeros,options_.order)'; %simulate series
stochastic_steady_state=IRF_no_shock_mat(1+burnin,:); % stochastic_steady_state/EMAS is any of the final points after burnin

shock_mat = zeros(burnin+IRF_periods,M_.exo_nbr);
shock_mat(1+burnin,strmatch('eps_pi',M_.exo_names,'exact'))= 1; %1
IRF_mat = simult_(oo_.dr.ys,oo_.dr,shock_mat,options_.order)';

IRF_mat_percent_from_SSS_eta_ii= (IRF_mat(1+burnin+1:1+burnin+IRF_periods,:)-IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,:))./repmat(stochastic_steady_state,IRF_periods,1); %only valid for variables not yet logged

tmp1        		= sqrt(4*(  max(IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(22)),1e-16) ) );
tmp2                = sqrt(4*(stochastic_steady_state(:,reorderidx(22))));
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(22))	            = log(tmp1./tmp2);

tmp1        		= sqrt(4*(  max(IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(23)),1e-16) ) );
tmp2                = sqrt(4*(stochastic_steady_state(:,reorderidx(23))));
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(23))	            = log(tmp1./tmp2);

%for logged variables
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(14)) = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(14)) -IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(14)));
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(15)) = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(15)) -IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(15)));
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(16)) = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(16)) -IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(16)));
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(17)) = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(17)) -IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(17)));
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(18)) = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(18)) -IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(18)));
%IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(19)) = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(19)) -IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(19)));
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(24)) = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(24)) -IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(24)));
IRF_mat_percent_from_SSS_eta_ii(:,reorderidx(25)) = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(25)) -IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,reorderidx(25)));