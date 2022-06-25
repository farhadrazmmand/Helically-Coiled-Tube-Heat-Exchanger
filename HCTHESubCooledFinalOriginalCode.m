clc
clear all
format short g
%% Start Code
%Input Parameters
%Gas Inlet Parameters (for now "Helium")
%Geometrical and Thermophysical Parameters about Tube and Helical based on Xiao's Research (Original1)
D_t_ex=0.0185;%Outer Diameter of Tube (m) (Variable)
r_t_ex=D_t_ex/2;
D_t_in=0.0135;%Inner Diameter of Tube (m) (Variable)
r_t_in=D_t_in/2;
D_h=0.203;%Coil Diameter (m)
H=0.135;%Helical Pitch (m) (Variable)
pi=3.14;
teta=atan((H/2)/(2*pi*((r_t_ex+r_t_in)/2)));%Helical Angle
k_t=207;%Tube Thermal Conductivity
L_helical=(((H.^2)+((pi*D_h).^2)).^0.5);%Helical Length for every loop
Ls=L_helical/4;
mdat=0.1163;%Mass Flow Rate for Fluid
mdat_g=0.14333;%Mass Flow Rate for Gas
D_shell=D_h+((D_t_ex-D_t_in)/2);
A_g=pi*((D_shell.^2)-(D_h.^2));
deltaPsat=1;
%Nano Particle Parameter
fi=input('Volume fraction of nano particle:');%Volume Fraction of Nano Particle
ro_p=19300;%Density of Nano Particle
ro_f_ref=999;%Fluid Density in the 20 degree celcius
d_p=67.5e-9;%Diameter of Nano Particle
k_p=310;%Thermal Conductivity of Nano Particle
C_p_p=129;%Specific Heat for Nano Paricle
n=3;%empirical shape factor of Nano Particle
Av=6.022*(10e23);%Avogadro
M=0.0180153;%Mole Mass of base Fluid
Anp=-4.34023;
Bnp=-105.60799;
Cnp=527.45631;
Dnp=-151.74505;
Enp=-903.41949;
Fnp=2814.30560;
alphanp=-0.00355;
bettanp=0.03564;
gammanp=-0.10898;
sigmanp=0.08845;
epsilonnp=0.02737;
omeganp=-0.04025;
%General Parameters
sigma=0.000000056;%Stefan–Boltzmann constant (sigma)=5.67*10^-8
T_g_in1=769.64;
%% Start Calculations
for N=1:20
    %Start of Section 1 (0-Ls)
    disp ('Number of Loop:');
    disp (N);
    e111=5;
    e11=5;%input('Number of radial element for cross section in Single-Phase Mode of First Section from Loop No.1:');
    e1=5;%input('Number of element in Single-Phase Mode of First Section from Loop No.1:');
    Length1=Ls/e1;
    l_shell=Length1*cos(teta);
    T_t_ex1=zeros(e1,1);
    T_t_in1=zeros(e1,1);
    T_out1=zeros(e111,e11,e1);
    T_g_out=zeros(e1,1);
    if N==1
        T_in1=367.49;%input('Value of Fluid Temperature at Start Process:');
    else
        T_in1=T3;
    end
    for i=1:e1
        i=i;
        r_t_in=D_t_in/2;
        r_0=r_t_in/e11;
        r_in=r_0;
        for j=1:e11
            for k=1:1:e111
                if k==1
                    Fi=0;
                elseif k==2
                    Fi=90;
                elseif k==3
                    Fi=180;
                elseif k==4
                    Fi=270;
                elseif k==5
                    Fi=360;
                end
                if N==2 || N>2
                    if i==1
                        T_in1=T3(1,1,:);
                    end
                end
                mu_g1=(4*(10.^(-5)))*exp(4170/T_g_in1);%Gas Dynamic Viscosity
                C_p_g1=976.78+(1.0637*T_g_in1);%Gas Specific Heat
                k_g1=0.43482+(5*(10.^(-4))*T_g_in1);%Gas Thermal Conductivity
                ro_g1=2729.3-0.73*T_g_in1;%Gas Density
                v_g1=mdat_g/(A_g*ro_g1);%Gas Velocity
                Pr_g1=(C_p_g1*mu_g1)/k_g1;%Prandtl Number for Gas
                Re_g1=(ro_g1*v_g1*(D_t_ex/cos(teta)))/mu_g1;%Reynolds Number for Gas
                Nu_D_g1=0.3+(((0.62*(Re_g1.^0.5)*(Pr_g1.^(1/3)))/((1+((0.4/Pr_g1).^(2/3))).^0.25))*((1+((Re_g1/282000).^(5/8))).^0.8));%Nusselt Number for Gas
                h_g1=Nu_D_g1*k_g1*cos(teta)/D_t_ex;%Heat Transfer Coefficient for Gas
                mu_f1=((2.1897e-11)*(T_in1.^4))-((3.055e-8)*(T_in1.^3))+((1.6028e-5)*(T_in1.^2))-(0.0037524*T_in1)+0.33158;%Fluid Dynamic Viscosity
                C_p_f1=((1.1105e-5)*(T_in1.^3))-(0.0031078*(T_in1.^2))-(1.478*T_in1)+4631.9;%Fluid Specific Heat
                k_f1=((1.5362e-8)*(T_in1.^3))-((2.261e-05)*(T_in1.^2))+(0.010879*T_in1)-1.0294;%Fluid Thermal Conductivity
                ro_f1=((-1.5629e-5)*(T_in1.^3))+(0.011778*(T_in1.^2))-(3.0726*T_in1)+1227.8;%Fluid Density
                ro_v=0.260;%Fluid Density
                d_f1=0.1*((6*M)/(Av*pi*ro_f_ref));
                mu_nf1=mu_f1/(1-(34.87*((d_p/d_f1).^(-0.3))*(fi.^1.03)));%Dynamic Viscosity Equation for Nano Fluid
                k_nf1=k_f1*((k_p+((n-1)*k_f1)-((n-1)*fi*(k_f1-k_p)))/(k_p+((n-1)*k_f1)+(fi*(k_f1-k_p))));%Thermal Conductivity Equation for Nano Fluid
                A11=h_g1;% A(number coefficient)(number section)
                r=(r_t_ex+r_t_in)/2;
                A21=k_t/(r*log(r_t_ex/r_t_in));
                g=9.81;
                ro_nf1=(ro_p*fi)+(ro_f1*(1-fi));%Dendity Equation for Nano Fluid
                P_g1=ro_nf1*g*Length1*cos(teta);%Fluid Pressure
                A_il=pi*((2*r_in).^2);
                u_m1=mdat/(A_il*ro_nf1);%Fluid Velocity
                C_p_nf1=((fi*ro_p*C_p_p)+((1-fi)*ro_f1*C_p_f1))/ro_nf1;%Specific Heat Equation for Nano Fluid
                Pr_nf1=(C_p_nf1*mu_nf1)/k_nf1;%Prandtl Number Equation for Nano Fluid
                Re_nf1=(ro_nf1*u_m1*(2*r_in))/(mu_nf1);%Reynolds Number Equation for Nano Fluid
                if Re_nf1<2300
                    Nu_nf1=4.36;%Nusselt Number
                else
                    if P_g1<0
                        P1=-P_g1;
                    else
                        P1=P_g1;
                    end
                    f1=(2*P1*D_t_in)/(Length1*ro_nf1*(u_m1.^2));%Friction Factor Equation for Nano Fluid
                    Nu_nf1=((f1/8)*(Re_nf1-1000)*Pr_nf1)/(1+(12.7*((f1/8).^0.5))*((Pr_nf1.^(2/3))-1));%Nusselt Number Equation for Nano Fluid
                end
                Nu_loc1=Nu_nf1*(((-2.331e-05)*(Fi.^2))+((8.424e-03)*Fi)+0.4576);%Local Nusselt Number
                if fi==0
                    hlv_nf1=2317000;%Enthalpy
                else
                    C01=(Anp*(fi.^5))+(Bnp*(fi.^4))+(Cnp*(fi.^3))+(Dnp*(fi.^2))+(Enp*fi)+Fnp;
                    C11=(alphanp*(fi.^5))+(bettanp*(fi.^4))+(gammanp*(fi.^3))+(sigmanp*(fi.^2))+(epsilonnp*fi)+omeganp;
                    hlv_nf1=C01*(P1.^C11);%Enthalpy
                end
                h_nf1=Nu_loc1*k_nf1/(2*r_in);
                G1=0.00122*(((k_nf1.^0.79)*(C_p_nf1.^0.45)*(ro_f1.^0.49))/((sigma.^0.5)*(mu_nf1.^0.29)*(hlv_nf1.^0.24)*(ro_v.^0.24)))*(P1.^0.75);
                S_chen1=1/(1+((2.53*(10.^-6))*(Re_nf1.^1.17)));
                if i==1
                    h_nb1=G1*((373.15-T_in1).^0.24);
                else
                    h_nb1=G1*((T_t_in1-T_out1(1,1,i)).^0.24);
                end
                A31=h_nf1;
                A311=h_nb1*S_chen1*pi*D_t_in;
                A41=(mdat_g*C_p_g1)/(l_shell*pi*(D_shell-D_h));
                A51=(mdat*C_p_nf1)/(Length1*pi*(2*r_in));
                b11=((A11)/2)*T_g_in1;
                b21=0;
                b31=-373.15*A311;
                b41=A41*T_g_in1;
                b51=(-A51)*T_in1;
                %Unknowns: q, T_g_out, T_t_ex, T_t_in, T_out
                Y1=[1,(A11)/2,-A11,0,0%Eq1
                    1,0,-A21,A21,0%Eq2
                    1,0,0,-(A31+A311),A31%Eq3
                    1,A41,0,0,0%Eq4
                    1,0,0,0,-A51];%Eq5
                B1=[b11;b21;b31;b41;b51];
                X1=Y1\B1;
                T_out1(k,j,i)=X1(5,:);
                TEMP11(k,j,i)=X1(5,:);
                TF1=TEMP11;
                Nusselt11(k,j,i)=Nu_loc1;
                Nu1=Nusselt11;
                Friction_Factor(j,i)=f1;
                FF1=Friction_Factor;
                Velocity_r(j,i)=u_m1;
                V_r1=Velocity_r;
                Reynolds11(j,i)=Re_nf1;
                Re1=Reynolds11;
                if Re1(j,:)<0
                    break
                end
                HT1(k,j,i)=h_nf1;
                HT11=HT1;
            end
            r_in=r_in+r_0;
        end
        T_in1=T_out1(1,1,i);
        T_t_in1=X1(4,:);
        TEMP21(i)=T_t_in1;
        TIS1=TEMP21;
        T_g_out1(i)=X1(2,:);
        T_g_in1=T_g_out1(i);
        TEMP31(i)=X1(2,:);
        TG1=TEMP31;
        if TF1(k,j,i)>373.15
            break
        end
    end
    disp 'Distribution of Reynolds Number at Single-Phase Mode of First Section:';
    disp (Re1);
    Re1=Re1(:,i)
    disp 'Distribution of Internal Surface Temperature at Single-Phase Mode of First Section:';
    disp (TIS1);
    TIS1=TIS1(i)
    disp 'Distribution of Nusselt Number at Single-Phase Mode of First Section:';
    disp (Nu1);
    Nu1=Nu1(:,:,i)
    disp 'Distribution of Friction Factor at Single-Phase Mode of First Section:';
    disp (FF1);
    disp 'Distribution of Radius Velocity at Single-Phase Mode of First Section:';
    disp (V_r1);
    disp 'Distribution of Fluid Temperature at Single-Phase Mode of First Section:';
    disp (TF1);
    T1=TF1(:,:,i)
    disp 'Distribution of Gas Temperature at Single-Phase Mode of First Section:';
    disp (TG1);
    T_g1=TEMP31(i);
    disp 'Distribution of Heat Transfer Coefficient at Single-Phase Mode of First Section:';
    disp (HT11);
    if TF1(k,j,i)>373.15
        break
    end
    %End of Section 1
    %% Start of Section 2 (Ls - 3*Ls) = Ls---->(2*Ls)
    e222=5;
    e22=5;%input('Number of radial element in Single-Phase Mode of Second Section from Loop No.1:');
    e2=5;%input('Number of element in Single-Phase Mode of Second Section from Loop No.1:');
    Length2=(2*Ls)/e2;
    l_shell=Length2*cos(teta);
    T_t_ex2=zeros(e2,1);
    T_t_in2=zeros(e2,1);
    T_out2=zeros(e222,e22,e2);
    T_g_out2=zeros(e2,1);
    T_g_in2=T_g1;
    for i=1:e2
        i=i;
        r_t_in=D_t_in/2;
        r_0=r_t_in/e22;
        r_in=r_0;
        for j=1:e22
            for k=1:1:e111
                if k==1
                    Fi=0;
                elseif k==2
                    Fi=90;
                elseif k==3
                    Fi=180;
                elseif k==4
                    Fi=270;
                elseif k==5
                    Fi=360;
                end
                if i==1
                    T_in2=T1(1,1,:);
                end
                mu_g2=(4*(10.^(-5)))*exp(4170/T_g_in2);%Gas Dynamic Viscosity
                C_p_g2=976.78+(1.0637*T_g_in2);%Gas Specific Heat
                k_g2=0.43482+(5*(10.^(-4))*T_g_in2);%Gas Thermal Conductivity
                ro_g2=2729.3-0.73*T_g_in2;%Gas Density
                v_g2=mdat_g/(A_g*ro_g2);
                Pr_g2=(C_p_g2*mu_g2)/k_g2;
                Re_g2=(ro_g2*v_g2*(D_t_ex/cos(teta)))/mu_g2;
                Nu_D_g2=0.3+(((0.62*(Re_g2.^0.5)*(Pr_g2.^(1/3)))/((1+((0.4/Pr_g2).^(2/3))).^0.25))*((1+((Re_g2/282000).^(5/8))).^0.8));
                h_g2=Nu_D_g2*k_g2*cos(teta)/D_t_ex;
                mu_f2=((2.1897e-11)*(T_in2.^4))-((3.055e-8)*(T_in2.^3))+((1.6028e-5)*(T_in2.^2))-(0.0037524*T_in2)+0.33158;
                C_p_f2=((1.1105e-5)*(T_in2.^3))-(0.0031078*(T_in2.^2))-(1.478*T_in2)+4631.9;
                k_f2=((1.5362e-8)*(T_in2.^3))-((2.261e-05)*(T_in2.^2))+(0.010879*T_in2)-1.0294;
                ro_f2=((-1.5629e-5)*(T_in2.^3))+(0.011778*(T_in2.^2))-(3.0726*T_in2)+1227.8;
                d_f2=0.1*((6*M)/(Av*pi*ro_f_ref));
                mu_nf2=mu_f2/(1-(34.87*((d_p/d_f2).^(-0.3))*(fi.^1.03)));
                k_nf2=k_f2*((k_p+((n-1)*k_f2)-((n-1)*fi*(k_f2-k_p)))/(k_p+((n-1)*k_f2)+(fi*(k_f2-k_p))));
                A12=h_g2;% A(number coefficient)(number section)
                r=(r_t_ex+r_t_in)/2;
                A22=k_t/(r*log(r_t_ex/r_t_in));
                g=-9.81;
                ro_nf2=(ro_p*fi)+(ro_f2*(1-fi));
                P_g2=ro_nf2*g*Length2*cos(teta);
                A_il=pi*((2*r_in).^2);
                u_m2=mdat/(A_il*ro_nf2);
                C_p_nf2=((fi*ro_p*C_p_p)+((1-fi)*ro_f2*C_p_f2))/ro_nf2;
                Pr_nf2=(C_p_nf2*mu_nf2)/k_nf2;
                Re_nf2=(ro_nf2*u_m2*D_t_in)/(mu_nf2);
                if Re_nf2<2300
                    Nu_nf2=4.36;
                else
                    if P_g2<0
                        P2=-P_g2;
                    else
                        P2=P_g2;
                    end
                    f2=(2*P2*D_t_in)/(Length2*ro_nf2*(u_m2.^2));
                    Nu_nf2=((f2/8)*(Re_nf2-1000)*Pr_nf2)/(1+(12.7*((f2/8).^0.5))*((Pr_nf2.^(2/3))-1));
                end
                Nu_loc2=Nu_nf2*(((-2.331e-05)*(Fi.^2))+((8.424e-03)*Fi)+0.4576);
                if fi==0
                    hlv_nf2=2317000;
                else
                    C02=(Anp*(fi.^5))+(Bnp*(fi.^4))+(Cnp*(fi.^3))+(Dnp*(fi.^2))+(Enp*fi)+Fnp;
                    C12=(alphanp*(fi.^5))+(bettanp*(fi.^4))+(gammanp*(fi.^3))+(sigmanp*(fi.^2))+(epsilonnp*fi)+omeganp;
                    hlv_nf2=C02*(P2.^C12);
                end
                h_nf2=Nu_loc2*k_nf2/(2*r_in);
                G2=0.00122*(((k_nf2.^0.79)*(C_p_nf2.^0.45)*(ro_f2.^0.49))/((sigma.^0.5)*(mu_nf2.^0.29)*(hlv_nf2.^0.24)*(ro_v.^0.24)))*(P2.^0.75);
                S_chen2=1/(1+((2.53*(10.^-6))*(Re_nf2.^1.17)));
                if i==1
                    h_nb2=G2*((373.15-T_in2).^0.24);
                else
                    h_nb2=G2*((T_t_in2-T_out2(1,1,i)).^0.24);
                end
                A32=h_nf2;
                A321=h_nb2*S_chen2*pi*D_t_in;
                A42=(mdat_g*C_p_g2)/(l_shell*pi*(D_shell-D_h));
                A52=(mdat*C_p_nf2)/(Length2*pi*(2*r_in));
                b12=((A12)/2)*T_g_in2;
                b22=0;
                b32=-373.15*A321;
                b42=A42*T_g_in2;
                b52=(-A52)*T_in2;
                %Unknowns: q, T_g_out, T_t_ex, T_t_in, T_out
                Y2=[1,(-A12)/2,A12,0,0%Eq1
                    1,0,-A22,A22,0%Eq2
                    1,0,0,-(A32+A321),A32%Eq3
                    1,A42,0,0,0%Eq4
                    1,0,0,0,-A52];%Eq5
                B2=[b12;b22;b32;b42;b52];
                X2=Y2\B2;
                T_out2(k,j,i)=X2(5,:);
                TEMP12(k,j,i)=X2(5,:);
                TF2=TEMP12;
                Nusselt22(k,j,i)=Nu_loc2;
                Nu2=Nusselt22;
                Friction_Factor(j,i)=f2;
                FF2=Friction_Factor;
                Velocity_r(j,i)=u_m2;
                V_r2=Velocity_r;
                Reynolds22(j,i)=Re_nf2;
                Re2=Reynolds22;
                if Re2(j,:)<0
                    break
                end
                HT2(k,j,i)=h_nf2;
                HT22=HT2;
            end
            r_in=r_in+r_0;
        end
        T_in2=T_out2(1,1,i);
        T_t_in2=X2(4,:);
        TEMP22(i)=T_t_in2;
        TIS2=TEMP22;
        T_g_out2(i)=X2(2,:);
        T_g_in2=T_g_out2(i);
        TEMP32(i)=X2(2,:);
        TG2=TEMP32;
        if TF2(k,j,i)>373.15
            break
        end
    end
    disp 'Distribution of Reynolds Number at Single-Phase Mode of Second Section:';
    disp (Re2);
    Re2=Re2(:,i)
    disp 'Distribution of Internal Surface Temperature at Single-Phase Mode of Second Section:';
    disp (TIS2);
    TIS2=TIS2(i)
    disp 'Distribution of Nusselt Number at Single-Phase Mode of Second Section:';
    disp (Nu2);
    Nu2=Nu2(:,:,i)
    disp 'Distribution of Friction Factor at Single-Phase Mode of Second Section:';
    disp (FF2);
    disp 'Distribution of Radius Velocity at Single-Phase Mode of Second Section:';
    disp (V_r2);
    disp 'Distribution of Fluid Temperature at Single-Phase Mode of Second Section:';
    disp (TF2);
    T2=TF2(:,:,i)
    disp 'Distribution of Gas Temperature at Single-Phase Mode of Second Section:';
    disp (TG2);
    T_g2=TEMP32(i);
    disp 'Distribution of Heat Transfer Coefficient at Single-Phase Mode of Second Section:';
    disp (HT22);
    if TF2(k,j,i)>373.15
        break
    end
    %End of Section 2
    %% Start of Section 3 (3*Ls - Ls) = Ls
    e333=5;
    e33=5;%input('Number of radial element in Single-Phase Mode of Third Section from Loop No.1:');
    e3=5;%input('Number of element in Single-Phase Mode of Third Section from Loop No.1:');
    Length3=Ls/e3;
    l_shell=Length3*cos(teta);
    T_t_ex3=zeros(e3,1);
    T_t_in3=zeros(e3,1);
    T_out3=zeros(e333,e33,e3);
    T_g_out3=zeros(e3,1);
    T_g_in3=T_g2;
    for i=1:e3
        i=i;
        r_t_in=D_t_in/2;
        r_0=r_t_in/e33;
        r_in=r_0;
        for j=1:e33
            for k=1:1:e111
                if k==1
                    Fi=0;
                elseif k==2
                    Fi=90;
                elseif k==3
                    Fi=180;
                elseif k==4
                    Fi=270;
                elseif k==5
                    Fi=360;
                end
                if i==1
                    T_in3=T2(1,1,:);
                end
                mu_g3=(4*(10.^(-5)))*exp(4170/T_g_in3);%Gas Dynamic Viscosity
                C_p_g3=976.78+(1.0637*T_g_in3);%Gas Specific Heat
                k_g3=0.43482+(5*(10.^(-4))*T_g_in3);%Gas Thermal Conductivity
                ro_g3=2729.3-0.73*T_g_in3;%Gas Density
                v_g3=mdat_g/(A_g*ro_g3);
                Pr_g3=(C_p_g3*mu_g3)/k_g3;
                Re_g3=(ro_g3*v_g3*(D_t_ex/cos(teta)))/mu_g3;
                Nu_D_g3=0.3+(((0.62*(Re_g3.^0.5)*(Pr_g3.^(1/3)))/((1+((0.4/Pr_g3).^(2/3))).^0.25))*((1+((Re_g3/282000).^(5/8))).^0.8));
                h_g3=Nu_D_g3*k_g3*cos(teta)/D_t_ex;
                mu_f3=((2.1897e-11)*(T_in3.^4))-((3.055e-8)*(T_in3.^3))+((1.6028e-5)*(T_in3.^2))-(0.0037524*T_in3)+0.33158;
                C_p_f3=((1.1105e-5)*(T_in3.^3))-(0.0031078*(T_in3.^2))-(1.478*T_in3)+4631.9;
                k_f3=((1.5362e-8)*(T_in3.^3))-((2.261e-05)*(T_in3.^2))+(0.010879*T_in3)-1.0294;
                ro_f3=((-1.5629e-5)*(T_in3.^3))+(0.011778*(T_in3.^2))-(3.0726*T_in3)+1227.8;
                d_f3=0.1*((6*M)/(Av*pi*ro_f_ref));
                mu_nf3=mu_f3/(1-(34.87*((d_p/d_f3).^(-0.3))*(fi.^1.03)));
                k_nf3=k_f3*((k_p+((n-1)*k_f3)-((n-1)*fi*(k_f3-k_p)))/(k_p+((n-1)*k_f3)+(fi*(k_f3-k_p))));
                A13=h_g3;% A(number coefficient)(number section)
                r=(r_t_ex+r_t_in)/2;
                A23=k_t/(r*log(r_t_ex/r_t_in));
                g=9.81;
                ro_nf3=(ro_p*fi)+(ro_f3*(1-fi));
                P_g3=ro_nf3*g*Length3*cos(teta);
                A_il=pi*((2*r_in).^2);
                u_m3=mdat/(A_il*ro_nf3);
                C_p_nf3=((fi*ro_p*C_p_p)+((1-fi)*ro_f3*C_p_f3))/ro_nf3;
                Pr_nf3=(C_p_nf3*mu_nf3)/k_nf3;
                Re_nf3=(ro_nf3*u_m3*D_t_in)/(mu_nf3);
                if Re_nf3<2300
                    Nu_nf3=4.36;
                else
                    if P_g3<0
                        P3=-P_g3;
                    else
                        P3=P_g3;
                    end
                    f3=(2*P3*D_t_in)/(Length3*ro_nf3*(u_m3.^2));
                    Nu_nf3=((f3/8)*(Re_nf3-1000)*Pr_nf3)/(1+(12.7*((f3/8).^0.5))*((Pr_nf3.^(2/3))-1));
                end
                Nu_loc3=Nu_nf3*(((-2.331e-05)*(Fi.^2))+((8.424e-03)*Fi)+0.4576);
                if fi==0
                    hlv_nf3=2317000;
                else
                    C03=(Anp*(fi.^5))+(Bnp*(fi.^4))+(Cnp*(fi.^3))+(Dnp*(fi.^2))+(Enp*fi)+Fnp;
                    C13=(alphanp*(fi.^5))+(bettanp*(fi.^4))+(gammanp*(fi.^3))+(sigmanp*(fi.^2))+(epsilonnp*fi)+omeganp;
                    hlv_nf3=C03*(P3.^C13);
                end
                h_nf3=Nu_loc3*k_nf3/(2*r_in);
                G3=0.00122*(((k_nf3.^0.79)*(C_p_nf3.^0.45)*(ro_f3.^0.49))/((sigma.^0.5)*(mu_nf3.^0.29)*(hlv_nf3.^0.24)*(ro_v.^0.24)))*(P3.^0.75);
                S_chen3=1/(1+((2.53*(10.^-6))*(Re_nf3.^1.17)));
                if i==1
                    h_nb3=G3*((373.15-T_in3).^0.24);
                else
                    h_nb3=G3*((T_t_in3-T_out3(1,1,i)).^0.24);
                end
                A33=h_nf3;
                A331=h_nb3*S_chen3*pi*D_t_in;
                A43=(mdat_g*C_p_g3)/(l_shell*pi*(D_shell-D_h));
                A53=(mdat*C_p_nf3)/(Length3*pi*(2*r_in));
                b13=((A13)/2)*T_g_in3;
                b23=0;
                b33=-373.15*A331;
                b43=A43*T_g_in3;
                b53=(-A53)*T_in3;
                %Unknowns: q, T_g_out, T_t_ex, T_t_in, T_out
                Y3=[1,(-A13)/2,A13,0,0%Eq1
                    1,0,-A23,A23,0%Eq2
                    1,0,0,-(A33+A331),A33%Eq3
                    1,A43,0,0,0%Eq4
                    1,0,0,0,-A53];%Eq5
                B3=[b13;b23;b33;b43;b53];
                X3=Y3\B3;
                T_out3(k,j,i)=X3(5,:);
                TEMP13(k,j,i)=X3(5,:);
                TF3=TEMP13;
                Nusselt33(k,j,i)=Nu_loc3;
                Nu3=Nusselt33;
                Friction_Factor(j,i)=f3;
                FF3=Friction_Factor;
                Velocity_r(j,i)=u_m3;
                V_r3=Velocity_r;
                Reynolds33(j,i)=Re_nf3;
                Re3=Reynolds33;
                if Re3(j,:)<0
                    break
                end
                HT3(k,j,i)=h_nf3;
                HT33=HT3;
            end
            r_in=r_in+r_0;
        end
        T_in3=T_out3(1,1,i);
        T_t_in3=X3(4,:);
        TEMP23(i)=T_t_in3;
        TIS3=TEMP23;
        T_g_out3(i)=X3(2,:);
        T_g_in3=T_g_out3(i);
        TEMP33(i)=X3(2,:);
        TG3=TEMP33;
        if TF3(k,j,i)>373.15
            break
        end
    end
    disp 'Distribution of Reynolds Number at Single-Phase Mode of Third Section:';
    disp (Re3);
    Re3=Re3(:,i)
    disp 'Distribution of Internal Surface Temperature at Single-Phase Mode of Third Section:';
    disp (TIS3);
    TIS3=TIS3(i)
    disp 'Distribution of Nusselt Number at Single-Phase Mode of Third Section:';
    disp (Nu3);
    Nu3=Nu3(:,:,i)
    disp 'Distribution of Friction Factor at Single-Phase Mode of Third Section:';
    disp (FF3);
    disp 'Distribution of Radius Velocity at Single-Phase Mode of Third Section:';
    disp (V_r3);
    disp 'Distribution of Fluid Temperature at Single-Phase Mode of Third Section:';
    disp (TF3);
    T3=TF3(:,:,i)
    disp 'Distribution of Gas Temperature at Single-Phase Mode of Third Section:';
    disp (TG3);
    T_g3=TEMP33(i);
    disp 'Distribution of Heat Transfer Coefficient at Single-Phase Mode of Third Section:';
    disp (HT33);
    T_g_in1=T_g3;
    T_t_in1=T_t_in3;
    if TF3(k,j,i)>373.15
        break
    end
end