TAN=10.7; % mM
Ccl=TAN;
pH_space=5.3:0.1:11;
vf=1.5*10^(-6);% m/s
lmh=vf*3.6*10^6; % L/m2/hr
pKa=9.25;
pKb=4.76;

% Transport coefficient
%Kf_na=0.2; %0.41  Value in Wang Li paper
%Kf_cl=0.2; % 0.41
%Kf_am=0.2; % am - NH4+;  amm - NH3
%Kf_amm=0.2; % Adjusted manually
Kdna=0.2;%Kf_na;%
Kdcl=0.2;%Kf_cl;%
Kdam=0.066;%Kf_am;
Kdamm=0.371;%Kf_amm;%
Kcna=0.2;%Kf_na;
Kccl=0.2;%Kf_cl;
Kcam=1.3;%Kf_am;
Kcamm=1.32;%Kf_amm;

Dna=1.33*10^-9; % Na+ diffusivity (m2/s)
Dcl=2.03*10^-9; % Cl diffusivity (m2/s)
Dam=1.954e-9; % NH4 diffusivity (m2/s)
Damm=3.35*10^-9; % NH3 diffusivity (m2/s)
Dh=8.24*10^-9; % H+ diffusivity (m2/s)
Doh=4.51*10^-9; % OH- diffusivity (m2/s)

p1a=0.09; % 0.07 (1-ldna)^2; % Na steric partitioning, fitting paprameter
p1b=0.078; % 0.07 (1-ldcl)^2; % Cl steric paritioning, fitting parameter
p1c=0.056; % NH4 steric paritioning, fitting parameter
p1d=0.01; % NH3 steric paritioning, fitting parameter
p2a=1;
p2b=1;
p2c=1;
p2d=1; % dielectric partitioning

%Hydrodynamics to calculate k
D=nthroot(Dna*Dcl*Dam,3); % mean diffusion (m2/s)
Re=312.3342; % Reynolds number calculated from Excel
Sc=450.03; % Schimit number
Sh=0.16*(Re^0.605)*(Sc^0.42); % Sherwood number
viscosity=1.23*10^(-6)*exp(0.124*TAN*0.001+6.59);
k=Sh*D/0.004; % mass transfer coefficient (m/s)

% Membrane
e=0.1; % effective porosity, typically 0.05, set to 0.1 here due to low pressure membrane
Lm=150*10^-9; %150 membrane thickness (m)
K_COOH1=10^(-5.23); %Parameters for membrane charge calculation
K_COOH2=10^(-8.97);
K_NH3=10^(-4.74);
X_NH2=36; % mM
X_COOH1=82; % mM
X_COOH2=350; % mM


X=zeros(1,length(pH_space));
pH_perm=zeros(1,length(pH_space));
TAN_perm=zeros(1,length(pH_space));
Rej=zeros(1,length(pH_space));
Rej_am=zeros(1,length(pH_space));
Rej_amm=zeros(1,length(pH_space));
Rej_na=zeros(1,length(pH_space));
NH3_perm=zeros(1,length(pH_space));
NH4_perm=zeros(1,length(pH_space));
P=zeros(1,length(pH_space));
CW=zeros(1,length(pH_space));
Ksol=zeros(1,length(pH_space)); % m/s
Don0=zeros(1,length(pH_space)); % feed side Donnan
Don1=zeros(1,length(pH_space)); % permeate side Donnan
Na_perm=zeros(1,length(pH_space)); % Na conc in perm
SP1=zeros(1,length(pH_space));  % TAN permeatability
SP2=zeros(1,length(pH_space));  % NH4 permeatability
SP3=zeros(1,length(pH_space));  % NH3 permeatability
flux_am=zeros(1,length(pH_space));
flux_amm=zeros(1,length(pH_space));
A=zeros(1,length(pH_space));
adv=zeros(1,length(pH_space));
diff=zeros(1,length(pH_space));
ele=zeros(1,length(pH_space));
adv_amm=zeros(1,length(pH_space));
diff_amm=zeros(1,length(pH_space));
NH3_partition=zeros(1,length(pH_space));

for j=1:length(pH_space)
    pH=pH_space(j);
    pOH=14-pH;
    H=10^-pH; %proton mol/L
    OH=10^-pOH; % mol/L
    CH=1000*10^-pH; % proton (mM=mol/m3)
    COH=1000*10^-pOH; % OH (mM)
    % Memrbane charge
    X(j)=X_NH2/(1+K_NH3/H)-X_COOH1/(1+H/K_COOH1)-X_COOH2/(1+H/K_COOH2);
    Camm=10^(-pKa)/(10^(-pH)+10^(-pKa))*TAN; % Free ammonia
    Cam=TAN-Camm; % Ammonium
    Cna=Camm+COH-1000*10^(-8.5); % Na
    NH3_partition(j)=Camm/TAN;

    % Concentration polarization
    Cw=zeros(4,200);
    Cw(1,1)=Cna;
    Cw(2,1)=Cam;
    Cw(3,1)=Ccl;
    Cw(4,1)=Camm;% initial guess/same as feed

    % Transport for ions
      for i=1:size(Cw,2)
        syms d0
        eq=Cw(1,i)*p1a*p2a*exp(-d0)+Cw(2,i)*p1c*p2c*exp(-d0)-Cw(3,i)*p1b*p2b*exp(d0)+X(j);
        [phiL]=vpasolve(eq,d0);
        CnaL=Cw(1,i)*p1a*p2a*exp(-phiL); %Na donnon partition
        CamL=Cw(2,i)*p1c*p2c*exp(-phiL); %Ammonia donnon
        CclL=Cw(3,i)*p1b*p2b*exp(phiL);
        %d0=x; % dimensionless
        % Entrance
        %Cna0=Cw(i)*p1a*p2a*exp(-d0); % Na at x=0
        %Ccl0=Cw(i)*p1b*p2b*exp(d0); % Cl at x=0
        syms a b am c d % a-Na conc b-Cl conc, an-ammonia conc, c-potential diff, d-donnan diff
        eq1=-Dna*Kdna*e*(a-CnaL)/Lm-Dna*Kdna*e*0.5*(CnaL+a)*(c-phiL)/Lm+Kcna*0.5*(a+CnaL)*vf-vf*a/(p1a*p2a*exp(-d));
        eq2=-Dam*Kdam*e*(am-CamL)/Lm-Dam*Kdam*e*0.5*(CamL+am)*(c-phiL)/Lm+Kcam*0.5*(am+CamL)*vf-vf*am/(p1c*p2c*exp(-d)); 
        eq3=Dcl*Kdcl*e*(b-CclL)/Lm+(-1)*Dcl*Kdcl*e*0.5*(b+CclL)*(c-phiL)/Lm-Kccl*0.5*(b+CclL)*vf+vf*b/(p1b*p2b*exp(d));
        eq4=a+am-b+X(j);
        eq5=a/(p1a*p2a*exp(-d))+am/(p1c*p2c*exp(-d))-b/(p1b*p2b*exp(d));
        [aa,bb,NH4,cc,dd]=vpasolve(eq1,eq2,eq3,eq4,eq5,a,b,am,c,d);
        CnaR=aa;
        CamR=NH4;
        CclR=bb;
        phi_p=cc-phiL; % potential diff in pore dimensionless
        d1=dd;
        cp_na=CnaR/(p1a*p2a*exp(-d1));
        cp_am=CamR/(p1c*p2c*exp(-d1));
        cp_cl=CclR/(p1b*p2b*exp(d1));
        
        % Concentration polarization
        C_W_na=(Cna-cp_na)*exp(vf/k)+cp_na; % membran surface concentration (mol/m3)
        C_W_am=(Cam-cp_am)*exp(vf/k)+cp_am;
        C_W_cl=C_W_am+C_W_na;
        if abs(Cw(1,i)-C_W_na)/C_W_na<1E-4
          Cw(1,i+1)=Cw(1,i);
        else
          Cw(1,i+1)=0.5*(C_W_na+Cw(1,i));
        end

        if abs(Cw(2,i)-C_W_am)/C_W_am<1E-4
          Cw(2,i+1)=Cw(2,i);
        else
          Cw(2,i+1)=0.5*(C_W_am+Cw(2,i));
        end
        
        if abs(Cw(3,i)-C_W_cl)/C_W_cl<1E-4
          Cw(3,i+1)=Cw(3,i);
        else
          Cw(3,i+1)=0.5*(C_W_cl+Cw(3,i));
        end
        % stop loop when all ions' polarization concentrations fit
        if abs(Cw(1,i)-C_W_na)/C_W_na<1E-4 && abs(Cw(2,i)-C_W_am)/C_W_am<1E-4 && abs(Cw(3,i)-C_W_cl)/C_W_cl<1E-4
          break;
        end
        Na_perm(j)=cp_na;
        %disp(i)
        
      end

      %Transport for NH3
      for iter=1:size(Cw,2)
          CammL=Cw(4,iter)*p1d*p2d; % NH3 at x=0
          syms f 
          eq6=-Damm*Kdamm*e*(f-CammL)/Lm+Kcamm*0.5*(f+CammL)*vf-vf*f/(p1d*p2d);
          [ff]=vpasolve(eq6,f);
          CammR=ff;
          cp_amm=CammR/(p1d*p2d);
          C_W_amm=(Camm-cp_amm)*exp(vf/k)+cp_amm; % membran surface concentration (mol/m3)
          if abs(Cw(4,iter)-C_W_amm)/C_W_amm<1E-4
            break
          else
            Cw(4,iter+1)=0.5*(C_W_amm+Cw(4,iter));
          end
          %disp(iter);
          
      end

     % Hydrolysis of the transported NH3
     % NH3 --> NH4 + OH  pKb=4.76
     syms x
     equation=(OH+x)*(cp_am*0.001+x)-10^(-pKb)*(cp_amm*0.001-x); % Convert concentration of NH34 into mol/L
     solution=solve(equation, x);
     x_solve=double(solution(2,1));

     NH3_final=cp_amm-x_solve*1000;
     NH4_final=cp_am+x_solve*1000;
     OH_final=OH+x_solve;
     pOH_final=-log10(OH_final);
     pH_perm(j)=14-pOH_final;
     TAN_perm(j)=NH3_final+NH4_final;
     NH3_perm(j)=NH3_final;
     NH4_perm(j)=NH4_final;

     Rej(j)=1-(NH3_final+NH4_final)/TAN;
     Rej_am(j)=1-cp_am/Cam;
     Rej_amm(j)=1-cp_amm/Camm;
     Rej_na(j)=1-cp_na/Cna;

     SP1(j)=(lmh*(cp_am+cp_amm))/(C_W_amm+C_W_am-cp_am-cp_amm);
     SP2(j)=(lmh*cp_am)/(C_W_am-cp_am);
     SP3(j)=(lmh*cp_amm)/(C_W_amm-cp_amm);

     flux_am(j)=vf*cp_am*3600*1000; % mol/m2/s --> mmol/m2/h
     flux_amm(j)=vf*cp_amm*3600*1000; % mol/m2/s--> mmol/m2/h

     adv(j)=Kcam*0.5*(CamR+CamL)*vf;
     diff(j)=-Dam*Kdam*e*(CamR-CamL)/Lm;
     ele(j)=-Dam*Kdam*e*0.5*(CamL+CamR)*(cc-phiL)/Lm;

     adv_amm(j)=Kcamm*0.5*(CammR+CammL)*vf;
     diff_amm(j)=-Damm*Kdamm*e*(CammR-CammL)/Lm;

     %Don0(j)=double(phiL);
     %Don1(j)=double(d1);

     % Pressure
    c1=0.5*(CamL+CamR);
    c2=0.5*(CclL+CclR);
    c3=0.5*(CnaL+CnaR);
    c4=0.5*(CammL+CammR);
    f1=Kcam/(Kdam*Dam); % NH4 friction with fluid (s/m2)
    f2=Kccl/(Kdcl*Dcl); % Cl friction with fluid (s/m2)
    f3=Kcna/(Kdna*Dna);
    f4=Kcamm/(Kdamm*Damm);
    f=1.75*10^15; %9.5 7.75 fluid membrane friction (mol*s/m5)
    p=-f*vf*Lm-c1*f1*(1-Kcam)*vf*Lm-c2*f2*(1-Kccl)*vf*Lm-c3*f3*(1-Kcna)*vf*Lm-c4*f4*(1-Kcamm)*vf*Lm-e*Kcam*(CamR-CamL)-e*Kccl*(CclR-CclL)-e*Kcna*(CnaR-CnaL)-e*Kcamm*(CammR-CammL)+1/3*(Kcam+Kccl+Kcna)*X(j)*e*phi_p;
    ph=-p+(C_W_am+C_W_cl+C_W_na+C_W_amm-cp_am-cp_cl-cp_na-cp_amm); % mole/m3  Osmotic pressure
    Ph=ph*1*(8.3144*298)/1000*0.01; % bar
    P(j)=Ph;
    A(j)=lmh/(Ph-8.3144*298*(C_W_am+C_W_cl+C_W_na+C_W_amm-cp_am-cp_cl-cp_na-cp_amm)/1000*0.01);
disp(j);
end