function out = eval_model(theta)
% theta = [p2a, p2b, p2c, p2d, ...
%          r_na, r_cl, r_am, r_amm, ...   % r_* = Kc_*/Kd_*
%          e, Lm_nm, k                    % porosity, thickness[nm], mass-transfer
%         ];


% --------- Parameters to be assessed ----------
p1a = theta(1);  p1b = theta(2);  p1c = theta(3);  p1d = theta(4);
Kcna = theta(5); Kccl = theta(6); Kcam = theta(7); Kcamm = theta(8);
Kdna=theta(12); Kdcl=theta(13); Kdam=theta(14); Kdamm=theta(15); 
e    = theta(9);
Lm   = theta(10) * 1e-9;   % nm -> m
k    = theta(11);          % m/s


% --------- Other fixed parameters ----------
TAN=10.7; % mM
Ccl=TAN;
pH_space=7;%5.3:0.1:11;
vf=1.5e-6;               % m/s
lmh=vf*3.6e6;            % L/m2/h
pKa=9.25; pKb=4.76;

Kf_na=0.2; Kf_cl=0.2; Kf_am=0.2; Kf_amm=0.2;


Dna=1.33e-9; Dcl=2.03e-9; Dam=1.954e-9; Damm=3.35e-9;
Dh=8.24e-9; Doh=4.51e-9;

p2a=1; p2b=1; p2c=1; p2d=1;

K_COOH1=10^(-5.23); K_COOH2=10^(-8.97); K_NH3=10^(-4.74);
X_NH2=36; X_COOH1=82; X_COOH2=350; % mM

% --------- main loop ----------

pH_perm=zeros(1,length(pH_space));
TAN_perm=zeros(1,length(pH_space));
Rej=zeros(1,length(pH_space));

for j=1:length(pH_space)
    pH=pH_space(j);
    pOH=14-pH;
    H=10^-pH; OH=10^-pOH;
    CH=1000*10^-pH; COH=1000*10^-pOH;

    % Membrane charge
    X = X_NH2/(1+K_NH3/H) - X_COOH1/(1+H/K_COOH1) - X_COOH2/(1+H/K_COOH2);

    % TAN speciation
    Camm=10^(-pKa)/(10^(-pH)+10^(-pKa))*TAN; % NH3
    Cam=TAN-Camm;                            % NH4+
    Cna=Camm+COH-1000*10^(-8.5);             % electroneutrality
    Ccl=TAN;

    % ---- CP ----
    Cw=zeros(4,200); Cw(1,1)=Cna; Cw(2,1)=Cam; Cw(3,1)=Ccl; Cw(4,1)=Camm;

    Na_perm=NaN; cp_am=NaN; cp_cl=NaN; cp_amm=NaN;
    CamL=NaN; CamR=NaN; CclL=NaN; CclR=NaN; CnaL=NaN; CnaR=NaN;
    phiL=NaN; cc=NaN; d1=NaN; CammL=NaN; CammR=NaN; 

    for i=1:size(Cw,2)
        syms d0
        eq=Cw(1,i)*p1a*p2a*exp(-d0)+Cw(2,i)*p1c*p2c*exp(-d0)-Cw(3,i)*p1b*p2b*exp(d0)+X;
        phiL=vpasolve(eq,d0);

        CnaL=Cw(1,i)*p1a*p2a*exp(-phiL);
        CamL=Cw(2,i)*p1c*p2c*exp(-phiL);
        CclL=Cw(3,i)*p1b*p2b*exp(phiL);

        syms a b am c d
        eq1=-Dna*Kdna*e*(a-CnaL)/Lm - Dna*Kdna*e*0.5*(CnaL+a)*(c-phiL)/Lm + Kcna*0.5*(a+CnaL)*vf - vf*a/(p1a*p2a*exp(-d));
        eq2=-Dam*Kdam*e*(am-CamL)/Lm - Dam*Kdam*e*0.5*(CamL+am)*(c-phiL)/Lm + Kcam*0.5*(am+CamL)*vf - vf*am/(p1c*p2c*exp(-d));
        eq3= Dcl*Kdcl*e*(b-CclL)/Lm +(-1)*Dcl*Kdcl*e*0.5*(b+CclL)*(c-phiL)/Lm - Kccl*0.5*(b+CclL)*vf + vf*b/(p1b*p2b*exp(d));
        eq4=a+am-b+X;
        eq5=a/(p1a*p2a*exp(-d)) + am/(p1c*p2c*exp(-d)) - b/(p1b*p2b*exp(d));
        [aa,bb,NH4,cc,dd]=vpasolve(eq1,eq2,eq3,eq4,eq5,a,b,am,c,d);

        CnaR=aa; CamR=NH4; CclR=bb; phi_p=cc-phiL; d1=dd;

        cp_na=CnaR/(p1a*p2a*exp(-d1));
        cp_am=CamR /(p1c*p2c*exp(-d1));
        cp_cl=CclR /(p1b*p2b*exp( d1));

        C_W_na=(Cna-cp_na)*exp(vf/k)+cp_na;
        C_W_am=(Cam-cp_am)*exp(vf/k)+cp_am;
        C_W_cl=C_W_am+C_W_na;

        if abs(Cw(1,i)-C_W_na)/C_W_na<1E-4, Cw(1,i+1)=Cw(1,i); else, Cw(1,i+1)=0.5*(C_W_na+Cw(1,i)); end
        if abs(Cw(2,i)-C_W_am)/C_W_am<1E-4, Cw(2,i+1)=Cw(2,i); else, Cw(2,i+1)=0.5*(C_W_am+Cw(2,i)); end
        if abs(Cw(3,i)-C_W_cl)/C_W_cl<1E-4, Cw(3,i+1)=Cw(3,i); else, Cw(3,i+1)=0.5*(C_W_cl+Cw(3,i)); end
        if abs(Cw(1,i)-C_W_na)/C_W_na<1E-4 && abs(Cw(2,i)-C_W_am)/C_W_am<1E-4 && abs(Cw(3,i)-C_W_cl)/C_W_cl<1E-4
            Na_perm=cp_na; break;
        end
    end

    % NH3 transport
    for iter=1:size(Cw,2)
        CammL=Cw(4,iter)*p1d*p2d;
        syms f
        eq6=-Damm*Kdamm*e*(f-CammL)/Lm + Kcamm*0.5*(f+CammL)*vf - vf*f/(p1d*p2d);
        ff=vpasolve(eq6,f);
        CammR=ff;
        cp_amm=CammR/(p1d*p2d);

        C_W_amm=(Camm-cp_amm)*exp(vf/k)+cp_amm;
        if abs(Cw(4,iter)-C_W_amm)/C_W_amm<1E-4, break; else, Cw(4,iter+1)=0.5*(C_W_amm+Cw(4,iter)); end
    end

    % NH3 to NH4+ + OH-
    syms x
    equation=(10^-pOH + x)*(cp_am*1e-3 + x) - 10^(-pKb)*(cp_amm*1e-3 - x);
    sol=solve(equation,x);
    x_solve=double(sol(2,1));
    NH3_final=cp_amm - x_solve*1000;
    NH4_final=cp_am + x_solve*1000;

    TAN_perm(j)=NH3_final + NH4_final;
    Rej(j)=1 - TAN_perm(j)/TAN;
    pH_perm(j)=14 - (-log10(10^-pOH + x_solve));
end

% --------- Define SA output ----------
out.R_mean = mean(Rej);                         
out.R_min  = min(Rej);                          
out.pHperm_mean = mean(pH_perm);               

out.Rej_curve = Rej(:);
end
