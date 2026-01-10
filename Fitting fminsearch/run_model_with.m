function [Rej_out, A_out, pHperm_out, SP_out] = run_model_with(p1a,p1b,p1c,p1d,f, pH_space)

% ----------------- (constants & baseline from your script) -----------------
TAN=10.7; Ccl=TAN;
vf=1.5e-6; lmh=vf*3.6e6;
pKa=9.25; pKb=4.76;

Kdna=0.2; Kdcl=0.2; Kdam=0.066; Kdamm=0.371;
Kcna=0.2; Kccl=0.2; Kcam=1.3;  Kcamm=1.32;

Dna=1.33e-9; Dcl=2.03e-9; Dam=1.954e-9; Damm=3.35e-9;

p2a=1; p2b=1; p2c=1; p2d=1;

% Hydrodynamics (k)
D = nthroot(Dna*Dcl*Dam,3);
Re=312.3342; Sc=450.03;
Sh=0.16*(Re^0.605)*(Sc^0.42);
viscosity=1.23e-6*exp(0.124*TAN*0.001+6.59); %#ok<NASGU>
k=Sh*D/0.004;

% Membrane
e=0.1; Lm=150e-9;
K_COOH1=10^(-5.23); K_COOH2=10^(-8.97); K_NH3=10^(-4.74);
X_NH2=36; X_COOH1=82; X_COOH2=350;

% ----------------- prealloc -----------------
npH = numel(pH_space);
Rej_out    = zeros(1,npH);
A_out      = zeros(1,npH);
pHperm_out = zeros(1,npH);
SP_out     = zeros(1,npH);   % ★ 新增：存 SP1(j)

for j=1:npH
    pH = pH_space(j); pOH=14-pH;
    H  = 10^-pH;  OH=10^-pOH;

    % membrane charge
    Xchg = X_NH2/(1+K_NH3/H) - X_COOH1/(1+H/K_COOH1) - X_COOH2/(1+H/K_COOH2);

    % TAN speciation
    Camm = 10^(-pKa)/(10^(-pH)+10^(-pKa))*TAN; % NH3
    Cam  = TAN - Camm;                          % NH4+
    Cna  = Camm + 1000*10^-pOH - 1000*10^(-8.5);
    Ccl  = TAN;

    % --- Concentration polarization loop ---
    Cw = zeros(4,200);
    Cw(1,1)=Cna; Cw(2,1)=Cam; Cw(3,1)=Ccl; Cw(4,1)=Camm;

    cp_am = NaN; cp_amm = NaN; cp_na = NaN; C_W_am=NaN; C_W_amm=NaN;
    CamL=NaN; CamR=NaN; CclL=NaN; CclR=NaN; CnaL=NaN; CnaR=NaN; cc=NaN; phiL=NaN; d1=NaN;

    try
        % ------- NH4+, Na+, Cl-  -------
        for i=1:size(Cw,2)
            syms d0
            eq=Cw(1,i)*p1a*p2a*exp(-d0) + Cw(2,i)*p1c*p2c*exp(-d0) ...
               - Cw(3,i)*p1b*p2b*exp(d0) + Xchg;
            phiL=vpasolve(eq,d0);

            CnaL=Cw(1,i)*p1a*p2a*exp(-phiL);
            CamL=Cw(2,i)*p1c*p2c*exp(-phiL);
            CclL=Cw(3,i)*p1b*p2b*exp( phiL);

            syms a b am c d
            eq1=-Dna*Kdna*e*(a-CnaL)/Lm - Dna*Kdna*e*0.5*(CnaL+a)*(c-phiL)/Lm ...
                + Kcna*0.5*(a+CnaL)*vf - vf*a/(p1a*p2a*exp(-d));
            eq2=-Dam*Kdam*e*(am-CamL)/Lm - Dam*Kdam*e*0.5*(CamL+am)*(c-phiL)/Lm ...
                + Kcam*0.5*(am+CamL)*vf - vf*am/(p1c*p2c*exp(-d));
            eq3= Dcl*Kdcl*e*(b-CclL)/Lm +(-1)*Dcl*Kdcl*e*0.5*(b+CclL)*(c-phiL)/Lm ...
                - Kccl*0.5*(b+CclL)*vf + vf*b/(p1b*p2b*exp(d));
            eq4=a+am-b+Xchg;
            eq5=a/(p1a*p2a*exp(-d)) + am/(p1c*p2c*exp(-d)) - b/(p1b*p2b*exp(d));
            [aa,bb,NH4,cc,dd]=vpasolve(eq1,eq2,eq3,eq4,eq5,a,b,am,c,d);

            CnaR=aa; CamR=NH4; CclR=bb; d1=dd;

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
                break;
            end
        end

        % ------- NH3 transport -------
        for iter=1:size(Cw,2)
            CammL=Cw(4,iter)*p1d*p2d;
            syms fNH3
            eq6=-Damm*Kdamm*e*(fNH3-CammL)/Lm + Kcamm*0.5*(fNH3+CammL)*vf ...
                - vf*fNH3/(p1d*p2d);
            ff=vpasolve(eq6,fNH3);
            CammR=ff;
            cp_amm=CammR/(p1d*p2d);
            C_W_amm=(Cw(4,iter)-cp_amm)*exp(vf/k)+cp_amm;
            if abs(Cw(4,iter)-C_W_amm)/C_W_amm<1E-4
                break;
            else
                Cw(4,iter+1)=0.5*(C_W_amm+Cw(4,iter));
            end
        end

        % -------  SP1(j) -------
        SP_out(j) = ( lmh * (cp_am + cp_amm) ) / ...
                    ( C_W_amm + C_W_am - cp_am - cp_amm );

        % ------- Hydrolysis to get permeate composition and pH -------
        syms x
        equation=(10^-pOH + x)*(cp_am*1e-3 + x) - 10^(-pKb)*(cp_amm*1e-3 - x);
        sol=solve(equation,x);
        x_solve=double(sol(2,1));
        NH3_final=cp_amm - x_solve*1000;
        NH4_final=cp_am + x_solve*1000;
        TAN_perm = NH3_final + NH4_final;

        Rej_out(j)    = 1 - TAN_perm/TAN;  
        pHperm_out(j) = 14 - (-log10(10^-pOH + x_solve));

        % ------- Pressure/perm. A -------
        c1=0.5*(CamL+CamR); c2=0.5*(CclL+CclR); c3=0.5*(CnaL+CnaR); c4=0.5*(CammL+CammR);
        f1=Kcam/(Kdam*Dam); f2=Kccl/(Kdcl*Dcl); f3=Kcna/(Kdna*Dna); f4=Kcamm/(Kdamm*Damm);
        p = -f*vf*Lm - c1*f1*(1-Kcam)*vf*Lm - c2*f2*(1-Kccl)*vf*Lm ...
            - c3*f3*(1-Kcna)*vf*Lm - c4*f4*(1-Kcamm)*vf*Lm ...
            - e*Kcam*(CamR-CamL) - e*Kccl*(CclR-CclL) - e*Kcna*(CnaR-CnaL) - e*Kcamm*(CammR-CammL) ...
            + 1/3*(Kcam+Kccl+Kcna)*Xchg*e*(cc-phiL);

        
        ph = -p + (C_W_am + (C_W_am+C_W_na) + C_W_na + C_W_amm ...
                   - cp_am - cp_na - cp_amm);
        Ph = ph * (8.3144*298)/1000*0.01; % bar
        A_out(j) = lmh/(Ph - (8.3144*298)*((C_W_am+(C_W_am+C_W_na)+C_W_na+C_W_amm) ...
                   - (cp_am+cp_na+cp_amm))/1000*0.01);

    catch
        % In case vpasolve fails: return NaN
        Rej_out(j)    = NaN;
        A_out(j)      = NaN;
        pHperm_out(j) = NaN;
        SP_out(j)     = NaN;
    end
end


Rej_out    = nan2pen(Rej_out,    1e6);
A_out      = nan2pen(A_out,      1e6);
pHperm_out = nan2pen(pHperm_out, 1e6);
SP_out     = nan2pen(SP_out,     1e6);
end

function v = nan2pen(v,pen)
    v(~isfinite(v)) = pen;
end
