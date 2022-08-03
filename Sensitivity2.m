function [Pao2]=Sensitivity2(Cv,Rp,Emax,HR,v)

dt = 0.5e-3 ; % [sec]
N_cycle =floor(60/(HR*dt)); % number of steps per heart cycle
Heart_cycles = 20 ; % total heart cycles

% given parameters
V0 = 15 ; % [ml] V0 for Plv calculation
% Vascular constants:
Ra = 0.1; % arterial resistance
Rv = 0.01; % venous filling resistance
Ca = 2.0; % arterial compliance
Ed=10/(120-V0);%diastolic elasticity is constant (from tirgul 4)

N_Systole = round((1/3)*N_cycle); % number of points per ventricle Systole
En(1:N_Systole)=0.5*(1+sin(2*pi*(1:N_Systole)/(N_Systole)-pi/2)); %normalized systolic E
En(N_Systole+1:N_cycle)=0; % diastole, elasticity=0
E=max(Ed,En*Emax) ; % combining diastolic and systolic elasticity functions


 

% Initiate variables:
%Volume [ml]
Vlv_s(1,1:N_cycle) = 120*v; % left ventricle
Va_s(1,1:N_cycle) = 270*v; % arteries
Vv_s(1,1:N_cycle) = 2700*v; % veins
%Pressure [mmHg]
Plv(1:N_cycle) = 0; % left ventricle
Pa(1:N_cycle) = 70; % arterial capacitor
Pv(1:N_cycle) = 9; % venous filling
Pao(1:N_cycle) = 70; % aorta
%Flow [ml/sec]
Qlv(1:N_cycle) = 0; % left ventricle (outflow)
Qp(1:N_cycle) = 0; % peripheral resistance
Qv(1:N_cycle) = 0; % ventricle filling (inflow)

Vlv2=  zeros(1,Heart_cycles*N_cycle); 
Va12=  zeros(1,Heart_cycles*N_cycle); 
Vv2=  zeros(1,Heart_cycles*N_cycle); 
Plv2= zeros(1,Heart_cycles*N_cycle);
Pa2=  zeros(1,Heart_cycles*N_cycle); 
Pv2= zeros(1,Heart_cycles*N_cycle); 
Pao2= zeros(1,Heart_cycles*N_cycle); 
Qv2= zeros(1,Heart_cycles*N_cycle); 
Qlv2= zeros(1,Heart_cycles*N_cycle); 
Qp2= zeros(1,Heart_cycles*N_cycle); 

%calculating all variables for each cycle at N points:
for CycleIdx=1:Heart_cycles
 for StepInCycle=2:N_cycle
      %Volume [ml] …
      Vlv_s(StepInCycle)=Vlv_s(StepInCycle-1)+(Qv(StepInCycle-1)-Qlv(StepInCycle-1))*dt;
      Va_s(StepInCycle)=Va_s(StepInCycle-1)+(Qlv(StepInCycle-1)-Qp(StepInCycle-1))*dt;
      Vv_s(StepInCycle)=Vv_s(StepInCycle-1)+(Qp(StepInCycle-1)-Qv(StepInCycle-1))*dt;
      
      %Pressure [mmHg] …
      Plv(StepInCycle)=E(StepInCycle)*(Vlv_s(StepInCycle)-V0);
      Pa(StepInCycle)=(1/Ca)*Va_s(StepInCycle);
      Pv(StepInCycle)=(1/Cv)*Vv_s(StepInCycle);
      
      if Plv(StepInCycle)>Pa(StepInCycle)
          Pao(StepInCycle)=Plv(StepInCycle);
      else
          Pao(StepInCycle)=Pa(StepInCycle);
      end
      
       %Flow [ml/sec] …
       Qv(StepInCycle)=max(0,(Pv(StepInCycle)-Plv(StepInCycle))/Rv);
       Qlv(StepInCycle)=max(0,(Pao(StepInCycle)-Pa(StepInCycle))/Ra);
       Qp(StepInCycle)=(Pa(StepInCycle)-Pv(StepInCycle))/Rp;
       

 end
 
  % Saving variables from cycle to continuous variables: …
 a=(CycleIdx-1)*(N_cycle)+1;
 b=CycleIdx*N_cycle;
 Plv2(a:b)=Plv;
 Pa2(a:b)=Pa;
 Pao2(a:b)=Pao;
 Pv2(a:b)=Pv;
 Vlv2(a:b)=Vlv_s;
 Val2(a:b)=Va_s;
 Qp2(a:b)=Qp;
 Qlv2(a:b)=Qlv;
 Vv2(a:b)=Vv_s;
 
  % Update the initial variables before the next cycle: …
 Plv(1)=Plv(N_cycle);
 Pa(1)=Pa(N_cycle);
 Pao(1)=Pao(N_cycle);
 Pv(1)=Pv(N_cycle);
 Vlv_s(1)=Vlv_s(N_cycle);
 Va_s(1)=Va_s(N_cycle);
 Qp(1)=Qp(N_cycle);
 Qlv(1)=Qlv(N_cycle);
 Vv_s(1)=Vv_s(N_cycle);
 
end
 
 
end

