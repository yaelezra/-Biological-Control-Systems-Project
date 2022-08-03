clear all
close all
clc
%% Flags:
Plot_flag = 1; % 0 = off , 1 = on
%% system parameters:
%part1
% Time parameters
HR = 74 ; % [BPM] % 60 + sum of last digits from all members
%HR = 80 ; % [BPM] – for 2nd part only
dt = 0.5e-3 ; % [sec]
Heart_cycles = 20 ; % total heart cycles
N_cycle =floor(60/(HR*dt)); % number of steps per heart cycle
T=N_cycle*dt;
% Heart parameters
V0 = 15 ; % [ml] V0 for Plv calculation
Emax = 2.0 ; % contractility
N_Systole = round((1/3)*N_cycle); % number of points per ventricle Systole
Ed=10/(120-V0);%diastolic elasticity is constant (from tirgul 4)
En(1:N_Systole)=0.5*(1+sin(2*pi*(1:N_Systole)/(N_Systole)-pi/2)); %normalized systolic E

En(N_Systole+1:N_cycle)=0; % diastole, elasticity=0
E = max(Ed,En*Emax) ; % combining diastolic and systolic elasticity functions
% Vascular constants:
Ra = 0.1; % arterial resistance
Rp = 1.0; % peripheral resistance
Rv = 0.01; % venous filling resistance
Ca = 2.0; % arterial compliance
Cv = 300.0; % venous compliance
% Initiate variables:
%Volume [ml]
Vlv(1) = 120; % left ventricle
Va(1) = 270; % arteries
Vv(1) = 2700; % veins
%Pressure [mmHg]
Plv(1) = 0; % left ventricle
Pa(1) = 70; % arterial capacitor
Pv(1) = 9; % venous filling
Pao(1) = 70; % aorta
%Flow [ml/sec]
Qlv(1) = 0; % left ventricle (outflow)
Qp(1) = 0; % peripheral resistance
Qv(1) = 0; % ventricle filling (inflow)
% initialization of continuous variables (Memory management):

%% Main Program

 
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
   
 for StepInCycle = 2 : N_cycle
      %Volume [ml] …
      Vlv(StepInCycle)=Vlv(StepInCycle-1)+(Qv(StepInCycle-1)-Qlv(StepInCycle-1))*dt;
      Va(StepInCycle)=Va(StepInCycle-1)+(Qlv(StepInCycle-1)-Qp(StepInCycle-1))*dt;
      Vv(StepInCycle)=Vv(StepInCycle-1)+(Qp(StepInCycle-1)-Qv(StepInCycle-1))*dt;
      
      
      %Pressure [mmHg] …
      Plv(StepInCycle)=E(StepInCycle)*(Vlv(StepInCycle)-V0);
      Pa(StepInCycle)=(1/Ca)*Va(StepInCycle);
      Pv(StepInCycle)=(1/Cv)*Vv(StepInCycle);
      
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
 Vlv2(a:b)=Vlv;
 Va2(a:b)=Va;
 Qp2(a:b)=Qp;
 Qlv2(a:b)=Qlv;
 Vv2(a:b)=Vv;
 
  % Update the initial variables before the next cycle: …
 Plv(1)=Plv(N_cycle);
 Pa(1)=Pa(N_cycle);
 Pao(1)=Pao(N_cycle);
 Pv(1)=Pv(N_cycle);
 Vlv(1)=Vlv(N_cycle);
 Va(1)=Va(N_cycle);
 Qp(1)=Qp(N_cycle);
 Qlv(1)=Qlv(N_cycle);
 Vv(1)=Vv(N_cycle);
 end

 

%% Plot
if Plot_flag  % Plotting commands: … ; end
    
    % plotting maximal aortic pressure vs cycles 1.1
    pks =findpeaks(Pao2); 
    figure; stem(1:20,pks(2:end),'b','filled'); % 1st peak isn't relevant since it's not the maximum, so we go from 2 till end
    ylabel('Aortic Pressure[mmHg]');
    xlabel('N [cycles]');
    title('Maximal Aortic Pressure per Heart Cycle');
    
    % plotting 4 types of pressure vs time 1.2
    time=dt*1000*(28000:30500); % time of 19th cycle
    MinPeakDistance=0.4;
    figure; subplot(2,1,1); plot(time,Pao2(28000:30500),'b'); hold on; plot(time,Plv2(28000:30500),'g'); hold on; plot(time,Pa2(28000:30500),'r'); hold on; plot(time,Pv2(28000:30500),'m'); 
    findpeaks(Pao2(28000:30500),time,'MinPeakDistance',MinPeakDistance); 
    legend('Pao','Plv','Pa','Pv');
    ylabel('pressure[mmHg]');
    xlabel('time [msec]');
    title('4 Pressure Parameters for One Cycle');
    subplot(2,1,2); plot(time,Vlv2(28000:30500),'b'); 
    ylabel('Left Ventricle Volume [mL]');
    xlabel('time [msec]');
    title('Volume of Left Ventricle');
    
   %pressure and volume for each given time 
   time_e=dt*1000*(1:N_cycle);
   te= [0,30,130,260 ];
   Ve=[114.6,114.6,66.06,65.85];
   Pe=[9.515,59.73,80.25,5.03];
   El=zeros(1,4);
   for i=1:4
      El(i)=(Pe(i))/(Ve(i)-V0);
   end
   
   figure;plot(te,El);hold on; plot(time_e,E);
   legend('elasticity according PV loop', 'original elasticity');
   xlabel('time[msec]');ylabel('elasticity[mmHg/mL]');
   figure;plot(te,El);
   xlabel('time[msec]');ylabel('elasticity[mmHg/mL]');



end

%part2


Cv = 300.0;
Rp = 1.0;
Emax = 2.0 ;
HR = 80;
%Volume [ml]
Vlv(1) = 120; % left ventricle
Va(1) = 270; % arteries
Vv(1) = 2700; % veins

vec=0.5:0.05:2.5;

Cv_vec=vec.*Cv;
Rp_vec = vec.*Rp;
Emax_vec = vec.*Emax ;
HR_vec = vec.*HR;

Pao_Cv=zeros(size(Cv_vec));
Pao_Rp=zeros(size(Rp_vec));
Pao_Emax=zeros(size(Emax_vec));
Pao_HR=zeros(size(HR_vec));
Pao_BV=zeros(size(vec));

for i=1:size(vec,2)
Pao_Cv(i)=mean(Sensitivity(Cv_vec(i),Rp,Emax,HR));
Pao_Rp(i)=mean(Sensitivity(Cv,Rp_vec(i),Emax,HR));
Pao_Emax(i)=mean(Sensitivity(Cv,Rp,Emax_vec(i),HR));
Pao_HR(i)=mean(Sensitivity(Cv,Rp,Emax,HR_vec(i)));
Pao_BV(i)=mean(Sensitivity2(Cv,Rp,Emax,HR,vec(i)));

end

figure;
plot (vec,Pao_Cv,'m');
hold on
plot (vec,Pao_Rp,'r');
hold on
plot (vec,Pao_Emax,'g');
hold on
plot(vec,Pao_HR,'b');
hold on
plot(vec,Pao_BV,'k');





