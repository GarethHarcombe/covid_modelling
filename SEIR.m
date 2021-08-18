function SEIR

% test a range of R0 values
R0 = [1.1:0.5:3.0];

% Number of hospital beds
HospCapacity = 1000;

% Days infectious for
InfDays = 10;

% infection fataility rate
IFR = 0.01;

% Hospitalisation rate
HospRate = 0.2;

% Time span
T = [0:2000];

% Run for a few values of R0
for I = 1:length(R0)
    % initial susceptible population
    S0 = 5e6;
    
    % key parameter beta
    beta = R0(I)/InfDays;
    
    % initial condition 100 infectious, no recovereds
    y0 = [1,0,100/S0,0];
    
    % Run the simulation for a very long time and return values of SIR each day
    % Sometimes a stiff system so use 23t
    % y = [S, E, I, R];
    opt = odeset('Events', @susceptableNegative);
    [t,y] = ode23t(@(t,y)RhsSEIR(t,y,beta, (5e6/365)/5e6.*(t>=365)), T, y0, opt);
    
    % When was the peak infection day and how big was the peak
    [InfPeak(I), PeakDay(I)] = max(y(:,3));
    
    % How many new cases each day (S(t) - S(t-1))
    NewCases(:,I) = (y(1:end-1,1) - y(2:end,1))*S0;
    
    
    TotalSusceptable(:,I) = y(1:end,1)*S0;
    TotalExposed(:,I) = y(1:end,2)*S0;
    % How many active cases are there each day
    TotalInfected(:,I) = y(1:end,3)*S0;
    TotalRecovered(:,I) = y(1:end,4)*S0;
    
    % How many people were infected when it's all over
    FinalInfected(I,1) = y(end,4)*S0;
    
    % How many people died?
    FinalDeaths(I,1) = y(end,4)*S0*IFR;
    
    % How many people in hospital each day?
    InHospital(:,I) = y(1:end,3)*S0*HospRate;
    
    % How many days was the hospital system overwhelmed?
    OverHospCap(I,1) = sum(InHospital(:,I)>HospCapacity);
end

HerdImmunity = S0*(1-1./R0);

% table(R0, FinalDeaths, FinalInfected, OverHospCap)
% 
% 
% figure(1)
% subplot(2,1,1)
% plot(T(2:end),NewCases)
% ylabel('Daily New Cases')
% xlabel('Time (days)')
% xlim([0 400])
% 
% subplot(2,1,2)
% plot(TotalInfected)
% ylabel('Infected')
% xlabel('Time (days)')
% xlim([0 400])
% legend(num2str(R0))
% 
% figure(2)
% subplot(3,1,1)
% plot(R0,FinalInfected,R0,HerdImmunity)
% ylabel('Total infected')
% xlabel('R_0')
% 
% subplot(3,1,2)
% plot(R0,PeakDay)
% ylabel('Peak day')
% xlabel('R_0')
% 
% subplot(3,1,3)
% plot(R0,OverHospCap)
% ylabel('Days Health services overwhelmed')
% xlabel('R_0')
% 
% figure(3)
% xrange = [0 1000];
% subplot(4,1,1)
% plot(TotalSusceptable)
% ylabel('Susceptable')
% xlabel('Time (days)')
% xlim(xrange)
% legend(num2str(R0))
% 
% subplot(4,1,2)
% plot(TotalExposed)
% ylabel('Exposed')
% xlabel('Time (days)')
% xlim(xrange)
% 
% subplot(4,1,3)
% plot(TotalInfected)
% ylabel('Infected')
% xlabel('Time (days)')
% xlim(xrange)
% 
% subplot(4,1,4)
% plot(TotalRecovered)
% ylabel('Recovered')
% xlabel('Time (days)')
% xlim(xrange)

figure(4)
S0 = 5e6;
lockdown = 0.9;
normal = 1.3;
y0 = [1,0,10/S0,0];
vaccineIntroduced = 280;

opt = odeset('Events', @susceptableNegative);
[t,y1] = ode23t(@(t,y)RhsSEIR(t,y,(normal/InfDays) + (lockdown/InfDays-1.2/InfDays).*(t>=30) + (normal/InfDays-0.9/InfDays).*(t>=200), 0), T, y0, opt);
[t,y2] = ode23t(@(t,y)RhsSEIR(t,y,(normal/InfDays) + (lockdown/InfDays-1.2/InfDays).*(t>=30) + (normal/InfDays-0.9/InfDays).*(t>=200), 1/850.*(t>=vaccineIntroduced)), T, y0, opt);
[t,y3] = ode23t(@(t,y)RhsSEIR(t,y,(normal/InfDays) + (lockdown/InfDays-1.2/InfDays).*(t>=30) + (normal/InfDays-0.9/InfDays).*(t>=200), 1/480.*(t>=vaccineIntroduced)), T, y0, opt);

sumCases_y1 = 0;
sumCases_y2 = 0;
sumCases_y3 = 0;
for i = T + 1
   if i > 1
       sumCases_y1 = sumCases_y1 + (1/4).*(y1(i,2)).*5e6;
       totalCases_y1(i) = sumCases_y1;

       sumCases_y2 = sumCases_y2 + (1/4).*(y2(i,2)).*5e6;
       totalCases_y2(i) = sumCases_y2;
       
       sumCases_y3 = sumCases_y3 + (1/4).*(y3(i,2)).*5e6;
       totalCases_y3(i) = sumCases_y3;
   end
end

plot(t, S0.*y1(:,3),t, S0.*y2(:,3),t, S0.*y3(:,3)) %, t, HospLine
xline(vaccineIntroduced, '--black', 'Vaccine Introduced')
ylabel('Active Cases')
xlabel('Time (days)')
legend('No Vaccination', 'Partial Vaccination', 'Full Vaccination') %, 'Hospital Capacity'

figure(7)
plot(t, HospRate.*S0.*y1(:,3),t, HospRate.*S0.*y2(:,3),t, HospRate.*S0.*y3(:,3)) %, t, HospLine
xline(vaccineIntroduced, '--black', 'Vaccine Introduced')
yline(HospCapacity, '--black', 'Hospital Capacity')
ylabel('Hospitalised Cases')
xlabel('Time (days)')
legend('No Vaccination', 'Partial Vaccination', 'Full Vaccination')

figure(5)
plot(t, totalCases_y2, t, totalCases_y3, t, totalCases_y1)
xline(vaccineIntroduced, '--black', 'Vaccine Introduced')
ylabel('Total Cases')
xlabel('Time (days)')
legend('No Vaccination', 'Partial Vaccination', 'Full Vaccination')

disp([num2str(totalCases_y1(end)*IFR),'  ', num2str(totalCases_y2(end)*IFR),'  ', num2str(totalCases_y3(end)*IFR)])
disp([num2str(totalCases_y1(end)),'  ', num2str(totalCases_y2(end)),'  ', num2str(totalCases_y3(end))])
disp([num2str(sum(y1(:,3).*S0>HospCapacity/HospRate)),'  ', num2str(sum(y2(:,3).*S0>HospCapacity/HospRate)),'  ', num2str(sum(y3(:,3).*S0>HospCapacity/HospRate))])

figure(6)
T = [0:280];
[t,y1] = ode23t(@(t,y)RhsSEIR(t,y,(normal/InfDays) + (lockdown/InfDays-1.2/InfDays).*(t>=30) + (normal/InfDays-0.9/InfDays).*(t>=200), 0), T, y0, opt);
[t,y2] = ode23t(@(t,y)RhsSEIR(t,y,(normal/InfDays) + (lockdown/InfDays-1.2/InfDays).*(t>=30) + (normal/InfDays-0.9/InfDays).*(t>=200), 1/850.*(t>=vaccineIntroduced)), T, y0, opt);
[t,y3] = ode23t(@(t,y)RhsSEIR(t,y,(normal/InfDays) + (lockdown/InfDays-1.2/InfDays).*(t>=30) + (normal/InfDays-0.9/InfDays).*(t>=200), 1/480.*(t>=vaccineIntroduced)), T, y0, opt);

plot(t, S0.*y1(:,3),t, S0.*y2(:,3),t, S0.*y3(:,3))
xline(vaccineIntroduced, '--black', 'Vaccine Released')
xline(200, '--black', 'Lockdown Eased')
ylabel('Total Cases')
xlabel('Time (days)')
end


function dy = RhsSEIR(t, y, beta, f)

% expected time to no longer be infectious is 10 days
gamma = 1/10; % 1/3.6 for SEIR
% assuming birth rate and death rate are equal, birth rate is 18.5 per 1000

sigma = 1/4;
a = 1/10;
v = 0.6;
omega = 1/365;


S = y(1);
E = y(2);
I = y(3);
R = y(4);

dy(1,1) = -beta.*S*I - v*f + omega*R;
dy(2,1) = beta*S*I - sigma*E;
dy(3,1) = sigma*E-a*I;
dy(4,1) = a*I + v*f - omega*R;

end

function [value, isterminal, direction] = susceptableNegative(t, y)
value = y(1);
isterminal = 1;   % Stop the integration
direction = 0;
end
