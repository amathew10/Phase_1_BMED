%Akash Mathew
%ECE 4784
%Problem 1 Steady state neuron

clear all
close all

t_sim = 100; %simulation time in milliseconds
dt = 1/t_sim; %time step
t = 0:dt:t_sim; %time array 

I(1:length(t)) = 0; %injected current of 0 amps for 100 ms

%conductances
gbar_K = 36;
gbar_Na = 120;
gbar_L = .3;

%gate voltages
E_K = -12;
E_Na = 115;
E_L = 10.6;

%membrane rest voltage and capacitance
V_rest = -70;
C_m = 1;


V_m = 0; %intial voltage to calculate initial value of m,n,h

a_m = .1*( (25-V_m) / (exp((25-V_m)/10)-1)); %alpha_m
b_m = 4*exp(-V_m/18);%beta_m

a_n = .01 * ( (10-V_m) / (exp((10-V_m)/10)-1)); %alpha_n
b_n = .125*exp(-V_m/80); %beta_n

a_h = .07*exp(-V_m/20); %alpha_h
b_h = 1/(exp((30-V_m)/10)+1); %beta_h

%initial values of m,n and h
m = a_m/(a_m+b_m); %m_0 
n = a_n/(a_n+b_n); %n_0
h = a_h/(a_h+b_h); %h_0


%calculation at each time step
for i=1:length(t)-1
    %Gating variables
    a_m(i) = .1*( (25-V_m(i)) / (exp((25-V_m(i))/10)-1) );
    b_m(i) = 4*exp(-V_m(i)/18);
    
    a_n(i) = .01 * ( (10-V_m(i)) / (exp((10-V_m(i))/10)-1) );
    b_n(i) = .125*exp(-V_m(i)/80);
    
    a_h(i) = .07*exp(-V_m(i)/20);
    b_h(i) = 1/(exp((30-V_m(i))/10)+1);
    
    %Currents
    I_Na = (m(i)^3) * gbar_Na * h(i) * (V_m(i)-E_Na); 
    I_K = (n(i)^4) * gbar_K * (V_m(i)-E_K); 
    I_L = gbar_L *(V_m(i)-E_L); 
    I_ion = I(i) - I_K - I_Na - I_L; 
    
    V_m(i+1) = V_m(i) + dt*I_ion/C_m;
    m(i+1) = m(i) + dt*(a_m(i) *(1-m(i)) - b_m(i) * m(i));
    n(i+1) = n(i) + dt*(a_n(i) *(1-n(i)) - b_n(i) * n(i)); 
    h(i+1) = h(i) + dt*(a_h(i) *(1-h(i)) - b_h(i) * h(i));
end

V_m = V_m + V_rest;

figure
Vm = plot(t,V_m);
legend(Vm,'Voltage')
ylabel('Voltage (mV)')
xlabel('Time (ms)')
title('Membrane Potential')


figure
gK = plot(t,gbar_K*n.^4);
hold on
gNa = plot(t,gbar_Na*(m.^3).*h,'r');
legend([gK, gNa], 'gK', 'gNa')
ylabel('Conductance (mS/cm^2)')
xlabel('Time (ms)')
title('gK and gNa')










