clc
clf
clear

global FA0 pA pB k1 Ea Kc2 T1 T2 Rgas T0 sum_theta_cp delta_Hrx 

% ideal gas constant
Rgas = 8.314;


FA0     = 965/3.6;     % CO
FB0     = 8021/3.6;    % H2O
FC0     = 3383/3.6;    % CO2
FD0     = 13633/3.6;   % H2
FME0    = 72/3.6;      % CH4
FNI0    = 4714/3.6;    % N2
FAR0    = 56/3.6;      % Ar

FT0 = FA0+FB0+FC0+FD0+FME0+FNI0+FAR0;
f0=[FA0;FB0;FC0;FD0;FME0;FNI0;FAR0];


T0 = 450;         %K
P0 = 3700000;     %Pa

X=0;
% FA= FA0*(1-X);
% FB= FB0-FA0*(1-X);
% FC= FC0+FA0*X;
% FD= FD0+FA0*X;

pA =((1-X)*FA0)/FT0*P0;
pB =(FB0-X*FA0)/FT0*P0;

% Define parameters for rate constant
k1 = 25.8E3;
Ea = 16E3;
T1= 484;

% Define parameters for equilibrium constant
Kc2 = 1.767*10^-2;
T2 = 484; %???

% Define heat of reaction
delta_Hrx = -41100; %Jmol-1


% Define heat capacities (constant pressure)
% CP_A = 6.60*T0+0.0006*T0^2;               % CO
% CP_B = 8.22*T0+0.000075*T0^2+0.000000447*T0^3; % H2O
% CP_C = 10.34*T0+0.00137*T0^2+195500/T0;   % CO2
% CP_D = 6.62*T0+0.000405*T0^2;               % H2
% CP_ME = 5.34*T0+0.00575*T0^2 ;               % CH4
% CP_NI = 6.5*T0+0.0005*T0^2;              % N2
% CP_AR = 4.97*T0;                       % Ar

% CP_A = 6.60+0.00120*T0;               % CO
% CP_B = 8.22+0.00015*T0+0.00000134*T0^2; % H2O
% CP_C = 10.34+0.00274*T0-195500/T0^2;   % CO2
% CP_D = 6.62+0.00081*T0;               % H2
% CP_ME = 5.34+0.0115*T0 ;               % CH4
% CP_NI = 6.50+0.00100*T0;              % N2
% CP_AR = 4.97;                       % Ar

CP_A = 6.60+0.00120*T0;               % CO
CP_B = 8.22+0.00015*T0+0.00000134*T0^2; % H2O
CP_C = 10.34+0.00274*T0-195500/T0^2;   % CO2
CP_D = 6.62+0.00081*T0;               % H2
CP_ME = 5.34+0.0115*T0 ;               % CH4
CP_NI = 6.50+0.00100*T0;              % N2
CP_AR = 4.97;                       % Ar

theta_A = FA0/FT0;
theta_B = FB0/FT0;
theta_C = FC0/FT0;
theta_D = FD0/FT0;
theta_ME = FME0/FT0;
theta_NI = FNI0/FT0;
theta_AR = FAR0/FT0;

sum_theta_cp=(theta_A*CP_A+theta_B*CP_B+theta_C*CP_C+theta_D*CP_D+theta_ME*CP_ME+theta_NI*CP_NI+theta_AR*CP_AR)*4.18;

% Define the initial conditions
X0 = 0;
% Define the volume for which the approximate solution should be
% determined
W0 = 0;
W_max = 5000;
W_span = [W0,W_max];
% Now solve the ODE
[W,X] = ode45('WGS_ode',W_span,X0);

% Re-calculate the temperature now that we know the conversion
T = T0 - delta_Hrx*X/sum_theta_cp;
% Re-calculate the rate of reaction 
k = k1 * exp(-Ea/Rgas*(1/T1-1./T));
% equilibrium constant
Kc = Kc2*exp(delta_Hrx/Rgas*(1/T2 - 1./T));
% Kc = Kc2*exp(4400*(1/T2 - 1./T));

% rate of reaction
rate = (k*pA*pB)/(1+127*pA*pB+26*pA);

% Calculate equilibrium conversion
Xequi = Kc./(1+Kc);

% plot the solution
% figure(1)
plot(W,X, W,Xequi,'r-.')
xlabel('W in g')
ylabel('X')
title('Adiabatic Water-Gas Shift, conversion')
legend('X', 'X_{equilibrium}')
% figure(2)
% plot(W,T)
% xlabel('V in m^3')
% ylabel('T in K')
% title('Adiabatic isomerisation of n-butane, temperature')

% figure(3)
% plot(W,rate)
% xlabel('V in m^3')
% ylabel('rate in mol / (m^3 h)')
% title('Adiabatic isomerisation of n-butane, rate of reaction')
