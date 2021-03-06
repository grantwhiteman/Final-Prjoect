close all; clear all; clc;

T = linspace(500,1000); % degrees K
t = T/1000;
% Setup equations for each species
% First we enter in the parameters and compute the enthalpy and entropy for each species.

%Hydrogen
%link

% T = 298-1000K valid temperature range
A =  33.066178;
B = -11.363417;
C =  11.432816;
D = -2.772874;
E = -0.158558;
F = -9.980797;
G =  172.707974;
H =  0.0;

Hf_29815_H2 = 0.0; % kJ/mol
S_29815_H2 = 130.68; % J/mol/K

dH_H2 = A*t + B*t.^2/2 + C*t.^3/3 + D*t.^4/4 - E./t + F - H;
S_H2 = (A*log(t) + B*t + C*t.^2/2 + D*t.^3/3 - E./(2*t.^2) + G);
%H2O
% link Note these parameters limit the temperature range we can examine, as these parameters are not valid below 500K. There is another set of parameters for lower temperatures, but we do not consider them here.

% 500-1700 K valid temperature range
A =   30.09200;
B =   6.832514;
C =   6.793435;
D =  -2.534480;
E =   0.082139;
F =  -250.8810;
G =   223.3967;
H =  -241.8264;

Hf_29815_H2O = -241.83; %this is Hf.
S_29815_H2O = 188.84;

dH_H2O = A*t + B*t.^2/2 + C*t.^3/3 + D*t.^4/4 - E./t + F - H;
S_H2O = (A*log(t) + B*t + C*t.^2/2 + D*t.^3/3 - E./(2*t.^2) + G);
%CO
% link

% 298. - 1300K valid temperature range
A =   25.56759;
B =   6.096130;
C =   4.054656;
D =  -2.671301;
E =   0.131021;
F =  -118.0089;
G =   227.3665;
H = -110.5271;

Hf_29815_CO = -110.53; %this is Hf kJ/mol.
S_29815_CO = 197.66;

dH_CO = A*t + B*t.^2/2 + C*t.^3/3 + D*t.^4/4 - E./t + F - H;
S_CO = (A*log(t) + B*t + C*t.^2/2 + D*t.^3/3 - E./(2*t.^2) + G);
%CO2
%link

% 298. - 1200.K valid temperature range
A =   24.99735;
B =   55.18696;
C =  -33.69137;
D =   7.948387;
E =  -0.136638;
F =  -403.6075;
G =   228.2431;
H =  -393.5224;

Hf_29815_CO2 = -393.51; %this is Hf.
S_29815_CO2 = 213.79;

dH_CO2 = A*t + B*t.^2/2 + C*t.^3/3 + D*t.^4/4 - E./t + F - H;
S_CO2 = (A*log(t) + B*t + C*t.^2/2 + D*t.^3/3 - E./(2*t.^2) + G);
% Standard state heat of reaction
% We compute the enthalpy and free energy of reaction at 298.15 K for the following reaction $CO + H2O \rightleftharpoons H2 + CO2$.

Hrxn_29815 = Hf_29815_CO2 + Hf_29815_H2 - Hf_29815_CO - Hf_29815_H2O;
Srxn_29815 = S_29815_CO2 + S_29815_H2 - S_29815_CO - S_29815_H2O;
Grxn_29815 = Hrxn_29815 - 298.15*(Srxn_29815)/1000;

sprintf('deltaH = %1.2f',Hrxn_29815)
sprintf('deltaG = %1.2f',Grxn_29815)
%ans = deltaH = -41.15
%ans =deltaG = -28.62

%Correcting $\Delta H$ and $\Delta G$
%we have to correct for temperature change away from standard state. We only correct the enthalpy for this temperature change. The correction looks like this:

%$$ \Delta H_{rxn}(T) = \Delta H_{rxn}(T_{ref}) + \sum_i \nu_i
%(H_i(T)-H_i(T_{ref}))$$

%Where $\nu_i$ are the stoichiometric coefficients of each species, with appropriate sign for reactants and products, and $(H_i(T)-H_i(T_{ref})$ is precisely what is calculated for each species with the equations

Hrxn = Hrxn_29815 + dH_CO2 + dH_H2 - dH_CO - dH_H2O;
%The entropy is on an absolute scale, so we directly calculate entropy at each temperature. Recall that H is in kJ/mol and S is in J/mol/K, so we divide S by 1000 to make the units match.

Grxn = Hrxn - T.*(S_CO2 + S_H2 - S_CO - S_H2O)/1000;
% Plot how the $\Delta G$ varies with temperature
% over this temperature range the reaction is exothermic, although near 1000K it is just barely exothermic. At higher temperatures we expect the reaction to become endothermic.

figure; hold all
plot(T,Grxn)
plot(T,Hrxn)
xlabel('Temperature (K)')
ylabel('(kJ/mol)')
legend('\Delta G_{rxn}', '\Delta H_{rxn}', 'location','best')

R = 8.314e-3; %kJ/mol/K
K = exp(-Grxn/R./T);

figure
plot(T,K)
xlim([500 1000])
xlabel('Temperature (K)')
ylabel('Equilibrium constant')
%A = CO;
%B = H2O;
%C = H2;
%D = CO2;
Pa0 = 5; Pb0 = 5; Pc0 = 0; Pd0 = 0;  % pressure in atm
R = 0.082;
Temperature = 1000;

% we can estimate the equilibrium like this. We could also calculate it
% using the equations above, but we would have to evaluate each term. Above
% we simple computed a vector of enthalpies, entropies, etc...
K_Temperature = interp1(T,K,Temperature);

% If we let X be fractional conversion then we have $C_A = C_{A0}(1-X)$,
% $C_B = C_{B0}-C_{A0}X$, $C_C = C_{C0}+C_{A0}X$, and $C_D =
% C_{D0}+C_{A0}X$. We also have $K(T) = (C_C C_D)/(C_A C_B)$, which finally
% reduces to $0 = K(T) - Xeq^2/(1-Xeq)^2$ under these conditions.

f = @(X) K_Temperature - X^2/(1-X)^2;

Xeq = fzero(f,[1e-3 0.999]);
sprintf('The equilibrium conversion for these feed conditions is: %1.2f',Xeq)

P_CO = Pa0*(1-Xeq);
P_H2O = Pa0*(1-Xeq);
P_H2 = Pa0*Xeq;
P_CO2 = Pa0*Xeq;

K_Temperature=(P_CO2*P_H2)/(P_CO*P_H2O);