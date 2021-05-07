
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