function dX_dW = WGS_ode(W,X)
% Reaction Engineering, Lecture 6
% Function the right hand side of the ODE and all explicit equations that
% need to be solved simutaneously for the adiabatic PFR problem discussed
% in the lecture.
% File written by Tina D?ren
% Last modified 9/11/2014
%
% components present in solution
% A - n-butane
% B - i-butane
% C - i-pentane
%
% Variable dictionary
% FA0         global   molar flow rate of n-butane at reator inlet (kmol/h)
% CA0         global   concentration of n-butane at reactor inlet (kmol/m3)
% T0          global   temperature at reactor inlet (K)
% k1          global   rate constant at temperature T1 (h-1)
% Ea          global   activation energy (J mol-1 K-1)
% T1          global   reference temperature for reaction constant k1 (K)
% Kc2         global   equilibrium constant at temperature T2 (-)
% T2          global   reference temperature for equilibrium constant Kc2 (K)
% sum_theta_cp global  sum of theta_i times CP (J/mol)
% delta_Hrx   global   heat of rection (J/mol)
% R_gas       global   universal gas constant
% V           input    volume / independent variable (m3)
% X           input    vector containing the conversion of A (-)
% dX_dV       output   vector containing the gradient 
% rA          local   rate of reaction of n-butane (kmol/m3/h)

global FA0 pA pB k1 Kc2 T1 T2 Ea Rgas T0 sum_theta_cp delta_Hrx 
% FA=O(1);
% FB=O(2);
% FC=O(3);
% FD=O(4);
% X=O(5);

% Temperature - equation
T = T0 - delta_Hrx*X/sum_theta_cp;
% rate constant - equation
k = k1 * exp(-Ea/Rgas*(1/T1-1./T));
% equilibrium constant - equation 
Kc = Kc2*exp(delta_Hrx/Rgas*(1/T2 - 1/T));
% Kc=Kc2*exp(4400/1/T2-1/T);
% rate of reaction - equation 
rA = ((k*pA*pB)/(1+127*pA*pB+26*pA))*(1-(1+1./Kc).*X);
% rA = (k*pA*pB)/(1+127*pA*pB+26*pA)
dX_dW = -rA /FA0;

% n=nomel(O);
% dO_dX=zeroes(n,1);

