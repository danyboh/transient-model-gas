function Mjt = met_nulp (p,T,xa,xy,Ro_st)
% Format: [Mjt] = met_nulp (p,T,xa,xy,Ro_st);
% Mjt - Joule-Thomson coefficient, K/MPa;
% p - absolute pressure, MPa;
% T - absolute temperature, K;
% xa - molar fraction of Azot;
% xy - molar fraction of CO2;
% Ro_st - density for standard conditions, kg/m.kub.
Pst = 0.101325; % MPa
Tst = 293.15; % K
% Знаходим псевдокритичну температуру Tpk
Tpk = 88.25*(0.9915 + 1.759*Ro_st - xy - 1.681*xa);
% Знаходим псевдокритичну густину Ropk
Ropk = 163.5*(Ro_st/0.6682)^0.6 + 62.62*xa + 163.359*xy;
tau = T/Tpk;
K = FGerg91(p,T,xa,xy,Ro_st);
Ro = p/Pst*Tst/T*Ro_st/K;
Om = Ro/Ropk;
% Cклад газу, за яким проводилася апроксимація
% x =[ 94.9 2.6 0.4 0.1 0.1 1.2 0.6 0.1];
% Ro_st = 0.70536
a=[-67.802 446.76 -1102.4 1209.1 -498.87
 143.8 -965.14 2433.4 -2736.8 1162.1
 -113.04 775.41 -2001.9 2307.9 -1003.15
 22.27 -158.22 420 -490.79 209.56
 14.178 -95.204 241.87 -276.03 118.25
 2.7783 -20.996 62.704 -90.622 56.442 ];
for i = 1 : 6,
 K_r(i) = polyval(a(i,:),tau);
end
Mjt = polyval(K_r,Om); 