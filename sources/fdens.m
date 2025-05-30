function [DRo,Rom,z,A1,A2,A3,A4]=fdens(p,T,Rom);
% Auxiliary function for calculation the compressibility factor and molar
% density by the equation VNIC SMV.
% Format: [DRo,Rom,z,A1,A2,A3,A4]=fdens(p,T,Rom)
% z - compressibility factor;
% A1,A2,A3 - parameters for calculation the isoentrope factor (Kapa);
% A4 - A4=dz/dTp, derivative for calculation drossel coefficient;
% p - absolute pressure, MPa;
% T - temperature, K;
% Rom - molar density, kmol/m.kub.
global Ropk Tpk c R
 r=10; % Parameters
 S=[7 6 6 5 5 4 3 3 3 2]; % for state equation
Rop=Rom/Ropk;
 Tp=T/Tpk;
 z=1; A1=0; A2=0; A3=0; A4=0;
 for k=1:r
 for l=0:S(k)
 z=z+c(k,l+1)*Rop^k/Tp^l;
 A1=A1+(k+1)*c(k,l+1)*Rop^k/Tp^l;
 A2=A2-(l-1)*c(k,l+1)*Rop^k/Tp^l;
 A3=A3-l*(l-1)/k*c(k,l+1)*Rop^k/Tp^l;
 A4=A4-l*c(k,l+1)*Rop^k/Tp^(l+1)*1/Tpk;
 end
 end
 DRo=(1e3*p-R*T*z*Rom)/(R*T*(1+A1));
 Rom=Rom+DRo;