function [Kgerg,Zgerg,zc]=FGerg91(p,T,xa,xy,Roc);
% Calculation of the compressibility coefficient by the
% modified equation GERG-91 mod. (GOST 30319.2)
% T = 250 - 340 K, P = 0.1 - 12 MPa
%
% Format: [Kgerg,Zgerg,zc]=fgerg91(p,T,xa,xy,Roc);
%
% Kgerg - compressibility coefficient (Kgerg=Zgerg/zc);
% Zgerg - compressibility factor;
% zc - compressibility factor for standard conditions;
% p - absolute pressure, MPa;
% T - absolute temperature, K;
% xa - molar fraction of Azot;
% xy - molar fraction of CO2;
% Roc - density for standard conditions, kg/m.kub.
pc=0.101325; Tc=293.15;
xe=1-xa-xy;
zc=1-(.0741*Roc-.006-.063*xa-.0575*xy)^2;
Me=(24.05525*zc*Roc-28.0135*xa-44.01*xy)/xe;
H=128.64+47.479*Me;
B1=-.425468+2.865e-3*T-4.62073e-6*(T^2)+(8.77118e-4-5.56281e-6*T+8.81514e9*(T^2))*H+(-8.24747e-7+4.31436e-9*T-6.08319e-12*T^2)*H^2;
B2=-.1446+7.4091e-4*T-9.1195e-7*(T^2);
B23=-.339693+1.61176e-3*T-2.04429e-6*(T^2);
B3=-.86834+4.0376e-3*T-5.1657e-6*(T^2);
C1=-.302488+1.95861e-3*T-3.16302e-6*(T^2)+(6.46422e-4-4.22876e-6*T+6.88157e9*(T^2))*H+(-3.32805e-7+2.2316e-9*T-3.67713e-12*(T^2))*H^2;
C2=7.8498e-3-3.9895e-5*T+6.1187e-8*(T^2);
C3=2.0513e-3+3.4888e-5*T-8.3703e-8*(T^2);
C223=5.52066e-3-1.68609e-5*T+1.57169e-8*(T^2);
C233=3.58783e-3+8.06674e-6*T-3.25798e-8*(T^2);
B_=.72+1.875e-5*((320-T)^2);
C_=.92+.0013*(T-270);
Bm=xe^2*B1+xe*xa*B_*(B1+B2)-1.73*xe*xy*(B1*B3)^.5+xa^2*B2+2*xa*xy*B23+xy^2*B3;
Cm=xe^3*C1+3*xe^2*xa*C_*(C1^2*C2)^(1/3)+2.76*xe^2*xy*(C1^2*C3)^(1/3)+3*xe*xa^2*C_*(C1*C2^2)^(1/3)+6.6*xe*xa*xy*(C1*C2*C3)^(1/3)+2.76*xe*xy^2*(C1*C3^2)^(1/3)+xa^3*C2+3*xa^2*xy*C223+3*xa*xy^2*C233+xy^3*C3;
b=1000*p/(2.7715*T);
C0=b^2*Cm;
B0=b*Bm;
A1=1+B0;
A0=1+1.5*(B0+C0);
A22=A0-(A0^2-A1^3)^.5;
znak=1;
if A22<0
 znak=-1;
end
A2=znak*(abs(A22))^(1/3);
Zgerg=(1+A2+A1/A2)/3;
zc=1-(0.0741*Roc-0.006-0.063*xa-0.0575*xy)^2;
Kgerg=Zgerg/zc;