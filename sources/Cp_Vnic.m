function [Cp_kg,Ro,Mm] = Cp_Vnic(p,T,x);
% Calculation of the isobaric
% by the equation VNIC SMV (GOST 30319.3).
% T = 240 - 480 K, P = 0.1 - 12 MPa
% Auxiliery files: DAT_VNIC, FVNIC, FDENS, DAT_CP0, CALKCP0.
% Format: [Cp_kg,Ro,Mm] = Cp_Vnic(p,T,x);
% p - absolute pressure, MPa;
% T - temperature, K;
% x - vector-row the mole fraction of 18(AGA-orientation),
% or 8(RD & VNIC-orientation) components.
% x = [ Metan Etan Propan n-Butan i-Butan Azot CO2 H2S ]
% Cp_kg - heat capacity, kDg/(kg*K);
% Ro - density for work conditions, kg/m3
% Mm - molar mass, kg/kmol
 R=8.31451; % kDg/(kMol*K)
 [Kvnic,z,A1,A2,A3,A4,Mm,Rom]=Fvnic(p,T,x);
 Ro = Mm * Rom ; % kg/m3
 Cpom=calkcpo(x,T); % kDg/(Kmol'*K)
 Cv=R*(Cvom/R+A3);
 Cp=R*(Cv/R+(1+A2)^2/(1+A1));
 Cp_kg = Cp/Mm;
