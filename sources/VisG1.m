function Mj=VisG1(p,T,xa,xy,Roc);
% Calculation of the dynamic viscosity by metods
% GOST 30319.1-96. Delta=6% in respect of equation VNIC SMV.
% p=0.1 - 12 MPa
% Format: Mj=VisG1(p,T,xa,xy,Roc)
% Mj - dynamic viscosity, mkPa*s
% p - absolute pressure, MPa;
% T - temperature, K;
% xa - mole fraction of N2;
% xy - mole fraction of CO2;
% Roc - density (st.cond.), kg/m.kub.
Ppk=2.9585*(1.608-0.05994*Roc+xy-0.392*xa);
Tpk=88.25*(0.9915+1.759*Roc-xy-1.681*xa);
Pp=p/Ppk;
Tp=T/Tpk;
Cmj=1+Pp^2/30/(Tp-1);
Mjt=3.24*(T^0.5+1.37-9.09*Roc^0.125)/(Roc^0.5+2.08-1.5*(xa+xy)); % mkPa*s
Mj=Mjt*Cmj;