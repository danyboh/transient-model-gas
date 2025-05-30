function [Kvnic,z,A1,A2,A3,A4,Mm,Rom_r]=Fvnic(p,T,x);
% Calculation of the compressibility coefficient by the equation VNIC SMV
% (GOST 30319.2). T = 250 - 340 K, P= 0.1 - 12 MPa.
% Auxiliary files: fdens.m, dat_vnic.m
% Format: [Kvnic,z,A1,A2,A3,A4]=Fvnic(p,T,x)
% Kvnic - compressibility coefficient;
% z - compressibility factor;
% A1,A2,A3 - parameters for calculation the isoentrope factor (Kapa);
% A4 - A4=dz/dTp, derivative for calculation drossel coefficient;
% Mm - molar mass, kg/kmol;
% Rom_r - molar density, kmol/m.kub.
% p - absolute pressure, MPa;
% T - temperature, K;
% x - vector-row the mole fraction of 18(AGA-orientation),
% or 20(KievTDC-orientation),
% or 8(RD & VNIC-orientation) components.
% x = [ Metan Etan Propan n-Butan i-Butan Azot CO2 H2S ]
if length(x) == 20
 Xaga=zeros(1,18);
 Xaga(1:5)=x(1:5); Xaga(6:7)=x(11:12);Xaga(8)=x(15);Xaga(9:13)=x(6:10);
 Xaga(13)=Xaga(13)+x(19)+x(20); Xaga(10)=Xaga(10)+x(18);
 Xaga(14)=x(16)+x(13); Xaga(15)=x(14); Xaga(16)=0; Xaga(17)=x(17); Xaga(18)=0;
  x=Xaga;
 %sum(x), pause
end
if length(x)==18
 Xvnic=zeros(1,8);
 Xvnic=x(1:8); Xvnic(4)=Xvnic(4)+sum(x(9:13)); % Correspondly to
 Xvnic(6)=Xvnic(6)+sum(x(14:18)); % GOST 30319.2-96
 x=Xvnic;
 %sum(x), pause
end
% 1. DATA CHAPTER
 global Ropk Tpk c R % For function Fdens.m
dat_vnic % Load the datas from files Dat_vnic.m to workspase
 % 2. CALCULATION CHAPTER
 % 2.1. Definition the composition and state of the natural gaz
 n=8; % Number the components of natural gaz
 % 2.2. Calculation.
 R=8.31451; % kDg/(kMol*K)
 Oij=zeros(n);
 Vk=zeros(n);
 Tkij=zeros(n);
 for i=1:n
 for j=1:n
 Oij(i,j)=(O(i)*M(i)/Rok(i)+O(j)*M(j)/Rok(j))/(M(i)/Rok(i)+M(j)/Rok(j));
 Vk(i,j)=(1-lambda(i,j))*(((M(i)/Rok(i))^(1/3)+(M(j)/Rok(j))^(1/3))/2)^3;
 Tkij(i,j)=(1-ksi(i,j))*(Tk(i)*Tk(j))^0.5;
 end
 end
 Xi=[x;x;x;x;x;x;x;x];
 Xj=Xi';
 Ropk=1/sum(sum(Xi.*Xj.*Vk));
 Om=Ropk*sum(sum(Xi.*Xj.*Vk.*Oij));
 Tkm=sum(sum(Xi.*Xj.*Vk.*Tkij.^2));
 Tpk=(Tkm*Ropk)^0.5;

 Ppk=1e-3*R*Ropk*Tpk*(0.28707-0.05559*Om);
 % -----------------------------------------------------------------------
 % Робочі умови
 Pp=p/Ppk;
 Rom=9e3*p/(R*T*(1.1*Pp+0.7));
 c=Akl+Bkl*Om;
 [DRo,Rom,z]=fdens(p,T,Rom); iter=1;
 iter = 1;
max_iter = 50; 
 while abs(DRo/Rom)>=1e-6 && iter < max_iter
 [DRo,Rom,z,A1,A2,A3,A4]=fdens(p,T,Rom);
 iter=iter+1;
 end
 if iter > max_iter
    warning('Fvnic did not converge after %d iterations at p = %.2f MPa, T = %.2f K', max_iter, p, T);
end
 Rom_r = Rom;
 % Стандартні умови
 Pc=0.101325; Tc=293.15;
 Pp=Pc/Ppk;
 Rom=9e3*Pc/(R*Tc*(1.1*Pp+0.7));
 c=Akl+Bkl*Om;
max_iter = 50;
 [DRo,Rom,zc]=fdens(Pc,Tc,Rom); iter=1;
 while abs(DRo/Rom)>=1e-6 && iter < max_iter
 [DRo,Rom,zc]=fdens(Pc,Tc,Rom);
 iter=iter+1;
 end
 if iter == max_iter
    warning('Fvnic: No convergence at standard conditions (Pc=%.3f MPa, Tc=%.1f K)', Pc, Tc);
end
 %Kc=[1 2 3 4 4]; % Number of Carbon-atoms in moleculs of components (C<k>H<2k+2>)
 %zc=1-(0.0458*sum(Kc.*x(1:5))-0.0022+0.0195*x(6)+0.075*x(7))^2; % GOST 30319.1-96
 Kvnic=z/zc;
 Mm=sum(x.*M);
