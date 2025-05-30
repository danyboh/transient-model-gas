function yp=Rozrax_Cp(u,y)
% x = [ Metan Etan Propan n-Butan i-Butan Azot CO2 H2S ]
 x = [ 93.635 3.075 0.881 0.141 0.170 1.181 0.917 0]/100;
% Ro_st = 0.7186 >граф.
Qc=25000/3600; % m^3/s
Roc = 0.7186;
qm=Qc*Roc; % kg/s
g=9.81;
dely=0; % m
L=10000;
Dv=0.4; %m
Dz=0.43; %m
p = y(1); % Pa
pm = p/1e6; % MPa
T = y(2); % K
[Cp_kg,Ro,Mm] = Cp_Vnic(pm,T,x); % [ kDg/(kg*K), kg/m3, kg/kmol ]
Cp_kg = Cp_kg*1000; % Dg/(kg*K)
xa = x(6)/100; xy = x(7)/100;
Kjt = met_nulp (pm,T,xa,xy,Roc); %Koeficient Dgoulya Tomsona, K/MPa
Kjt = Kjt / 1e6; % K/Pa
Mj=VisG1(pm,T,xa,xy,Roc); % mkPa*s
[Kgerg,z,zc]=FGerg91(pm,T,xa,xy,Roc);
F=(pi*Dv.^2)/4; % Ploshca poperechnogo peretinu truboprovodu, m2
R=8314.472; % Gazova stala, Dg/(kMol*K)
kt=1.75; %koeficient teploperedachi, Vt/(m2*K)
Re = 10000;
lyam = 0.316/(Re).^(1/4); %koeficient gidravlichnogo opory
Tg = 273.15+14; %Temperatura grunty
yp=[ -((Mm*g*dely*p)/(z*T*R*L)+(qm.^2*lyam*z*R*T)/(2*Dv*F.^2*Mm*p)), -(kt*pi*Dz)/(qm*Cp_kg)*(T-Tg) - (g*dely)/(L*qm*Cp_kg) - Kjt*((Mm*g*dely*p)/(z*T*R*L)+(qm.^2*lyam*z*R*T)/(2*Dv*F.^2*Mm*p))];
u0=0; uf=10000;
y0=[3e6 293.15];
[u,y]=ode45('Rozrax_Cp',u0,uf,y0);
figure (1); plot(u,y(:,1));grid;title('Po Tisky');
figure (2); plot(u,y(:,2));grid;title('Po Temperature');