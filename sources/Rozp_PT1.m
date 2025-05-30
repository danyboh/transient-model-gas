u0=0; uf=10000;y0=[3e6 293.15];
[u,y]=ode45('Rozrax_Cp',u0,uf,y0);
figure (1); plot(u,y(:,1)); grid; title('Po Tisky'); figure (2); plot(u,y(:,2));grid;
title('Po Temperature'); 