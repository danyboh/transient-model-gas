function RHS = compute_rhs(U, p, dx)
% Обчислення правої частини для SSPRK2 (метод Rusanov)

Nx = size(U, 2);
RHS = zeros(size(U));

% Консервативні змінні
rho = U(1,:);
rhou = U(2,:);
rhoe = U(3,:);

u = rhou ./ rho;
e = rhoe ./ rho;

% Потоки
F1 = rhou;
F2 = rhou .* u + p;
F3 = (rhoe + p) .* u;

% Дискретизація методом Rusanov
for i = 2:Nx-1
    % ліві й праві потоки
    UL = U(:,i-1); UR = U(:,i);
    uL = UL(2)/UL(1); uR = UR(2)/UR(1);
    eL = UL(3)/UL(1); eR = UR(3)/UR(1);
    TL = (eL - 0.5*uL^2) / (1.5 * 518.3);
    TR = (eR - 0.5*uR^2) / (1.5 * 518.3);
    pL = compute_pressure(UL, TL);
    pR = compute_pressure(UR, TR);
    
    FL = [UL(2); UL(2)^2/UL(1) + pL; (UL(3) + pL)*uL];
    FR = [UR(2); UR(2)^2/UR(1) + pR; (UR(3) + pR)*uR];

    a = max(abs(uL) + sqrt(1.3*pL/UL(1)), abs(uR) + sqrt(1.3*pR/UR(1)));

    F_half = 0.5*(FL + FR) - 0.5*a*(UR - UL);

    RHS(:,i) = -(F_half - ([UL(2); UL(2)^2/UL(1) + pL; (UL(3)+pL)*uL])) / dx;
end

% Граничні умови: відображення
RHS(:,1) = RHS(:,2);
RHS(:,Nx) = RHS(:,Nx-1);

end
