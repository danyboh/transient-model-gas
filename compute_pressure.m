function p = compute_pressure(U, T)
% Обчислення тиску з використанням p = Z*rho*R*T

R = 518.3;
Z = 1; % у спрощеному вигляді
rho = U(1,:);
p = Z .* rho .* R .* T;

end
