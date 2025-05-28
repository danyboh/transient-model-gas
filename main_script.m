% === ЄДИНИЙ СКРИПТ ДЛЯ МОДЕЛЮВАННЯ НЕСТАЦІОНАРНОГО РУХУ ГАЗУ ===

clc; clear;

% === Параметри моделі ===
L = 1000;      % довжина трубопроводу [м]
Nx = 1000;     % кількість комірок
T_end = 5;     % час моделювання [с]
CFL = 0.5;     % число Куранта
x = linspace(0, L, Nx); % координатна сітка

% === Початкові умови ===
p0 = 5.5e6;        % тиск на вході [Па]
pL = 5.0e6;        % тиск на виході [Па]
T0 = 288.15;       % температура [K]
R = 518.3;         % газова стала для природного газу
x_mol = [0.94 0.02 0.01 0.005 0.005 0.01 0.008 0.002];  % склад газу

% Обчислення Z, Mm, Cv
A = 0.005; B = 0.001;
Z = 1 + A * (p0 / 1e6) - B * T0;
M = [16.043, 30.07, 44.097, 58.12, 58.12, 28.01, 44.01, 34.08];
Mm = sum(x_mol .* M);
Cv = 1.5 * 8.31451 / Mm * 1000;  % [Дж/(кг*К)]

% Розподіли по трубопроводу
p_profile = linspace(p0, pL, Nx);
rho_profile = p_profile ./ (Z * R * T0);
u_profile = linspace(1, 3, Nx);
e_profile = Cv * T0 + 0.5 * u_profile.^2;

U = zeros(3, Nx);
U(1, :) = rho_profile;
U(2, :) = rho_profile .* u_profile;
U(3, :) = rho_profile .* e_profile;

% === Числове розв'язання (SSPRK2) ===
dx = x(2) - x(1);
t = 0;
while t < T_end
    rho = U(1, :);
    u = U(2, :) ./ rho;
    e = U(3, :) ./ rho;
    T = (e - 0.5 * u.^2) / Cv;
    p = rho .* R .* T;
    c = sqrt(1.3 * p ./ rho);
    umax = max(abs(u) + c);
    dt = CFL * dx / umax;
    if t + dt > T_end
        dt = T_end - t;
    end

    F = [rho .* u;
         rho .* u.^2 + p;
         (U(3, :)+p) .* u];
    F_plus = [F(:, 2:end), F(:, end)];
    F_minus = [F(:, 1), F(:, 1:end-1)];
    U_plus = [U(:, 2:end), U(:, end)];
    U_minus = [U(:, 1), U(:, 1:end-1)];
    lambda = max(abs(u) + c);
    Flux = 0.5 * (F_plus + F) - 0.5 * lambda * (U_plus - U);
    RHS1 = -(Flux - [Flux(:,1), Flux(:,1:end-1)]) / dx;

    U1 = U + dt * RHS1;

    rho = U1(1, :);
    u = U1(2, :) ./ rho;
    e = U1(3, :) ./ rho;
    T = (e - 0.5 * u.^2) / Cv;
    p = rho .* R .* T;
    F = [rho .* u;
         rho .* u.^2 + p;
         (U1(3, :)+p) .* u];
    F_plus = [F(:, 2:end), F(:, end)];
    F_minus = [F(:, 1), F(:, 1:end-1)];
    U_plus = [U1(:, 2:end), U1(:, end)];
    U_minus = [U1(:, 1), U1(:, 1:end-1)];
    lambda = max(abs(u) + c);
    Flux = 0.5 * (F_plus + F) - 0.5 * lambda * (U_plus - U1);
    RHS2 = -(Flux - [Flux(:,1), Flux(:,1:end-1)]) / dx;

    U = 0.5 * (U + U1 + dt * RHS2);
    t = t + dt;
end

% === Побудова графіків ===
rho = U(1,:);
u = U(2,:) ./ rho;
e = U(3,:) ./ rho;
T = (e - 0.5 * u.^2) / Cv;

invalid_T = T <= 0 | isnan(T) | ~isreal(T);
T(invalid_T) = 273.15;
p = rho .* R .* T;
p(~isfinite(p)) = 1e5;

T_C = T - 273.15;  % Температура в Цельсіях

positions = [100 600 560 420;
             700 600 560 420;
             100 100 560 420;
             700 100 560 420];

f1 = figure; set(f1, 'Position', positions(1,:));
plot(x, p / 1e5); ylabel('Тиск [бар]'); xlabel('Відстань [м]'); grid on; title('Розподіл тиску');

f2 = figure; set(f2, 'Position', positions(2,:));
plot(x, T_C); ylabel('Температура [°C]'); xlabel('Відстань [м]'); grid on; title('Розподіл температури');

f3 = figure; set(f3, 'Position', positions(3,:));
plot(x, u); ylabel('Швидкість [м/с]'); xlabel('Відстань [м]'); grid on; title('Розподіл швидкості');

f4 = figure; set(f4, 'Position', positions(4,:));
plot(x, rho); ylabel('Густина [кг/м^3]'); xlabel('Відстань [м]'); grid on; title('Розподіл густини');

valid = isfinite(p) & isfinite(T) & T > 0;
K_c = trapz(x(valid), p(valid) ./ (Z * R * T(valid)));
fprintf('Запас газу в трубопроводі (K_c): %.3f\n', K_c);