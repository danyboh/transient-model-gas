% characteristic_method.m
% Транзієнтний рух газу в магістральному трубопроводі методом характеристик
% На основі [11] та попередніх специфікацій

close all; clear all; clc

%% 1. Параметри труби та газу
D     = 0.5;                   % діаметр труби, м
A     = pi*D^2/4;              % площа перерізу, м^2
L     = 1000;                  % довжина труби, м
R     = 518;                   % газова стала, Дж/(кг·K)
T     = 288;                   % температура, K
M     = 16e-3;                 % молярна маса, кг/моль
alpha = 0.02;                  % коефіцієнт Z(p), 1/МПа
beta  = 1.0;                   % зсув Z(p)
nu    = 1.2e-5;                % кінематична в’язкість, м^2/с
lambda_f = @(Re) 0.3164*Re.^(-0.25);  % формула Дарсі–Вайсбаха

%% 2. Просторова та часова сітка
Nx    = 101;                   
x     = linspace(0, L, Nx)';   % вектор x як стовпець
dx    = L/(Nx-1);

% обчислення dt за умовою CFL
p_max = 6;                     
c_max = sqrt(R*T*(alpha*p_max+beta)/M);
dt    = 0.5 * dx / c_max;      % безпечний крок за часом
t_end = 50;                    % кінцевий час моделювання, с
Nt    = ceil(t_end/dt);
t     = (0:Nt-1)*dt;           % вектор часу

%% 3. Ініціалізація полів
p   = zeros(Nx, Nt);   % тиск, МПа
m   = zeros(Nx, Nt);   % масовий потік, кг/с
Z   = zeros(Nx, Nt);   % фактор стисливості
rho = zeros(Nx, Nt);   % густина, кг/м^3

% граничні та початкові умови
p_in   = 6;    % вхідний тиск, МПа
p_out  = 3;    % вихідний тиск, МПа
m0     = 100;  % початковий масовий потік, кг/с

% початковий розподіл тиску лінійно між inlet і outlet
p(:,1)   = p_in - (p_in - p_out)*(x/L);
m(:,1)   = m0;
Z(:,1)   = alpha*p(:,1) + beta;
rho(:,1) = p(:,1).*M./(Z(:,1)*R*T);

%% 4. Основний цикл часу (Method of Characteristics)
for n = 1:Nt-1
    % 4.1 локальні швидкості характеристик
    c       = sqrt(R*T.*Z(:,n)./M);
    lambda1 =  c;  % напрям C+
    lambda2 = -c;  % напрям C-
    
    % 4.2 втрати на тертя
    Re = abs(m(:,n))./(rho(:,n)*A) * D ./ nu;
    fr = lambda_f(Re);
    J  = fr .* R*T .* Z(:,n) .* m(:,n) .* abs(m(:,n)) ./ (2*D*A*M*p(:,n));
    
    % 4.3 foot points
    xp = x - lambda1*dt;
    xm = x - lambda2*dt;
    xp = real(xp);  xm = real(xm);
    xp = min(max(xp, x(1)), x(end));
    xm = min(max(xm, x(1)), x(end));
    
    % 4.4 інтерполяція m та Z у точках xp, xm
    m_p = interp1(x, m(:,n), xp, 'linear');
    m_m = interp1(x, m(:,n), xm, 'linear');
    Z_p = interp1(x, Z(:,n), xp, 'linear');
    Z_m = interp1(x, Z(:,n), xm, 'linear');
    
    % 4.5 оновлення m та Z
    m(:,n+1) = 0.5*(m_p + m_m) - J*dt;
    Z(:,n+1) = 0.5*(Z_p + Z_m);
    
    % 4.6 граничні умови
    m(end,   n+1) = m0;        
    p(1,     n+1) = p_in;      
    Z(1,     n+1) = alpha*p_in + beta;
    
    % 4.7 обчислення p та rho
    p(:, n+1)   = m(:,n+1) .* R*T ./ (A .* Z(:,n+1));
    rho(:, n+1) = p(:,n+1) .* M ./ (Z(:,n+1)*R*T);
end

%% 5. Підготовка даних для графіків
p_initial    = p(:,1);
time_vec     = t/3600;                       % години
p_inlet_vec  = p(1, :);
q_std_vec    = m(1, :) ./ rho(1, :) * 3600;  % м^3/год
p_outlet_vec = p(end, :);
q_out_vec    = m(end, :) ./ rho(end, :) * 3600;
Z_vec        = Z(end, :);
p_surface    = p;                            % Nx?Nt

%% === Plots ===
figure('Position', [100, 600, 700, 420]);
plot(x, p_initial, 'LineWidth', 1.2); grid on;
xlabel("L, m"); ylabel("P, MPa");
title("Початковий розподіл тиску");

figure('Position', [900, 600, 700, 420]);
subplot(2,1,1);
plot(time_vec, p_inlet_vec, 'LineWidth', 1.2); grid on;
ylabel("P, MPa");
title("Вхідний тиск");

subplot(2,1,2);
plot(time_vec, q_std_vec, 'LineWidth', 1.2); grid on;
xlabel("t, hours");
ylabel("Q, m^3/hours");
title("Вхідна об’ємна витрата");

figure('Position', [100, 100, 700, 600]);
subplot(3,1,1);
plot(time_vec, p_outlet_vec, 'LineWidth', 1.2); grid on;
ylabel("P, MPa");
title("Вихідний тиск");

subplot(3,1,2);
plot(time_vec, q_out_vec, 'LineWidth', 1.2); grid on;
ylabel("Q, m^3/hours");
title("Вихідна об’ємна витрата");

subplot(3,1,3);
plot(time_vec, Z_vec, 'LineWidth', 1.2); grid on;
xlabel("t, hours");
ylabel("Z");
title("Фактор стисливості на виході");

% 3D поверхня тиску p(x,t)
[XGrid, TGrid] = meshgrid(x, t/3600);   % TGrid у годинах
Psurf = p.';                           % Nt?Nx

figure('Position', [900, 100, 800, 500]);
surf(XGrid, TGrid, Psurf, 'EdgeColor', 'none');
xlabel("L, m"); ylabel("t, hours"); zlabel("P, MPa");
title("Поверхня зміни тиску газу по довжині газопроводу в часі при нестаціонарному режимі руху газу");
view(45, 30);
colorbar;
shading interp;
