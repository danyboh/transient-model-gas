% characteristic_method.m
% Транзієнтний рух газу в магістральному трубопроводі методом характеристик
% На основі [11] та попередніх специфікацій

close all; clear all; clc

%% 1. Параметри труби та газу
D     = 0.5;                    % діаметр труби, м
A     = pi*D^2/4;               % площа перерізу, м^2
L     = 1000;                   % довжина труби, м
R     = 518;                    % газова стала, Дж/(кг·K)
T     = 288;                    % температура, K
M     = 16e-3;                  % молярна маса, кг/моль
alpha = 0.02;                   % коефіцієнт Z(p), 1/МПа
beta  = 1.0;                    % зсув Z(p)
nu    = 1.2e-5;                 % кінематична в’язкість, м^2/с
lambda_f = @(Re) 0.3164 .* Re.^(-0.25);  % формула Дарсі–Вайсбаха

%% 2. Просторова та часова сітка
Nx     = 101;
x      = linspace(0, L, Nx)';    % вектор x як стовпець
dx     = L/(Nx-1);

% CFL-умова: крок за часом у секундах
p_max  = 6;                     
c_max  = sqrt(R*T*(alpha*p_max+beta)/M);
dt_sec = 0.5 * dx ./ c_max;      % крок за часом, с

% задаємо часовий діапазон в годинах
t_end_h = 0.02;                  % кінець моделювання, год
dt_h    = dt_sec/3600;           % крок у годинах
Nt      = ceil(t_end_h/dt_h);
time_h  = (0:Nt-1) * dt_h;        % вектор часу, год

%% 3. Ініціалізація полів
p   = zeros(Nx, Nt);    % тиск, МПа
m   = zeros(Nx, Nt);    % масовий потік, кг/с
Z   = zeros(Nx, Nt);    % фактор стисливості
rho = zeros(Nx, Nt);    % густина, кг/м^3

% граничні та початкові умови
p_in   = 6;    % вхідний тиск, МПа
p_out  = 3;    % вихідний тиск, МПа
m0     = 100;  % початковий масовий потік, кг/с

% початковий розподіл тиску лінійно між inlet і outlet
p(:,1)   = p_in - (p_in - p_out) .* (x/L);
m(:,1)   = m0;
Z(:,1)   = alpha .* p(:,1) + beta;
rho(:,1) = p(:,1) .* M ./ (Z(:,1) * R * T);

%% 4. Основний цикл часу (Method of Characteristics)
for n = 1:Nt-1
    % 4.1 локальні швидкості характеристик з захистом під коренем
    arg = (R * T .* Z(:,n) ./ M);
    arg(arg < 0) = 0;  % уникнути негативного під коренем
    c       = sqrt(arg);
    lambda1 =  c;
    lambda2 = -c;
    
    % 4.2 втрати на тертя з елементною арифметикою
    Re = abs(m(:,n)) ./ (rho(:,n) .* A) .* D ./ nu;
    Re(Re <= 0) = eps;  % уникнути нульових чи від’ємних
    fr = lambda_f(Re);
    J  = fr .* R .* T .* Z(:,n) .* m(:,n) .* abs(m(:,n)) ./ (2 .* D .* A .* M .* p(:,n));
    
    % 4.3 foot points (використовуємо dt_sec)
    xp = x - lambda1 * dt_sec;
    xm = x - lambda2 * dt_sec;
    xp = min(max(xp, x(1)), x(end));
    xm = min(max(xm, x(1)), x(end));
    
    % 4.4 інтерполяція m та Z
    m_p = interp1(x, m(:,n), xp, 'linear');
    m_m = interp1(x, m(:,n), xm, 'linear');
    Z_p = interp1(x, Z(:,n), xp, 'linear');
    Z_m = interp1(x, Z(:,n), xm, 'linear');
    
    % 4.5 оновлення m та Z
    m(:,n+1) = 0.5 .* (m_p + m_m) - J .* dt_sec;
    Z(:,n+1) = 0.5 .* (Z_p + Z_m);
    
    % 4.6 граничні умови
    m(end,   n+1) = m0;
    p(1,     n+1) = p_in;
    Z(1,     n+1) = alpha * p_in + beta;
    
    % 4.7 обчислення p та rho
    p(:, n+1)   = m(:,n+1) .* R .* T ./ (A .* Z(:,n+1));
    rho(:, n+1) = p(:,n+1) .* M ./ (Z(:,n+1) * R * T);
end

%% 5. Підготовка даних для графіків
p_initial    = p(:,1);
p_inlet_vec  = p(1, :);
q_std_vec    = m(1, :) ./ rho(1, :) * 3600;
p_outlet_vec = p(end, :);
q_out_vec    = m(end, :) ./ rho(end, :) * 3600;
Z_vec        = Z(end, :);
Psurf        = p.';  % Nt?Nx

%% === Plots ===
figure('Position', [100, 600, 700, 420]);
plot(x, p_initial, 'LineWidth', 1.2); grid on;
xlabel("L, m"); ylabel("P, MPa");
title("Початковий розподіл тиску");

figure('Position', [900, 600, 700, 420]);
subplot(2,1,1);
plot(time_h, p_inlet_vec, 'LineWidth', 1.2); grid on;
ylabel("P, MPa"); title("Вхідний тиск");
subplot(2,1,2);
plot(time_h, q_std_vec, 'LineWidth', 1.2); grid on;
xlabel("t, hours"); ylabel("Q, m^3/hours");
title("Вхідна об’ємна витрата");

figure('Position', [100, 100, 700, 600]);
subplot(3,1,1);
plot(time_h, p_outlet_vec, 'LineWidth', 1.2); grid on;
ylabel("P, MPa"); title("Вихідний тиск");
subplot(3,1,2);
plot(time_h, q_out_vec, 'LineWidth', 1.2); grid on;
ylabel("Q, m^3/hours"); title("Вихідна об’ємна витрата");
subplot(3,1,3);
plot(time_h, Z_vec, 'LineWidth', 1.2); grid on;
xlabel("t, hours"); ylabel("Z");
title("Фактор стисливості на виході");

[XGrid, TGrid] = meshgrid(x, time_h);
figure('Position', [900, 100, 800, 500]);
surf(XGrid, TGrid, Psurf, 'EdgeColor', 'none');
xlabel("L, m"); ylabel("t, hours"); zlabel("P, MPa");
title("Поверхня зміни тиску газу по довжині газопроводу в часі при нестаціонарному режимі руху газу");
view(45, 30); colorbar; shading interp;
