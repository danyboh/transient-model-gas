% === ЄДИНИЙ СКРИПТ ДЛЯ МОДЕЛЮВАННЯ НЕСТАЦІОНАРНОГО РУХУ ГАЗУ ===

close all;
clc; clear;

% === Параметри моделі ===
L = 1000; Nx = 1000; T_end = 5; CFL = 0.5;
x = linspace(0, L, Nx);

% === Початкові умови ===
p0 = 5.5e6; pL = 5.0e6; T0 = 288.15;
R = 518.3; x_mol = [0.94 0.02 0.01 0.005 0.005 0.01 0.008 0.002];

% === Властивості газу ===
A = 0.005; B = 0.001;
Z = 1 + A * (p0 / 1e6) - B * T0;
M = [16.043, 30.07, 44.097, 58.12, 58.12, 28.01, 44.01, 34.08];
Mm = sum(x_mol .* M);
Cv = 1.5 * 8.31451 / Mm * 1000;

% === Початкові розподіли ===
p_profile = linspace(p0, pL, Nx);
rho_profile = p_profile ./ (Z * R * T0);
u_profile = linspace(1, 3, Nx);
e_profile = Cv * T0 + 0.5 * u_profile.^2;

U = zeros(3, Nx);
U(1,:) = rho_profile;
U(2,:) = rho_profile .* u_profile;
U(3,:) = rho_profile .* e_profile;

p_initial = p_profile / 1e6;

% === Ініціалізація ===
dx = x(2) - x(1);
t = 0; nt = 0;
time_vec = [];
p_inlet_vec = [];
q_std_vec = [];
p_outlet_vec = [];
q_out_vec = [];
Z_vec = [];
p_surface = [];

while t < T_end
    rho = U(1,:);
    u = U(2,:) ./ rho;
    e = U(3,:) ./ rho;
    T = (e - 0.5 * u.^2) / Cv;
    p = rho .* R .* T;
    c = sqrt(1.3 * p ./ rho);
    umax = max(abs(u) + c);
    dt = CFL * dx / umax;
    if t + dt > T_end
        dt = T_end - t;
    end

    nt = nt + 1;
    time_vec(nt) = t / 3600;
    p_inlet_vec(nt) = p(1) / 1e6;
    q_std_vec(nt) = u(1) * rho(1) * (T0 / T(1)) / (Z * R) * 3600;
    p_outlet_vec(nt) = p(end) / 1e6;
    q_out_vec(nt) = u(end) * rho(end) * (T0 / T(end)) / (Z * R) * 3600;
    Z_vec(nt) = 1 + A * (p(end) / 1e6) - B * T(end);
    p_surface(nt,:) = p / 1e6;

    F = [rho .* u;
         rho .* u.^2 + p;
         (U(3,:)+p) .* u];
    F_plus = [F(:,2:end), F(:,end)];
    U_plus = [U(:,2:end), U(:,end)];
    lambda = max(abs(u) + c);
    Flux = 0.5 * (F_plus + F) - 0.5 * lambda * (U_plus - U);
    RHS1 = -(Flux - [Flux(:,1), Flux(:,1:end-1)]) / dx;

    U1 = U + dt * RHS1;

    rho = U1(1,:);
    u = U1(2,:) ./ rho;
    e = U1(3,:) ./ rho;
    T = (e - 0.5 * u.^2) / Cv;
    p = rho .* R .* T;
    F = [rho .* u;
         rho .* u.^2 + p;
         (U1(3,:)+p) .* u];
    F_plus = [F(:,2:end), F(:,end)];
    U_plus = [U1(:,2:end), U1(:,end)];
    lambda = max(abs(u) + c);
    Flux = 0.5 * (F_plus + F) - 0.5 * lambda * (U_plus - U1);
    RHS2 = -(Flux - [Flux(:,1), Flux(:,1:end-1)]) / dx;

    U = 0.5 * (U + U1 + dt * RHS2);
    t = t + dt;
end

% === Графік №1 ===
figure('Position', [100, 600, 700, 420]);
plot(x, p_initial, 'LineWidth', 1.6);
grid on;
xlabel("Відстань по трубопроводу [м]");
ylabel("Тиск [МПа]");
title("Розподіл тиску природного газу в початковий момент часу");

% === Графік №2 ===
figure('Position', [900, 600, 700, 420]);
subplot(2,1,1);
plot(time_vec, p_inlet_vec, 'LineWidth', 1.5);
ylabel("Тиск на вході [МПа]");
grid on; title("Зміна тиску на вході в часі");

subplot(2,1,2);
plot(time_vec, q_std_vec, 'LineWidth', 1.5);
xlabel("Час [год]"); ylabel("Витрата [м^3/год]");
grid on; title("Зміна об''ємної витрати (зведеної до стандартних умов)");

% === Графік №3 ===
figure('Position', [100, 100, 700, 600]);
subplot(3,1,1);
plot(time_vec, p_outlet_vec, 'LineWidth', 1.5);
ylabel("Тиск на виході [МПа]");
grid on; title("Зміна тиску на виході в часі");

subplot(3,1,2);
plot(time_vec, q_out_vec, 'LineWidth', 1.5);
ylabel("Витрата [м^3/год]");
grid on; title("Зміна витрати на виході в часі");

subplot(3,1,3);
plot(time_vec, Z_vec, 'LineWidth', 1.5);
xlabel("Час [год]"); ylabel("Коеф. стисливості Z");
grid on; title("Зміна коефіцієнта стисливості в часі");

% === Графік №4: Поверхня тиску ===
[TimeGrid, XGrid] = meshgrid(time_vec, x);
figure('Position', [900, 100, 800, 500]);
surf(XGrid, TimeGrid, p_surface', 'EdgeColor', 'none');
xlabel("Відстань [м]"); ylabel("Час [год]"); zlabel("Тиск [МПа]");
title("Поверхня зміни тиску газу по довжині газопроводу в часі");
view(45, 30); colorbar; shading interp;
