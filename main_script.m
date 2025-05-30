
% === SINGLE SCRIPT FOR NONSTATIONARY GAS FLOW MODELING WITH LOOPED Fvnic + SAFETY CHECKS ===

close all;
clc; clear;

L = 10; Nx = 2; T_end = 5; CFL = 0.5;
x = linspace(0, L, Nx);

p0 = 5.5e6; pL = 5.0e6; T0 = 288.15;
R = 518.3; x_mol = [0.94 0.02 0.01 0.005 0.005 0.01 0.008 0.002];
M = [16.043, 30.07, 44.097, 58.12, 58.12, 28.01, 44.01, 34.08];
Mm = sum(x_mol .* M);

p_profile = linspace(p0, pL, Nx);
T_init = T0 * ones(1, Nx);
Z = zeros(1, Nx);
for i = 1:Nx
    [~, Z(i), ~, ~, ~, ~, ~] = Fvnic(p_profile(i), T_init(i), x_mol);
end
rho_profile = p_profile ./ (Z .* R .* T_init);
u_profile = linspace(0.1, 0.5, Nx);
[Cp_init, ~, ~] = Cp_Vnic(p_profile, T0, x_mol);
Cv_init = real(Cp_init - R);
Cv_init(Cv_init < 10) = 10;
e_profile = Cv_init .* T0 + 0.5 * u_profile.^2;

U = zeros(3, Nx);
U(1,:) = rho_profile;
U(2,:) = rho_profile .* u_profile;
U(3,:) = rho_profile .* e_profile;

p_initial = p_profile / 1e6;

dx = x(2) - x(1);
t = 0; nt = 0;
time_vec = []; p_inlet_vec = []; q_std_vec = [];
p_outlet_vec = []; q_out_vec = []; Z_vec = []; p_surface = [];

while t < T_end
    rho = U(1,:); u = real(U(2,:) ./ rho); e = real(U(3,:) ./ rho);
    T = real(e - 0.5 * u.^2);
    T(T < 100) = 100;
    Cp = zeros(1, Nx); Cv = zeros(1, Nx); Z = zeros(1, Nx);
    for i = 1:Nx
        [Cp(i), ~, ~] = Cp_Vnic(rho(i) * R * T(i), T(i), x_mol);
        [~, Z(i), ~, ~, ~, ~, ~] = Fvnic(rho(i) * R * T(i), T(i), x_mol);
    end
    Cv = real(Cp - R); Cv(Cv < 10) = 10;
    T = real(e - 0.5 * u.^2) ./ Cv;
    T(T < 100) = 100;
    p = real(rho .* R .* T .* Z);
    c = sqrt(abs(1.3 * p ./ rho));
    umax = max(abs(u) + c);
    dt = CFL * dx / umax;
    if t + dt > T_end; dt = T_end - t; end

    nt = nt + 1;
    time_vec(nt) = real(t / 3600);
    p_inlet_vec(nt) = real(p(1) / 1e6);
    q_std_vec(nt) = real(u(1) * rho(1) * (T0 / T(1)) / (Z(1) * R) * 3600);
    p_outlet_vec(nt) = real(p(end) / 1e6);
    q_out_vec(nt) = real(u(end) * rho(end) * (T0 / T(end)) / (Z(end) * R) * 3600);
    Z_vec(nt) = real(Z(end));
    p_surface(nt,:) = real(p / 1e6);

    F = [rho .* u;
         rho .* u.^2 + p;
         (U(3,:)+p) .* u];
    F_plus = [F(:,2:end), F(:,end)];
    U_plus = [U(:,2:end), U(:,end)];
    lambda = max(abs(u) + c);
    Flux = 0.5 * (F_plus + F) - 0.5 * lambda * (U_plus - U);
    RHS1 = -(Flux - [Flux(:,1), Flux(:,1:end-1)]) / dx;

    U1 = U + dt * RHS1;

    rho = U1(1,:); u = real(U1(2,:) ./ rho); e = real(U1(3,:) ./ rho);
    T = real(e - 0.5 * u.^2);
    T(T < 100) = 100;
    for i = 1:Nx
        [Cp(i), ~, ~] = Cp_Vnic(rho(i) * R * T(i), T(i), x_mol);
        [~, Z(i), ~, ~, ~, ~, ~] = Fvnic(rho(i) * R * T(i), T(i), x_mol);
    end
    Cv = real(Cp - R); Cv(Cv < 10) = 10;
    T = real(e - 0.5 * u.^2) ./ Cv;
    T(T < 100) = 100;
    p = real(rho .* R .* T .* Z);
    F = [rho .* u;
         rho .* u.^2 + p;
         (U1(3,:)+p) .* u];
    F_plus = [F(:,2:end), F(:,end)];
    U_plus = [U1(:,2:end), U1(:,end)];
    lambda = max(abs(u) + c);
    Flux = 0.5 * (F_plus + F) - 0.5 * lambda * (U_plus - U1);
    RHS2 = -(Flux - [Flux(:,1), Flux(:,1:end-1)]) / dx;

    U = 0.5 * (U + U1 + dt * RHS2);
    t = real(t + dt);
end

figure('Position', [100, 600, 700, 420]);
plot(x, p_initial, 'LineWidth', 1.6);
grid on;
xlabel("L, m");
ylabel("P, MPa");

figure('Position', [900, 600, 700, 420]);
subplot(2,1,1);
plot(time_vec, p_inlet_vec, 'LineWidth', 1.5);
ylabel("P, MPa");
grid on;

subplot(2,1,2);
plot(time_vec, q_std_vec, 'LineWidth', 1.5);
xlabel("t, hours"); ylabel("Q, m^3/hours");
grid on;

figure('Position', [100, 100, 700, 600]);
subplot(3,1,1);
plot(time_vec, p_outlet_vec, 'LineWidth', 1.5);
ylabel("P, MPa");
grid on;

subplot(3,1,2);
plot(time_vec, q_out_vec, 'LineWidth', 1.5);
ylabel("Q, m^3/hours");
grid on;

subplot(3,1,3);
plot(time_vec, Z_vec, 'LineWidth', 1.5);
xlabel("t, hours"); ylabel("Z");
grid on;

[TimeGrid, XGrid] = meshgrid(time_vec, x);
figure('Position', [900, 100, 800, 500]);
surf(XGrid, TimeGrid, p_surface', 'EdgeColor', 'none');
xlabel("L, m"); ylabel("t, hours"); zlabel("P, MPa");
view(45, 30); colorbar; shading interp;
