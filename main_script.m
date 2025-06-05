% === main_script_ideal.m ===

close all;
clc; clear;

% === Parameters ===
L = 2000; Nx = 500; T_end = 300; CFL = 0.5;
x = linspace(0, L, Nx);
dx = x(2) - x(1);

% === Constants (ideal gas, methane-based) ===
R = 518.3;        % [J/(kg·K)]
Cv = 35;          % [J/(kg·K)]
Cp = Cv + R;      % [J/(kg·K)]
Z = 1;            % ideal gas

% === Initial conditions ===
p0 = 5.5e6; pL = 5.0e6; T0 = 288.15;
p_profile = linspace(p0, pL, Nx);
T_init = T0 * ones(1, Nx);
rho_profile = p_profile ./ (Z * R * T_init);
u_profile = linspace(0.1, 0.5, Nx);
e_profile = Cv * T_init + 0.5 * u_profile.^2;

U = zeros(3, Nx);
U(1,:) = rho_profile;
U(2,:) = rho_profile .* u_profile;
U(3,:) = rho_profile .* e_profile;

p_initial = p_profile / 1e6;
t = 0; nt = 0;
time_vec = []; p_inlet_vec = []; q_std_vec = [];
p_outlet_vec = []; q_out_vec = []; Z_vec = []; p_surface = [];

while t < T_end
    rho = U(1,:);
    u = U(2,:) ./ rho;
    e = U(3,:) ./ rho;
    T = (e - 0.5 * u.^2) / Cv;
    p = rho .* R .* T;

    c = sqrt(1.3 * p ./ rho);
    umax = max(abs(u) + c);
    dt = CFL * dx / umax;
    if t + dt > T_end; dt = T_end - t; end

    nt = nt + 1;
    time_vec(nt) = t / 3600;
    p_inlet_vec(nt) = p(1) / 1e6;
    q_std_vec(nt) = u(1) * rho(1) * (T0 / T(1)) / (Z * R) * 3600;
    p_outlet_vec(nt) = p(end) / 1e6;
    q_out_vec(nt) = u(end) * rho(end) * (T0 / T(end)) / (Z * R) * 3600;
    Z_vec(nt) = Z;
    p_surface(nt,:) = p / 1e6;

    F = [rho .* u;
         rho .* u.^2 + p;
         (U(3,:) + p) .* u];

    F_plus = [F(:,2:end), F(:,end)];
    U_plus = [U(:,2:end), U(:,end)];
    lambda = max(abs(u) + c);
    Flux = 0.5 * (F_plus + F) - 0.5 * lambda * (U_plus - U);
    RHS1 = -(Flux - [Flux(:,1), Flux(:,1:end-1)]) / dx;
    U1 = U + dt * RHS1;

    rho = U1(1,:); u = U1(2,:) ./ rho; e = U1(3,:) ./ rho;
    T = (e - 0.5 * u.^2) / Cv; p = rho .* R .* T;
    F = [rho .* u;
         rho .* u.^2 + p;
         (U1(3,:) + p) .* u];
    F_plus = [F(:,2:end), F(:,end)];
    U_plus = [U1(:,2:end), U1(:,end)];
    lambda = max(abs(u) + c);
    Flux = 0.5 * (F_plus + F) - 0.5 * lambda * (U_plus - U1);
    RHS2 = -(Flux - [Flux(:,1), Flux(:,1:end-1)]) / dx;

    U = 0.5 * (U + U1 + dt * RHS2);
    t = t + dt;
end

% === Plots ===
figure('Position', [100, 600, 700, 420]);
plot(x, p_initial); grid on; xlabel("L, m"); ylabel("P, MPa");

figure('Position', [900, 600, 700, 420]);
subplot(2,1,1); plot(time_vec, p_inlet_vec); ylabel("P, MPa"); grid on;
subplot(2,1,2); plot(time_vec, q_std_vec); xlabel("t, hours"); ylabel("Q, m^3/hours"); grid on;

figure('Position', [100, 100, 700, 600]);
subplot(3,1,1); plot(time_vec, p_outlet_vec); ylabel("P, MPa"); grid on;
subplot(3,1,2); plot(time_vec, q_out_vec); ylabel("Q, m^3/hours"); grid on;
subplot(3,1,3); plot(time_vec, Z_vec); xlabel("t, hours"); ylabel("Z"); grid on;

[TimeGrid, XGrid] = meshgrid(time_vec, x);
figure('Position', [900, 100, 800, 500]);
surf(XGrid, TimeGrid, p_surface', 'EdgeColor', 'none');
xlabel("L, m"); ylabel("t, hours"); zlabel("P, MPa");
view(45, 30); colorbar; shading interp;
