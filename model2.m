% characteristic_method_perturbation.m
% Transient gas flow in a pipeline with short-term perturbation
% Method of Characteristics, modified for 1-minute simulation

close all; clear all; clc

%% 1. Pipeline and gas parameters
D     = 0.5;                    % pipe diameter, m
A     = pi*D^2/4;               % cross-sectional area, m^2
L     = 1000;                   % pipeline length, m
R     = 518;                    % gas constant, J/(kg·K)
T     = 288;                    % temperature, K
M     = 16e-3;                  % molar mass, kg/mol
alpha = 0.02e-6;                % Z(p) coefficient, 1/Pa (converted from 1/MPa)
beta  = 1.0;                    % Z(p) offset
nu    = 1.2e-5;                 % kinematic viscosity, m^2/s
lambda_f = @(Re) 0.3164 .* Re.^(-0.25);  % Darcy–Weisbach friction

%% 2. Spatial and temporal grid (1 minute = 1/60 hours)
Nx      = 201;  % Increased resolution
x       = linspace(0, L, Nx)';    % x coordinates as column vector
dx      = L/(Nx-1);

% CFL condition: time step in seconds - use much stricter condition
p_max   = 6.5e6;                  % max pressure in Pa (including pulse)
c_max   = sqrt(R*T/M);            % ideal gas sound speed
dt_sec  = 0.05 * dx/c_max;        % much stricter CFL condition (0.05 instead of 0.1)

% Shorten simulation to focus on wave transit (about 0.4 seconds)
t_end_h = 0.4/3600;  % 0.4 seconds in hours
dt_h    = dt_sec/3600;           % time step, hours
Nt      = ceil(t_end_h/dt_h);
time_h  = (0:Nt-1) * dt_h;       % time vector, hours

%% 3. Inlet perturbation – Gaussian pressure pulse (in Pa)
p0       = 6e6;                  % base pressure, Pa
p_out    = 3e6;                  % outlet pressure, Pa
A_p      = 0.5e6;                % pulse amplitude, Pa
t_center = 0.2/3600;             % pulse center at 0.2 seconds
sigma    = 0.1/3600;             % narrower pulse width
p_inlet_vec = p0 + A_p * exp(-((time_h - t_center).^2)/(2*sigma^2));

%% 4. Field initialization (all in SI units)
p   = zeros(Nx, Nt);    % pressure, Pa
m   = zeros(Nx, Nt);    % mass flow, kg/s
Z   = zeros(Nx, Nt);    % compressibility factor
rho = zeros(Nx, Nt);    % density, kg/m^3

m0     = 4.44;          % initial mass flow, kg/s (for 25,000 m³/h)

% initialize with pressure gradient from inlet to outlet
p(:,1)   = p0 - (p0 - p_out) * (x/L);
m(:,1)   = m0;
Z(:,1)   = alpha * p(:,1) + beta;
rho(:,1) = p(:,1) .* M ./ (Z(:,1) * R * T);

% Verify initial conditions are consistent
c_init = sqrt(R * T / M);  % ideal gas sound speed

%% 5. Time-marching loop (Method of Characteristics) - FIXED FOOT POINT HANDLING
% Add debugging for interior nodes
mid_node = round(Nx/2);
quarter_node = round(Nx/4);

for n = 1:Nt-1
    % 5.1 MOC with ideal gas and friction
    % Use constant sound speed for ideal gas
    c = sqrt(R * T / M);  % constant sound speed for ideal gas
    
    % 5.2 Compute velocity and Reynolds number for friction
    V = m(:,n) ./ (rho(:,n) .* A);  % velocity in m/s
    Re = abs(V) * D / nu;  % Reynolds number
    Re = max(Re, 2300);    % minimum for turbulent flow
    lambda = lambda_f(Re); % friction factor
    
    % 5.3 Foot points for characteristics
    xp = x - c * dt_sec;  % C+ characteristic
    xm = x + c * dt_sec;  % C- characteristic
    
    % 5.4 Handle out-of-bounds foot points properly
    % Initialize arrays for pressure and mass flow at foot points
    p_p = zeros(size(x));
    p_m = zeros(size(x));
    m_p = zeros(size(x));
    m_m = zeros(size(x));
    
    for i = 1:Nx
        % Handle C+ characteristic foot point
        if xp(i) < x(1)
            % Use inlet boundary condition (time-varying pulse)
            p_p(i) = p_inlet_vec(n);
            m_p(i) = m(1,n);
        elseif xp(i) > x(end)
            % Use outlet boundary value
            p_p(i) = p(end,n);
            m_p(i) = m(end,n);
        else
            % Foot point is within domain - interpolate
            p_p(i) = interp1(x, p(:,n), xp(i), 'linear');
            m_p(i) = interp1(x, m(:,n), xp(i), 'linear');
        end
        
        % Handle C- characteristic foot point
        if xm(i) < x(1)
            % Use inlet boundary condition (time-varying pulse)
            p_m(i) = p_inlet_vec(n);
            m_m(i) = m(1,n);
        elseif xm(i) > x(end)
            % Use outlet boundary value
            p_m(i) = p(end,n);
            m_m(i) = m(end,n);
        else
            % Foot point is within domain - interpolate
            p_m(i) = interp1(x, p(:,n), xm(i), 'linear');
            m_m(i) = interp1(x, m(:,n), xm(i), 'linear');
        end
    end
    
    % 5.5 Update interior points using characteristic equations with friction
    % For compressible flow: dp/dt ± ρc*dm/dt = friction terms along characteristics
    for i = 2:Nx-1
        % C+ characteristic: dp/dt + ρc*dm/dt = -friction
        % C- characteristic: dp/dt - ρc*dm/dt = -friction
        rho_c = rho(i,n) * c;
        
        % Friction term: λ|V|V/(2D) * ρc * dt
        V_i = V(i);
        friction_term = lambda(i) * abs(V_i) * V_i / (2 * D) * rho_c * dt_sec;
        
        % CORRECTED SIGNS: Use characteristic equations with proper signs
        % When inlet pressure increases, interior pressure should increase
        p_new = 0.5 * (p_p(i) + p_m(i)) + 0.5 * rho_c * (m_p(i) - m_m(i)) / A - friction_term;
        m_new = 0.5 * (m_p(i) + m_m(i)) + 0.5 * A * (p_p(i) - p_m(i)) / rho_c - friction_term * A / rho_c;
        
        p(i,n+1) = p_new;
        m(i,n+1) = m_new;
    end
    
    % 5.6 Apply boundary conditions AFTER interior update
    % Inlet: prescribed pressure
    p(1,n+1) = p_inlet_vec(n+1);
    
    % Outlet: constant mass flow (but let pressure evolve naturally)
    m(end,n+1) = m0;
    % Enforce outlet pressure to initial profile
    p(end, n+1) = p0 - (p0 - p_out);
    % Don't overwrite outlet pressure - let the wave propagate there
    
    % 5.7 compute density (ideal gas law)
    rho(:,n+1) = p(:,n+1) .* M ./ (R * T);
    
    % 5.8 bounds checking (prevent negative/zero pressure)
    p(:,n+1) = max(p(:,n+1), 1e5);   % Minimum 0.1 MPa
    p(:,n+1) = min(p(:,n+1), 1e7);   % Maximum 10 MPa
    m(:,n+1) = max(m(:,n+1), -1000); % Allow some reverse flow
    m(:,n+1) = min(m(:,n+1), 1000);  % Maximum mass flow
end

%% 6. Prepare and visualize results (convert to engineering units)
P0      = p(:,1) / 1e6;         % MPa
pin     = p(1,:) / 1e6;         % MPa
pout    = p(end,:) / 1e6;       % MPa

% Calculate volumetric flow rates in m³/s, then convert to m³/h
qstd    = m(1,:) ./ rho(1,:) * 3600;  % m³/h
qout    = m(end,:) ./ rho(end,:) * 3600;  % m³/h

Zout    = Z(end,:);

% Remove any NaN/Inf values before plotting
Psurf = (p.' / 1e6);  % MPa
Psurf(isnan(Psurf) | isinf(Psurf)) = 0;  % Replace with zeros

% Debug: Check wave propagation
fprintf('\nSimulation results:\n');
fprintf('Max inlet pressure: %.2f MPa\n', max(pin));
fprintf('Min inlet pressure: %.2f MPa\n', min(pin));
fprintf('Max outlet pressure: %.2f MPa\n', max(pout));
fprintf('Min outlet pressure: %.2f MPa\n', min(pout));
fprintf('Max inlet flow: %.0f m³/h\n', max(qstd));
fprintf('Min inlet flow: %.0f m³/h\n', min(qstd));

% Check for any NaN or Inf values
if any(isnan(p(:))) || any(isinf(p(:)))
    fprintf('WARNING: NaN or Inf values detected in pressure!\n');
end
if any(isnan(m(:))) || any(isinf(m(:)))
    fprintf('WARNING: NaN or Inf values detected in mass flow!\n');
end
if any(isnan(Z(:))) || any(isinf(Z(:)))
    fprintf('WARNING: NaN or Inf values detected in Z-factor!\n');
end

% Verify units and wave propagation
fprintf('\nUnit verification:\n');
fprintf('Inlet pressure range: %.1f to %.1f Pa\n', min(p(1,:)), max(p(1,:)));
fprintf('Inlet mass flow range: %.1f to %.1f kg/s\n', min(m(1,:)), max(m(1,:)));
fprintf('Expected wave speed: %.0f m/s\n', c_init);
fprintf('Expected transit time: %.2f s\n', L/c_init);
fprintf('Simulation duration: %.2f s\n', time_h(end)*3600);

% Check wave propagation at different points
quarter_node = round(Nx/4);
three_quarter_node = round(3*Nx/4);
fprintf('\nWave propagation check:\n');
fprintf('Quarter-pipe (x=%.0f m): P range %.2f-%.2f MPa\n', ...
        x(quarter_node), min(p(quarter_node,:))/1e6, max(p(quarter_node,:))/1e6);
fprintf('Mid-pipe (x=%.0f m): P range %.2f-%.2f MPa\n', ...
        x(mid_node), min(p(mid_node,:))/1e6, max(p(mid_node,:))/1e6);
fprintf('Three-quarter (x=%.0f m): P range %.2f-%.2f MPa\n', ...
        x(three_quarter_node), min(p(three_quarter_node,:))/1e6, max(p(three_quarter_node,:))/1e6);

% initial pressure distribution
figure('Position',[100,600,700,420]);
plot(x, P0, 'LineWidth',1.2); grid on;
xlabel('L, m'); ylabel('P, MPa');
title('Initial Pressure Distribution');

% inlet pressure and flow
figure('Position',[900,600,700,420]);
subplot(2,1,1);
plot(time_h*3600, pin, 'LineWidth',1.2); grid on;
ylabel('P, MPa'); title('Inlet Pressure');
xlabel('t, s');
subplot(2,1,2);
plot(time_h*3600, qstd, 'LineWidth',1.2); grid on;
xlabel('t, s'); ylabel('Q, m^3/h'); title('Inlet Flow');

% outlet pressure, flow, and Z
figure('Position',[100,100,700,600]);
subplot(3,1,1);
plot(time_h*3600, pout, 'LineWidth',1.2); grid on;
ylabel('P, MPa'); title('Outlet Pressure');
subplot(3,1,2);
plot(time_h*3600, qout, 'LineWidth',1.2); grid on;
ylabel('Q, m^3/h'); title('Outlet Flow');
subplot(3,1,3);
plot(time_h*3600, Zout, 'LineWidth',1.2); grid on;
xlabel('t, s'); ylabel('Z'); title('Outlet Compressibility Factor');

% 3D pressure surface
[Xg, Tg] = meshgrid(x, time_h*3600);
figure('Position',[900,100,800,500]);
surf(Xg, Tg, Psurf, 'EdgeColor', 'none');
shading interp; colorbar;
xlabel('L, m'); ylabel('t, s'); zlabel('P, MPa');
title('Pressure Surface over Length and Time');
view(45,30);

% Add interior node tracking plot
figure('Position',[100,50,800,400]);
subplot(2,1,1);
plot(time_h*3600, p(mid_node,:)/1e6, 'LineWidth',1.2); grid on;
ylabel('P, MPa'); title(sprintf('Mid-pipe Pressure (x = %.1f m)', x(mid_node)));
xlabel('t, s');
subplot(2,1,2);
plot(time_h*3600, m(mid_node,:), 'LineWidth',1.2); grid on;
xlabel('t, s'); ylabel('m, kg/s'); title(sprintf('Mid-pipe Mass Flow (x = %.1f m)', x(mid_node)));

% Add wave propagation comparison plot
figure('Position',[100,200,800,500]);
plot(time_h*3600, p(quarter_node,:)/1e6, 'LineWidth',1.5, 'DisplayName', sprintf('x=%.0f m', x(quarter_node))); hold on;
plot(time_h*3600, p(mid_node,:)/1e6, 'LineWidth',1.5, 'DisplayName', sprintf('x=%.0f m', x(mid_node)));
plot(time_h*3600, p(three_quarter_node,:)/1e6, 'LineWidth',1.5, 'DisplayName', sprintf('x=%.0f m', x(three_quarter_node)));
plot(time_h*3600, pin, '--', 'LineWidth',1.2, 'DisplayName', 'Inlet');
grid on; legend('Location', 'best');
xlabel('t, s'); ylabel('P, MPa');
title('Wave Propagation Comparison');
