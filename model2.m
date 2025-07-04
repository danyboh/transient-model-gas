% characteristic_method_perturbation.m
% Transient gas flow in a pipeline with short-term perturbation
% Method of Characteristics, modified for 1-minute simulation
% CORRECTED: Physically realistic inlet boundary condition

close all; clear all; clc

%% 1. Pipeline and gas parameters
D     = 0.5;                    % pipe diameter, m
A     = pi*D^2/4;               % cross-sectional area, m^2
L     = 1000;                   % pipeline length, m
R     = 518;                    % specific gas constant, J/(kg·K)
T     = 288;                    % temperature, K
M     = 16e-3;                  % molar mass, kg/mol (for reference only)
% FIXED: Use proper Z-factor for natural gas (should be close to 1, not linear with pressure)
Z_base = 0.85;                  % base compressibility factor for natural gas
nu    = 1.2e-5;                 % kinematic viscosity, m^2/s
lambda_f = @(Re) 0.3164 .* Re.^(-0.25);  % Darcy–Weisbach friction

%% 2. Source and consumer parameters (NEW: Physically realistic boundary models)
% Source operates as: P_source - P_inlet = R_source * Q_inlet^2 + L_source * dQ/dt
% This represents a compressor or pressure regulator with internal resistance
p_source_base = 6e6;            % base source pressure, Pa
R_source = 5e2;                 % source resistance coefficient, Pa·s²/m⁶ (reduced for more dynamics)
L_source = 2e2;                 % source inductance coefficient, Pa·s²/m³ (increased for inertia)
C_source = 1e-8;                % source capacitance (pressure buffer), m³/Pa

% Consumer/outlet parameters
% Consumer operates as: P_outlet = P_consumer - R_consumer * Q_outlet^2
p_consumer_base = 3e6;          % base consumer pressure demand, Pa
R_consumer = 2e2;               % consumer resistance coefficient, Pa·s²/m⁶

%% 3. Spatial and temporal grid 
Nx      = 101;   % Finer grid for better wave resolution (increased from 51)
x       = linspace(0, L, Nx)';    % x coordinates as column vector
dx      = L/(Nx-1);

% Adaptive time step for long-term simulation
p_max   = 6.5e6;                  % max pressure in Pa (including pulse)
c_max   = sqrt(R*T*Z_base/M);     % sound speed with proper Z-factor

% Use adaptive time step for better wave resolution
% For transient phase (first 5 seconds): small time step
% For steady phase (after 5 seconds): larger time step
dt_sec_transient = 0.05 * dx/c_max;  % smaller time step for first 5 seconds
dt_sec_steady = 0.5;                 % smaller steady time step for better dynamics
dt_sec = dt_sec_transient;          % start with small time step

% Full-scale simulation for 30 minutes
t_end_h = 30/60;  % 30 minutes in hours

% Create time vector with adaptive sampling
% Dense sampling for first 5 seconds, sparse for rest
t_dense = 0:dt_sec_transient/3600:5/3600;  % first 5 seconds
t_sparse = 5/3600:dt_sec_steady/3600:t_end_h;  % rest of simulation
time_h = [t_dense, t_sparse(2:end)];  % combine, avoid duplicate at 5s
Nt = length(time_h);

fprintf('Total time steps: %d\n', Nt);
fprintf('Expected simulation time: ~%.1f minutes\n', Nt/1000);

%% 4. Source operational scenario – Extended pressure variation for 30-minute simulation
% Instead of rigid pressure, define source pressure potential that interacts with pipeline
p_source_vec = zeros(size(time_h));

for i = 1:length(time_h)
    t_min = time_h(i) * 60;  % time in minutes
    
    if t_min < 2  % Initial startup (0-2 min)
        p_source_vec(i) = p_source_base + 0.5e6 * (1 - exp(-t_min/0.5));  % gradual pressure rise
    elseif t_min < 5  % Stabilization (2-5 min)
        p_source_vec(i) = p_source_base + 0.5e6;  % high pressure
    elseif t_min < 10  % Load increase (5-10 min)
        p_source_vec(i) = p_source_base + 0.5e6 + 0.3e6 * sin(2*pi*(t_min-5)/5);  % oscillating load
    elseif t_min < 15  % Normal operation (10-15 min)
        p_source_vec(i) = p_source_base + 0.2e6;  % moderate pressure
    elseif t_min < 20  % Demand spike (15-20 min)
        p_source_vec(i) = p_source_base + 0.8e6 * exp(-((t_min-17.5)/2)^2);  % Gaussian spike
    elseif t_min < 25  % Gradual reduction (20-25 min)
        p_source_vec(i) = p_source_base + 0.3e6 * (1 - (t_min-20)/5);  % linear decrease
    else  % Final phase (25-30 min)
        p_source_vec(i) = p_source_base + 0.1e6 * cos(2*pi*(t_min-25)/5);  % small oscillations
    end
end

% Ensure no negative pressures
p_source_vec = max(p_source_vec, 2e6);

%% 5. Field initialization (all in SI units)
p   = zeros(Nx, Nt);    % pressure, Pa
m   = zeros(Nx, Nt);    % mass flow, kg/s
Z   = zeros(Nx, Nt);    % compressibility factor
rho = zeros(Nx, Nt);    % density, kg/m^3

% Source state variables
p_source_buffer = zeros(1, Nt);  % source pressure buffer
m_source_prev = 0;               % previous source mass flow

p_out = p_consumer_base;  % outlet pressure equals consumer base pressure, Pa
m0     = 4.44;          % initial mass flow, kg/s (for 25,000 m³/h)

% initialize with pressure gradient from inlet to outlet
p(:,1)   = p_source_base - (p_source_base - p_out) * (x/L);
m(:,1)   = m0;
% FIXED: Use proper Z-factor that doesn't go to zero
Z(:,1)   = Z_base * ones(Nx, 1);  % constant Z-factor
% CORRECTED: Use proper density formula with specific gas constant
rho(:,1) = p(:,1) ./ (Z(:,1) * R * T);  % kg/m³

% Initialize source buffer
p_source_buffer(1) = p_source_vec(1);

% Verify initial conditions are consistent
c_init = sqrt(R * T * Z_base);  % sound speed with proper Z-factor, m/s

%% 6. Time-marching loop (Method of Characteristics) - ADAPTIVE TIME STEP
fprintf('Starting 30-minute simulation with realistic source model...\n');
tic;  % start timer

for n = 1:Nt-1
    % Progress indicator
    if mod(n, round(Nt/10)) == 0
        fprintf('Progress: %.0f%% (%.1f min simulated)\n', n/Nt*100, time_h(n)*60);
    end
    
    % Adaptive time step: small for transient, large for steady
    current_time_sec = time_h(n) * 3600;
    if current_time_sec < 5.0  % first 5 seconds: small time step
        dt_sec = dt_sec_transient;
    else  % after 5 seconds: larger time step
        dt_sec = dt_sec_steady;
    end
    
    % 6.1 MOC with proper gas dynamics
    % Use sound speed with proper Z-factor
    c = sqrt(R * T * Z_base);  % constant sound speed, m/s
    
    % 6.2 Compute velocity and Reynolds number for friction
    V = m(:,n) ./ (rho(:,n) .* A);  % velocity in m/s
    Re = abs(V) * D / nu;  % Reynolds number
    Re = max(Re, 2300);    % minimum for turbulent flow
    lambda = lambda_f(Re); % friction factor
    
    % 6.3 Foot points for characteristics
    xp = x - c * dt_sec;  % C+ characteristic
    xm = x + c * dt_sec;  % C- characteristic
    
    % 6.4 Handle out-of-bounds foot points properly
    % Initialize arrays for pressure and mass flow at foot points
    p_p = zeros(size(x));
    p_m = zeros(size(x));
    m_p = zeros(size(x));
    m_m = zeros(size(x));
    
    for i = 1:Nx
        % Handle C+ characteristic foot point
        if xp(i) < x(1)
            % Use inlet boundary values
            p_p(i) = p(1,n);
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
            % Use inlet boundary values
            p_m(i) = p(1,n);
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
    
    % 6.5 Update interior points using characteristic equations with friction
    for i = 2:Nx-1
        rho_c = rho(i,n) * c;
        
        % Friction term: λ|V|V/(2D) * ρc * dt
        V_i = V(i);
        friction_term = lambda(i) * abs(V_i) * V_i / (2 * D) * rho_c * dt_sec;
        
        % Characteristic equations
        p_new = 0.5 * (p_p(i) + p_m(i)) + 0.5 * rho_c * (m_p(i) - m_m(i)) / A - friction_term;
        m_new = 0.5 * (m_p(i) + m_m(i)) + 0.5 * A * (p_p(i) - p_m(i)) / rho_c - friction_term * A / rho_c;
        
        p(i,n+1) = p_new;
        m(i,n+1) = m_new;
    end
    
    % 6.6 CORRECTED: Physically realistic inlet boundary condition
    % Source model: P_source - P_inlet = R_source * Q_inlet^2 + L_source * dQ/dt
    % Combined with C- characteristic from pipeline
    
    % Get characteristic information from pipeline
     if xm(1) <= x(end) && xm(1) >= x(1)
         % Interpolate from interior
        m_char = interp1(x, m(:,n), xm(1), 'linear');
        p_char = interp1(x, p(:,n), xm(1), 'linear');
     else
         % Use current values
        m_char = m(1,n);
        p_char = p(1,n);
    end
    
         % Source pressure buffer dynamics (represents source capacitance)
     % Make buffer more responsive to changes
     tau_source = 5.0;  % time constant in seconds (reduced for faster response)
     p_source_buffer(n+1) = p_source_buffer(n) + (p_source_vec(n+1) - p_source_buffer(n)) * dt_sec / tau_source;
    
         % CORRECTED: Work with mass flow rate for better consistency
     m_inlet_current = m(1,n);  % kg/s
     
     % Source resistance relationship: P_source_buffer - P_inlet = R_source * Q_inlet^2 + L_source * dm/dt
     % Convert mass flow to volumetric for resistance calculation
     Q_inlet_vol = m_inlet_current / rho(1,n);  % m³/s
     
     % Calculate mass flow rate change (more consistent)
     dm_dt = (m_inlet_current - m_source_prev) / dt_sec;  % kg/s²
     m_source_prev = m_inlet_current;
    
    % Combine source equation with C- characteristic
     rho_c_inlet = rho(1,n) * c;
    
    % From C- characteristic: p_inlet = p_char + (m_char - m_inlet) * rho_c / A
    % From source: p_inlet = p_source_buffer - R_source * Q_inlet^2 - L_source * dQ/dt
    % Solve iteratively or use approximation
    
         % Simplified approach: estimate inlet pressure from source, then adjust mass flow
     % Use volumetric flow for resistance and mass flow change for inductance
     p_inlet_source = p_source_buffer(n+1) - R_source * Q_inlet_vol^2 - L_source * dm_dt / rho(1,n);
    
         % Combine with characteristic using weighted average
     alpha_source = 0.5;  % weight for source influence (reduced for more pipeline interaction)
     p_inlet_combined = alpha_source * p_inlet_source + (1-alpha_source) * p_char;
    
    % Calculate mass flow from C- characteristic
    m_inlet_new = m_char + A * (p_inlet_combined - p_char) / rho_c_inlet;
    
    % Apply bounds and smoothing
    p(1,n+1) = max(min(p_inlet_combined, 1e7), 1e5);  % reasonable bounds
    m(1,n+1) = max(min(m_inlet_new, 20), 0.1);        % reasonable bounds
    
         % 6.7 IMPROVED: Outlet boundary condition - allow wave propagation with soft pressure control
     % Use C+ characteristic to calculate outlet conditions
     if xp(end) >= x(1) && xp(end) <= x(end)
         % Interpolate from interior using C+ characteristic
         m_outlet_char = interp1(x, m(:,n), xp(end), 'linear');
         p_outlet_char = interp1(x, p(:,n), xp(end), 'linear');
         
         % Calculate outlet conditions using C+ characteristic
         rho_c_outlet = rho(end,n) * c;
         
                  % CORRECTED: Consumer boundary condition with resistance characteristic
          % Consumer model: P_outlet = P_consumer_base - R_consumer * Q_outlet^2
          
          % Calculate outlet volumetric flow
          Q_outlet_vol = m_outlet_char / rho(end,n);  % m³/s
          
          % Consumer pressure with flow-dependent drop
          p_consumer_actual = p_consumer_base - R_consumer * Q_outlet_vol^2;
          
          % Combine consumer demand with wave propagation
          alpha_consumer = 0.3;  % weight for consumer influence (allows more wave dynamics)
          p_outlet_combined = alpha_consumer * p_consumer_actual + (1-alpha_consumer) * p_outlet_char;
          
          p(end,n+1) = p_outlet_combined;
          
          % Calculate mass flow from combined pressure and characteristics
          m(end,n+1) = m_outlet_char + A * (p_outlet_char - p(end,n+1)) / rho_c_outlet;
         
     else
                   % Before wave arrival: use consumer characteristic
          if n == 1
              p(end,n+1) = p_consumer_base;
         m(end,n+1) = m0;
          else
              % Use consumer model with current flow
              Q_outlet_current = m(end,n) / rho(end,n);
              p_consumer_current = p_consumer_base - R_consumer * Q_outlet_current^2;
              p(end,n+1) = p_consumer_current;
              m(end,n+1) = m(end,n);
          end
     end
    
    % 6.8 Update Z-factor and density
    Z(:,n+1) = Z_base * ones(Nx, 1);  % constant Z-factor
    % CORRECTED: Use proper density formula with specific gas constant
    rho(:,n+1) = p(:,n+1) ./ (Z(:,n+1) * R * T);  % kg/m³
    
    % 6.9 Final bounds checking (prevent negative/zero pressure)
     p(:,n+1) = max(p(:,n+1), 1e5);   % Minimum 0.1 MPa
     p(:,n+1) = min(p(:,n+1), 1e7);   % Maximum 10 MPa
    m(:,n+1) = max(m(:,n+1), 0.1);   % Minimum positive flow
     m(:,n+1) = min(m(:,n+1), 20);    % Reasonable maximum mass flow
     
     % Check for numerical instability
     if any(isnan(p(:,n+1))) || any(isinf(p(:,n+1)))
         fprintf('Warning: NaN/Inf detected in pressure at time step %d\n', n);
         break;
     end
     if any(isnan(m(:,n+1))) || any(isinf(m(:,n+1)))
         fprintf('Warning: NaN/Inf detected in mass flow at time step %d\n', n);
         break;
     end
end

elapsed_time = toc;
fprintf('Simulation completed in %.1f seconds (%.1f minutes)\n', elapsed_time, elapsed_time/60);

%% 7. RTTM leak detection calculations
% Calculate volumetric flow rates
Q_inlet = m(1,:) ./ rho(1,:) * 3600;      % m³/h at inlet
Q_outlet = m(end,:) ./ rho(end,:) * 3600;  % m³/h at outlet

% Calculate gas inventory in pipeline
K_c = zeros(1, Nt);
for n = 1:Nt
    rho_current = rho(:,n);
    K_c(n) = trapz(x, rho_current) * A;  % kg total gas mass in pipeline
end

% Calculate derivative of gas inventory dK_c/dt
dK_c_dt = zeros(1, Nt);
for n = 2:Nt-1
    dt_actual = (time_h(n+1) - time_h(n-1)) * 3600;
    dK_c_dt(n) = (K_c(n+1) - K_c(n-1)) / dt_actual;
end
dK_c_dt(1) = (K_c(2) - K_c(1)) / ((time_h(2) - time_h(1)) * 3600);
dK_c_dt(end) = (K_c(end) - K_c(end-1)) / ((time_h(end) - time_h(end-1)) * 3600);

% Convert to volumetric units
rho_avg = mean(rho(:,:), 1);
dK_c_dt_vol = dK_c_dt ./ rho_avg * 3600;  % m³/h

% Calculate instantaneous volumetric imbalance
d_imbalance = Q_inlet - Q_outlet - dK_c_dt_vol;  % m³/h

% Wave propagation analysis
c_sound = sqrt(R * T * Z_base);  % m/s
t_theory = x / c_sound;
wave_arrival_times = zeros(Nx, 1);
wave_detected = false(Nx, 1);
pressure_change_threshold = 5000;  % Pa

for i = 1:Nx
    for n = 2:Nt
        if abs(p(i,n) - p(i,n-1)) > pressure_change_threshold && ~wave_detected(i)
            wave_arrival_times(i) = time_h(n) * 3600;
            wave_detected(i) = true;
            break;
        end
    end
end

% Simplified leak detection logic for no-leak model
mean_imbalance = mean(d_imbalance);
std_imbalance = std(d_imbalance);

% Use realistic leak threshold - 5% of average flow or 3-sigma rule, whichever is larger
flow_based_threshold = 0.05 * mean(Q_inlet);  % 5% of average inlet flow
statistical_threshold = 3 * std_imbalance;    % 3-sigma statistical threshold
leak_threshold = max(flow_based_threshold, statistical_threshold);

% Exclude initial transient period (first 5 minutes) from leak detection
transient_mask = time_h * 60 > 5;  % after 5 minutes
steady_imbalance = d_imbalance(transient_mask);

% Simple leak detection: persistent imbalance above threshold
if ~isempty(steady_imbalance)
    persistent_leak = abs(steady_imbalance) > leak_threshold;
    total_leak_time = sum(persistent_leak) / length(persistent_leak) * 100;
    leak_confirmed = total_leak_time > 10;  % leak confirmed if present >10% of steady time
else
    leak_confirmed = false;
    total_leak_time = 0;
end

% 8. Physics validation
fprintf('\n=== LEAK DETECTION ANALYSIS (RTTM Method) ===\n');
fprintf('Average volumetric imbalance: %.2f m³/h\n', mean_imbalance);
fprintf('Standard deviation of imbalance: %.2f m³/h\n', std_imbalance);
fprintf('Leak threshold: %.2f m³/h\n', leak_threshold);
fprintf('Time with imbalance > threshold: %.1f%% of steady operation\n', total_leak_time);

% Wave validation
fprintf('\n=== WAVE PROPAGATION VALIDATION ===\n');
fprintf('Theoretical wave speed: %.0f m/s\n', c_sound);
fprintf('Expected transit time: %.2f s\n', L/c_sound);
detected_waves = sum(wave_detected);
fprintf('Waves detected at %d/%d monitoring points\n', detected_waves, Nx);

% Final leak assessment
if leak_confirmed
    fprintf('\n*** LEAK DETECTED: Persistent imbalance exceeds threshold ***\n');
else
    fprintf('\n*** NO LEAK DETECTED: System operating normally ***\n');
end

% Remove any NaN/Inf values before plotting
d_imbalance(isnan(d_imbalance) | isinf(d_imbalance)) = 0;

%% LEAK DETECTION VISUALIZATION (RTTM Method)

% DIAGNOSTIC Figure: K_c and Pressure at Boundaries
figure('Position',[50,800,900,400]);
subplot(1,2,1);
plot(time_h*60, K_c, 'LineWidth',2, 'Color', 'blue');
grid on;
xlabel('t, min'); ylabel('K_c, kg');
title('Gas Inventory K_c(t)');

subplot(1,2,2);
plot(time_h*60, p(1,:)/1e6, 'LineWidth',1.5, 'DisplayName', 'Inlet'); hold on;
plot(time_h*60, p(end,:)/1e6, 'LineWidth',1.5, 'DisplayName', 'Outlet');
grid on; legend('Location', 'best');
xlabel('t, min'); ylabel('P, MPa');
title('Pressure at Boundaries');

% Figure 1: Primary flow rates and dK_c/dt
figure('Position',[100,700,900,500]);
subplot(2,1,1);
plot(time_h*60, Q_inlet, 'LineWidth',1.5, 'DisplayName', 'Q_{inlet}(t)'); 
hold on;
plot(time_h*60, Q_outlet, 'LineWidth',1.5, 'DisplayName', 'Q_{outlet}(t)');
grid on; legend('Location', 'best');
xlabel('t, min'); ylabel('Q, m³/h');
title('Volumetric Flow Rates Q_{inlet}(t) and Q_{outlet}(t)');

subplot(2,1,2);
plot(time_h*60, dK_c_dt_vol, 'LineWidth',1.5, 'Color', 'blue');
grid on;
xlabel('t, min'); ylabel('dK_c/dt, m³/h');
title('Gas Inventory Change Rate dK_c/dt');

% Figure 2: Imbalance analysis
figure('Position',[100,100,900,500]);
subplot(2,1,1);
plot(time_h*60, d_imbalance, 'LineWidth',1.5, 'Color', 'black');
hold on;
xlim_vals = xlim;
plot(xlim_vals, [leak_threshold leak_threshold], 'r--', 'LineWidth',2, 'DisplayName', 'Leak Threshold');
plot(xlim_vals, [-leak_threshold -leak_threshold], 'r--', 'LineWidth',2);
grid on; legend('Location', 'best');
xlabel('t, min'); ylabel('d(t), m³/h');
title('Instantaneous Volumetric Imbalance d(t) = Q_{in} - Q_{out} - dK_c/dt');

subplot(2,1,2);
plot(time_h*60, Q_inlet - Q_outlet, 'LineWidth',1.5, 'DisplayName', 'Q_{in} - Q_{out}');
hold on;
plot(time_h*60, dK_c_dt_vol, 'LineWidth',1.5, 'DisplayName', 'dK_c/dt');
grid on; legend('Location', 'best');
xlabel('t, min'); ylabel('Flow, m³/h');
title('Balance Components: Flow Difference vs Inventory Change');

% Figure 3: Pressure distribution p(x,t)
figure('Position',[1000,700,800,500]);
[Xg, Tg] = meshgrid(x, time_h*60);
Psurf = (p.' / 1e6);  % MPa
surf(Xg, Tg, Psurf, 'EdgeColor', 'none');
shading interp; colorbar;
xlabel('L, m'); ylabel('t, min'); zlabel('P, MPa');
title('Pressure Distribution p(x,t) - for leak localization');
view(45,30);

% Figure 4: Wave validation and pressure profiles (removed pressure change rates)
figure('Position',[1000,100,900,500]);
subplot(2,1,1);
detected_points = find(wave_detected);
if ~isempty(detected_points)
    plot(x(detected_points), wave_arrival_times(detected_points), 'ro', 'MarkerSize', 8, 'DisplayName', 'Detected waves');
    hold on;
    plot(x, t_theory, 'b-', 'LineWidth',2, 'DisplayName', 'Theoretical arrival');
    grid on; legend('Location', 'best');
    xlabel('L, m'); ylabel('t, s');
    title('Wave Arrival Time Validation');
end

subplot(2,1,2);
quarter_node = round(Nx/4);
mid_node = round(Nx/2);
three_quarter_node = round(3*Nx/4);

plot(time_h*60, p(1,:)/1e6, 'LineWidth',1.5, 'DisplayName', 'Inlet (x=0)'); hold on;
plot(time_h*60, p(quarter_node,:)/1e6, 'LineWidth',1.5, 'DisplayName', sprintf('x=%.0f m', x(quarter_node)));
plot(time_h*60, p(mid_node,:)/1e6, 'LineWidth',1.5, 'DisplayName', sprintf('x=%.0f m', x(mid_node)));
plot(time_h*60, p(three_quarter_node,:)/1e6, 'LineWidth',1.5, 'DisplayName', sprintf('x=%.0f m', x(three_quarter_node)));
plot(time_h*60, p(end,:)/1e6, 'LineWidth',1.5, 'DisplayName', 'Outlet (x=L)');
grid on; legend('Location', 'best');
xlabel('t, min'); ylabel('P, MPa');
title('Pressure Wave Profiles (for Negative Wave Method)');

fprintf('\n=== VISUALIZATION COMPLETE ===\n');
fprintf('Generated 4 figures for RTTM leak detection analysis\n');

