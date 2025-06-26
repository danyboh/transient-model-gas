% characteristic_method.m
% ����糺����� ��� ���� � ������������� ����������� ������� �������������
% �� ����� [11] �� ��������� ������������

close all; clear all; clc

%% 1. ��������� ����� �� ����
D     = 0.5;                   % ������ �����, �
A     = pi*D^2/4;              % ����� �������, �^2
L     = 1000;                  % ������� �����, �
R     = 518;                   % ������ �����, ��/(��K)
T     = 288;                   % �����������, K
M     = 16e-3;                 % ������� ����, ��/����
alpha = 0.02;                  % ���������� Z(p), 1/���
beta  = 1.0;                   % ���� Z(p)
nu    = 1.2e-5;                % ���������� �������, �^2/�
lambda_f = @(Re) 0.3164*Re.^(-0.25);  % ������� ������������

%% 2. ���������� �� ������ ����
Nx    = 101;                   
x     = linspace(0, L, Nx)';   % ������ x �� ��������
dx    = L/(Nx-1);

% ���������� dt �� ������ CFL
p_max = 6;                     
c_max = sqrt(R*T*(alpha*p_max+beta)/M);
dt    = 0.5 * dx / c_max;      % ��������� ���� �� �����
t_end = 50;                    % ������� ��� �����������, �
Nt    = ceil(t_end/dt);
t     = (0:Nt-1)*dt;           % ������ ����

%% 3. ����������� ����
p   = zeros(Nx, Nt);   % ����, ���
m   = zeros(Nx, Nt);   % ������� ����, ��/�
Z   = zeros(Nx, Nt);   % ������ ����������
rho = zeros(Nx, Nt);   % �������, ��/�^3

% ������� �� �������� �����
p_in   = 6;    % ������� ����, ���
p_out  = 3;    % �������� ����, ���
m0     = 100;  % ���������� ������� ����, ��/�

% ���������� ������� ����� ����� �� inlet � outlet
p(:,1)   = p_in - (p_in - p_out)*(x/L);
m(:,1)   = m0;
Z(:,1)   = alpha*p(:,1) + beta;
rho(:,1) = p(:,1).*M./(Z(:,1)*R*T);

%% 4. �������� ���� ���� (Method of Characteristics)
for n = 1:Nt-1
    % 4.1 ������� �������� �������������
    c       = sqrt(R*T.*Z(:,n)./M);
    lambda1 =  c;  % ������ C+
    lambda2 = -c;  % ������ C-
    
    % 4.2 ������ �� �����
    Re = abs(m(:,n))./(rho(:,n)*A) * D ./ nu;
    fr = lambda_f(Re);
    J  = fr .* R*T .* Z(:,n) .* m(:,n) .* abs(m(:,n)) ./ (2*D*A*M*p(:,n));
    
    % 4.3 foot points
    xp = x - lambda1*dt;
    xm = x - lambda2*dt;
    xp = real(xp);  xm = real(xm);
    xp = min(max(xp, x(1)), x(end));
    xm = min(max(xm, x(1)), x(end));
    
    % 4.4 ������������ m �� Z � ������ xp, xm
    m_p = interp1(x, m(:,n), xp, 'linear');
    m_m = interp1(x, m(:,n), xm, 'linear');
    Z_p = interp1(x, Z(:,n), xp, 'linear');
    Z_m = interp1(x, Z(:,n), xm, 'linear');
    
    % 4.5 ��������� m �� Z
    m(:,n+1) = 0.5*(m_p + m_m) - J*dt;
    Z(:,n+1) = 0.5*(Z_p + Z_m);
    
    % 4.6 ������� �����
    m(end,   n+1) = m0;        
    p(1,     n+1) = p_in;      
    Z(1,     n+1) = alpha*p_in + beta;
    
    % 4.7 ���������� p �� rho
    p(:, n+1)   = m(:,n+1) .* R*T ./ (A .* Z(:,n+1));
    rho(:, n+1) = p(:,n+1) .* M ./ (Z(:,n+1)*R*T);
end

%% 5. ϳ�������� ����� ��� �������
p_initial    = p(:,1);
time_vec     = t/3600;                       % ������
p_inlet_vec  = p(1, :);
q_std_vec    = m(1, :) ./ rho(1, :) * 3600;  % �^3/���
p_outlet_vec = p(end, :);
q_out_vec    = m(end, :) ./ rho(end, :) * 3600;
Z_vec        = Z(end, :);
p_surface    = p;                            % Nx?Nt

%% === Plots ===
figure('Position', [100, 600, 700, 420]);
plot(x, p_initial, 'LineWidth', 1.2); grid on;
xlabel("L, m"); ylabel("P, MPa");
title("���������� ������� �����");

figure('Position', [900, 600, 700, 420]);
subplot(2,1,1);
plot(time_vec, p_inlet_vec, 'LineWidth', 1.2); grid on;
ylabel("P, MPa");
title("������� ����");

subplot(2,1,2);
plot(time_vec, q_std_vec, 'LineWidth', 1.2); grid on;
xlabel("t, hours");
ylabel("Q, m^3/hours");
title("������ �ᒺ��� �������");

figure('Position', [100, 100, 700, 600]);
subplot(3,1,1);
plot(time_vec, p_outlet_vec, 'LineWidth', 1.2); grid on;
ylabel("P, MPa");
title("�������� ����");

subplot(3,1,2);
plot(time_vec, q_out_vec, 'LineWidth', 1.2); grid on;
ylabel("Q, m^3/hours");
title("������� �ᒺ��� �������");

subplot(3,1,3);
plot(time_vec, Z_vec, 'LineWidth', 1.2); grid on;
xlabel("t, hours");
ylabel("Z");
title("������ ���������� �� �����");

% 3D �������� ����� p(x,t)
[XGrid, TGrid] = meshgrid(x, t/3600);   % TGrid � �������
Psurf = p.';                           % Nt?Nx

figure('Position', [900, 100, 800, 500]);
surf(XGrid, TGrid, Psurf, 'EdgeColor', 'none');
xlabel("L, m"); ylabel("t, hours"); zlabel("P, MPa");
title("�������� ���� ����� ���� �� ������ ����������� � ��� ��� ��������������� ����� ���� ����");
view(45, 30);
colorbar;
shading interp;
