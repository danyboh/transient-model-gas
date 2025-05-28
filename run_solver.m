function [U, x] = run_solver()
%    ��� ������� ������ �������� ������

% ���������
L = 1000;      % ������� ������������ [�]
Nx = 200;      % ������� ������
T_end = 5;     % ��� ����������� [�]
CFL = 0.5;     % ����� �������
[U, x] = init_conditions([], Nx);
dx = x(2) - x(1);
t = 0;

% ����� ��� ������ ���������� (�������)
store_steps = 100;
U_store = zeros(3, Nx, store_steps);
k = 1;

% ������� ����
while t < T_end
    [u, e] = recover_primitives(U);
    T = e ./ (1.5 * 518.3); % �������� ����� Cv = const
    p = compute_pressure(U, T);
    c = sqrt(1.3 * p ./ U(1,:));
    umax = max(abs(u) + c);
    dt = CFL * dx / umax;
    if t + dt > T_end, dt = T_end - t; end

    % SSPRK2
    RHS1 = compute_rhs(U, p, dx);
    U1 = U + dt * RHS1;
    RHS2 = compute_rhs(U1, p, dx);
    U = 0.5 * (U + U1 + dt * RHS2);

    if mod(k, store_steps/10) == 0
        U_store(:, :, round(k/(store_steps/10))) = U;
    end

    t = t + dt;
    k = k + 1;
end

plot_results(U, x);

end
