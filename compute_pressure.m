function p = compute_pressure(U, T)
% ���������� ����� � ������������� p = Z*rho*R*T

R = 518.3;
Z = 1; % � ���������� ������
rho = U(1,:);
p = Z .* rho .* R .* T;

end
