function [u, e] = recover_primitives(U)
% ³��������� �������� u �� ������ ����㳿 e � ������� U

rho = U(1,:);
rhou = U(2,:);
rhoe = U(3,:);

u = rhou ./ rho;
e = rhoe ./ rho;

end
