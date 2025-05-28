function [u, e] = recover_primitives(U)
% Відновлення швидкості u та питомої енергії e з вектора U

rho = U(1,:);
rhou = U(2,:);
rhoe = U(3,:);

u = rhou ./ rho;
e = rhoe ./ rho;

end
