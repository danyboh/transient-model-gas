function plot_results(U, x)
% �������� ������� ���������� ����������� � �������

rho = U(1,:);
u = U(2,:) ./ rho;
e = U(3,:) ./ rho;
T = (e - 0.5 * u.^2) / (1.5 * 518.3);

% ������ �� ���������� ��� ���������� ����������
invalid_T = T <= 0 | isnan(T) | ~isreal(T);
T(invalid_T) = 273.15;

% ����� ������� p ����� ��� ������� �����
p = rho .* 518.3 .* T;
p(~isfinite(p)) = 1e5;  % ������� NaN/Inf ����� �� ��������� ��������

fprintf('min(p): %.3e, min(T): %.3f, min(rho): %.3f\n', min(p), min(T), min(rho));

% ������������ ������� �� �����
positions = [100 600 560 420; 700 600 560 420; 100 100 560 420; 700 100 560 420; 1300 600 560 420; 1300 100 560 420; 100 900 560 420];
i_fig = 1;

f = figure; set(f, 'Position', positions(i_fig,:)); i_fig = i_fig + 1;
plot(x, p/1e5); ylabel('���� [���]'); xlabel('³������ [�]'); grid on; title('������� �����');

f = figure; set(f, 'Position', positions(i_fig,:)); i_fig = i_fig + 1;
plot(x, T); ylabel('����������� [K]'); xlabel('³������ [�]'); grid on; title('������� �����������');

f = figure; set(f, 'Position', positions(i_fig,:)); i_fig = i_fig + 1;
plot(x, u); ylabel('�������� [�/�]'); xlabel('³������ [�]'); grid on; title('������� ��������');

Z = 1; R = 518.3;
valid = isfinite(p) & isfinite(T) & T > 0;
K_c = trapz(x(valid), p(valid) ./ (Z * R * T(valid)));

f = figure; set(f, 'Position', positions(i_fig,:)); i_fig = i_fig + 1;
plot(x, rho); ylabel('������� [��/�^3]'); xlabel('³������ [�]'); grid on; title('������� �������');

fprintf('����� ���� � ����������� (K_c): %.3f\n', K_c);

if exist('U_store', 'var')
    Nt = size(U_store, 3);
    K_vec = zeros(1, Nt);
    for j = 1:Nt
        Uj = U_store(:,:,j);
        rhoj = Uj(1,:);
        uj = Uj(2,:) ./ rhoj;
        ej = Uj(3,:) ./ rhoj;
        Tj = (ej - 0.5 * uj.^2) / (1.5 * R);
        Tj(Tj <= 0 | isnan(Tj) | ~isreal(Tj)) = 273.15;
        pj = rhoj .* R .* Tj;
        K_vec(j) = trapz(x, pj ./ (Z * R * Tj));
    end

    t_vec = linspace(0, 5, Nt);
    dK_dt = [0, diff(K_vec)./diff(t_vec)];

    f = figure; set(f, 'Position', positions(i_fig,:)); i_fig = i_fig + 1;
    plot(t_vec, K_vec); xlabel('��� [�]'); ylabel('K_c'); grid on; title('����� ���� � ���');

    f = figure; set(f, 'Position', positions(i_fig,:));
    plot(t_vec, dK_dt); xlabel('��� [�]'); ylabel('dK_c/dt'); grid on; title('������� ������ ����');
end

end