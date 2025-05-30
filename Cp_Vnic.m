
function [Cp, rho, Mm] = Cp_Vnic(p, T, x)
% Psevdofunktsiia dlia rozrakhunku teploemnosti pryrodnoho hazu vidpovidno do HOST 30319.3-96
% Na tsʹomu etapi povertaiemo sproshchene Cp yak funktsiiu temperatury

% Vkhid:
% p  - tysʹk [Pa]
% T  - temperatura [K]
% x  - vektor molʹnykh chastok komponentiv hazu

% Vykhid:
% Cp - izobarna teploemnistʹ [J/kg·K]
% rho - hustyna (ne vykorystovuietʹsia na tsʹomu etapi)
% Mm - moliarna masa

% Moliarni masy komponentiv
M = [16.043, 30.07, 44.097, 58.12, 58.12, 28.01, 44.01, 34.08];  % [g/mol]
Mm = sum(x .* M); % serednia moliarna masa, g/mol
Mm = Mm / 1000;   % [kg/mol]

% Hazova stala
R = 8.31451;  % J/(mol·K)

% Sproshchena modelʹ Cp(T)
Cp_mol = 35 + 0.1 * (T - 273); % [J/(mol·K)] — nablyzheno dlia demonstratsii
Cp = Cp_mol / Mm; % [J/(kg·K)]

rho = p ./ (R / Mm * T); % Idealʹne nablyzhennia hustyny
end
