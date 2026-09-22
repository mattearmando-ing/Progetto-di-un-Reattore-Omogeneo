clear; close all; clc
% Intervallo del rapporto di moderazione
M = 1:1:200;
% === Costanti nucleari (JEFF-3.3) ===
% Sezioni d’urto in cm^2 (1 barn = 1e-24 cm^2)
sigma_f_U235 = 585e-24;
sigma_g_U235 = 99e-24;
sigma_a_U235 = sigma_f_U235 + sigma_g_U235;
sigma_s_U235 = 10e-24;
sigma_a_U238 = 2.85e-24;
sigma_s_U238 = 8.1e-24;
sigma_a_H = 0.33e-24;
sigma_s_H = 20e-24;
sigma_a_S = 0.53e-24; % assorbimento
sigma_s_S = 2.6e-24; % scattering
xi_H = 1.000; % idrogeno
xi_O = 0.122; % ossigeno
xi_U235= 0.0084; % U-235
xi_U238 = 0.02; % valore tipico
% Costanti empiriche JEFF-3.3 per U238 in acqua leggera
a = 0.123; % coefficiente empirico integrale risonanza
c = 0.5; % esponente empirico integrale risonanza
% === Densità atomiche ===
NA = 6.022e23;
enr = [0.02, 0.05, 0.1, 0.1975]; % arricchimento in U-235 (19.75%)
MW_U235 = 235;
MW_U238 = 238;
MM_H2O = 18.015; % g/mol
rho_H2O = 0.9778;
figure(1);
yline(1, '--r', 'k_{\infty} = 1');
hold on
for i=1:length(enr)
 % n_H [atomi/cm³]
 n_H2O = rho_H2O / MM_H2O * NA;
 n_H = 2 *n_H2O;
 % n_U [atomi/cm³]
 n_F = n_H./M;
 % n_U235 [atomi/cm³]
 n_U235 = n_F/8 * enr(i); %perchè UO2SO4

 % n_U238 [atomi/cm³]
 n_U238 = n_F/8 * (1 - enr(i));

 % rho [gUO2SO4/L]
 MM_UO2SO4 = 235 * enr(i) + 235 * (1 -enr(i)) + 32.065 + 6 * 15.999;
 rho = n_F .* MM_UO2SO4 ./ NA *1000;
 % === Eta ===
 nu = 2.43;
 eta = nu * sigma_f_U235*n_U235 /
(1.4*(sigma_a_U235*n_U235+sigma_a_U238*n_U238));
 % -presente fattore correttivo

 % === Epsilon ===
 epsilon = (1 + 0.690 * n_U238 ./ n_H) ...
 ./ (1 + 0.563 * n_U238 ./ n_H);

 % === p ===
 Sigma_p = sigma_s_U235 .* n_U235 + sigma_s_U238 .* n_U238 + sigma_s_H .* n_H;
 I = a .* (Sigma_p ./ n_U238).^c;
 xi_bar = (xi_H .* sigma_s_H .* M + xi_U235 .* sigma_s_U235) ...
 ./ (sigma_s_H .* M + sigma_s_U235);
 %p = exp(- n_U238 .* I ./ (xi_bar .* Sigma_p));
 p =
(sigma_s_H.*M)./(sigma_f_U235./(1+1/enr(i))+sigma_a_U235./(1+1/enr(i))+sigma_a_U2
38./(enr(i)+1)+sigma_s_H.*M);
 % === f ===
 f =
(sigma_a_U235./(1+1/enr(i))+sigma_a_U238./(enr(i)+1))./(sigma_a_U235./(1+1/enr(i))+
sigma_a_U238./(enr(i)+1)+sigma_a_H.*M);

 % === k_inf ===
 k_inf = eta .* epsilon .* p .* f;

 if i==4
 plot(M, k_inf, 'k', 'LineWidth', 2);
 hold on
 else
 plot(M, k_inf,'--', 'LineWidth', 1);
 hold on
 xlabel('Rapporto di moderazione \theta (H/U)');
 ylabel('k_{\infty}');
 title('k_{\infty}(M)');
 grid on;
 end
end
% legend('Arricchimento U')
h = legend('k_{\infty}','2%','5%','10%','19,75%');
title(h, 'Arricchimento uranio (%)','FontSize', 14);
set(h, 'FontSize', 12);
legend('Location','southwest')
set(gca, 'FontSize', 14);
% === Grafici ===
% Fattori
figure;
subplot(2,2,1)
plot(M, eta * ones(size(M)), 'b', 'LineWidth', 2);
title('\eta (riproduzione termica)'); xlabel('\theta (H/U)'); grid on; ylim([0 3]);
set(gca, 'FontSize', 18);
subplot(2,2,2)
plot(M, epsilon, 'g', 'LineWidth', 2);
title('\epsilon (fissione veloce)'); xlabel('\theta (H/U)'); grid on;
set(gca, 'FontSize', 18);
subplot(2,2,3)
plot(M, p, 'r', 'LineWidth', 2);
title('p (non assorbimento in moderazione)'); xlabel('\theta (H/U)'); grid on; ylim([0 1]);
set(gca, 'FontSize', 18);
subplot(2,2,4)
plot(M, f, 'm', 'LineWidth', 2);
title('f (utilizzazione termica)'); xlabel('\theta (H/U)'); grid on; ylim([0 1]);
set(gca, 'FontSize', 18);
% k_inf finale
figure;
plot(M, k_inf, 'b', 'LineWidth', 2);
ylim([0 1.5])
hold on
xlabel('Rapporto di moderazione \theta (H/U)');
ylabel('k_{\infty}');
title('k_{\infty}(M)');
yline(1, '--k', 'k_{\infty} = 1');
grid on;
[max_value, index] = max(k_inf);
fprintf('k_inf_max = %.5f \n', max_value);
fprintf('Theta_max = %.2f \n', index);
plot(index,max_value, 'ro', 'MarkerFaceColor', 'r','MarkerSize',12); % Punto rosso
text(65, 0.8, sprintf(' M = %.2f\n k_{\\infty} = %.5f', index, max_value), ...
 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 18);
xline(index, '--r');
plot(50,k_inf(50), 'ro', 'MarkerFaceColor', 'g','MarkerSize',12); % Punto verde
text(25, 1.15, sprintf(' M = %.2f\n k_{\\infty} = %.5f', 50, k_inf(50)), ...
 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 18);
xline(50, '--g');
set(gca, 'FontSize', 18);
sigma_tot=n_U235*sigma_a_U235+n_H*sigma_a_H+n_U238*sigma_a_U238;
%% quantità veleno
% Dati geometrici
H = 33.6778; % cm
R = 16.8389; % cm
% kinf
kinf = k_inf(50);
M_ottimale = 50;
n_F = n_H / M_ottimale;
enr = 0.1975;
n_U235_0 = (n_F/8) * enr;
n_U238_0 = (n_F/8) * (1 - enr);
% Calcolo Sigma_a0 (senza veleno)
Sigma_a0 = sigma_a_U235*n_U235_0 + sigma_a_U238*n_U238_0 + sigma_a_H*n_H; %
cm^-1
% Coefficiente di diffusione neutronica D
% Ipotizzo un valore tipico realistico per reattore omogeneo ad acqua leggera
D = 0.8; % cm (valore tipico realistico, generalmente tra 0.5 e 1 cm)
% Buckling geometrico B^2 per cilindro
B2 = (pi/H)^2 + (2.4048/R)^2; % cm^-2
% Calcolo densità atomica boro (B-10)
sigma_a_B10 = 3837e-24; % cm² (B-10)
N_B = (D*B2/(kinf-1) - Sigma_a0) / sigma_a_B10; % atomi/cm³
% Concentrazione acido borico H3BO3 (B-10 è circa il 20% del boro naurale, H3BO3
contiene il 17.5% di boro naturale)
ppm_H3BO3 = (N_B/6.022e23)*10*1e3/(0.2*0.175); %(g/L)
% Output risultati
fprintf('Concentrazione di acido borico necessaria = %.2f g/L \n', ppm_H3BO3);
%% Andamento U-235
% === Parametri di input ===
V = pi * R^2 * H; % Volume cilindro in cm^3
% Potenza
P = 2e5; % Watt (200 kW)
% Energia per fissione
E_fiss = 200 * 1.602e-13; % J (200 MeV)
% Calcolo del flusso medio usando distribuzione armonica fondamentale
phi = P / (E_fiss * sigma_f_U235 * n_U235_0 * V); % [neutroni/cm^2/s]
phi_0 = 2 *phi;
% Tempo di simulazione (in secondi)
t_max = 5e8; % ~30 anni
N = 1000;
t = linspace(0, t_max, N);
% Evoluzione di N_F(t)
N_F = n_U235_0 * exp(-sigma_a_U235 * phi * t);
% Calcolo di k_inf(t) = costante * N_F(t)
% ipotizziamo k_inf(t) proporzionale alla densità fissile (lineare)
k_inf_t = kinf * (N_F / n_U235_0);
% Trova tempo in cui k_inf scende a 1
idx = find(k_inf_t <= 1, 1);
if isempty(idx)
 fprintf('k_inf non raggiunge mai 1 nel tempo simulato.\n');
 t_critico = NaN;
else
 t_critico = t(idx);
 fprintf('Tempo in cui k_inf scende a 1: %.2f giorni\n', t_critico / 86400);
end
% === Plot ===
figure;
subplot(2,1,1);
plot(t/86400, N_F, 'LineWidth', 2);
xlabel('Tempo [giorni]'); ylabel('N_F [atomi/cm^3]');
title('Evoluzione della densità del combustibile'); grid on;
subplot(2,1,2);
plot(t/86400, k_inf_t, 'LineWidth', 2);
xlabel('Tempo [giorni]'); ylabel('k_{\infty}(t)');
title('Andamento di k_{\infty} nel tempo'); grid on;
yline(1,'--r','k_{\infty}=1');
xlim([0 100])
%% distribuzione flusso
% Parametri
alpha01 = 2.4048; % Primo zero della Bessel J0 % Flusso massimo
arbitrario
N = 200; % Numero di punti
theta = linspace(0, 2*pi, N);
r = linspace(0, R, N);
[Theta, Rad] = meshgrid(theta, r);
% Coordinate cartesiane per il disco
X = Rad .* cos(Theta);
Y = Rad .* sin(Theta);
% Scegli piano z (max flusso a z=0)
z_plot = 0;
Phi_disk = phi_0 * besselj(0, alpha01 * Rad / R) * cos(pi * z_plot / H);
% Distribuzione radiale lungo il diametro (r da -R a +R)
r_line = linspace(-R, R, N);
Phi_line = phi_0 * besselj(0, alpha01 * abs(r_line) / R); % abs per simmetria
% === Subplot ===
figure;
% 1. Plot su disco
subplot(1,2,1);
surf(X, Y, Phi_disk, 'EdgeColor', 'none');
xlabel('x [cm]'); ylabel('y [cm]'); zlabel('\phi(r,0)');
title('Distribuzione radiale (disco)');
colorbar;
axis equal; view(2); % Vista dall'alto
colormap('hot'); % Cambia qui la scala colori (vedi sotto per alternative)
set(gca, 'FontSize', 12);
% 2. Plot cartesiano
subplot(1,2,2);
plot(r_line, Phi_line, 'r', 'LineWidth', 4);
xlabel('r [cm]');
ylabel('\phi(r,0)');
title('Sezione radiale lungo il diametro');
grid on;
xlim([-R, R]);
sgtitle('Flusso neutronico armonica fondamentale (cilindro)');
set(gca, 'FontSize', 12);
%% radionuclidi modo 1
% === Parametri della pastiglia ===
V_p = 1; % volume pastiglia [cm^3]
rho_Mo = 10.2; % densita' molibdeno [g/cm^3]
M_Mo = 95.94; % massa molare molibdeno [g/mol]
frac_98 = 0.24; % frazione isotopica 98Mo
% Calcolo densita' atomica iniziale di 98Mo nella pastiglia
n_Mo_total = (rho_Mo * V_p / M_Mo) * NA; % atomi totali di Mo nella pastiglia
n_Mo98_0 = frac_98 * n_Mo_total / V_p; % [atomi/cm^3] di Mo-98
% === Parametri di reazione ===
sigma_c_Mo = 0.13e-24; % sezione d'urto (n,gamma) [cm^2]
% === Dati decadimento Mo-99 ===
T12_Mo99 = 65.94 * 3600; % emivita Mo-99 in secondi (65,94 h)
lambda_Mo99 = log(2) / T12_Mo99; % costante decadimento [1/s]
% === Evoluzione temporale ===
t_max = 3e6; % secondi (~115 giorni)
N = 100000;
t = linspace(0, t_max, N);
% === Produzione Mo-99 con decadimento ===
a = sigma_c_Mo * phi_0 + lambda_Mo99;
n_Mo99_1d = n_Mo98_0 * sigma_c_Mo * phi_0 / a .* (1 - exp(-a * t));
% === Calcolo della saturazione ===
N_Mo99_sat = (sigma_c_Mo * n_Mo98_0 * phi_0) / (sigma_c_Mo * phi_0 +
lambda_Mo99);
% === Plot ===
figure;
plot(t/86400, n_Mo99_1d * V_p, 'LineWidth', 2); hold on;
yline(N_Mo99_sat * V_p, '--r', 'Saturazione', 'LineWidth', 2);
xlabel('Tempo [giorni]');
ylabel('Atomi di ^{99}Mo nella pastiglia');
title('Produzione di ^{99}Mo per attivazione di ^{98}Mo (con decadimento)');
grid on;
legend('Produzione', 'Saturazione');
legend('Location','best')
hold off;
%% secondo modo
% === Parametri del combustibile ===
% U-235 puro al 100%
n_U235_100 = (n_F/8);
yield_Mo99 = 0.06; % resa in Mo-99 per fissione (6%)
T12_Mo99 = 65.94 * 3600; % semivita Mo-99 in secondi (65,94 h)
lambda_Mo99 = log(2) / T12_Mo99; % costante di decadimento [1/s]
% === Evoluzione temporale ===
t_max = 5e8; % secondi (~115 giorni)
N = 100000;
t = linspace(0, t_max, N);
N_U235 = n_U235_100 * exp(-sigma_a_U235 * phi * t);
prod_sorgente = yield_Mo99 * sigma_f_U235 * phi * n_U235_100;
a = sigma_a_U235 * phi + lambda_Mo99;
n_Mo99_2d = prod_sorgente/a * (1 - exp(-a*t));
% === Rateo di produzione ===
dn_Mo99_dt = yield_Mo99 * sigma_f_U235 * phi * N_U235;
% === Plot ===
figure;
subplot(2,1,1);
plot(t/86400, n_Mo99_2d, 'LineWidth', 2);
xlabel('Tempo [giorni]');
ylabel('Atomi di ^{99}Mo');
title('Produzione cumulativa di ^{99}Mo da fissione di ^{235}U (con decadimento)');
grid on;
xlim([0 25])
% === Calcolo saturazione ===
N_Mo99_sat_fiss = (yield_Mo99 * sigma_f_U235 * n_U235_100 * phi) / (sigma_a_U235
* phi + lambda_Mo99);
% === Plot esempio (aggiungi a un grafico già esistente, ad es. n_Mo99_2d) ===
hold on
yline(N_Mo99_sat_fiss, '--r', 'Saturazione fissione', 'LineWidth', 2);
legend('Produzione cumulativa', 'Saturazione');
legend('Location','best')
hold off
subplot(2,1,2);
plot(t/86400, dn_Mo99_dt, 'LineWidth', 2);
xlabel('Tempo [giorni]');
ylabel('Rateo di produzione [atomi/s/cm^3]');
title('Rateo istantaneo di produzione di ^{99}Mo');
grid on;
%% confonto
figure;
plot(t/86400, n_Mo99_1d * V_p, 'b-', 'LineWidth', 2); hold on;
plot(t/86400, n_Mo99_2d, 'r', 'LineWidth', 2);
xlabel('Tempo [giorni]');
ylabel('Atomi di ^{99}Mo prodotti');
legend('Attivazione ^{98}Mo (pastiglia)', 'Fissione ^{235}U (combustibile)', 'Location',
'Best');
title('Confronto produzione cumulativa di ^{99}Mo');
grid on;
