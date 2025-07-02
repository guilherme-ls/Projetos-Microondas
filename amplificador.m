%% Script para Projeto de Amplificador de Baixo Ruído (LNA)
% Este script segue a metodologia do RF Toolbox da MathWorks para projetar
% um LNA. Ele determina os pontos de casamento e, em seguida, calcula as
% dimensões físicas para as redes de casamento em microfita.
% A simulação do circuito completo foi removida para focar no design.
%
% NOTA: Os parâmetros de ruído são assumidos para fins de demonstração.

clear;
clc;
close all;

%% 1. Definições Iniciais do Circuito e do Dispositivo
Z0 = 50; % Impedância do sistema (Ohms)
f0 = 10e9; % Frequência central de operação (Hz)
c = physconst('LightSpeed'); % Velocidade da luz no vácuo

% Definição dos materiais para a linha de microfita
met = metal('Name', 'Cobre', 'Conductivity', 5.8e7, ...      
       'Thickness', 17e-6);
diel = dielectric('Name','RT5880','EpsilonR', 2.2, ...
       'LossTangent', 0.0009, 'Thickness', 0.787e-3);
er = 2.2;      % Epsilon Relativo do dielétrico (RT5880)
h = 0.787e-3;  % Altura do dielétrico em metros (RT5880)

% Parâmetros-S do Transistor em 10 GHz
s = sparameters('ATF-13336.s2p');
freqs = s.Frequencies;
[~, idx] = min(abs(freqs - f0));  % índice mais próximo

S11 = s.Parameters(1,1,idx);
S21 = s.Parameters(2,1,idx);
S12 = s.Parameters(1,2,idx);
S22 = s.Parameters(2,2,idx);
S_params = [S11 S12; S21 S22];

% Parâmetros de Ruído em 10 GHz
F_min_dB = 1.3;
Gamma_opt = 0.47 * exp(1i*deg2rad(-96));
Rn = 0.5 * Z0;

% Cria o objeto amplificador para análise
lna_ckt = read(rfckt.amplifier, 'ATF-13336.s2p');
figure;
plot(lna_ckt,'Gmag','Ga','Gt','dB');
figure;
plot(lna_ckt,'Fmin','NF','dB')
figure;
hsm = smithplot;
circle(lna_ckt,f0,'Stab','In','Stab','Out','Ga',8:1:11, 'NF',1.4:0.2:2.6,hsm);
legend('Location','SouthEast')

%% 2. Análise de Estabilidade e Plot dos Círculos de Operação
fprintf('--- ANÁLISE INICIAL A %.1f GHz ---\n', f0/1e9);

delta = S11*S22 - S12*S21;
K = (1 - abs(S11)^2 - abs(S22)^2 + abs(delta)^2) / (2*abs(S12*S21));

fprintf('K = %.4f e delta = %.4f\n', K, delta);
if K > 1 && abs(delta) < 1
    disp('Transistor absolutamente estável.');
else
    disp('Transistor potencialmente instável. Use rede estabilizadora.');
end

% Testando se pode ser tratado como unilateral
U = abs(S11) * abs(S12) * abs(S21) * abs(S22) / ((1 - abs(S11)^2) * (1 - abs(S22)^2));
ULB = 1 / (1 + U^2);
UUB = 1 / (1 - U^2);
fprintf('%.4f < Gt - Gtu < %.4f\n', 10*log10(ULB), 10*log10(UUB));

% Objetivo: Figura de Ruído (NF) <= 2 dB com o maior ganho possível.
F_design_dB = 1.5;
F_design = 10^(F_design_dB/10);
F_min = 10^(F_min_dB/10);

% Calcula os parâmetros do círculo de ruído
N_param = (F_design - F_min) * abs(1 + Gamma_opt)^2 / (4 * Rn / Z0);
Cf = Gamma_opt / (N_param + 1);
Rf = (1 / (N_param + 1)) * sqrt(N_param * (N_param + 1 - abs(Gamma_opt)^2));

% Calcula os parâmetros do círculo de ganho
x_values = 0:0.0001:(2*pi);
Gamma_s_values = Cf + Rf * exp(1i * x_values);
Gain_values = 10*log10((1 - abs(Gamma_s_values).^2)./abs(1 - S11*Gamma_s_values).^2);
[maxGs, idx] = max(Gain_values);
Gamma_s = Gamma_s_values(idx);

% Parametriza círculo Gs
gs = (1 - abs(Gamma_s)^2) / abs(1 - S11 * Gamma_s)^2 * (1 - abs(S11)^2);
Cs = gs * conj(S11) / (1 - (1 - gs) * abs(S11)^2);
Rs = sqrt(1 - gs) * (1 - abs(S11)^2) / (1 - (1 - gs) * abs(S11)^2);

% Calcula o Gamma_l necessário para casamento conjugado na saída
Gamma_out = S22 + (S12*S21*Gamma_s)/(1-S11*Gamma_s);
Gamma_l = conj(Gamma_out);

%% 3. Análise de Ganho e Ruído para o Ponto de Projeto
fprintf('\n--- PONTOS DE PROJETO ESCOLHIDOS ---\n');
fprintf('Coeficiente de Reflexão da Fonte (Γs) = %.3f ∠ %.3f°\n', abs(Gamma_s), rad2deg(angle(Gamma_s)));
fprintf('Coeficiente de Reflexão da Carga (Γl) = %.3f ∠ %.3f°\n', abs(Gamma_l), rad2deg(angle(Gamma_l)));

% Figura de Ruído para o Gamma_s escolhido (deve ser F_design_dB)
NF_design_dB = 10*log10(F_min + (4*Rn/Z0)*abs(Gamma_s-Gamma_opt)^2 / ((1-abs(Gamma_s)^2)*abs(1+Gamma_opt)^2) );

% Ganho Disponível para o Gamma_s escolhido (Fórmula Bilateral Exata)
Ga_design_val = (abs(S21)^2 * (1 - abs(Gamma_s)^2)) / (abs(1 - S11*Gamma_s)^2 * (1 - abs(Gamma_out)^2));
Ga_design_dB = 10*log10(Ga_design_val);

% Ganho do Transistor (G0) e Ganhos das Redes (unilateral)
Gs_u = 10*log10((1 - abs(Gamma_s)^2)/abs(1 - S11*Gamma_s)^2);
G0 = 10*log10(abs(S21)^2);
Gl_u = 10*log10(1 / (1 - abs(S22)^2));
GT_u = Gs_u + G0 + Gl_u;

% Ganho de Transdução Bilateral (Exato)
GT = 10*log10((abs(S21)^2 * (1-abs(Gamma_s)^2) * (1-abs(Gamma_l)^2)) / ...
     (abs((1-S11*Gamma_s)*(1-S22*Gamma_l) - S12*S21*Gamma_l*Gamma_s)^2));
GTu = 10*log10((abs(S21)^2 * (1-abs(Gamma_s)^2) * (1-abs(Gamma_l)^2)) / ...
     (abs(1-S11*Gamma_s)^2 * abs(1-S22*Gamma_l)^2));

fprintf('\n--- ANÁLISE DE GANHO E RUÍDO ---\n');
fprintf('Figura de Ruído Esperada: %.3f dB\n', NF_design_dB);
fprintf('Ganho Disponível (Ga): %.3f dB\n', Ga_design_dB);
fprintf('Ganho da Rede de Entrada (Gs, unilateral): %.3f dB\n', Gs_u);
fprintf('Ganho do Transistor (G0): %.3f dB\n', G0);
fprintf('Ganho da Rede de Saída (Gl, unilateral): %.3f dB\n', Gl_u);
fprintf('Ganho Total (GT, unilateral): %.3f dB\n', GT_u);
fprintf('Ganho Total de Transdução (GT, bilateral exato): %.3f dB\n', GT);


%% 4. Plot dos Círculos de Projeto e Pontos de Casamento
figure;
hsm = smithplot;
title('Círculos de Projeto de LNA e Pontos de Casamento');
hold on;

% Plota o círculo de Ganho Disponível de projeto
%circle(lna_ckt, f0, 'Ga', Ga_design_dB);
% Plota o círculo de Figura de Ruído de projeto
circle(lna_ckt, f0, 'NF', F_design_dB);
xpoints = 0:0.001:(2*pi);
Gs_points = Cs + Rs * exp(1i * xpoints);
plot(real(Gs_points), imag(Gs_points), 'LineStyle','-','Color','r','DisplayName', 'Gs(Freq=10[GHz])');

% Plota os pontos de casamento
plot(real(Gamma_s), imag(Gamma_s), '*', 'MarkerSize', 8, 'Color', 'g', 'DisplayName', 'Γ_s');
plot(real(Gamma_l), imag(Gamma_l), '*', 'MarkerSize', 8, 'Color', 'm', 'DisplayName', 'Γ_l');
legend('show', 'Location', 'northeast');
hold off;

%% 5. Projeto das Redes de Casamento em Microfita
fprintf('\n--- PROJETO DAS REDES DE CASAMENTO EM MICROFITA ---\n');

% Impedâncias alvo para as redes de casamento
% A rede de entrada deve apresentar Zs ao transistor, vindo de Z0.
% A rede de saída deve apresentar Zl ao transistor, vindo de Z0.
Zs = Z0 * (1+Gamma_s)/(1-Gamma_s);
Zl = Z0 * (1+Gamma_l)/(1-Gamma_l);
fprintf('Impedância Requerida na Entrada do Transistor (Zs) = %.3f + %.3fj Ω\n', real(Zs), imag(Zs));
fprintf('Impedância Requerida na Saída do Transistor (Zl) = %.3f + %.3fj Ω\n', real(Zl), imag(Zl));

% Projeto da rede de ENTRADA (transforma Z0 em Zs)
[d_in, l_in, w_in] = design_matcher(Z0, conj(Zs), f0, met, diel);
% Projeto da rede de SAÍDA (transforma Z0 em Zl)
[d_out, l_out, w_out] = design_matcher(Z0, conj(Zl), f0, met, diel);

fprintf('\n-- Rede de Entrada (50Ω -> Zs) --\n');
fprintf('Largura da linha (W): %.4f mm\n', w_in*1000);
fprintf('Distância do transistor ao stub (d): %.4f mm\n', d_in*1000);
fprintf('Comprimento do stub em aberto (l): %.4f mm\n', l_in*1000);

fprintf('\n-- Rede de Saída (50Ω -> Zl) --\n');
fprintf('Largura da linha (W): %.4f mm\n', w_out*1000);
fprintf('Distância do transistor ao stub (d): %.4f mm\n', d_out*1000);
fprintf('Comprimento do stub em aberto (l): %.4f mm\n', l_out*1000);

epsilon_eff = (diel.EpsilonR + 1)/2 + ((diel.EpsilonR - 1)/2) .* (1 ./ sqrt(1 + 12 * diel.Thickness ./ w_in));
d_common = physconst('LightSpeed') / (2 * f0 * sqrt(epsilon_eff));
common_line = rfckt.microstrip("EpsilonR", diel.EpsilonR, "LossTangent", diel.LossTangent, ...
    "Height",diel.Thickness,"Thickness",met.Thickness, ...
    "SigmaCond",met.Conductivity,"LineLength",d_common,"Width",w_in);

input_match = rfckt.microstrip("EpsilonR", diel.EpsilonR, "LossTangent", diel.LossTangent, ...
    "Height",diel.Thickness,"Thickness",met.Thickness,"StubMode","Shunt","Termination","Open", ...
    "SigmaCond",met.Conductivity,"LineLength",l_in,"Width",w_in);
output_match = rfckt.microstrip("EpsilonR", diel.EpsilonR, "LossTangent", diel.LossTangent, ...
    "Height",diel.Thickness,"Thickness",met.Thickness,"StubMode","Shunt","Termination","Open", ...
    "SigmaCond",met.Conductivity,"LineLength",l_out,"Width",w_out);
out_line = rfckt.microstrip("EpsilonR", diel.EpsilonR, "LossTangent", diel.LossTangent, ...
    "Height",diel.Thickness,"Thickness",met.Thickness, ...
    "SigmaCond",met.Conductivity,"LineLength",d_out,"Width",w_in);
in_line = rfckt.microstrip("EpsilonR", diel.EpsilonR, "LossTangent", diel.LossTangent, ...
    "Height",diel.Thickness,"Thickness",met.Thickness, ...
    "SigmaCond",met.Conductivity,"LineLength",d_in,"Width",w_in);
LNA = rfckt.cascade('ckts', {common_line, input_match, in_line, lna_ckt, out_line, output_match, common_line});

% Alternativo

[Xl, Bl] = calc_reactances(Z0, conj(Zl));
Ll = Xl(1) / (2 * pi * f0);
Cl = Bl(1) / (2 * pi * f0);
%Ls = 0.42 * 1e-9;
%Cs = 0.5818 * 1e-12;
[Xs, Bs] = calc_reactances(Z0, conj(Zs));
Ls = Xs(1) / (2 * pi * f0);
Cs = Bs(1) / (2 * pi * f0);
%Ls = 0.418 * 1e-9;
%Cs = 0.525 * 1e-12;
fprintf('Capacitor s: %.4f pF\n', Cs*1e12);
fprintf('Indutor s: %.4f nH\n', Ls*1e9);
fprintf('Capacitor l: %.4f pF\n', Cl*1e12);
fprintf('Indutor l: %.4f nH\n', Ll*1e9);

%input_match = rfckt.cascade('Ckts',{rfckt.shuntrlc('C',Cs),rfckt.seriesrlc('L',Ls)});
%output_match = rfckt.cascade('Ckts',{rfckt.seriesrlc('L',Ll),rfckt.shuntrlc('C',Cl)});
%LNA = rfckt.cascade('Ckts', {input_match,lna_ckt,output_match});

% Fim alternativo

freq = 7e9:10e6:12e9;
analyze(LNA,freq,Z0,Z0,Z0);
figure;
plot(LNA,'Ga','Gt','dB');
figure;
plot(LNA,'NF','dB')
figure;
s = sparameters(LNA);
rfplot(s,1,1);
figure;
rfplot(s,2,1)

% Calcula largura de banda e Q
s21 = rfparam(s, 2, 1);         % extrai S21
mag = 20*log10(abs(s21));       % converte para dB

% Limiar para banda passante (ex: -3 dB)
fprintf('Ganho S21: %.2 dB\n', max(mag));
threshold_pass = -3 + max(mag);
in_band = find(mag > threshold_pass);
f_start = freq(in_band(1));
f_stop  = freq(in_band(end));
BW_pass = f_stop - f_start;
fprintf('Banda útil: %.2f MHz (%.2f GHz – %.2f GHz)\n', ...
    BW_pass/1e6, f_start/1e9, f_stop/1e9);
Q = f0 / BW_pass;
fprintf('Fator de Qualidade (Q): %.2f\n', Q);

%% Funções Auxiliares para o Projeto
% Esta função projeta uma rede de casamento com shunt stub
function [d, l, w] = design_matcher(Z0, Z_target, f0, met, diel)

    RL = real(Z_target);
    XL = imag(Z_target);

    %xpoints = 0:0.001:(2*pi);
    %Gs_points = Cs + Rs * exp(1i * xpoints);
    %plot(real(Gs_points), imag(Gs_points), 'LineStyle','-','Color','r');


    line = microstripLine("Conductor",met,"Substrate",diel,"Height",diel.Thickness);
    line = design(line,f0);
    w = line.Width;
    epsilon_eff = (diel.EpsilonR + 1)/2 + ((diel.EpsilonR - 1)/2) .* (1 ./ sqrt(1 + 12 * diel.Thickness ./ w));

    lambda = physconst('LightSpeed') / (f0 * sqrt(epsilon_eff));
    t = (XL - sqrt(RL * ((Z0 - RL)^2 + XL^2) / Z0 )) / (RL - Z0);
    if t >= 0
        d = 1 / (2 * pi) * atan(t) * lambda;
    else
        d = 1 / (2 * pi) * (pi + atan(t)) * lambda;
    end

    B = (RL^2 * t - (Z0 - XL * t) * (XL + Z0 * t)) / (Z0 * (RL^2 + (XL + Z0 * t)^2));
    l = - 1 / (2 * pi) * atan(B * Z0) * lambda;
    if l < 0
        l = l + lambda/2;
    end
end

function [X, B] = calc_reactances(Z0, ZL)
    RL = real(ZL);
    XL = imag(ZL);

    % Calcula X (reatância da linha principal)
    disc_X = RL*(Z0 - RL);
    X_pos =  sqrt(disc_X) - XL;
    X_neg = -sqrt(disc_X) - XL;

    % Calcula B (susceptância do stub)
    disc_B = (Z0 - RL) / RL;
    B_pos =  sqrt(disc_B) / Z0;
    B_neg = -sqrt(disc_B) / Z0;

    % Retorna todas as soluções possíveis
    X = [X_pos, X_neg];
    B = [B_pos, B_neg];
end
