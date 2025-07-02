% =====================
% PARÂMETROS DE ENTRADA
% =====================
n = 3;                  % Ordem do filtro
fc = 10e9;              % Frequência central (Hz)
bw = 2e8;             % Banda passante (Hz)
Z0 = 50;               % Impedância do sistema (Ohms)
f0 = sqrt((fc + bw/2) * (fc - bw/2));
met = metal('Name', 'Cobre', 'Conductivity', 5.8e7, ...      
       'Thickness', 17e-6);
ripple_db = 0.5;        % Ripple na banda passante (dB)
diel = dielectric('Name','RT5880','EpsilonR', 2.2, ...
       'LossTangent', 0.0009, 'Thickness', 0.787e-3);
fbw = bw/f0;
c = physconst('LightSpeed');
%len = c / (2 * f0 * sqrt(er));

% =====================
% CALCULO LARGURA DA MICROFITA
% =====================

function W = microstrip_width(Z0, er, d)
    % Constantes intermediárias
    A = (Z0 / 60) * sqrt((er + 1)/2) + ((er - 1)/(er + 1)) * (0.23 + 0.11/er);
    B = (377 * pi) / (2 * Z0 * sqrt(er));
    
    % Determina W/d
    if (Z0 <= 120 * pi / sqrt(er + 1))  % W/d <= 2 regime
        Wd = (8 * exp(A)) / (exp(2*A) - 2);
    else  % W/d > 2 regime
        Wd = (2/pi) * ( B - 1 - log(2*B - 1) + ...
            ((er - 1)/(2*er)) * (log(B - 1) + 0.39 - 0.61/er) );
    end
    
    % Calcula largura W
    W = Wd * d;
end

function l = microstrip_length(f0, diel, w)
    eeff = (diel.EpsilonR + 1)/2 + ((diel.EpsilonR - 1)/2) .* (1 ./ sqrt(1 + 12 * diel.Thickness ./ w));
    l = physconst('LightSpeed') / (4 * f0 * sqrt(eeff));
end

Winicial = microstrip_width(Z0, diel.EpsilonR, diel.Thickness);
Linicial = microstrip_length(f0, diel, Winicial);

% =====================
% TABELA g DO POZAR (ripple 0.5 dB)
% =====================
g_values = containers.Map('KeyType','double','ValueType','any');
g_values(3) = [1.5963, 1.0967, 1.5963, 1.0000];
g_values(5) = [1.7058, 1.2296, 2.5408, 1.2296, 1.7058, 1.0000];

g = g_values(n);

% =====================
% CÁLCULO DOS ACOPLAMENTOS J E IMPEDÂNCIAS
% =====================
J = zeros(1, n + 1); % n estágios = n+1 acoplamentos
Zeven = zeros(1, n + 1);
Zodd = zeros(1, n + 1);
W = zeros(1, n + 1);
s = zeros(1, n + 1);


J(1) = sqrt( (pi * fbw) / (2 * g(1)) ) / Z0;
for k = 2:n
    J(k) = ((pi * fbw) / (2 * sqrt(g(k-1) * g(k)))) / Z0;
end
J(n+1) = sqrt( (pi * fbw) / (2 * g(n) * g(n+1)) ) / Z0;

for k = 1:(n+1)
    Zeven(k) = Z0 * (1 + J(k) * Z0 + (J(k) * Z0)^2);
    Zodd(k)  = Z0 * (1 - J(k) * Z0 + (J(k) * Z0)^2);
    
    % Busca dimensões físicas usando função auxiliar
    %[W(k), s(k)] = find_microstrip_dims(Zeven(k), Zodd(k), f0, diel, met);
    %opts = optimset('Display','off','TolX',1e-6);
    %[x_opt, ~] = fminsearch(@(x) optimfun(x, Zeven(k), Zodd(k), diel), [30, 25], opts);
    %W_opt = x_opt(1);
    %s_opt = x_opt(2);

    % NEW FUNC
    line = coupledMicrostripLine("Conductor",met,"Substrate",diel,"Height", diel.Thickness);
    line = design(line, f0, "Z0e", Zeven(k), "Z0o", Zodd(k));
    W(k) = line.Width;
    s(k) = line.Spacing;
end

% =====================
% CÁLCULO DOS COMPRIMENTOS
% =====================

L = zeros(1, n+1);
for k = 1:(n+1)
    L(k) = microstrip_length(f0, diel, W(k));
end

% Linha da entrada
line = microstripLine("Conductor",met,"Substrate",diel,"Height",diel.Thickness);
line = design(line, f0, "Z0",Z0);
Winicial = line.Width;
%Linicial = line.Length;

% =====================
% CONSTRUÇÃO DO MODELO DO FILTRO
% =====================
boardsize = 2 * (W(1) + Winicial);
for i = 1:(n+1)
    boardsize = boardsize + W(i) + s(i);
end

filter = filterCoupledLine("CoupledLineSpacing",s,"CoupledLineWidth",W, ...
    "FilterOrder",n,"Substrate", diel, "PortLineWidth", Winicial, ...
    "Conductor", met, "CoupledLineLength", L, "PortLineLength", Linicial, ...
    "GroundPlaneWidth",boardsize, "Height", diel.Thickness);

disp("Jn")
disp(J .* 1000)
disp("Zeven:")
disp(Zeven)
disp("Zodd:")
disp(Zodd)
disp("W (mm):")
disp(W.*1000)
disp("s (mm):")
disp(s.*1000)
disp("L (mm):")
disp(filter.CoupledLineLength .* 1000)
disp("W e L inicial (extra)")
disp(Winicial * 1000)
disp(Linicial * 1000)

figure;
show(filter)
exportgraphics(gcf, 'filter.png', 'Resolution', 300);

% =====================
% SIMULAÇÃO E PLOT
% =====================
freq = linspace(f0 - 500e7, f0 + 500e7, 10001);
S = sparameters(filter, freq, 'Behavioral',1);

figure;
rfplot(S)
title('Parâmetros S do Filtro');
xlabel('Frequência (GHz)');
ylabel('Magnitude (dB)');
exportgraphics(gcf, 'S_filter.png', 'Resolution', 300);

%rfb = rfbudget(filter,f0,-30,Bwpass);

figure;
rfplot(S, 1, 1);
xlabel('Frequência (GHz)');
ylabel('Magnitude (dB)');
exportgraphics(gcf, 'S11_filter.png', 'Resolution', 300);

figure;
rfplot(S, 2, 1);
xlabel('Frequência (GHz)');
ylabel('Magnitude (dB)');
exportgraphics(gcf, 'S21_filter.png', 'Resolution', 300);

% =====================
% BANDA DE PASSAGEM E BANDA DE REJEIÇÃO
% =====================

s21 = rfparam(S, 2, 1);         % extrai S21
mag = 20*log10(abs(s21));       % converte para dB

% Limiar para banda passante (ex: -3 dB)
threshold_pass = -3;
in_band = find(mag > threshold_pass);

% Frequências de início e fim
f_start = freq(in_band(1));
f_stop  = freq(in_band(end));

% Banda de passagem
BW_pass = f_stop - f_start;
fprintf('Banda de passagem: %.2f MHz (%.4f GHz – %.4f GHz)\n', ...
    BW_pass/1e6, f_start/1e9, f_stop/1e9);

% Rejeição
threshold_stop = -60;
stopband = freq(mag < threshold_stop);  % pontos da rejeição

% Rejeição abaixo de x GHz
stop_left = stopband(stopband < f_start);
fprintf('Rejeição à esquerda até %.4f GHz\n', max(stop_left)/1e9);

stop_right = stopband(stopband > f_stop);
fprintf('Rejeição à direita começa em %.4f GHz\n', min(stop_right)/1e9);


