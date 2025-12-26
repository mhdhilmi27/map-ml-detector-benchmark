% Created by Muhammad Hilmi

clear; clc; close all;

%% Prior Distribution: symbol probability

% for random probability generator uncomment this
%p = rand(1,4);
%p = p/norm(p,1);

p = [1/15 6/15 2/15 6/15];

%% Messages: Generate symbols following the symbol probability
Nmes    = 10000;
symbols = [0 1 2 3];

cs = cumsum(p); 

m = zeros(1, Nmes);
for i = 1:Nmes
    rand_msg = rand;

    if rand_msg < cs(1)
        m(i) = 0;
    elseif rand_msg < cs(2)
        m(i) = 1;
    elseif rand_msg < cs(3)
        m(i) = 2;
    elseif rand_msg < cs(4)
        m(i) = 3;
    end
end

% symbol distribution check
check_prob(m, symbols, p)

%% Sampled Signal Alternatives
Ts = 1;
k = 10;
fc = k/((2*Ts));
Nsamples = 100;

t = linspace(0, Ts-(Ts/Nsamples), Nsamples);

% rectangular pulse g(t) = 1/Ts
g_t = (t >= 0 & t < Ts) * (1/Ts);

% Basis Function
phi1 = g_t .* sqrt(2*Ts) .* cos(2*pi*fc*t);
phi2 = -g_t .* sqrt(2*Ts) .* sin(2*pi*fc*t);

% Signal space
Al = [1/sqrt(2*Ts), -1/sqrt(2*Ts), -1/sqrt(2*Ts), 1/sqrt(2*Ts)];
Bl = [1/sqrt(2*Ts), 1/sqrt(2*Ts), -1/sqrt(2*Ts), -1/sqrt(2*Ts)];

% Construct s_l
s_l = Al(:)*phi1 + Bl(:)*phi2;

% Plot s_l for l=0,1,2,3 (Signal Alternatives)
plot_signal_alt(t, Ts, s_l)

%% Sampled Signal

% mapping message with signal alternatives 
signal_tx = zeros(1, Nmes*Nsamples);
idx = 1;
for i = m
    temp = s_l(i+1,:);
    signal_tx(idx:idx+Nsamples-1) = temp;
    idx = idx + Nsamples;
end

% time axis
t_signal_tx = linspace(0, (Nmes*Ts) - (Ts/Nsamples), Nmes*Nsamples);

plot_signal(t_signal_tx, signal_tx, Nmes, Ts, 's(t)', 'Transmit signal s(t)')

%% signal in orthonormal basis (for visualization purpose)
Ts = 1;

A = [ 1 -1 -1  1] / sqrt(2*Ts);
B = [ 1  1 -1 -1] / sqrt(2*Ts);

figure; hold on; grid on; axis equal;
scatter(A, B, 120, 'filled');

labels = {'s0','s1','s2','s3'};
for k = 1:4
    text(A(k)+0.05, B(k)+0.05, labels{k}, 'FontSize', 12);
end

xlabel('A (along \phi_1)');
ylabel('B (along \phi_2)');
title('4-QAM in Orthonormal Signal Space');

xline(0); yline(0);

xlim([-1 1]);
ylim([-1 1]);

%% Gaussian channel: Adding Noise
N0 = 0.5; %linear value
noise = randn(1, Nmes*Nsamples) * sqrt(Nsamples*N0/2);

signal_rx = signal_tx + noise;
t_signal_rx = t_signal_tx;

plot_signal(t_signal_rx, noise, Nmes, Ts, 'Noise', 'Noise')
plot_signal(t_signal_rx, signal_rx, Nmes, Ts, 's(t)', 'Received signal s(t)')


%% Projection Onto Basis Vectors

r1 = zeros(1, Nmes);
r2 = zeros(1, Nmes);
for k_msg = 1:Nmes
    idx = ((k_msg-1)*Nsamples) + (1:Nsamples);
    r1(k_msg) = (1/Nsamples) * sum(phi1 .* (signal_rx(idx)));
    r2(k_msg) = (1/Nsamples) * sum(phi2 .* (signal_rx(idx)));
end

%% MAP and ML Receiver
% MAP Receiver
mrec_map = (r1.'*Al) + (r2.'*Bl) + (N0/2)*log(p);
[~, highest_idx] = max(mrec_map, [], 2);
symbol_rx = highest_idx-1;
mrec_map = symbol_rx.';

% ML Receiver
mrec_ml = r1.'*Al + r2.'*Bl;
[~, highest_idx] = max(mrec_ml, [], 2);
symbol_rx = highest_idx-1;
mrec_ml = symbol_rx.';

%% Decision Region MAP vs ML
phi1_min = min(Al) - 0.5;
phi1_max = max(Al) + 0.5;
phi2_min = min(Bl) - 0.5;
phi2_max = max(Bl) + 0.5;

phi1_vals = linspace(phi1_min, phi1_max, 200);
phi2_vals = linspace(phi2_min, phi2_max, 200);
[phi1_grid, phi2_grid] = meshgrid(phi1_vals, phi2_vals);

% Treat grid points as "virtual" received projections
r1g = phi1_grid(:);
r2g = phi2_grid(:);      

% MAP decision region
scores_map_grid = r1g*Al + r2g*Bl + (N0/2)*log(p);  
[~, idx_map_grid] = max(scores_map_grid, [], 2);     
regions_map = reshape(idx_map_grid, size(phi1_grid));

% ML decision on the grid
scores_ml_grid = r1g*Al + r2g*Bl;             
[~, idx_ml_grid] = max(scores_ml_grid, [], 2);
regions_ml = reshape(idx_ml_grid, size(phi1_grid));

%% MAP figure: Received vs Sent
figure;

% MAP Sent
subplot(1,2,1);
scatter(r1, r2, 10, m, 'filled'); hold on;
contour(phi1_grid, phi2_grid, regions_map, 'w', 'LineWidth', 1.8);
colormap(jet(4)); clim([0 3]);
axis equal;
xlim([phi1_min phi1_max]); ylim([phi2_min phi2_max]);
xlabel('\phi_1'); ylabel('\phi_2');
title('Messages Received');

% MAP Received
subplot(1,2,2);
scatter(r1, r2, 10, mrec_map, 'filled'); hold on;
contour(phi1_grid, phi2_grid, regions_map, 'w', 'LineWidth', 1.8);
colormap(jet(4)); clim([0 3]);
axis equal;
xlim([phi1_min phi1_max]); ylim([phi2_min phi2_max]);
xlabel('\phi_1'); ylabel('\phi_2');
title('MAP Receiver Decision');

sgtitle('MAP Receiver');

%% ML figure: Received vs Sent
figure;

% ML Sent
subplot(1,2,1);
scatter(r1, r2, 10, m, 'filled'); hold on;
contour(phi1_grid, phi2_grid, regions_ml, 'w', 'LineWidth', 1.8);
colormap(jet(4)); clim([0 3]);
axis equal;
xlim([phi1_min phi1_max]); ylim([phi2_min phi2_max]);
xlabel('\phi_1'); ylabel('\phi_2');
title('Messages Received');

% ML Received
subplot(1,2,2);
scatter(r1, r2, 10, mrec_ml, 'filled'); hold on;
contour(phi1_grid, phi2_grid, regions_ml, 'w', 'LineWidth', 1.8);
colormap(jet(4)); clim([0 3]);
axis equal;
xlim([phi1_min phi1_max]); ylim([phi2_min phi2_max]);
xlabel('\phi_1'); ylabel('\phi_2');
title('ML Receiver Decision');

sgtitle('ML Receiver');

%% Symbols Error Rate
wrong_map = find((mrec_map - m) ~= 0);
sym_error_map = length(wrong_map)/ Nmes;

wrong_ml = find((mrec_ml - m) ~= 0);
sym_error_ml = length(wrong_ml)/ Nmes;

fprintf('Symbol Error Probability MAP : %.4f\n', sym_error_map);
fprintf('Symbol Error Probability ML : %.4f\n', sym_error_ml);

%% Bit Error Rate Investigation
first_mapping = [0 0; 0 1; 1 0; 1 1;];
second_mapping = [1 1; 0 0; 0 1; 1 0;];
Nbps = 2;

m_bit = zeros(1, (length(m)*Nbps));
idx = 1;
for i = 1:length(m)
    sym = m(i);
    temp = first_mapping(sym+1,:);
    m_bit(idx:idx+Nbps-1) = temp;
    idx = idx + Nbps;
end

mrec_1_map_bit = zeros(1, (length(m)*Nbps));
idx = 1;
for i = 1:length(m)
    sym = mrec_map(i);
    temp = first_mapping(sym+1,:);
    mrec_1_map_bit(idx:idx+Nbps-1) = temp;
    idx = idx + Nbps;
end

mrec_2_map_bit = zeros(1, (length(m)*Nbps));
idx = 1;
for i = 1:length(m)
    sym = mrec_map(i);
    temp = second_mapping(sym+1,:);
    mrec_2_map_bit(idx:idx+Nbps-1) = temp;
    idx = idx + Nbps;
end

% Bit Error Rate
wrong_map_1 = find((mrec_1_map_bit - m_bit) ~= 0);
bit_error_map_1 = length(wrong_map_1)/ length(m_bit);

wrong_map_2 = find((mrec_2_map_bit - m_bit) ~= 0);
bit_error_map_2 = length(wrong_map_2)/ length(m_bit);

fprintf('First Mapping Bit Error Probability : %.4f\n', bit_error_map_1);
fprintf('Second Mapping Bit Error Probability : %.4f\n', bit_error_map_2);

%% Helper Function
function check_prob(m, symbols, p)
    Nmes = length(m);
    
    fprintf('Symbol   Count   Generated P   Theory P\n');
    fprintf('---------------------------------------\n');
    
    for i = 1:length(symbols)
        s = symbols(i);
        count_s = sum(m == s);      
        prob_s  = count_s / Nmes;  
        prob_theory = p(i);
    
        fprintf('%d        %d       %.4f      %.4f\n', s, count_s, prob_s, prob_theory);
    end

    fprintf('---------------------------------------\n');
end

function plot_signal_alt(t, Ts, s_l)
    figure;
    titles = ["s_0(t)";"s_1(t)";"s_2(t)";"s_3(t)"];
    
    for i = 1:4
        subplot(4,1,i);
        plot(t, s_l(i,:), 'LineWidth',1.3); grid on;
        xlim([0 Ts]);
        xlabel('t/T_s'); ylabel('s(t)');
        title(titles{i});
    end
    sgtitle('Signal Alternatives');
end

function plot_signal(t, signal, Nmes, Ts, y_label, fig_title)
    figure; 
    plot(t, signal, 'LineWidth', 1.1); 
    grid on;
    xlim([0, Nmes*Ts]);
    xlabel('t/Ts');
    ylabel(y_label);
    title(fig_title);
end


