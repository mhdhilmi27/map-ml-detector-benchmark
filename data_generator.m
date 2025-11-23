% Created by Muhammad Hilmi

clear; clc; close all;

list_N0 = linspace(0.05, 0.95, 19);

for N0 = list_N0
    %% Prior Distribution: symbol probability
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
    
    %% Gaussian channel: Adding Noise                   
    noise = randn(1, Nmes*Nsamples) * sqrt(Nsamples*N0/2);
    
    signal_rx = signal_tx + noise;
    t_signal_rx = t_signal_tx;
    
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
    
    %% write data
    data_N0 = N0 * ones(1, length(m));
    data = [m; r1; r2; mrec_map; mrec_ml; data_N0];
    data = data.';
    
    data(:,2) = round(data(:,2), 3);  
    data(:,3) = round(data(:,3), 3);   
    
    path = 'training_data';
    filename = sprintf('%s/dataset_n0_%.2f.csv', path, N0);
    T = array2table(data, 'VariableNames', ...
        {'m','r1','r2','mrec_map','mrec_ml','N0'});
    writetable(T, filename);
end