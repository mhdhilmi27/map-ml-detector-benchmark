%{
plot_signal(t_signal_tx, signal_tx, Nmes, Ts, 's(t)', 'Transmit signal s(t)')
plot_signal(t_signal_rx, noise, Nmes, Ts, 'Noise', 'Noise')
plot_signal(t_signal_rx, signal_rx, Nmes, Ts, 's(t)', 'Received signal s(t)')
check_prob(m, symbols, p);

figure; 
plot(t_signal_rx, noise, 'LineWidth', 1.1); 
grid on;
xlim([0, Nmes*Ts]);
xlabel('t/Ts'); ylabel('Noise'); title('Noise');

figure; 
plot(t_signal_rx, signal_rx, 'LineWidth', 1.1); 
grid on;
xlim([0, Nmes*Ts]);
xlabel('t/Ts'); ylabel('s(t)'); title('Receive signal s(t)');


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

%{ 
function for checking the prior distibution of generated m message
(can be ignored)
%}

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
end


%is_first_matrix = 1;
%for i = m
%   if is_first_matrix == 1
%       signal_tx = s_l(i+1,:);
%       is_first_matrix = 0;
%   elseif is_first_matrix == 0
%       temp = s_l(i+1,:);
%       signal_tx = [signal_tx, temp];
%   end
%end

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

figure; 
plot(t_signal_tx, signal_tx, 'LineWidth', 1.1); 
grid on;
xlim([0, Nmes*Ts]);
xlabel('t/Ts'); ylabel('s(t)'); title('Transmit signal s(t)');
%}

list_N0 = linspace(0.05, 0.95, 19);
for i = list_N0
    disp(i);
end