clc
clearvars -except data_det T Fs
close all

num_switch = 1;



%% data
%%% loading data
if exist('data_det') ~=1
    load('E:\oneDrive_WPI\OneDrive - Worcester Polytechnic Institute (wpi.edu)\Reza\Data\emery_data\eeganes07laplac250_detrend_all.mat');
%     load('D:\OneDrive - Worcester Polytechnic Institute (wpi.edu)\Reza\Data\emery_data\eeganes07laplac250_detrend_all.mat');

end

% load('E:\oneDrive_WPI\OneDrive - Worcester Polytechnic Institute (wpi.edu)\-0\switching_data\ali_2.mat');
% 
% data_det = strct_switch_data.data_det; 
% t_sec = strct_switch_data.t_sec;
% Fs = strct_switch_data.Fs;
% Ts = strct_switch_data.Ts;
% 
% % data_reza = data_reza.data;
% data_ = data_det.';
% l_d = size(data_det, 2);
% 
% t_s = t_sec;
% t_min = t_s/60;
% ch_n = size(data_det, 1);

% data_det =data_;
method_GC = 'PNAS';

% this method affetcs only on angle_eig_vec_fun
% method_angle = 1;
method_angle = 2;

t_min = T./60;

ch_n = 5;

%% config 
%%% analyze configuration
config.ch_num = ch_n;                                       % number of channels
config.Fs = Fs;                                         % sampling rate of data
config.sample_r = Fs;                                   % size of slepian window
config.win_sec = ch_n*2;
config.win_length = config.win_sec*config.sample_r;      % (sample) length of window that we want to calc Global Coherence on them
% config.seg_num = config.win_sec;                                      % (number) number of segments that we want to devide a window into them  (number)
config.seg_num = 1;
config.seg_length = config.win_length./config.seg_num;   % (sample)
config.f_l = 12;                                          % lower bound of desired freq interval
config.f_u = config.f_l + 0;                            % upper bound of desired freq interval
config.over_lap = .0;                                   % (of 1) length of overlap window based win length
config.method_angle = 2;
config.num_slepian = 5;
config.halfbandwidth = 2*config.num_slepian-1;
config.t_min = t_min;

%%% data extraction
m=size(data_det);
data_chs = zeros(m(1) , config.ch_num);

%%% Time Plot
plt_time = 0;
while plt_time == 1
    plt_time = 0;
    ch_plt_vec = [6 40];
    fnt_size = 25;
    
    for i=1 : config.ch_num
%         ch_plt = ch_plt_vec(i);
        figure('units','normalized','outerposition',[0 0 1 1]),
        plot(t_min, data_det(:,i), 'Color', [27, 27, 28]./255)
        
        ylim([-3000 3000])
        
        set(gca,'FontSize', fnt_size)
        
% 
        set(gca, 'YTIck', [-500:500:500])
        set(gca, 'YAxisLocation', 'right')
        
        


        xlabel('Time (mins)')
        ylabel('Amplitude')
        xlim([0 t_min(end)])
        ylim([-800 800])
        
        set(gca, 'XTIck', [20:40:140])
        
        
        set(gcf, 'PaperPosition', [0 0 6 4]);
%         print('time_40','-dpng','-r600')
        
        
        
    end
    close all
end


%% extracting data
switch ch_n
    case 32
%         ch_rand = [1 2 3 4 5 6 7 11 12 13 14  19 24  29 35 38 40 42 44 48 49 50 51 52 53 54 58 59 60 62 63 64];
        ch_rand = [1 2 3 4 5 6 7 11 12 13 14 17 18 19 24 25 29 33 34 35 36 37 38 40 41 42 43 48 49 54 55 60];
%         ch_rand = [18 17 48 54 24 12 11 6 5 2 36 49 60 55 3 7 8 34 39 44 50 61 56 45 40 10 16 30 27 46 52 28];
%         ch_rand = sort(ch_rand);
    case 20
        ch_rand = [1 2  4 5  11 12  14  18  25  33  35 36 37  41 42 43 48  54 55 60];

    
    case 10 
%         ch_rand = [17 12 11 6 36 3 10 27 46 28];
%         ch_rand = [17 48 12 11 55 7 27 46 28 64];
%         ch_rand = [18, 41, 19, 6, 2, 36,43, 34, 40, 16];
%         ch_rand = [18, 41, 6, 2, 36 , 34, 40, 16, 22, 46];
        ch_rand = [18, 17, 41, 6, 2, 36, 33, 53, 59, 28];
        
    case 8 
%         ch_rand = [17 48 12 11 27 46 28 64];
        ch_rand = [17, 12, 41, 6, 2, 1, 22, 53];
        
    case 5
%         ch_rand = [6, 36, 23, 53, 58];
        ch_rand = [6, 2, 36, 34, 10];
        
    case 3
        ch_rand = [6, 36, 2];
        
    case 64
        ch_rand = 1: config.ch_num;
end

% data_chs = data_det; 
for i=1: length(ch_rand)
    data_chs(:, i) = data_det(:,ch_rand(i));
end
% plt_time_data(data_chs, t_min)
config.num_ch = ch_rand;
% config.num_ch = 19;
% normalized time data
% max_val = max(abs(data_chs(:)));
% data_chs = data_chs./max_val;

 %% multitaper & ordinary spectrogram
plt_mltp = 0;
while plt_mltp==1
    plt_mtlp =0;
    % input structure of multitaper spectrogram
    cfg_mtp_spect = struct;
    
    cfg_mtp_spect.T = T;               % time indices
    cfg_mtp_spect.ch_des_num = length(ch_rand);     % desired channel for plotting spectrogram
    cfg_mtp_spect.f_des_l = 1;         % (freq(Hz)) the lower range for plotting spectrogram
    cfg_mtp_spect.f_des_u = 29;         % (freq(Hz)) the lower range for plotting spectrogram
    cfg_mtp_spect.mtp_NW = 3.5;         % NW (halfbandwidth) value for calculating multitaper
    cfg_mtp_spect.mtp_win_length = 60;  % (sec) win length for calculating multitaper
    cfg_mtp_spect.mtp_over_lap = .3;   % (of 1) length of overlap window based win length
    
    for i =1: length(ch_rand)
        cfg_mtp_spect.ch_des_num = ch_rand(i);       % desired channel for plotting spectrogram
        
        % this fuction plot spectrogram based multitaper approach for desired channel
        multiTaper_spectrogram(data_det , cfg_mtp_spect)
    end
    % this fuction plot spectrogram based multitaper approach for desired channel
    
end



%% calc GC and Igen info
% [GC , sorted_eig_info] = calc_GC_steps(data_chs , config);
[GC_overal , sorted_eig_info, fft_seg_mean_removed] = calc_GC_steps_rsf(data_chs , config);

% str_save = sprintf('anesthesia_ch%d_win%d_overLap%.1f_fL%d_fU%d_numSlep%d.mat', config.ch_num, config.win_sec, config.over_lap, config.f_l, config.f_u, config.num_slepian);
% save(str_save, 'sorted_eig_info')

% we can load instead of calculating again GC adn EIGEN info if we saved it
% for desired freq before

% struct_sub6_gc = [];

% load('anesthesia_ch32_win8_overLap0_fL11_fU11_numSlep5')
% load('ch32_win32_seg32_frL25_frU25_overLap0')

%% (ANALYZING RESULTS)

win_sec_GCoh_empircal = config.win_sec;
win_sec_observation = config.win_sec / config.seg_num;
fprintf('\n Number of channel: %d' , ch_n);
fprintf(' \n Length of Window for calc GCOh empirical: %d sec', win_sec_GCoh_empircal)
fprintf(' \n Length of subWindow for calc Cov mat of GCOh empirical (Window for SS-Coh): %d sec', win_sec_observation)
fprintf('\n------------------------------\n')

fnt_size = 20;

%%% Global Coherence PLOTTING
t_end = floor(t_min(end));

if num_switch ==3
    swtich_vec = [1, 2, 3, 2, 1];
    t_min_switch = [0, 45, 55, 95, 105, t_min(end)];
    
elseif num_switch ==2
    swtich_vec = [1, 2, 1];
    t_min_switch = [0, 45, 105, t_min(end)];

elseif num_switch ==1
    swtich_vec = [1];
    t_min_switch = [0, t_min(end)];
end
plot_GC_fun_switch(GC_overal , config, fnt_size, t_min_switch);

% close all;

% save('E:\oneDrive_WPI\OneDrive - Worcester Polytechnic Institute (wpi.edu)\-0\figure_1\gc_empr.mat', 'GC_overal', 'config');

%%% calculating angle between largest corresponding eigenvectors over time (1d)
% angle_eig_vec_fun(sorted_eig_info , config, fnt_size);
% close all;

%% reza for whole cross spectral/Yk-observation- (PAPER EMBC)

initial_paper_switch(fft_seg_mean_removed, config, GC_overal, swtich_vec, t_min_switch, num_switch)
k=1
%% sum eigenvalues
% [max_eig_val , sum_eig_val] = max_sum_eig_value(sorted_eig_info , config );
% close all

% ang_diff_phase_ev(sorted_eig_info , config, fnt_size)
% close all 

% trace_ev(sorted_eig_info , config, fnt_size )
% close all


%% nested functions
function plt_time_data(data, t_min)

my = size(data);
x = t_min;

max_val = max(data(:));
y_lim = [-1 max_val+1];

tick_num_step = 3;
% y_ticks = [0 4];
%% setup values
left_pos = .15; 
bot_pos = .15;
width_val = .7;
height_val = .20;


bot_diff = .22;

line_width = .2;
fnt_size= 14;

%% plots 
h =figure;
ax = subplot(3,1,3);
plot(x, data(:, 1), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 0*bot_diff) width_val height_val];
ylabel('Ch 1')
xlabel('Time (min)')
ylim([-600 600])
xlim([0 t_min(end)])

ax.FontSize = fnt_size;
% xtick
x_num_step = tick_num_step;
x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested fun
ax.XTick = x_tick_values;

% ytick
% y_num_step = tick_num_step-1;
% y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
% ax.YTick = y_tick_values;
ax.YTick = [-300 300];

ax.XAxis.FontSize = fnt_size + 3;
ax.XLabel.FontSize = fnt_size + 3;


ax = subplot(3,1,2);
plot(x, data(:, 4), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 1*bot_diff) width_val height_val];
ylabel('Ch 4')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ytick
% y_num_step = tick_num_step-1;
% y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
% ax.YTick = y_tick_values;

ylim([-600 600])
xlim([0 t_min(end)])
ax.FontSize = fnt_size;
ax.YTick = [-300 300];

ax = subplot(3,1,1);
plot(x, data(:, 3), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 2*bot_diff) width_val height_val];
ylabel('Ch 3')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ytick
y_num_step = tick_num_step-1;
y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
ax.YTick = y_tick_values;

% ylim(y_lim)
xlim([0 t_min(end)])
ax.FontSize = fnt_size;

ax = subplot(5,1,2);
plot(x, data(:, 4), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 3*bot_diff) width_val height_val];
ylabel('Ch 4')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ytick
y_num_step = tick_num_step-1;
y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
ax.YTick = y_tick_values;

% ylim(y_lim)
xlim([0 t_min(end)])
ax.FontSize = fnt_size;

ax = subplot(5,1,1);
plot(x, data(:, 5), 'LineWidth', line_width)
ax.Position = [left_pos (bot_pos+ 4*bot_diff) width_val height_val];
ylabel('Ch 5')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
% ytick
y_num_step = tick_num_step-1;
y_tick_values = fun_finding_tick(ax, y_num_step, 'y'); % nested fun
ax.YTick = y_tick_values;

% ylim(y_lim)
xlim([0 t_min(end)])
ax.FontSize = fnt_size;



% [left bottom width height]

end

% finding the proper val for ticks
function tick_val = fun_finding_tick(ax, num_step, axis)

switch axis
    case 'x'
        lim_val = ax.XLim;
        div_num = 2;
    case 'y'
        lim_val = ax.YLim;
        div_num = 3;
end

step_tick = ((lim_val(end) - lim_val(1))./num_step);
cnt = find_order(step_tick); % nested fun
if cnt> 1
    step_tick_n = round(step_tick, -(cnt-1));
    tick_val = round(lim_val(1) + (step_tick_n./div_num), -(cnt-1)) : step_tick_n: lim_val(end);
else
    step_tick_n = round(step_tick, cnt);
    tick_val = round(lim_val(1) + (step_tick_n./div_num), cnt) : step_tick_n: lim_val(end);
end


end

% finding the order of input value
function cnt = find_order(x)

if x> 1
    stp_con = 0;
    cnt = -1;
    while stp_con == 0
        cnt = cnt +1;
        a = floor(x./10^cnt);
        
        if a == 0
            stp_con=1;
        end
        
    end
else
    stp_con = 0;
    cnt = 0;
    while stp_con == 0
        cnt = cnt +1;
        a = floor(x./10^(-cnt));
        
        if a > 0
            stp_con=1;
        end
        
    end
end

end
