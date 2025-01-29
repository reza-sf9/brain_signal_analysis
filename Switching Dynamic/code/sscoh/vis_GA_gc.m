function []= vis_GA_gc(GC_val, GC_empirical, t_min, iter)

% GC_val(1, :) - middle part of GC
% GC_val(2, :) - upper bankd of GC
% GC_val(1, :) - lower band of GC
x_ind = 1:length(GC_val);
step_x = t_min(end)./length(x_ind);
x_ind_t = 0: step_x: t_min(end)-step_x;

x_num_step = 3;
fnt_size = 20;

step_num = 99;

gc_mean = GC_val(1, :);

gc_upper = GC_val(2, :);
gc_lower = GC_val(3, :); 

ind_larg = find(gc_lower>1);
gc_lower(ind_larg)=1;
ind_small = find(gc_lower<0);
gc_lower(ind_small)=0;

ind_larg = find(gc_upper>1);
gc_upper(ind_larg)=1;
ind_small = find(gc_upper<0);
gc_upper(ind_small)=0;

y_step = 0: 1/step_num: 1;

K = length(gc_upper);
GC_img = zeros(step_num+1, K);
for k=1:K
    
    gc_mean_temp = gc_mean(k)+.2;
    gc_upper_temp = gc_upper(k)+.2;
    gc_lower_temp = gc_lower(k)+.2;
    
    if gc_upper_temp<gc_lower_temp
        temp = gc_upper_temp;
        gc_upper_temp = gc_lower_temp;
        gc_lower_temp = temp;
    end
    
    temp_mean = abs(y_step-gc_mean_temp);
    ind_min = find(temp_mean==min(temp_mean));
    GC_img(ind_min, k)=1;
    
    % upper values 
    temp_upper =  10*(abs(y_step-gc_upper_temp));
    ind_upper = find(temp_upper==min(temp_upper));
    
    ind_diff_upper = abs(ind_upper-ind_min)+1;
    m_upper = (0-1)./(ind_diff_upper);
    x_dif_upper = 1:ind_diff_upper;
    y_dff_upper = 1+ m_upper.*x_dif_upper;
    
    GC_img(ind_min+1: ind_min+x_dif_upper(end), k) = y_dff_upper;
    
    % lower vlaues
    temp_lower =  10*(abs(y_step-gc_lower_temp));
    ind_lower = find(temp_lower==min(temp_lower));
    
    ind_diff_lower = 2*(abs(ind_min- ind_lower)+1);
    m_lower = (0-1)./(ind_diff_lower);
    x_dif_lower = 1:ind_diff_lower;
    y_dff_lower = 1+ m_lower.*x_dif_lower;
    
    GC_img(ind_min-1: -1:ind_min - x_dif_lower(end), k) = y_dff_lower;
    
    
    
end
figure, 
imagesc(x_ind_t, y_step, GC_img)
xlim([0 x_ind_t(end)])
xlabel('Time(min)')
set(gca,'YDir','normal')
colorbar


% colormap jet
col_ = 4;
if col_ == 1 %'orange'
    c1 = [255, 255, 255]./255;
    c2 = [255, 242, 237]./255;
    c3 = [255, 237, 230]./255;
    c4 = [252, 225, 215]./255;
    c5 = [255, 222, 209]./255;
    c6 = [255, 213, 196]./255;
    c7 = [255, 202, 181]./255;
    c8 = [255, 195, 171]./255;
    c9 = [255, 186, 158]./255;
    c10 = [255, 173, 140]./255;
    c11 = [255, 162, 125]./255;
    c12 = [255, 151, 110]./255;
    c13 = [252, 140, 96]./255;
    c14 = [252, 129, 81]./255;
    c15 = [252, 120, 68]./255;
    c16 = [252, 111, 56]./255;
    c17 = [252, 100, 40]./255;
    c18 = [252, 90, 25]./255;
    c19 = [255, 73, 0]./255;
elseif col_ == 2 % 'pink'
    c19 = [255, 0, 174]./255;
    c18 = [255, 10, 177]./255;
    c17 = [255, 18, 179]./255;
    c16 = [255, 48, 189]./255;
    c15 = [255, 61, 193]./255;
    c14 = [255, 77, 198]./255;
    c13 = [255, 84, 200]./255;
    c12 = [255, 97, 204]./255;
    c11 = [255, 105, 206]./255;
    c10 = [255, 117, 210]./255;
    c9 = [255, 130, 214]./255;
    c8 = [255, 140, 217]./255;
    c7 = [255, 145, 219]./255;
    c6 = [255, 153, 222]./255;
    c5 = [255, 161, 224]./255;
    c4 = [255, 168, 226]./255;
    c3 = [247, 203, 232]./255;
    c2 = [250, 227, 242]./255;
    c1 = [255, 255, 255]./255;
    
elseif col_== 3 %'cyan'
    c19 = [0, 242, 255]./255;
    c18 = [13, 243, 255]./255;
    c17 = [28, 244, 255]./255;
    c16 = [36, 244, 255]./255;
    c15 = [48, 245, 255]./255;
    c14 = [56, 245, 255]./255;
    c13 = [69, 246, 255]./255;
    c12 = [89, 247, 255]./255;
    c11 = [105, 248, 255]./255;
    c10 = [120, 249, 255]./255;
    c9 = [138, 250, 255]./255;
    c8 = [150, 251, 255]./255;
    c7 = [166, 252, 255]./255;
    c6 = [184, 253, 255]./255;
    c5 = [197, 252, 252]./255;
    c4 = [210, 252, 252]./255;
    c3 = [227, 255, 255]./255;
    c2 = [240, 252, 252]./255;
    c1 = [255, 255, 255]./255;

elseif col_== 4 % red'
    c19 = [245, 7, 31]./255;
    c18 = [255, 15, 39]./255;
    c17 = [250, 27, 50]./255;
    c16 = [252, 38, 60]./255;
    c15 = [252, 50, 71]./255;
    c14 = [252, 61, 81]./255;
    c13 = [255, 74, 93]./255;
    c12 = [250, 82, 100]./255;
    c11 = [255, 97, 114]./255;
    c10 = [252, 109, 124]./255;
    c9 = [252, 121, 135]./255;
    c8 = [255, 133, 146]./255;
    c7 = [255, 145, 157]./255;
    c6 = [255, 158, 169]./255;
    c5 = [252, 172, 180]./255;
    c4 = [252, 184, 191]./255;
    c3 = [252, 220, 223]./255;
    c2 = [252, 237, 239]./255;
    c1 = [255, 255, 255]./255;
end

colormap([c1; c2; c3; c4; c5; c6; c7; c8; c9; c10; c11; c12; c13; c14; c15; c16; c17; c18; c19;]);
ax = gca;
x_tick_values = fun_finding_tick(ax, x_num_step, 'x'); % nested function 
ax.XTick = x_tick_values;

ax.YTick = [.2 .5 .8];



hold on

% figure,
plt_gc_empr = 1; 
if plt_gc_empr==1
    x_ind = 1:length(GC_empirical);
    step_x = t_min(end)./length(x_ind);
    x_ind_t_empirical = 0: step_x: t_min(end)-step_x;
    
    filt_ = 1;
    if filt_ == 1
        sigma = .8; % pick sigma value for the gaussian
        gaussFilter = gausswin(6*sigma + 1)';
        gaussFilter = gaussFilter / sum(gaussFilter); % normalize
        % use the convolution:
        plt_GC_emp = conv(GC_empirical, gaussFilter, 'same');
        l_width = 1.5;
    else
        plt_GC_emp = GC_empirical;
        l_width = 0.1;
    end
    
    plot(x_ind_t_empirical, plt_GC_emp, 'LineWidth', l_width, 'color', 'k')
   
    %
%     ax = gca;
%     ax.FontSize = fnt_size;
%     ax.XTick = x_tick_values;
%     xlabel(ax,'Time (min)')
% 
%     ax.YTick = [.2 .5 .8];

end



str_save = sprintf('/Results/gc_val_Iter%d.png', iter);
saveas(gcf,[pwd str_save]);
end

% finding the proper val for ticks
function tick_val = fun_finding_tick(ax, num_step, axis)

switch axis
    case 'x'
        lim_val = ax.XLim;
    case 'y'
        lim_val = ax.YLim;
end

step_tick = ((lim_val(end) - lim_val(1))./num_step);
cnt = find_order(step_tick); % call find_order
if cnt> 1
    step_tick_n = round(step_tick, -(cnt-1));
    tick_val = round(lim_val(1) + (step_tick_n./2), -(cnt-1)) : step_tick_n: lim_val(end);
else
    step_tick_n = round(step_tick, cnt);
    tick_val = round(lim_val(1) + (step_tick_n./2), cnt) : step_tick_n: lim_val(end);
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
