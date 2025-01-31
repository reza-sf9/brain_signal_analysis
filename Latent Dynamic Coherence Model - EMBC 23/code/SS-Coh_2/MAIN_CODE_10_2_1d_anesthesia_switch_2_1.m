clc
clear
close all

dim = 1; % 1 or 2
ch_num = 5; % 5, 8, 10, 20
Iter = 11;
num_switch = 1;


use_bias = 0; % if 1, use all bias elements  for W
if use_bias==0
    use_bias_2_end = 1; % if 1 use all bias elements except first
end

using_simulation = 0; % 1: using sim, 0: using real data
if using_simulation==0

    coef_noise_W = 0.5;
    coef_noise_mu = 0;

    coef_w = 1;

    coef_length_sim = 1;

    new_val = 1;
    new_W = new_val;
    new_SIM_data = new_val;

    fprintf('Using Real Data Set for Analysing SS-GCoh Method')
    fprintf('\n------------------------------\n')

elseif using_simulation==1
    coef_noise_W = 1;

    coef_noise_mu = 3;

    coef_w = .1;

    coef_length_sim = .1;

    new_val = 1;
    new_W = new_val;
    new_SIM_data = new_val;

    fprintf('Using Simulation Data Set for Analysing SS-GCoh Method')
    fprintf('\n------------------------------\n')
end


col_W_1 = 2: ch_num;
col_W_2 = 2:ch_num;
if use_bias ==1
    col_W_3 = 1: ch_num;

elseif use_bias_2_end==1
    col_W_3 = 2: ch_num;
end

QQ_sim = [15 0; 0 80].*10^-3;

%% PARAMETER UPDATE
ParamaUpdate.sigma_update = 0; % 1 update, 0 don't update
ParamaUpdate.mu_update = 0;
ParamaUpdate.L_update  = 0;
ParamaUpdate.W_update  = 1;

% based on switch we need to run the model K times and then connect them to
% each other

if num_switch==3
    structure_switch = [1, 2, 3, 2, 1];
elseif num_switch==2
    structure_switch = [1, 2, 1];
elseif num_switch ==1
    structure_switch = [1];
end

Num_seg = length(structure_switch);


%%% 
PARAM_q1= cell(1, Iter+1);
Bayes_Step_q1 = cell(1, Iter);

if num_switch==2
    PARAM_q2= cell(1, Iter+1);
    Bayes_Step_q2 = cell(1, Iter);
elseif num_switch==3

    PARAM_q2= cell(1, Iter+1);
    Bayes_Step_q2 = cell(1, Iter);

    PARAM_q3= cell(1, Iter+1);
    Bayes_Step_q3 = cell(1, Iter);
end

for q_k=1: num_switch


    % load inilial points
    [cfg_Init, Y_k]  = gc_load_yk_cfg(ch_num, num_switch);

    l_tot = length(Y_k);
    if q_k ==1
        gc_print_fun_1(cfg_Init, ch_num) % print some info
    end



    %% calc switch information
    empirical_switch_point = cfg_Init.empirical_switch_pionts;
    vec_1 = [empirical_switch_point, 0];
    vec_2 = [0, empirical_switch_point];

    vec_diff = abs(vec_1-vec_2);
    %     figure()
    %     stem(vec_diff)

    ind_jump = find(vec_diff>0);
    ind_jump(end) = ind_jump(end);
    mat_switch_ind = zeros(3, length(structure_switch));
    mat_switch_ind(1, :) = structure_switch;

    for ii=1: length(structure_switch)
        mat_switch_ind(2, ii) = ind_jump(ii);
        mat_switch_ind(3, ii) = ind_jump(ii+1)-1;
    end

    mat_switch_active_ind = mat_switch_ind;
    vec_1 = mat_switch_ind(1, :);
    ind_active = find(vec_1==q_k);
    vec_active = zeros(1, length(vec_1));
    for i_=1: length(ind_active)
        vec_active(1, ind_active(i_)) = 1;
    end
    mat_switch_active_ind(1, :) = vec_active;

    %% end of index


    %% segment each state and run one state model for it


    % generate new W
    if new_W==1
        W = gc_gen_W_init(coef_w, ch_num);
    end

    L_init = cfg_Init.L;
    L_init_switch = L_init{q_k};



    % print state model
    str_print = sprintf('State Space Model = %d/%d ', q_k, num_switch);
    disp(str_print)


    t_min = cfg_Init.cfgGen.t_min;




    y_end = length(Y_k);



    K = length(Y_k);         % number of all stepes
    KK = floor(coef_length_sim*K);
    KK = 1000;

    if new_SIM_data == 1
        [xk_sample_tot, Yk_sample_tot, cnfg_sim_tot] = gc_sim_data(dim, QQ_sim, L_init, KK, Y_k, ch_num, q_k);
        [xk_sample, Yk_sample, cnfg_sim] = gc_sim_data(dim, QQ_sim, L_init, KK, Y_k, ch_num, q_k);
        %         figure, plot(xk_sample)
    else
        str_load = sprintf('X_Y_sim%d_good.mat', dim);
        load(str_load)
        %     figure, plot(xk_sample)
    end


    %% NEW INITIAL PARAMS
    [init_Param, update_coef_W, W_sim] = gc_init_param(ParamaUpdate, cnfg_sim, coef_noise_mu, col_W_1, col_W_2, col_W_3, coef_noise_W, ch_num, dim);

    W_init = init_Param.W;

    %% init param
    init_Param.incongruent_vec = cnfg_sim.incongruent_vec;
    init_Param.mu_1_0 = [0; 0];
    init_Param.sigma_1_0 = 300*[1 0; 0 1]*10^-3;
    init_Param.F = cnfg_sim.F;
    init_Param.D = cnfg_sim.D;
    init_Param.dim = cnfg_sim.dim;
    init_Param.L = L_init_switch;

    %% initialize Param and Bayes container
    if num_switch==1
        if q_k==1
            PARAM_q1{1}= init_Param;
        end
    elseif num_switch ==2

        if q_k==1
            PARAM_q1{1}= init_Param;
        elseif q_k==2
            PARAM_q2{1}= init_Param;
        end

    elseif num_switch==3
        if q_k==1
            PARAM_q1{1}= init_Param;
        elseif q_k==2
            PARAM_q2{1}= init_Param;
        elseif q_k==3
            PARAM_q3{1}= init_Param;
        end
    end


    %% training model

    sv_update = [];

    step_interval = 1:length(Y_k);
    break_con = 0;
    thr_q = 10^-5;


    % ParamaUpdate.W_update  = 0;

    if using_simulation == 1
        Yk = Yk_sample;
    end

    % model training section
    for iter=1:Iter

        disp(iter)
        disp('-----------')

        if q_k==1
            param_curr = PARAM_q1{iter};
        elseif q_k==2
            param_curr = PARAM_q2{iter};
        elseif q_k==3
            param_curr = PARAM_q3{iter};
        end


        %%% filter - smoother
        Bayes = gc_filter_smoother_switch(Y_k, param_curr, mat_switch_active_ind);

        %%% update parameters
        [updated_param, error_tot] = gc_parameter_update_switch(param_curr, ParamaUpdate, Bayes, Y_k, iter, W_sim, update_coef_W, mat_switch_active_ind);


        if q_k==1
            Bayes_Step_q1{iter} = Bayes;
            PARAM_q1{iter+1} = updated_param;
        elseif q_k==2
            Bayes_Step_q2{iter} = Bayes;
            PARAM_q2{iter+1} = updated_param;
        elseif q_k==3
            Bayes_Step_q3{iter} = Bayes;
            PARAM_q3{iter+1} = updated_param;
        end

    end
end

iter_first = 1;
iter_last = Iter-1;

close all

% Plot State Process and GCoh 
for q_k=1: num_switch

    if q_k==1
        Bayes = Bayes_Step_q1;
        PARAM = PARAM_q1;
    elseif q_k==2
        Bayes = Bayes_Step_q2;
        PARAM = PARAM_q2;
     elseif q_k==3
        Bayes = Bayes_Step_q3;
        PARAM = PARAM_q3;
    end
    

    bayes_first = Bayes{iter_first};
    bayes_last = Bayes{iter_last};
    

    % finding the y range
    strct_min_max = fun_plt_find_yrange(bayes_first, bayes_last, xk_sample_tot, dim, using_simulation);

    %% plot sim
    % fun_plt_simulated_data(xk_sample, dim, strct_min_max)
    % fun_plot_abs_yk(Yk_sample, using_simulation, dim)


    %% smoother plot
%     vis_GA_x(bayes_first, xk_sample, dim, iter_first, strct_min_max, using_simulation, t_min, q_k);
%     vis_GA_x(bayes_last, xk_sample, dim, iter_last, strct_min_max, using_simulation, t_min, q_k);


    %% GC
    vis_GA_gc(bayes_first.GC_val, cfg_Init.GC_empirical, t_min, iter_first, q_k)
    vis_GA_gc(bayes_last.GC_val, cfg_Init.GC_empirical, t_min, iter_last, q_k)


end



BAYES_tot = [];
PARAM_tot = [];

if num_switch ==1 
    BAYES_tot.bayes_1 = Bayes_Step_q1;

    PARAM_tot.param_1 = PARAM_q1;
elseif num_switch ==2 
    BAYES_tot.bayes_1 = Bayes_Step_q1;
    BAYES_tot.bayes_2 = Bayes_Step_q2;
    
    PARAM_tot.param_1 = PARAM_q1;
    PARAM_tot.param_2 = PARAM_q2;

elseif num_switch ==3
    BAYES_tot.bayes_1 = Bayes_Step_q1;
    BAYES_tot.bayes_2 = Bayes_Step_q2;
    BAYES_tot.bayes_3 = Bayes_Step_q3;
    
    PARAM_tot.param_1 = PARAM_q1;
    PARAM_tot.param_2 = PARAM_q2;
    PARAM_tot.param_3 = PARAM_q3;
end


str_save = sprintf('LDCM_numSwitch_%d_numIter_%d_numCh_%d.mat', num_switch, Iter, ch_num);


LDCM_strct = [];
LDCM_strct.Bayes = BAYES_tot;
LDCM_strct.PARAM = PARAM_tot;
LDCM_strct.cfg_Init = cfg_Init;
save(str_save, 'LDCM_strct');






close all


%% plot error tot iteration
% vis_GA_error_TOT(abs_error_TOT_iter, dim, Iter)
% close all

%% W plot
% vis_GA_W(W_sim, W_iter, dim, using_simulation, use_bias)


%% mu plot
% vis_GA_mu(cnfg_sim.mu, MU_iter, dim, using_simulation)

%% cholseky plot

if using_simulation==1
    %% Calc Cholesky decomposition - plot
    vis_GA_cholesky(xk_sample, Yk, cnfg_sim, W_iter, MU_iter, dim, bayes_first, bayes_last, Iter)
    close all

elseif using_simulation==0

    w_type = 'init';
    x_tick_init = [-2, 0, 2];
    x_lim_init = [-4 4];
    y_tick_init = [-2, 0, 2];
    y_lim_init = [-4 4];
    W_init = squeeze(W_iter(1, :, :));
    x_inint = bayes_first.MU_smoother;
    mu_init = MU_iter(1, :).';
    L = cnfg_sim.L;
    vis_GA_cholesky_real(x_inint, Yk, L, W_init, mu_init, w_type, dim, x_tick_init, x_lim_init, y_tick_init, y_lim_init)
    close all


    w_type = 'final';
    x_tick_fin = [-2, 0, 2];
    x_lim_fin = [-4 4];
    y_tick_fin = [-2, 0, 2];
    y_lim_fin = [-4 4];
    x_final = bayes_last.MU_smoother;
    W_final = squeeze(W_iter(Iter, :, :));
    mu_final = MU_iter(Iter, :).';
    L = cnfg_sim.L;
    vis_GA_cholesky_real(x_final, Yk, L, W_final, mu_final, w_type, dim,  x_tick_fin, x_lim_fin, y_tick_fin, y_lim_fin)
    close all


end



%% nested functions

%% calculate expectation of D^-1
function [E_Sigma_W_x_inv] = gc_expectation_D_inv(MU_smoother, W, incongruent_vec, dim)

K = length(MU_smoother);
sigma_dim = length(W);
Sigma_inv = zeros(K, sigma_dim, sigma_dim);

for k=1 : K
    x_k = MU_smoother(:, k);

    Sigma_W_x_inv = zeros(sigma_dim);

    for j=1: sigma_dim
        w_m = W(j, 1);
        w_bias = W(1,2);
        lambda_m = exp(x_k(1,1) * w_m(1,1)+ w_bias);
        Sigma_W_x_inv(j, j) = inv(lambda_m);
    end

    Sigma_inv(k, :, :) = Sigma_W_x_inv;

end
E_Sigma_W_x_inv = Sigma_inv;
end



% load data
function [cfg_Init, Y_k]= gc_load_yk_cfg(ch_num, num_switch)

str_dict = '.\load_folder\';

switch ch_num
    case 32
        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr11_sr256_winSec8_segNum8_TRIAL', ch_num);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr11_sr256_winSec8_segNum8_TRIAL', ch_num);

        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec32_segNum32_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec32_segNum32_TRIAL', ch_num);

    case 20
%         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr12_sr256_winSec40_segNum40_TRIAL', ch_num);
%         str_name_2 = sprintf('Anethesia_cfg_Init_chNum20_fr12_sr256_winSec40_segNum40_TRIAL', ch_num);

        str_name_1 = sprintf('Anethesia_switch_%d_Yk_chNum%d_fr12_sr250_winSec%d_segNum%d_TRIAL', num_switch, ch_num, ch_num, ch_num);
        str_name_2 = sprintf('Anethesia_switch_%d_cfg_Init_chNum%d_fr12_sr250_winSec%d_segNum%d_TRIAL', num_switch, ch_num, ch_num, ch_num);

    case 10
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL_good', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL_good', ch_num);

        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL', ch_num);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec20_segNum20_TRIAL', ch_num);

    case 8
        str_name_1 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec16_segNum16_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec16_segNum16_TRIAL', ch_num);

    case 5
        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec8_segNum8_TRIAL', ch_num);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec8_segNum8_TRIAL', ch_num);

        % 'Anethesia_switch_Yk_chNum5_fr12_sr250_winSec10_segNum10_TRIAL.mat';
        % 'Anethesia_switch_cfg_Init_chNum5_fr12_sr250_winSec10_segNum10_TRIAL.mat';

        str_name_1 = sprintf('Anethesia_switch_%d_Yk_chNum%d_fr12_sr250_winSec10_segNum10_TRIAL', num_switch, ch_num);
        str_name_2 = sprintf('Anethesia_switch_%d_cfg_Init_chNum%d_fr12_sr250_winSec10_segNum10_TRIAL', num_switch, ch_num);

        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr12_sr256_winSec10_segNum10_TRIAL', ch_num);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr12_sr256_winSec10_segNum10_TRIAL', ch_num);

    case 3
        str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec15_segNum15_TRIAL', ch_num);
        str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec15_segNum15_TRIAL', ch_num);

        %         str_name_1 = sprintf('Anethesia_Yk_chNum%d_fr13_sr256_winSec10_segNum10_TRIAL', ch_sim);
        %         str_name_2 = sprintf('Anethesia_cfg_Init_chNum%d_fr13_sr256_winSec10_segNum10_TRIAL', ch_sim);
end

str_load_1 = sprintf('%s%s', str_dict, str_name_1);
str_load_2 = sprintf('%s%s', str_dict, str_name_2);


load(str_load_1);
load(str_load_2);
% yk = load(str_load_1);
% Y_k = cell2mat(struct2cell(yk));
% b = load(str_load_2);
% b2 = struct2cell(b);
end

% generate initial values for W
function [W] = gc_gen_W_init(coef_w, ch_num)
w_1 = coef_w*(0 + 1*rand(ch_num, 1));
w_1(1,1) = 1;
w_2 = coef_w*(0 +1*randn(ch_num, 1));
w_2(1,1) = 1;
w_3 = coef_w*randn(ch_num, 1);
w_3(1,1) = 0;

% if use_bias ==0 && use_bias_2_end==0
%     w_3 = 0;
% end

W = zeros(ch_num, 2);
W(:, 1) = w_1;
W(:, 2) = w_2;
W(:, 3) = w_3;
save('W.txt', 'W', '-ascii')

end

% generate simulation data
function [xk_sample, Yk_sample, cnfg_sim] = gc_sim_data(dim, QQ_sim, L_init, KK, Y_k, ch_num, switch_num)
%% data simulation
cnfg_sim.dim = dim;
cnfg_sim.x0 = [0;0];
cnfg_sim.F = eye(2);
cnfg_sim.D = [0;0];
cnfg_sim.Q = QQ_sim;
cnfg_sim.L = L_init;
cnfg_sim.incongruent_vec = ones(KK, 1);
cnfg_sim.W = load('W.txt');
muu = mean(Y_k).';
cnfg_sim.mu = muu(1:ch_num, 1);

% call function to simulate data
[xk_sample, Yk_sample] = SIM_samp_x_yk_switch(cnfg_sim, switch_num);


str_save = sprintf('X_Y_sim%d.mat', dim);
save(str_save, 'xk_sample', 'Yk_sample', 'cnfg_sim')
end

% print info
function []= gc_print_fun_1(cfg_Init, ch_num)

win_sec_GCoh_empircal = cfg_Init.cfgGen.win_sec;
win_sec_observation = cfg_Init.cfgGen.win_sec / cfg_Init.cfgGen.seg_num;
fprintf('\n Number of channel: %d' , ch_num);
fprintf(' \n Length of Window for calc GCOh empirical: %d sec', win_sec_GCoh_empircal)
fprintf(' \n Length of subWindow for calc Cov mat of GCOh empirical (Window for SS-Coh): %d sec', win_sec_observation)
fprintf('\n------------------------------\n')

end

% generate initial vlaues parameters
function [init_Param, update_coef_W, W_sim] = gc_init_param(ParamaUpdate, cnfg_sim, coef_noise_mu, col_W_1, col_W_2, col_W_3, coef_noise_W, ch_num, dim)
count = 0;
%% Q
fprintf('\n Updating Parameters: \n')
init_Param = [];
if ParamaUpdate.sigma_update == 1
    count = count + 1;
    fprintf('%d- Sigma \n', count)

    init_Param.Q = 100*[1 0; 0 1].*10^-3;

elseif ParamaUpdate.sigma_update == 0
    init_Param.Q = cnfg_sim.Q;

end

if dim==1
    init_Param.Q = init_Param.Q(1, 1);
end

%% mu
if ParamaUpdate.mu_update == 1
    count = count + 1;
    fprintf('%d- MU \n', count)

    init_Param.mu = cnfg_sim.mu - coef_noise_mu.*(0.5 + 1i*0.5);

elseif ParamaUpdate.mu_update == 0
    init_Param.mu = cnfg_sim.mu;
end

%% L
if ParamaUpdate.L_update == 1
    count = count + 1;
    fprintf('%d- L \n', count)

    load('L1');
    init_Param.L  = L1;

elseif ParamaUpdate.L_update == 0
    init_Param.L  = cnfg_sim.L;
end

%% W
if ParamaUpdate.W_update == 1
    count = count + 1;
    fprintf('%d- W\n', count)

    W_sim = cnfg_sim.W;
    mw = size(W_sim);
    update_coef_W = zeros(mw);
    if exist('col_W_1', 'var')
        update_coef_W(col_W_1, 1) =1;
    end
    if exist('col_W_2', 'var')
        update_coef_W(col_W_2, 2) =1;
    end
    if exist('col_W_3', 'var')
        update_coef_W(col_W_3, 3) =1;
    end

    str_noise = sprintf('noise_%d_3.mat', ch_num);

    if isfile(str_noise)

        %         mu = 0; sigma = 1;
        %         pdN = makedist('Normal',mu,sigma);
        %         pdHN = truncate(pdN,0,inf);
        %         step = 1/ch_num;
        %         X = (4:step:5-step)';
        %         Xt = X(X>0);
        %         a_1 = pdf(pdHN,Xt);
        %         a_2 = pdf(pdHN,Xt);
        %         a_3 = pdf(pdHN,Xt);
        %         nn = [a_1, a_2, a_3];


        noise = randn(ch_num, 3);
        noise_source = noise;
        %         noise_source = load(str_noise);
    else
        noise = randn(ch_num, 3);
        str_save = sprintf('noise_%d_3.mat', ch_num);
        save(str_save, 'noise')
    end

    %     noise_source = cell2mat(struct2cell(noise_source));
    W_noise = coef_noise_W.*noise_source.*update_coef_W;

    W_inint = W_sim + W_noise;

    if dim==1
        W_inint(:, 2) = [];
        update_coef_W(:, 2) = [];
        W_sim(:, 2) = [];
    end

    init_Param.W = W_inint;

elseif ParamaUpdate.W_update == 0
    init_Param.W = cnfg_sim.W;
end
fprintf('------------------------------\n')

end

