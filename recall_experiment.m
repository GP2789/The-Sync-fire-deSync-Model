function [ filename ] = recall_experiment( filename, n_patterns, n_trials, set_par_N, set_par_x  )
%% rename filename & create directories
filename = [filename '_' int2str(n_patterns) 'P_' int2str(n_trials) 'T'];
if(exist(filename, 'dir')==0); mkdir(filename); end
if(exist([filename '/Data'], 'dir')==0); mkdir([filename '/Data']); end

NC_stims_t_all = cell(n_patterns,1);
NC_stims_N_all = cell(n_patterns,1);
rand_phase = cell(n_patterns,1);
for b = 1:n_patterns
    NC_stims_t_all{b} = cell(n_trials,1);
    NC_stims_N_all{b} = cell(n_trials,1);
    rand_phase{b} = cell(n_trials,1);
    for t = 1:n_trials
        rand_phase{b}{t} = cell(2,1);
    end
end

%% MAIN SIMULATION
par_tmp = cell(n_patterns,1);
parfor b = 1 : n_patterns
    %% initialise parameters
    par = struct;
    par = set_parameters(par);
    par.trials = n_trials; par.B = n_patterns;
    for n = 1:length(set_par_N)
        if(ischar(set_par_x{n}))
            par.(set_par_N{n}) = {set_par_x{n}};
        else; par.(set_par_N{n}) = set_par_x{n};
        end
    end
    par.pre_stim_L = par.NC_pre_stim_L;
    if(par.B > 1); B_n = ['/B' int2str(b)]; else; B_n = []; end
    for t = 1:n_trials
        %% SIMULATE LEARNING
        tic
        fprintf('\nsimulating P (%.0f/%.0f) T (%.0f/%.0f) learning phase ...\n',b,par.B,t,par.trials)
        if(exist([filename '/Data' B_n '/T' int2str(t)], 'dir')==0); mkdir([filename '/Data' B_n '/T' int2str(t)]); end
        par.sim_type = par.sim_order{1};
        % random phase
        par = rmfield(par, 'NC_pre_stim_L'); par = rmfield(par, 'TC_pre_stim_L'); 
        par = set_parameters(par);
        r = round(rand()*13)*10; rand_phase{b}{t}{1} = r;
        par.NC_pre_stim_L = par.NC_pre_stim_L + r;
        par.TC_pre_stim_L = par.TC_pre_stim_L + r;
        par.sim_length = par.NC_pre_stim_L*2 + par.NC_stim_L;
        par.NC_stims_t = sort(round(rand(1, par.NC_n_stim)*par.NC_stim_L), 2);
        
        % create network
        [ I, weight_matrix, weight_matrix_STDP, tau_syn, wm_max, delay, par ] = create_network(par);
        NC_stims_t_all{b}{t} = par.NC_stims_t;
        NC_stims_N_all{b}{t} = par.NC_stims_N;
        % simulate
        [sim_data, ~, ~] = simulate_network(filename, 'off', ...
            par, I, weight_matrix, weight_matrix_STDP, tau_syn, wm_max, delay);
        % save
        parsave([filename '/Data' B_n '/T' int2str(t) '/' par.sim_order_n{1} '.mat'], sim_data);
        
        %% SIMULATE RECALL
        fprintf('\nsimulating P (%.0f/%.0f) T (%.0f/%.0f) recall phase ...\n',b,par.B,t,par.trials)
        if(exist([filename '/Data' B_n '/T' int2str(t)], 'dir')==0); mkdir([filename '/Data' B_n '/T' int2str(t)]); end
        par.sim_type = par.sim_order{2};
        % random phase
        par = rmfield(par, 'NC_pre_stim_L'); par = rmfield(par, 'TC_pre_stim_L'); 
        par = set_parameters(par);
        r = round(rand()*13)*10; rand_phase{b}{t}{2} = r;
        par.NC_pre_stim_L = par.NC_pre_stim_L + r;
        par.TC_pre_stim_L = par.TC_pre_stim_L + r;
        par.sim_length = par.NC_pre_stim_L*2 + par.NC_stim_L;
        
        % create netowk
        [ I, ~, ~, ~, ~, ~, par] = create_network(par);
        % carry forward weights from previous simulation
        weight_matrix = sim_data.weight_matrix;
        % extract IDs for neuronal groups
        TC_E_L1 = find(cellfun(@isempty,strfind(par.n_ID,'TC_E_L1'))==0);
        TC_E_L2 = find(cellfun(@isempty,strfind(par.n_ID,'TC_E_L2'))==0);
        BP_E = find(cellfun(@isempty,strfind(par.n_ID,'BP_E'))==0);
        NC_E = find(cellfun(@isempty,strfind(par.n_ID,'NC_E'))==0);
        % reverse neuromodulated pathways for recall
        % TC -> BP +
        weight_matrix(TC_E_L1, BP_E) = weight_matrix(TC_E_L1, BP_E) + sim_data.WM_p(TC_E_L1, BP_E) * par.TC_L1_BP_NM;
        wm_max(TC_E_L1, BP_E) = wm_max(TC_E_L1, BP_E) + par.TC_L1_BP_NM;
        weight_matrix(TC_E_L2, BP_E) = weight_matrix(TC_E_L2, BP_E) + sim_data.WM_p(TC_E_L2, BP_E) * par.TC_L2_BP_NM;
        wm_max(TC_E_L2, BP_E) = wm_max(TC_E_L2, BP_E) + par.TC_L2_BP_NM;
        % BP -> NC +
        weight_matrix(BP_E, NC_E) = weight_matrix(BP_E, NC_E) + sim_data.WM_p(BP_E, NC_E) * par.BP_NC_NM;
        wm_max(BP_E, NC_E) = wm_max(BP_E, NC_E) + par.BP_NC_NM;
        % NC -> BP -
        weight_matrix(NC_E, BP_E) = 0;
        
        % simulate
        [sim_data, ~, ~] = simulate_network(filename, 'off', ...
            par, I, weight_matrix, weight_matrix_STDP, tau_syn, wm_max, delay);
        % save
        parsave([filename '/Data' B_n '/T' int2str(t) '/' par.sim_order_n{2} '.mat'], sim_data);
        M = floor(toc/60);
        S = toc - M*60;
        fprintf('\nsimulation of P (%.0f/%.0f) T (%.0f/%.0f) took %.0fM %.0fS ...\n',b,par.B,t,par.trials, M,S)
    end
    par_tmp{b} = par;
end
par = par_tmp{1};
par.rand_phase = rand_phase;
par.NC_stims_t_all = NC_stims_t_all;
par.NC_stims_N_all = NC_stims_N_all;
save([filename '/Data/parameters.mat'],'par')
patterns_t = vertcat(NC_stims_t_all{:});
patterns_t = vertcat(patterns_t{:});
save([filename '/Data/patterns.mat'],'patterns_t')

end

function [] = parsave(filename, sim_data); save(filename, 'sim_data'); end


