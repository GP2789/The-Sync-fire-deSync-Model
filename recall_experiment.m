function [ filename ] = recall_experiment( filename, single_pattern, trials, timer )
% This function simulates the recall experiment simulated by the model. Use
% 'filename' to specify directory for data to be stored. Use
% 'single_pattern' as a 'yes' or 'no' answer to similate either a single
% trial with a random sequence (as in Figure 4), or to simulate all
% possible patterns dependent on the number of events and stimulation
% period (as in Figure 5). Use 'trials' to specify how many trials of each
% pattern to be simulated. 

%% DECLARATIONS
global weight_matrix; global p_w;
global weight_matrix_STDP; global wm_max; global tau_syn; global delay;

clearvars -global par;
set_parameters;
global par;
par.trials = trials;

%% determine how many combinations of events can fit into defined stimulation period
gap = 285; patterns_t = 0:gap:par.NC_stim_L; 
 
x = []; 
x_all = []; stop = 0;
while stop ~= 1
    x = zeros(1, length(patterns_t));
    y = 1:length(patterns_t);
    for i = 1:par.NC_n_stim
        i1 = ceil(rand()*length(y));
        x(y(i1)) = 1; y(i1) = [];
    end
    if(isempty(x_all) == 1 || ismember(length(patterns_t), sum(x_all == x,2)) ~= 1)
        x_all = [x_all; x]; c = 0;
    else; c = c + 1;
    end
    if(c == 1000); stop = 1; end
end

if(strcmp(single_pattern, 'yes') ~= 1) 
    % select all possible patterns in given window
    par.B = size(x_all, 1); x_all = sortrows(x_all);
else
    % or select random pattern for single trial
    x_all = x_all(ceil(rand()*size(x_all,1)), :); par.B = 1;
end
 
for i = 1:par.B; patterns{i} = find(x_all(i,:) == 1); end
if(par.B > 1); h3 = waitbar(0, 'Blink Experiment', 'Units', 'normalized', 'Position', [0.5 0.55 0.2 0.1]); end

%% rename filename & create directories
filename = [filename '_' int2str(par.B) 'P_' int2str(trials) 'T'];
if(exist(filename, 'dir')==0); mkdir(filename); end
if(exist([filename '/Data'], 'dir')==0); mkdir([filename '/Data']); end

par.rand_phase = zeros(length(par.sim_order), trials, par.B); par.pre_stim_L = par.NC_pre_stim_L;
%% MAIN SIMULATION
for b = 1:par.B
    if(par.B > 1); B_n = ['/B' int2str(b)]; else; B_n = []; end
    if(trials>1); h1 = waitbar(0, 'Trials', 'Units', 'normalized', 'Position', [0.5 0.4 0.2 0.1]); end
    for t = 1:trials
        if(length(par.sim_order)>1); h2 = waitbar(0, 'Recall Experiment', 'Units', 'normalized', 'Position', [0.5 0.25 0.2 0.1]); end
        for n = 1:length(par.sim_order)
            %% CREATE NETWORK
            if(exist([filename '/Data' B_n '/T' int2str(t)], 'dir')==0); mkdir([filename '/Data' B_n '/T' int2str(t)]); end
            par.sim_type = par.sim_order{n};
            if(length(par.sim_order)>1); waitbar(n/length(par.sim_order),h2); end

            % NEW RANDOM PHASE FOR EVENT STIMULATIONS
            % vary events to add stochasticity for every trial/pattern
            par = rmfield(par, 'NC_pre_stim_L'); par = rmfield(par, 'TC_pre_stim_L'); set_parameters;
            par.rand_phase(n, t, b) = round(rand()*13)*10; 
            par.NC_pre_stim_L = par.NC_pre_stim_L + par.rand_phase(n, t, b);
            par.TC_pre_stim_L = par.TC_pre_stim_L + par.rand_phase(n, t, b);
            par.sim_length = par.NC_pre_stim_L*2 + par.NC_stim_L;
            
            par.NC_stims_t = patterns_t( patterns{ b } ); 
            
            create_network();
            
            % RECORD EACH SIMULATED PATTERN FOR EVERY TRIAL
            if(sum(strcmp('NC', par.nG)) > 0)
                par.NC_stims_t_all{b}{t}{n} = par.NC_stims_t;
                par.NC_stims_N_all{b}{t}{n} = par.NC_stims_N;
            end

            % assign connection variables from previous part of simulation
            if(n>1); weight_matrix_STDP = stdp; wm_max = mx; delay = del; tau_syn = ts; end
            if(strcmp(par.sim_type, 'recall')==1)
               TC_E_L1 = find(cellfun(@isempty,strfind(par.n_ID,'TC_E_L1'))==0);
               TC_E_L2 = find(cellfun(@isempty,strfind(par.n_ID,'TC_E_L2'))==0);
               BP_E = find(cellfun(@isempty,strfind(par.n_ID,'BP_E'))==0);
               NC_E = find(cellfun(@isempty,strfind(par.n_ID,'NC_E'))==0);
               
               % APPLY PREVIOUS BINDING POOL LEARNING
               weight_matrix(BP_E, BP_E) = wm(BP_E, BP_E);
               % APPLY DIRECTIONAL NEUROMODULATORS
               % TC -> BP +
               weight_matrix(TC_E_L1, BP_E) = weight_matrix(TC_E_L1, BP_E) + p_w(TC_E_L1, BP_E) * par.TC_L1_BP_NM;
               wm_max(TC_E_L1, BP_E) = wm_max(TC_E_L1, BP_E) + par.TC_L1_BP_NM;
               weight_matrix(TC_E_L2, BP_E) = weight_matrix(TC_E_L2, BP_E) + p_w(TC_E_L2, BP_E) * par.TC_L2_BP_NM;
               wm_max(TC_E_L2, BP_E) = wm_max(TC_E_L2, BP_E) + par.TC_L2_BP_NM;
               % BP -> NC +
               weight_matrix(BP_E, NC_E) = weight_matrix(BP_E, NC_E) + p_w(BP_E, NC_E) * par.BP_NC_NM;
               wm_max(BP_E, NC_E) = wm_max(BP_E, NC_E) + par.BP_NC_NM;
               % set NC -> BP weights to zero during recall
               weight_matrix(NC_E, BP_E) = 0; 
            end

            %% SIMULATE NETWORK
            [sim_data, ~, stop] = simulate_network(filename, timer);
            if(stop ~= 1)
                % carry forward connection variables
                wm = weight_matrix; stdp = weight_matrix_STDP; mx = wm_max; del = delay; ts = tau_syn;
                if(par.B > 1); B_n = ['/B' int2str(b)]; else; B_n = []; end
                % SAVE DATA
                save([filename '/Data' B_n '/T' int2str(t) '/' par.sim_order_n{n} '.mat'],'sim_data', '-v7.3')
            end
        end
        if(trials>1); waitbar(t/trials,h1); end
        if(length(par.sim_order)>1); close(h2); end
    end
   if(par.B > 1); waitbar(b/par.B,h3); end
end
if(stop ~= 1)
    save([filename '/Data/parameters.mat'],'par')
    if(trials>1);close(h1);end
end
end


