function [  ] = preprocess_data( filename )
%% Load and initialise
if(exist([filename '/Processed Data'], 'dir')==0); mkdir([filename '/Processed Data']); end
load([filename '/Data/parameters' ]);
fprintf('\npreprocessing data ...\n'); fpr = false;
time = nan(par.B,par.trials);
T_BS = min(750, par.NC_pre_stim_L);
T_BS = T_BS : min(T_BS + 500, par.NC_pre_stim_L-50);
phase_BS = repmat(-pi:0.0492:pi, [1 3]);
% initialise LFP/TFA parameters
TFA.UBF = 1; TFA.OBF = 30; TFA.SR = 1000; TFA.Gamma = 0.5;  pow_s = [3 13];
%% loop through patterns
for b = 1:par.B
    if(par.B > 1); P_n = ['/B' int2str(b)]; else; P_n = [];  end
    for t = 1:par.trials % loop through trials
        tic
        if(exist([filename '/Processed Data/WM_data_P' int2str(b) '.mat'],'file') ~= 2 || ...
                exist([filename '/Processed Data/spike_data_P' int2str(b) '.mat'],'file') ~= 2 || ...
                exist([filename '/Processed Data/calcium_data_P' int2str(b) '.mat'],'file') ~= 2 || ...
                exist([filename '/Processed Data/freq_data_P' int2str(b) '.mat'],'file') ~= 2)
            for n = 1:length(par.sim_order_n) % loop through simulation stages
                %% LOAD & INITIALISE
                load([filename '/Data' P_n '/T' int2str(t) '/' par.sim_order_n{n}]);
                SD = sim_data.spike_detector; %V_m = sim_data.V_m;
                WM = sim_data.weight_matrix; p_w = sim_data.WM_p;
                TS = sim_data.TS;
                if(isfield(par,'pre_stim_L')); par.sim_length = par.pre_stim_L * 2 + par.NC_stim_L;
                else; par.sim_length = par.NC_pre_stim_L * 2 + par.NC_stim_L;
                end
                if(isfield(par,'rand_phase'))
                    r_phase = par.rand_phase{b}{t}{n};
                else; r_phase = 0;
                end
                sim_period = r_phase + 1 : par.sim_length + r_phase;
                SD(:,2) = SD(:,2) - r_phase;
                SD(SD(:,2)<1,:) = [];
                SD(SD(:,2)>par.sim_length,:) = [];
                %% EXTRACT MEAN WEIGHT MATRICES
                if(exist([filename '/Processed Data/WM_data_P' int2str(b) '.mat'],'file') ~= 2)
                    if(t==1) % INITIALISE
                        WM_all.([par.sim_order_n{n}]) = zeros(par.network_size, par.network_size);
                        WM_p_all.([par.sim_order_n{n}]) = zeros(par.network_size, par.network_size);
                        TS_all.([par.sim_order_n{n}]) = zeros(par.network_size, par.network_size);
                    end
                    % INCREMENT
                    WM_all.([par.sim_order_n{n}]) = WM_all.([par.sim_order_n{n}]) + WM;
                    WM_p_all.([par.sim_order_n{n}]) = WM_p_all.([par.sim_order_n{n}]) + p_w;
                    TS_all.([par.sim_order_n{n}]) = TS_all.([par.sim_order_n{n}]) + TS;
                    % AVERAGE
                    if(t==par.trials)
                        WM_all.([par.sim_order_n{n}]) = WM_all.([par.sim_order_n{n}]) / par.trials;
                        WM_p_all.([par.sim_order_n{n}]) = WM_p_all.([par.sim_order_n{n}]) / par.trials;
                        TS_all.([par.sim_order_n{n}]) = TS_all.([par.sim_order_n{n}]) / par.trials;
                        if(n==length(par.sim_order_n)); save([filename '/Processed Data/WM_data_P' int2str(b) '.mat'],...
                                'WM_all', 'WM_p_all', 'TS_all'); end
                    end
                end
                %% EXTRACT SPIKE TIMES
                if(exist([filename '/Processed Data/spike_data_P' int2str(b) '.mat'],'file') ~=2)
                    groups = unique(par.n_ID);
                    for g = 1:length(groups)
                        SPIKES.([par.sim_order_n{n}]).(['T' int2str(t)]).(groups{g}) = SD(ismember(SD(:,1),find(contains(par.n_ID,groups{g}))),:);
                    end
                    if(n==length(par.sim_order_n) && t == par.trials); save([filename '/Processed Data/spike_data_P' int2str(b) '.mat'], 'SPIKES'); end
                end
                %% EXTRACT CALCIUM AND WEIGHT CHANGE
                if(exist([filename '/Processed Data/calcium_data_P' int2str(b) '.mat'],'file') ~= 2)
                    if(n == 1 && ~isempty(sim_data.WM))
                        stims_N = par.NC_stims_N_all{b}{t};
                        stims_T = par.NC_stims_t_all{b}{t} + par.pre_stim_L; win = 250;
                        BP_E = find(contains(par.n_ID,'BP_E'));
                        TC_L1_E = find(contains(par.n_ID,'TC_E_L1')); 
                        per_G = numel(TC_L1_E) / par.n_TC_P(1);
                        TC_L1_E = mat2cell(TC_L1_E,ones(numel(TC_L1_E)/per_G,1)*per_G, 1);
                        TC_L2_E = find(contains(par.n_ID,'TC_E_L2'));
                        TC_L2_E = mat2cell(TC_L2_E,ones(numel(TC_L2_E)/per_G,1)*per_G, 1);
                        k = extract_IDs(size(WM), sim_data.WM, BP_E, BP_E);
                        p_w_all.BP_capacity = mean(sim_data.p_w(k, sim_period),1);
                        BP_G_ID = cell(length(stims_N), 1);
                        NC_G_ID = cell(length(stims_N), 1);
                        for i = 1:length(stims_N)
                            %% BP -> NC weights
                            NC_G_ID{i} = stims_N(i)-ceil(par.NC_stim_N_cos):stims_N(i)+ceil(par.NC_stim_N_cos);
                            k = extract_IDs(size(WM), sim_data.WM, BP_E, NC_G_ID{i});
                            p_w = sim_data.p_w(k, sim_period); p_w = p_w(:, stims_T(i):stims_T(i)+win);
                            p_w = p_w(:,2:end) - p_w(:,1:end-1); k = k(sum(p_w,2) > 0.5); clear p_w; 
                            calcium_all.BP_NC{i,1}(t,:) = sum(sim_data.calcium(k, sim_period),1)/numel(k);
                            p_w_all.BP_NC{i,1}(t,:) = sum(sim_data.p_w(k, sim_period),1)/numel(k);
                            %% BP <-> BP weights
                            [BP_i, ~] = ind2sub([1 1]*par.network_size,sim_data.WM(k));
                            BP_G_ID{i} = unique(BP_i); 
                            k = extract_IDs(size(WM), sim_data.WM, BP_G_ID{i}, BP_G_ID{i});
                            calcium_all.BP_BP{i,1}(t,:) = sum(sim_data.calcium(k, sim_period),1)/numel(k);
                            p_w_all.BP_BP{i,1}(t,:) = sum(sim_data.p_w(k, sim_period),1)/numel(k);
                            %% TC L1 -> BP weights
                            for j = 1:length(TC_L1_E)
                                k = extract_IDs(size(WM), sim_data.WM, TC_L1_E{j}, BP_G_ID{i});
                                p_w = sim_data.p_w(k, sim_period); p_w = p_w(:, stims_T(i):stims_T(i)+win);
                                p_w = p_w(:,2:end) - p_w(:,1:end-1); k = k(sum(p_w,2) > 0.25); clear p_w; 
                                calcium_all.TCL1_BP{i,j}(t,:) = sum(sim_data.calcium(k, sim_period),1)/numel(k);
                                p_w_all.TCL1_BP{i,j}(t,:) = sum(sim_data.p_w(k, sim_period),1)/numel(k);
                            end
                            %% TC L2 -> BP weights
                            for j = 1:length(TC_L2_E)
                                k = extract_IDs(size(WM), sim_data.WM, TC_L2_E{j}, BP_G_ID{i});
                                p_w = sim_data.p_w(k, sim_period); p_w = p_w(:, stims_T(i):stims_T(i)+win);
                                p_w = p_w(:,2:end) - p_w(:,1:end-1); k = k(sum(p_w,2) > 0.25); clear p_w; 
                                calcium_all.TCL2_BP{i,j}(t,:) = sum(sim_data.calcium(k, sim_period),1)/numel(k);
                                p_w_all.TCL2_BP{i,j}(t,:) = sum(sim_data.p_w(k, sim_period),1)/numel(k);
                            end
                        end
                        if(t == par.trials); save([filename '/Processed Data/calcium_data_P' int2str(b) '.mat'], ...
                                'calcium_all', 'p_w_all', 'BP_G_ID', 'NC_G_ID'); end
                    end
                end
                
                %% EXTRACT PHASE & POWER DATA
                if(exist([filename '/Processed Data/freq_data_P' int2str(b) '.mat'],'file') ~= 2)
                    %load([filename '/Processed Data/spike_data_P' int2str(b) '.mat']);
                    % create LFP
                    [LFP] = create_LFP(SPIKES.([par.sim_order_n{n}]).(['T' int2str(t)]).NC_E, 1, 1, pow_s, par.sim_length, 1);
                    % time-frequency analysis
                    [P, ~] = GaborFilter(LFP, TFA);
                    freq.NC = TFA.UBF : (TFA.OBF-TFA.UBF)/(size(P,1)-1) : TFA.OBF;
                    F = find(freq.NC>=pow_s(1),1,'first') : find(freq.NC<=pow_s(2),1,'last');
                    power = mean(abs(P(F,:))); power_BS = mean(power(T_BS));
                    power = ((power - power_BS) / power_BS) * 100;
                    phase = angle(mean(P(F,:)));
                    % find intrinsic frequency & pre-stim stationary period
                    dPh = phase(750:par.NC_pre_stim_L-250);
                    dPh = dPh(2:end) - dPh(1:end-1);
                    pi_cycle = round((2*pi)/mean(dPh(dPh > 0)));
                    phase_BS = repmat(-pi:(2*pi)/pi_cycle:pi, [1 3]);
                    % extract non-stationarity of signal
                    phase_NS = phase_compare( phase, phase_BS );
                    diff = length(phase) - length(phase_NS);
                    phase_NS = [ones(1, ceil(diff/2)) phase_NS ones(1, floor(diff/2))];
                    PHASE.([par.sim_order_n{n}]).NC(t,:) = phase;
                    PHASE_NS.([par.sim_order_n{n}]).NC(t,:) = phase_NS;
                    % extract power of signal
                    POWER.([par.sim_order_n{n}]).NC(t,:) = power;
                    % initialise
                    if(t==1); POW_all.([par.sim_order_n{n}]).NC = zeros(size(P)); end
                    % increment
                    LFP_all.([par.sim_order_n{n}]).NC(t,:) = LFP;
                    POW_all.([par.sim_order_n{n}]).NC = POW_all.([par.sim_order_n{n}]).NC + abs(P);
                    % average
                    if(t == par.trials)
                        POW_all.([par.sim_order_n{n}]).NC = POW_all.([par.sim_order_n{n}]).NC / par.trials;
                        if(n==length(par.sim_order_n)); save([filename '/Processed Data/freq_data_P' int2str(b) '.mat'], ...
                                'freq', 'LFP_all','POW_all', 'PHASE_NS', 'PHASE', 'POWER'); end
                    end
                end
            end
        end
        %% TIMING & PROGRESS
        B = dir([filename '/Processed Data']);
        B_rem = par.B - sum(contains({B(:).name}, 'calcium_data'));
        time(b,t) = toc;
        min_rem = (nanmean(nanmean(time))/60) * ((B_rem)*par.trials + (par.trials-t));
        H = floor(min_rem/60);
        M = floor(min_rem - H*60);
        per_done = ((par.B-B_rem)*par.trials + t)/(par.B*par.trials)*100;
        if(~fpr); fprintf('\b\b\b\b%3.0f%% time remaining %2.0fh %2.0fm\n', per_done, H, M); fpr = true;
        else; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f%% time remaining %2.0fh %2.0fm\n', per_done, H, M);
        end
    end
end
end

%% EXTRACT SYNAPSE ID FUNCTION
function [syn_ID] = extract_IDs(WM_size, STDP_on, pre, post)
x = zeros(WM_size);
x(pre,:) = x(pre,:)+1; % pre
x(:,post) = x(:,post)+1; % post
syn_ID = find(ismember(STDP_on, find(x==2))==1); % find synapse IDs
end

