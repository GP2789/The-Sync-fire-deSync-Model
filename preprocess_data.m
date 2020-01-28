function [  ] = preprocess_data( filename )
%% Load and initialise
if(exist([filename '/Processed Data'], 'dir')==0); mkdir([filename '/Processed Data']); end
load([filename '/Data/parameters' ]); 
SPIKES = [];
if(par.B*par.trials>1); h = waitbar(0, 'LOADING DATA', 'Units', 'normalized', 'Position', [0.5 0.55 0.2 0.1]); end

for b = 1:par.B % loop through patterns
    if(par.B > 1); P_n = ['/P' int2str(b)]; else; P_n = [];  end
    for t = 1:par.trials % loop through trials
        if(par.B*par.trials>1); waitbar(((b-1)*par.B+t)/(par.B*par.trials), h); end
        for n = 1:length(par.sim_order_n) % loop through simulation stages
            %% LOAD & INITIALISE
            disp('Loading data ...'); 
            load([filename '/Data' P_n '/T' int2str(t) '/' par.sim_order_n{n}]);
            
            SD = sim_data.spike_detector; V_m = sim_data.V_m; 
            WM = sim_data.weight_matrix; p_w = sim_data.WM_p;
            TS = sim_data.TS; 
            I = sim_data.I; 
            TFA.SR = 1000; TFA.Gamma = 0.5; 
            N = length(par.nG);
            par.sim_length = par.pre_stim_L * 2 + par.NC_stim_L;
            r_phase = (size(V_m,2) - (par.pre_stim_L * 2 + par.NC_stim_L))/2;
            sim_period = r_phase + 1 : par.sim_length + r_phase;
            c_L = 0; T = sim_period - r_phase;
            
            %% EXTRACT STDP WEIGHT CHANGE DATA
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
            end 
            
            %% EXTRACT SPIKE DATA
            n_G_ID = []; S_uK = []; S_nK = [];
            for g=1:N % loop through regions
                for L = 1:par.(['n_SG_' par.nG{g}]) % loop through layers
                    for g2 = 1:length(par.([par.nG{g} '_N'])) % loop through neuron type
                        %% extract layer spikes
                        spikes_L = SD(SD(:,1)<=par.(['n_' par.nG{g} '_' par.([par.nG{g} '_N']){g2}])(L) + c_L,:); 
                        spikes_L = spikes_L(spikes_L(:,1)>c_L,:); 
                        spikes_L = spikes_L(spikes_L(:,2) > min(sim_period),:); 
                        spikes_L = spikes_L(spikes_L(:,2) < max(sim_period),:);
                        spikes_L(:,2) = spikes_L(:,2) - r_phase; 
                        
                        if(strcmp(par.nG{g},'TC')) 
                            % cluster synfire spikes by neuron type, based on
                            % pre-determined synfire chain parameters
                           S_nK{g}{L}{g2} = par.n_TC_P(L);
                           n_L = par.(['n_' par.nG{g} '_' par.([par.nG{g} '_N']){g2}])(L) / par.n_TC_P(L);
                           for i = 1 : par.n_TC_P(L)
                               n_I = n_L*(i-1)+c_L+1 : n_L*i+c_L;
                               spikes_L2 = spikes_L(ismember(spikes_L(:,1), n_I), :);
                               n_G_ID{g}{L}{g2}{i} = n_I;
                               SPIKES{t}{n}{g}{L}{g2}{i} = spikes_L2;
                           end
                           
                        elseif(strcmp(par.nG{g},'BP') && strcmp(par.([par.nG{g} '_N']){g2}, 'E'))
                            % use kmeans to cluster binding pool excitatory spikes
                            if(n == 1) % just for encoding (presumed same n k at recall)
                                for K = fliplr(1:par.NC_n_stim)
                                    % kmeans into k clusters
                                    idx = kmeans(spikes_L, K);
                                    % count spikes in each cluster
                                    C = zeros(K,1);
                                    for i = 1 : K; C(i) = length(idx(idx == i)); end
                                    % if <=5 spikes in any 1 cluster then use a smaller k
                                    if(isempty(find(C <= 5))); break; end
                                end
                            else; idx = kmeans(spikes_L, K);
                            end
                            % store spikes into k clusters
                            for i = 1 : K
                                n_G_ID{g}{L}{g2}{i} = unique(spikes_L(idx == i, 1));
                                SPIKES{t}{n}{g}{L}{g2}{i} = spikes_L(idx == i, :);
                            end
                        else
                            % lump all NC spikes into 1 cluster.
                            S_nK{g}{L}{g2} = 1;
                            n_G_ID{g}{L}{g2}{1} = unique(spikes_L(:,1));
                            SPIKES{t}{n}{g}{L}{g2}{1} = spikes_L;
                        end
                        c_L = c_L + par.(['n_' par.nG{g} '_' par.([par.nG{g} '_N']){g2}])(L);
                    end
                end
            end
            c = 0;
            
            %% EXTRACT STDP / NEURON INTRINSIC DATA
            for g=1:N % loop through regions
                Sg = [];
                for L = 1:par.(['n_SG_' par.nG{g}]) % loop through layers
                    if(par.(['n_SG_' par.nG{g}])<=1); nL = ''; else; nL = ['_L' int2str(L)]; end
                    for g2=1:length(par.([par.nG{g} '_N'])) % loop through neuron type
                        for i=1:length(n_G_ID{g}{L}{g2}) % loop through clusters
                            %% STORE DATA FOR STDP CALCIUM CHANGES
                            % need to loop through region/layer/neuron type/cluster again to extract pairings
                            for g_2 = 1:N 
                                for L2 = 1:par.(['n_SG_' par.nG{g_2}])
                                    if(par.(['n_SG_' par.nG{g_2}])<=1); nL2 = ''; else; nL2 = ['_L' int2str(L2)]; end
                                    for g2_2=1:length(par.([par.nG{g_2} '_N']))
                                        for i2 = 1:length(n_G_ID{g_2}{L2}{g2_2})
                                            n_pre = n_G_ID{g}{L}{g2}{i}; % pre syn ID
                                            n_post = n_G_ID{g_2}{L2}{g2_2}{i2}; % post syn ID
                                            bx = 10; % filter data with 10ms bin width
                                            % identify synapses with STDP enabled using ID against stored weight matrix
                                            x = zeros(size(WM)); x(n_pre,:) = x(n_pre,:)+1; x(:,n_post) = x(:,n_post)+1;
                                            k = find(ismember(sim_data.WM, find(x==2))==1); 
                                            % initialise
                                            if(t==1) 
                                                calcium_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} = zeros(2, round(max(T)/bx));
                                                p_w_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} = zeros(2, round(max(T)/bx));
                                            end
                                            % if snapses with STDP exist in inter-region connection
                                            if(isempty(k) ~= 1) 
                                                C = sim_data.calcium(k,sim_period); % extract calcium over time
                                                P = sim_data.p_w(k,sim_period); % extract weights over time
                                                if(size(C,1) == 1); calc_max = C; else; calc_max = max(C); end
                                                P_M = mean(P,1); % mean of weights
                                                P_diff = [0 P_M(2:end) - P_M(1:end-1)]; % gradient change
                                                P_diff = P_diff / max([max(P_diff) abs(min(P_diff))]); % ^ normalised
                                                calc_M = mean(C,1) + std(C); % mean of calcium
                                                % define filtered data stores
                                                P_M_f = zeros(1, round(max(T) / bx)); P_diff_f = zeros(1, round(max(T) / bx));
                                                calc_M_f = zeros(1, round(max(T) / bx)); calc_max_f = zeros(1, round(max(T) / bx)); 
                                                T_f = bx/2:bx:max(T)-bx/2;
                                                for x_f = 1:length(P_diff_f) % downsample data
                                                   P_diff_f(x_f) = mean(P_diff((x_f-1)*bx+1:x_f*bx));
                                                   P_M_f(x_f) = mean(P_M((x_f-1)*bx+1:x_f*bx));
                                                   calc_M_f(x_f) = mean(calc_M((x_f-1)*bx+1:x_f*bx));
                                                   calc_max_f(x_f) = mean(calc_max((x_f-1)*bx+1:x_f*bx));
                                                end
                                                % increment 
                                                calcium_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} ...
                                                    = calcium_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} ...
                                                    + [calc_M_f; calc_max_f]; 
                                                p_w_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} = ...
                                                    p_w_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} ...
                                                    + [P_M_f; P_diff_f;];
                                                % average
                                                if(t==par.trials)
                                                    calcium_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} = ...
                                                        calcium_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} / par.trials;
                                                    p_w_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} = ...
                                                        p_w_all.(['B' int2str(b)]).([par.sim_order_n{n}]){g, g_2}{L, L2}{g2, g2_2}{i, i2} / par.trials;
                                                end
                                            end 
                                        end
                                    end
                                end
                            end
                            %% STORE DATA FOR INTRINSIC PROPERTIES / CURRENTS
                            % STORE VOLTAGE DATA
                            % initialise
                            if(t==1 || isfield(V_all.([par.sim_order_n{n}]).([par.nG{g}]),['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) == 0) % initialise 
                                V_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                    zeros(1, par.sim_length);
                            end
                            % incrememnt
                            V_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                V_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) + ...
                                V_m(n_G_ID{g}{L}{g2}{i},sim_period);
                            % average
                            if(t==par.trials)
                                V_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                    V_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) / par.trials;
                            end
                            
                            % STORE CURRENT DATA
                            % initialise
                            if(t==1 || isfield(Na_all.([par.sim_order_n{n}]).([par.nG{g}]),['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) == 0) % initialise 
                                if(strcmp(par.model_type, 'HH')==1)
                                    Na_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                        zeros(1, par.sim_length);
                                    K_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                        zeros(1, par.sim_length);
                                end
                                L_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                    zeros(1, par.sim_length);
                                SYN_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                    zeros(1, par.sim_length, length(par.uTS));
                            end
                            % incrememnt
                            if(strcmp(par.model_type, 'HH')==1)
                                Na_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                    Na_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) + ...
                                    mean(I.Na(n_G_ID{g}{L}{g2}{i}, sim_period),1);
                                K_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                    K_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) +  ...
                                    mean(I.K(n_G_ID{g}{L}{g2}{i}, sim_period),1);
                            end
                            L_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                L_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) +...
                                mean(I.L(n_G_ID{g}{L}{g2}{i}, sim_period),1);
                            SYN_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                SYN_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) + ...
                                mean(I.SYN(n_G_ID{g}{L}{g2}{i}, sim_period,:),1);
                            
                            % average
                            if(t==par.trials)
                                if(strcmp(par.model_type, 'HH')==1)
                                    Na_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                        Na_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) / par.trials;
                                    K_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                        K_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) / par.trials;
                                end
                                L_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                    L_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) / par.trials;
                                SYN_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) = ...
                                    SYN_all.([par.sim_order_n{n}]).([par.nG{g}]).(['L' int2str(L) '_G' int2str(i) '_' par.([par.nG{g} '_N']){g2}]) / par.trials;
                            end
                            
                            %% concatenate spikes
                            Sg = [Sg; SPIKES{t}{n}{g}{L}{g2}{i}];
                        end
                    end
                end
                %% pre process LFP and time frequency analysis for NC
                if(strcmp(par.nG{g},'NC')==1)
                    % initialise LFP/TFA parameters
                    pow_s = [3 13]; hann_on = 1; sample_on = 1;
                    TFA.UBF = 1; TFA.OBF = 30; 
                    pre_stim_L = min(par.pre_stim_L, 500);
                    % create LFP
                    if(isempty(Sg)~=1)
                        Sg(:,1) = Sg(:,1) - min(Sg(:,1)) + 1;
                        [LFP] = create_LFP(Sg, max(Sg(:,1)), hann_on, sample_on, pow_s, par.sim_length, 1);
                    else; LFP = zeros(1, par.sim_length);
                    end
                    % time-frequency analysis
                    [P, ~] = GaborFilter(LFP, TFA);
                    freq.([par.nG{g}]) = TFA.UBF : (TFA.OBF-TFA.UBF)/(size(P,1)-1) : TFA.OBF; 
                    F = min(find(freq.([par.nG{g}])>=pow_s(1))) : max(find(freq.([par.nG{g}])<=pow_s(2)));
                    phase = angle(mean(P(F,:)));
                    y_p = phase(max(pre_stim_L/2,par.pre_stim_L/2)-pre_stim_L/2+1:max(pre_stim_L/2,par.pre_stim_L/2)+pre_stim_L/2);
                    xy_d = phase_compare( phase, y_p );
                    PHASE.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]).(['T' int2str(t)]) = phase;
                    POWER.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]).(['T' int2str(t)]) = mean(abs(P(F,:)));
                    
                    % initialise 
                    if(t==1) 
                        LFP_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = zeros(1, par.sim_length);
                        TFA_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = zeros(size(P));
                        POW_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = zeros(size(P));
                    end
                    % increment
                    LFP_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = ...
                        LFP_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) + LFP;
                    TFA_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = ...
                        TFA_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) + P;
                    POW_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = ...
                        POW_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) + abs(P);
                    %TFA_trials.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]).(['T' int2str(t)]) = P; clear('P');
                    % average
                    if(t == par.trials) 
                        LFP_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = ...
                            LFP_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) / par.trials;
                        TFA_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = ...
                            TFA_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) / par.trials;
                        POW_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = ...
                            POW_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) / par.trials;
                        %PHASE_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)]) = ...
                            %angle(mean(TFA_all.([par.sim_order_n{n}]).([par.nG{g} nL]).(['B' int2str(b)])(F,:)));
                    end
                end
            end
        end
    end
    save([filename '/Processed Data/WM_data_P' int2str(b) '.mat'], 'WM_all', 'WM_p_all', 'TS_all', '-v7.3')
    save([filename '/Processed Data/spike_data_P' int2str(b) '.mat'], 'SPIKES', 'n_G_ID', '-v7.3')
    save([filename '/Processed Data/intrinsic_data_P' int2str(b) '.mat'], 'V_all', 'Na_all', 'K_all', 'L_all', 'SYN_all', '-v7.3')
end

save([filename '/Processed Data/calcium_data.mat'], 'calcium_all', 'p_w_all', 'T_f', '-v7.3')
save([filename '/Processed Data/freq_data.mat'], 'freq', 'T', 'LFP_all', 'TFA_all', 'POW_all', 'PHASE', 'POWER', '-v7.3')

if(par.B*par.trials>1); close(h); end

end

