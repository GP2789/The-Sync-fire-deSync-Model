function [  ] = content_specificity( filename )
% This will analysis similarity across and between patterns to reproduce
% Figure 5 of the paper.
%% LOAD FILES & INITIALISE
load([filename '/Data/parameters' ]);
load([filename '/Data/patterns' ]);
%load([filename '/pattern_data.mat']);
load([filename '/Processed Data/pattern_info.mat']);
conditions = {has_a_miss, double_bindings, all_else};
conditions_n = {'has_miss', 'double_binding', 'other_hits'};

par.sim_length = par.pre_stim_L * 2 + par.NC_stim_L;
if(exist([filename '/Processed Data/power_phase_data.mat'],'file')~=2)
    %% find ID of similar patterns for each pattern
    win = 50; % (lags of +/- 50 ms)
    same_trial = cell(size(patterns_t,1),1);
    patterns = patterns_t - patterns_t(:,1) + 1;
    for i = 1:size(patterns_t,1)
        pat_diff = patterns - patterns(i,:);
        same_trial{i} = find(all(abs(pat_diff) <= win,2));
        same_trial{i}(same_trial{i} == i) = [];
    end
    %% extract phase info by stimuli centred window
    des_pre_stim_sec = 3; pre_stim_L = des_pre_stim_sec*1000; win = 750;
    phase_DL = zeros(par.B, par.trials, par.pre_stim_L*2);
    phase_AL = zeros(par.B, par.trials, pre_stim_L*2);
    fprintf('importing phase & power data ... \n'); tic
    for b = 1:par.B % loop over patterns
        load([filename '/Processed Data/freq_data_P' int2str(b) '.mat']);
        if(b == 1)
            POW.DL = zeros(size(POW_all.DL.NC)); POW.AL = zeros(size(POW_all.AL.NC));
        end
        POW.DL = POW.DL + POW_all.DL.NC;
        POW.AL = POW.AL + POW_all.AL.NC;
        
        pre_stim_i = win;
        post_stim_i = par.pre_stim_L*2+par.NC_stim_L-win;
        for t = 1:par.trials
            % get pre-stim frequency
            dPh = PHASE.AL.NC(pre_stim_i:par.pre_stim_L-100);
            dPh = dPh(2:end) - dPh(1:end-1);
            % find how many ms per cycle
            pi_cycle = round((2*pi)/mean(dPh(dPh > 0)));
            % find how many cycles are needed for desired pre-stim period
            n_cycles = round((1000/pi_cycle) * (des_pre_stim_sec+1));
            % create elongated pre-stim period
            phase_BS = repmat(-pi:(2*pi)/pi_cycle:pi, [1 n_cycles]);
            % find phase of NC at join
            phase_at_zero = PHASE.AL.NC(pre_stim_i);
            x_i = find(phase_at_zero >= phase_BS,1,'last');
            phase_BS = phase_BS(1:x_i);
            phase_BS = phase_BS(length(phase_BS)-des_pre_stim_sec*1000+2+(par.pre_stim_L-win):end);
            % create elongated post-stim period
            phase_PS = repmat(-pi:(2*pi)/pi_cycle:pi, [1 n_cycles]);
            % find phase of NC at join
            phase_at_end = PHASE.AL.NC(post_stim_i);
            x_i = find(phase_at_end <= phase_PS,1,'first');
            phase_PS = phase_PS(x_i:end);
            phase_PS = phase_PS(1:des_pre_stim_sec*1000-(par.pre_stim_L-win));
            % conjoin elongated pre/post stim periods to trimmed data
            phase_NC = [phase_BS ... % pre stim
                PHASE.AL.NC(t, pre_stim_i:post_stim_i) ...% trimmed data
                phase_PS]; % post stim
            % recentre over time of first stimuli
            x1 = patterns_t(b,1) + 1;
            phase_AL(b, t, :) = phase_NC(t, x1 : x1 + pre_stim_L * 2 - 1);
            % just extract series aligned with first stimuli for DL
            phase_DL(b, t, :) = PHASE.DL.NC(t, x1 : x1 + par.pre_stim_L * 2 - 1);
        end
        fprintf('\b\b\b\b\b%3.0f%%\n', b/par.B*100)
    end
    POW.DL = POW.DL / par.B;
    POW.AL = POW.AL / par.B;
    save([filename '/Processed Data/power_phase_data.mat'],...
        'phase_AL','phase_DL', 'POW', 'freq', 'same_trial', 'patterns', 'pre_stim_L')
    fprintf('\b\b\b\b\bcompleted in %.0f minutes\n', toc/60)
else
    fprintf('loading phase & power data\n');
    load([filename '/Processed Data/power_phase_data.mat'])
end

%% PHASE ANALYSIS OVER/BETWEEN PATTERNS
fprintf('\ncalculating content similarity ...\n')
time = nan(par.B,1); fpr = false; n_trials = par.trials;
win = 400; % desired window around stimulus times
win = min(win, pre_stim_L*2 - (pre_stim_L+par.NC_stim_L));
patterns_i = 1:par.B;
for b = 1:par.B
    if(exist([filename '/Processed Data/similarity_data_P' int2str(b) '.mat'],'file')~=2)
        tic
        n_sim_patterns = length(same_trial{b});
        n_diff_patterns = par.B - n_sim_patterns - 1;
        sim_to_same = nan(n_sim_patterns, par.trials, pre_stim_L*2);
        sim_to_diff = nan(n_diff_patterns, par.trials, pre_stim_L*2);
        stim_L = par.pre_stim_L : min(par.pre_stim_L + patterns(b,end) + win, par.pre_stim_L*2);
        X_L = round(length(stim_L)/2):round((pre_stim_L*2)-length(stim_L)/2);
        % CALCULATE RSA SIMILARITY BETWEEN PATTERNS
        for t1 = 1:par.trials % loop through all trials
            %% TO SAME PATTERNS (TRIAL VARIATION)
            phase_DL_t = squeeze(phase_DL(b, t1, stim_L))';
            to_same_t = nan(n_sim_patterns, par.trials, pre_stim_L*2);
            parfor b2 = 1:n_sim_patterns % loop through all other trials
                for t2 = 1:n_trials
                    phase_AL_t = squeeze(phase_AL(same_trial{b}(b2), t2, :))';
                    to_same_t( b2, t2, X_L) = phase_compare( phase_AL_t, phase_DL_t );
                end
            end
            sim_to_same(:, t1, :) = squeeze(nanmean(to_same_t, 2)); clear to_same_t;
            
            %% TO DIFFERENT PATTERNS MEAN (PATTERN VARIATION)
            to_diff_t = nan(n_diff_patterns, par.trials, pre_stim_L*2);
            diff_trial = patterns_i(~ismember(patterns_i,same_trial{b}));
            diff_trial = diff_trial(diff_trial ~= b);
            parfor b2 = 1:n_diff_patterns %loop through all other patterns
                for t2 = 1:n_trials % loop through all trials
                    phase_AL_t = squeeze(phase_AL(diff_trial(b2), t2, :))';
                    to_diff_t( b2, t2, X_L ) = phase_compare( phase_AL_t, phase_DL_t );
                end
            end
            sim_to_diff(:, t1, :) = squeeze(nanmean(to_diff_t, 2)); clear to_diff_t;
        end
        if(n_sim_patterns ~= 0); sim_to_same = ...
                reshape(squeeze(nanmean(sim_to_same,2)),[n_sim_patterns pre_stim_L*2]);
        else; sim_to_same = nan(1, pre_stim_L*2);
        end
        if(n_diff_patterns ~= 0); sim_to_diff = ...
                reshape(squeeze(nanmean(sim_to_diff,2)),[n_diff_patterns pre_stim_L*2]);
        else; sim_to_diff = nan(1, pre_stim_L*2);
        end
        %% extract similarity based on conditions
        sim_to_same_cond = struct;
        sim_to_diff_cond = struct;
        content_spec_cond = struct;
        for c = 1:length(conditions_n)
            x_i = find(conditions{c});
            sim_to_same_cond.([conditions_n{c}]) = nanmean(sim_to_same(ismember(same_trial{b}, x_i), :), 1);
            not_sim_to_same = sim_to_same(~ismember(same_trial{b}, x_i), :);
            if(isnan(not_sim_to_same)); sim_to_diff_cond.([conditions_n{c}]) = nanmean(sim_to_diff,1);
            else; sim_to_diff_cond.([conditions_n{c}]) = nanmean([sim_to_diff; not_sim_to_same],1);
            end
            content_spec_cond.([conditions_n{c}]) = sim_to_same_cond.([conditions_n{c}]) - sim_to_diff_cond.([conditions_n{c}]);
        end
        sim_to_same = nanmean(sim_to_same, 1);
        sim_to_diff = nanmean(sim_to_diff, 1);
        content_spec = sim_to_same - sim_to_diff;
        % RECORD VARIABLES
        save([filename '/Processed Data/similarity_data_P' int2str(b) '.mat'], ...
            'sim_to_same', 'sim_to_diff', 'content_spec', ...
            'sim_to_same_cond', 'sim_to_diff_cond', 'content_spec_cond');
        
        B = dir([filename '/Processed Data']);
        B_rem = par.B - sum(contains({B(:).name}, 'similarity_data'));
        time(b) = toc;
        min_rem = (nanmean(time)/60) * B_rem;
        H = floor(min_rem/60);
        M = floor(min_rem - H*60);
        per_done = (par.B-B_rem)/par.B*100;
        if(~fpr); fprintf('\b\b\b\b%3.0f%% time remaining %2.0fh %2.0fm\n', per_done, H, M); fpr = true;
        else; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f%% time remaining %2.0fh %2.0fm\n', per_done, H, M);
        end
    end
end

if(exist([filename '/Processed Data/similarity_data_all.mat'],'file')~=3)
    fprintf('\ncollating data ... \n', b, par.B)
    for b = 1:par.B
        fprintf('\b\b\b\b\b%3.0f%%\n', b/par.B*100)
        load([filename '/Processed Data/similarity_data_P' int2str(b) '.mat'])
        if(b==1)
            content_spec_all = nan(par.B, size(sim_to_same,2));
            sim_to_same_all = nan(par.B, size(sim_to_same,2));
            sim_to_diff_all = nan(par.B, size(sim_to_same,2));
        end
        sim_to_same_all(b,:) = sim_to_same;
        sim_to_diff_all(b,:) = sim_to_diff;
        content_spec_all(b,:) = content_spec;
    end
    T = -size(content_spec_all,2)/2+1:size(content_spec_all,2)/2;
    T1 = find(T >= -1000,1,'first');
    T2 = find(T >= -250,1,'first');
    pre_s = mean(sim_to_same_all(:,T1:T2),2);
    sim_to_same_all = sim_to_same_all - pre_s;
    pre_s = mean(sim_to_diff_all(:,T1:T2),2);
    sim_to_diff_all = sim_to_diff_all - pre_s;

    save([filename '/Processed Data/similarity_data_all.mat'], ...
        'sim_to_same_all', 'content_spec_all', 'sim_to_diff_all', 'T');
end

end

