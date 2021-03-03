function [ ] = extract_pattern_info( filename )
%% initialise data
fprintf('extracting pattern info ...\n');
load([filename '\Data\patterns.mat'])
load([filename '\Data\parameters.mat'])
win = 500;
time = nan(par.B,1); fpr = false;
BP_replay = nan(par.B, par.NC_n_stim);
dW_BP = zeros(par.B, par.NC_n_stim);
BP_mp = zeros(par.B, par.NC_n_stim);
dW_NC = zeros(par.B, par.NC_n_stim);
dW_TC1 = zeros(par.B, par.NC_n_stim, par.n_TC_P(1));
dW_TC2 = zeros(par.B, par.NC_n_stim, par.n_TC_P(2));
dNS_DL = zeros(par.B, par.NC_n_stim, win);
dPC_DL = zeros(par.B, par.NC_n_stim, win);
dNS_AL = zeros(par.B, par.NC_n_stim, win);
dPC_AL = zeros(par.B, par.NC_n_stim, win);
W_cap = zeros(par.B, 1);
non_stationarity = zeros(par.B, 2);
BS_power_change = zeros(par.B, 2);
phase_at_stim = zeros(par.B, par.NC_n_stim);
pi_cycle = zeros(par.B, 1);
cont_spec = zeros(par.B, par.pre_stim_L);

%% load/extract pattern info
if(exist([filename '/pattern_data.mat'], 'file') ~=3)
    for b = 1:par.B
        tic; pattern = patterns_t(b,:);
        %% extract information about each stimulus
        load([filename '/Processed Data/calcium_data_P' int2str(b) '.mat'])
        load([filename '/Processed Data/freq_data_P' int2str(b) '.mat'])
        load([filename '/Processed Data/spike_data_P' int2str(b) '.mat'])
        %load([filename '/Processed Data/similarity_data_P' int2str(b) '.mat'])
        
        %% phase at stimulus presentation
        dPh = PHASE.DL.NC(750:par.NC_pre_stim_L-250);
        dPh = dPh(2:end) - dPh(1:end-1);
        pi_cycle(b) = round((2*pi)/mean(dPh(dPh > 0)));
        
        %% BP network capacity
        W_cap(b) = min(p_w_all.BP_capacity)*100;
        
        %% NC non-stationarity after learning
        for n = 1:2
            phase_NC = PHASE_NS.(par.sim_order_n{n}).NC(win+1:end-win);
            [~,i,~,p] = findpeaks(-phase_NC); i = i(p>=0.3); i_m = round(mean(i));
            if(isnan(i_m)); i_m = round((max(pattern) - min(pattern))/2)+par.pre_stim_L; end
            i_d = round((max(i)-min(i))/2); if(isempty(i_d)); i_d = 0; end; clear p i;
            non_stationarity(b, n) = mean(phase_NC(max(1,i_m-win-i_d+1):min(length(phase_NC),i_m+win+i_d)));

            %% NC power change from baseline after learning
            power_NC = POWER.(par.sim_order_n{n}).NC(win+1:end-win);
            [~,i,~,p] = findpeaks(-power_NC); i = i(p>=30); i_m = round(mean(i));
            if(isnan(i_m)); i_m = round((max(pattern) - min(pattern))/2)+par.pre_stim_L; end
            i_d = round((max(i)-min(i))/2); if(isempty(i_d)); i_d = 0; end; clear p i;
            BS_power_change(b, n) = mean(power_NC(max(i_m-win-i_d+1,1):min(length(phase_NC),i_m+win+i_d)));
        end
        
        %% weight change for each pattern
        for p = 1:length(pattern)
            i1 = pattern(p) + par.pre_stim_L + 1; 
            i2 = i1 + win - 1;
            
            % extract BP mid-point weight change
            BP_W_t = mean(p_w_all.BP_BP{p}(:,i1:i2),1);
            x1 = find(BP_W_t <= 0.95, 1, 'first') + pattern(p);
            if(isempty(x1)); x1 = NaN; end
            BP_mp(b, p) = x1;
            dW_BP(b, p) = max(BP_W_t) - min(BP_W_t);
            
            dW_NC(b, p) = max(p_w_all.BP_NC{p}(:,i1:i2));
            dW_TC1(b, p, :) = cellfun(@(x)max(x(:,i1:i2)),p_w_all.TCL1_BP(p,:));
            dW_TC2(b, p, :) = cellfun(@(x)max(x(:,i1:i2)),p_w_all.TCL2_BP(p,:));
                       
            dNS_DL(b, p, :) = mean(PHASE_NS.DL.NC(i1:i2),1);
            dPC_DL(b, p, :) = mean(POWER.DL.NC(i1:i2),1);
            dNS_AL(b, p, :) = mean(PHASE_NS.AL.NC(i1:i2),1);
            dPC_AL(b, p, :) = mean(POWER.AL.NC(i1:i2),1);
           
            % find spikes DL - AL difference
            spks = SPIKES.AL.T1.BP_E;
            spks = spks(ismember(spks(:,1), BP_G_ID{p}), :);
            spks = spks(spks(:,2) > i1 & spks(:,2) <= i1 + 250, 1);
            spks = unique(spks);
            BP_replay(b, p) = numel(spks) / numel(BP_G_ID{p}); 
        end
        
        %% timers
        clear calcium_all p_w_all PHASE POWER PHASE_NS POW_all LFP_all freq;
        time(b) = toc; min_rem = (nanmean(time)/60) * (par.B-b);
        H = floor(min_rem/60); M = floor(min_rem - H*60); per_done = b/par.B*100;
        if(~fpr); fprintf('\b\b\b\b%3.0f%% time remaining %2.0fh %2.0fm\n', per_done, H, M); fpr = true;
        else; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f%% time remaining %2.0fh %2.0fm\n', per_done, H, M);
        end
    end
    save([filename '/pattern_data.mat'], 'dW_BP', 'dW_NC', 'dW_TC1', 'dW_TC2', 'W_cap', ...
        'dPC_DL', 'dNS_DL', 'dPC_AL', 'dNS_AL', 'BP_replay',...
        'non_stationarity', 'BS_power_change', 'pi_cycle');
else; load([filename '/pattern_data.mat'])
end

%% find ID of similar patterns for each pattern
win = 50; % (lags of +/- 50 ms)
same_trial = cell(size(patterns_t,1),1);
patterns = patterns_t - patterns_t(:,1) + 1;
for i = 1:size(patterns_t,1)
    pat_diff = patterns - patterns(i,:);
    same_trial{i} = find(all(abs(pat_diff) <= win,2));
    same_trial{i}(same_trial{i} == i) = [];
end
n_trials = cellfun(@length,same_trial);

%% find chunking
blink_t = [];
for i = 1:par.NC_n_stim - 1
    blink_t(:,i) = patterns_t(:,i+1)-patterns_t(:,i);
end
blink_t = [inf(par.B,1) blink_t];
has_chunking = blink_t < 110;

%% classify hits
active_units = ~isnan(dW_BP) & ~isnan(dW_NC) & any(~isnan(dW_TC1),3) & any(~isnan(dW_TC2),3);
syn_misses = ~isnan(dW_BP) & all(isnan(dW_TC2),3);% | all(isnan(dW_TC2),3));
BP_misses = isnan(dW_BP);
is_hit = active_units;
% NaNs recorded at preprocessing indicate no active
% binding pool units were found for any given stimulus.

%% find which stimuli are in synfire chain 'gaps' of fast & slow hierarchies
load([filename '/Processed Data/spike_data_P1.mat']);
TC_SPKS = SPIKES.DL.T1.TC_E_L2; TC_LFP = zeros(1, par.sim_length);
for t = 1:par.sim_length; TC_LFP(t) = numel(TC_SPKS(TC_SPKS(:,2) == t, 2)); end

TC_LFP_t = filter(compute_psp(25,1:100),1,TC_LFP);
TC_LFP_t = filter(compute_psp(25,1:100),1,TC_LFP_t);
%TC_LFP = filter(compute_psp(15,1:100),1,TC_LFP);

[~,i,~,~] = findpeaks(-TC_LFP_t); clear TC_LFP_t;
TC_F = 1000 / mean(diff(i(1:7)));
TC_F = round(TC_F*[0.75 1.25]);
[TC_LFP] = create_LFP(TC_SPKS, 1, 1, TC_F, par.sim_length, 1);
TFA.UBF = 1; TFA.OBF = 30; TFA.SR = 1000; TFA.Gamma = 0.5; 
[P, ~] = GaborFilter(TC_LFP, TFA);
freq = TFA.UBF : (TFA.OBF-TFA.UBF)/(size(P,1)-1) : TFA.OBF;
F = find(freq >= TC_F(1),1,'first') : find(freq <= TC_F(2),1,'last');
TC_phase = angle(mean(P(F,:)));


%% find repetitions
is_repetition = false(par.B,par.NC_n_stim);
for b = 1:par.B
    is_rep = abs(diff(par.NC_stims_N_all{b}{1})) < ceil(par.NC_stim_N_cos);
    for p = 1:par.NC_n_stim-1
        if(is_rep(p)); is_repetition(b,p+1) = true; end
    end
end

%% calculate attentional blink
lag = 50;
AB_lags = lag:lag:par.NC_stim_L;
AB_lag_t = nan(par.B, par.NC_n_stim-1);
AB_hits = zeros(par.NC_n_stim, length(AB_lags));
AB_miss = zeros(par.NC_n_stim, length(AB_lags));
in_blink = false(par.B, par.NC_n_stim);
for i=1:par.B
    T_i1 = find(blink_t(i,2)<=AB_lags,1,'first');
    AB_lag_t(i,1) = AB_lags(T_i1);
    %if(~has_chunking(i,1) && ~is_repetition(i,1))
    if(is_hit(i,1))
        AB_hits(1,T_i1) = AB_hits(1,T_i1) + 1; % T1 hit
        if(is_hit(i,2))
            AB_hits(2,T_i1) = AB_hits(2,T_i1) + 1; % T2|T1 hit
            if(par.NC_n_stim > 2 && ~has_chunking(i,2))
                T_i2 = find(sum(blink_t(i,3))<=AB_lags,1,'first');
                AB_lag_t(i,2) = AB_lags(T_i2);
                if(is_hit(i,3)) % T2|T1 hit
                    AB_hits(3,T_i2) = AB_hits(3,T_i2) + 1; % T3|T2|T1 hit
                else; AB_miss(3,T_i2) = AB_miss(3,T_i2) + 1; in_blink(i,3) = true;  % T3|T2|T1 miss
                end
            end
        else; AB_miss(2,T_i1) = AB_miss(2,T_i1) + 1; in_blink(i,2) = true; % T2|T1 miss
        end
    else; AB_miss(1,T_i1) = AB_miss(1,T_i1) + 1; in_blink(i,1) = true;  % T1 miss
    end
    % end
end
AB_accuracy = AB_hits ./ (AB_hits+AB_miss) * 100;
%if(exist([filename '\Figures'],'dir')~=7); mkdir([filename '\Figures']); end
%saveas(fig,[filename '\Figures\pattern_info.jpg']); close(fig);

misses = sum(BP_misses | syn_misses,2);
has_a_miss = misses > 0;
TC_n = sum(~isnan(dW_TC2),3);
double_bindings = sum(TC_n > 1,2) > 0 & ~has_a_miss;
all_else = ~double_bindings & ~has_a_miss;
%%
% subplot(2,3,3); 
% NC_H = histc(phase_at_stim(~BP_misses(:,1) & ~in_blink(:,1),1), -pi:pi/4:pi);% -pi:pi/4:pi);
% NC_M = histc(phase_at_stim(BP_misses & ~in_blink), -pi:pi/4:pi);
% NC_M_p = NC_M ./ (NC_H + NC_M);
% NC_H_p = NC_H ./ (NC_H + NC_M); 
% bar(-pi:pi/4:pi, NC_H_p, 'facecolor', 'none'); hold on
% bar(-pi:pi/4:pi, NC_M_p, 'facecolor', 'r'); alpha(0.5); hold off;
%histogram(phase_at_stim(~BP_misses(:,1),1),-pi:pi/4:pi,'facecolor','none');hold on;
%histogram(phase_at_stim(BP_misses(:,1),1),-pi:pi/4:pi,'facecolor','r');hold off


%% NC ALPHA HIT/MISS
% NC_phase_H = phase_at_stim(is_hit & ~in_blink);
% NC_phase_M = phase_at_stim(~is_hit & ~in_blink);
% % t = -pi:pi/16:pi; x = cos(t);
% % 
% % fig = figure(); set(fig,'position',[0 0 1250 1000]);
% % subplot(2,2,1)
% polarhistogram(NC_phase_H,4,'normalization','pdf');
% title('hits'); set(gca,'fontsize',18)
% subplot(2,2,2)
% polarhistogram(NC_phase_M,4,'normalization','pdf');
% title('misses'); set(gca,'fontsize',18)
% 
% subplot(2,2,3)
% yyaxis left; histogram(NC_phase_H, [-pi:pi/4:pi],'normalization','pdf'); ylabel('norm. freq.');
% yyaxis right; plot(t, x, 'k', 'linewidth', 2); 
% xlim([-pi pi]); xticks([-pi 0 pi]); xticklabels({'-\pi', '0' '\pi'}); xlabel('phase');
% title(['hits n=' int2str(length(NC_phase_H))]); set(gca,'fontsize',18)
% ax = gca; ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = 'k';
% 
% subplot(2,2,4)
% yyaxis left; histogram(NC_phase_M, [-pi:pi/4:pi],'normalization','pdf'); 
% yyaxis right; plot(t, x, 'k', 'linewidth', 2); ylabel('amplitude');
% xlim([-pi pi]); xticks([-pi 0 pi]); xticklabels({'-\pi', '0' '\pi'}); xlabel('phase');
% title(['misses n=' int2str(length(NC_phase_M))]); set(gca,'fontsize',18)
% ax = gca; ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = 'k';
% 
% saveas(fig, [filename '\Figures\NC_alpha_misses.jpg']); close(fig);
%%

save([filename '\Processed Data\pattern_info.mat'],...
    'same_trial','n_trials','patterns','has_chunking','blink_t','is_hit','is_repetition',...
    'AB_accuracy', 'AB_hits', 'AB_miss', 'AB_lags', 'AB_lag_t', 'in_blink', ...
    'syn_misses', 'BP_misses', 'BP_mp', 'TC_SPKS', 'TC_LFP', 'TC_phase', ...
    'has_a_miss', 'all_else', 'double_bindings')


end

