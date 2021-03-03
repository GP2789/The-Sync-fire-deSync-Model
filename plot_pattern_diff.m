function [  ] = plot_pattern_diff( filename )
%% load and initialise
fprintf('extracting pattern info ...\n');
extract_patterns( filename );
win = 500;
load([filename '/Data/parameters' ]);
par.sim_length = par.pre_stim_L * 2 + par.NC_stim_L;
load([filename '/Processed Data/pattern_info.mat' ]);
load([filename '/Data/patterns.mat' ]);
T = 0:1:50;
PSP = (exp(1)*T./5).*exp(-T./5).*heaviside(T);
if(exist([filename '/Processed Data/pattern_sim_data.mat'],'file')~=3)
    dW_BP_NC = cell(par.B, par.NC_n_stim); dW_BP_BP = cell(par.B, par.NC_n_stim);
    dW_TCL1_BP = cell(par.B, par.NC_n_stim); dW_TCL2_BP = cell(par.B, par.NC_n_stim);
    dBL_act = cell(par.B, par.NC_n_stim); dAL_act = cell(par.B, par.NC_n_stim);
    dNC_phase = zeros(par.B, par.NC_n_stim);
    time = nan(par.B,1); fpr = false;
    for b = 1:par.B
        tic
        pattern = patterns_t(b,:);
        %% extract DL weight change
        load([filename '/Processed Data/calcium_data_P' int2str(b) '.mat'])
        for p = 1:length(pattern)
            i1 = find(0:par.sim_length <= pattern(p) + par.pre_stim_L, 1, 'last');
            i2 = i1 + win;
            dW_BP_NC{b,p} = sum(p_w_all.BP_NC{p}(:,i1:i2),1); 
            dW_BP_BP{b,p} = sum(p_w_all.BP_BP{p}(:,i1:i2),1);
            dW_TCL1_BP{b,p} = sum(p_w_all.TCL1_BP{p}(:,i1:i2),1);
            dW_TCL2_BP{b,p} = sum(p_w_all.TCL2_BP{p}(:,i1:i2),1);
        end
        clear calcium_all p_w_all;
        
        %% extract BL/AL NC spike activity
        load([filename '/Processed Data/spike_data_P' int2str(b) '.mat'])
        BL_spks = SPIKES.DL.T1.NC_E;
        AL_spks = SPIKES.AL.T1.NC_E;
        for p = 1:length(pattern)
            n_ID = par.NC_stims_N_all{b}{1}{1}(p);
            BL_t = pattern(p) + par.pre_stim_L;
            AL_t = pattern(p) + par.pre_stim_L; % allow lag at recall?
            BL_spks_t = BL_spks(BL_spks(:,1) == n_ID,2);
            AL_spks_t = AL_spks(AL_spks(:,1) == n_ID,2);
            BL_spks_t = BL_spks_t(BL_spks_t > BL_t & BL_spks_t < BL_t + win)-BL_t;
            AL_spks_t = AL_spks_t(AL_spks_t > AL_t & AL_spks_t < AL_t + win)-AL_t;
            dBL_act{b,p} = zeros(1, win); dBL_act{b,p}(BL_spks_t) = 1;
            dAL_act{b,p} = zeros(1, win); dAL_act{b,p}(AL_spks_t) = 1;
            dBL_act{b,p} = filter(PSP,1, dBL_act{b,p});
            dAL_act{b,p} = filter(PSP,1, dAL_act{b,p});
        end
        clear SPIKES BL_spks BL_spks_t AL_spks AL_spks_t n_ID;
        
        %% extract NC phase at time of NC stim
        load([filename '/Processed Data/freq_data_P' int2str(b) '.mat'])
        dNC_phase(b,:) = PHASE.DL.NC(:,pattern+par.NC_pre_stim_L);% + par.NC_stim_tau);
        clear freq LFP_all POW_all PHASE POWER;
        time(b) = toc;
        min_rem = (nanmean(time)/60) * (par.B-b);
        H = floor(min_rem/60);
        M = floor(min_rem - H*60);
        per_done = b/par.B*100;
        if(~fpr); fprintf('\b\b\b\b%3.0f%% time remaining %2.0fh %2.0fm\n', per_done, H, M); fpr = true;
        else; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f%% time remaining %2.0fh %2.0fm\n', per_done, H, M);
        end
    end
    
    save([filename '/Processed Data/pattern_sim_data.mat'], ...
        'dW_BP_NC','dW_BP_BP','dW_TCL1_BP','dW_TCL2_BP',...
        'dBL_act','dAL_act','dNC_phase')
    
else;load([filename '/Processed Data/pattern_sim_data.mat'])
end

%% get rid of in gap stimuli
dW_BP_NC = dW_BP_NC(~in_gap,:);
dW_BP_BP = dW_BP_BP(~in_gap,:);
dW_TCL1_BP = dW_TCL1_BP(~in_gap,:);
dW_TCL2_BP = dW_TCL2_BP(~in_gap,:);
dNC_phase = dNC_phase(~in_gap,:);
n_S = sum(~in_gap);
blink_t = blink_t(~in_gap,:);

%% ATTENTIONAL BLINK AND LEARNING
gap = 100;
T = gap:gap:par.NC_stim_L;
is_hit = false(n_S, size(patterns_t,2));
for b = 1:n_S
    BP_NC_h = vertcat(dW_BP_NC{b,:});
    BP_BP_h = vertcat(dW_BP_BP{b,:});
    TCL1_BP_h = vertcat(dW_TCL1_BP{b,:});
    TCL2_BP_h = vertcat(dW_TCL2_BP{b,:});
    BP_NC_h = (max(BP_NC_h,[],2) - min(BP_NC_h,[],2)) > 0.3;
    BP_BP_h = (max(BP_BP_h,[],2) - min(BP_BP_h,[],2)) > 0.3;
    TCL1_BP_h = (max(TCL1_BP_h,[],2) - min(TCL1_BP_h,[],2)) > 0.3;
    TCL2_BP_h = (max(TCL2_BP_h,[],2) - min(TCL2_BP_h,[],2)) > 0.3;
    is_hit(b,:) = all([BP_NC_h BP_BP_h TCL1_BP_h TCL2_BP_h],2);
end

T1_H = zeros(length(T),1); T1_M = zeros(length(T),1);
T2_H = zeros(length(T),1); T2_M = zeros(length(T),1);
T3_H = zeros(length(T),1); T3_M = zeros(length(T),1);
for i=1:n_S
    T_i1 = find(blink_t(i,1)<=T,1,'first');
    if(is_hit(i,1))
        T1_H(T_i1) = T1_H(T_i1) + 1; % T1 hit
        if(is_hit(i,2))
            T2_H(T_i1) = T2_H(T_i1) + 1; % T2|T1 hit
            T_i2 = find(sum(blink_t(i,:))<=T,1,'first');
            if(is_hit(i,3)) % T2|T1 hit
                T3_H(T_i2) = T3_H(T_i2) + 1; % T3|T2|T1 hit
            else; T3_M(T_i2) = T3_M(T_i2) + 1; % T3|T2|T1 miss
            end
        else; T2_M(T_i1) = T2_M(T_i1) + 1; % T2|T1 miss
        end
    else; T1_M(T_i1) = T1_M(T_i1) + 1; % T1 miss
    end
end
T1_per = (T1_H ./ (T1_H + T1_M)) * 100;
T2_per = (T2_H ./ (T2_H + T2_M)) * 100;
T3_per = (T3_H ./ (T3_H + T3_M)) * 100;
fig = figure(); set(fig,'position',[0 0 1500 400]);
subplot(1,3,1)
plot(T, T1_per,'b', 'linewidth', 2); hold on
plot(T, T2_per, 'r', 'linewidth', 2);
plot(T, T3_per, 'r--', 'linewidth', 2);
xlim([gap 1200])
title('accuracy by lag'); set(gca,'fontsize',14)
xlabel('lag (ms)'); ylabel('accuracy (%)')
legend('T1','T2|T1','T3|T2|T1', 'location','southeast')
subplot(1,3,2)
plot(T,T1_H/sum(T1_H),'b-', 'linewidth', 2); hold on
plot(T,T2_H/sum(T2_H),'r-', 'linewidth', 2);
plot(T,T3_H/sum(T3_H),'r--', 'linewidth', 2);
xlim([gap 1200])
title('hit proportion by lag'); set(gca,'fontsize',14)
xlabel('lag (ms)'); ylabel('norm. freq.')
legend('T1','T2|T1','T3|T2|T1', 'location','northeast')
subplot(1,3,3)
plot(T,T1_M/sum(T1_M),'b-', 'linewidth', 2); hold on
plot(T,T2_M/sum(T2_M),'r-', 'linewidth', 2);
plot(T,T3_M/sum(T3_M),'r--', 'linewidth', 2);
xlim([gap 1200])
title('miss proportion by lag'); set(gca,'fontsize',14)
xlabel('lag (ms)'); ylabel('norm. freq.')
legend('T1','T2|T1','T3|T2|T1', 'location','northeast')

saveas(fig, [filename '\Figures\attentional_blink.jpg']); close(fig);

%% NC PHASE & BLINK PAR
NC_phase_H = dNC_phase(is_hit(:,1),1);
NC_phase_H = NC_phase_H + pi;
NC_phase_M = dNC_phase(~is_hit(:,1),1);
NC_phase_M = NC_phase_M + pi;
t = 0:0.01:(2*pi); x = cos(t);

fig = figure(); set(fig,'position',[0 0 1250 1000]);
subplot(2,2,1)
polarhistogram(NC_phase_H,4,'normalization','pdf');
title('hits'); set(gca,'fontsize',18)
subplot(2,2,2)
polarhistogram(NC_phase_M,4,'normalization','pdf');
title('misses'); set(gca,'fontsize',18)

subplot(2,2,3)
yyaxis left; histogram(rad2deg(NC_phase_H), 12,'normalization','pdf'); ylabel('norm. freq.');
yyaxis right; plot(rad2deg(t), x, 'color', [1 0.45 0], 'linewidth', 2); ylabel('amplitude');
xlim([0 360]); xlabel('phase');
title(['hits n=' int2str(length(NC_phase_H))]); set(gca,'fontsize',18)
subplot(2,2,4)
yyaxis left; histogram(rad2deg(NC_phase_M), 12,'normalization','pdf'); ylabel('norm. freq.');
yyaxis right; plot(rad2deg(t), x, 'color', [1 0.45 0], 'linewidth', 2); ylabel('amplitude');
xlim([0 360]); xlabel('phase');
title(['misses n=' int2str(length(NC_phase_M))]); set(gca,'fontsize',18)

saveas(fig, [filename '\Figures\NC_alpha_misses.jpg']); close(fig);

end

