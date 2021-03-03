function [ ] = plot_content_spec( filename )
%% load data
load([filename '\Data\parameters.mat'])
load([filename '\Data\patterns.mat'])
load([filename '\Processed Data\pattern_info.mat'])
load([filename '/pattern_data.mat'])
if(exist([filename '/Figures'],'dir')~=7); mkdir([filename '/Figures']); end

 
%% 7A stim times distribution

hold on; box on; set(gca,'fontsize',12);
S1 = histcounts(patterns_t(:,1),0:100:par.NC_stim_L);
S2 = histcounts(patterns_t(:,2),0:100:par.NC_stim_L);
S3 = histcounts(patterns_t(:,3),0:100:par.NC_stim_L);
S_sum = S1 + S2 + S3; S1 = S1 ./ (S_sum);  S2 = S2 ./ (S_sum);  S3 = S3 ./ (S_sum);
b = bar([S1;S2;S3]','stacked','facecolor','flat', 'edgecolor','none','barwidth',1);
b(1).XData = 50:100:par.NC_stim_L;
b(2).XData = 50:100:par.NC_stim_L;
b(3).XData = 50:100:par.NC_stim_L;
legend('T1','T2','T3','location','north','orientation','horizontal')
ylabel('probability'); xlabel('time from onset (ms)')
title('presentation sequence');
xlim([0 par.NC_stim_L]); xticks([0:500:par.NC_stim_L])
ax = gca; ax.YAxis(1).Color = 'b';
saveas(fig, [filename '/Figures/7A_presentation-sequence']); close(fig);

%% 7B binding misses by lag
fig = figure(); hold on; box on; set(gca,'fontsize',12);
title('binding pool misses by lag');
yyaxis left
blink_h = histc(blink_t(:), 0:100:par.NC_stim_L) / par.B / (par.NC_n_stim-1);
b = bar( 50:100:par.NC_stim_L+50,blink_h, 'facecolor','flat','edgecolor','none','barwidth',1);
%b(1).XData = 50:100:par.NC_stim_L;
xlim([0 par.NC_stim_L]); ylim([0 0.2]); yticks([0 max(ylim)])
AB_t = 0:100:par.NC_stim_L; AB_m = nan(1, length(AB_t)-1); ylabel('lag probability')
blink_t(:,1) = blink_t(:,2);
for i = 2:length(AB_t)
    blink_i = blink_t < AB_t(i) & blink_t > AB_t(i-1);
    AB_m(i-1) = sum(BP_misses(blink_i));
end
AB_m = AB_m / par.B / par.NC_n_stim;
yyaxis right; plot(AB_t(2:end)-50, AB_m, 'r-s')
xlabel('lag between targets (ms)'); ylabel('miss probability')
ylim([0 0.05]); yticks([0 max(ylim)]);
ax = gca; ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = 'r';
saveas(fig, [filename '/Figures/7B_binding-pool-misses-by-lag']); close(fig);

%% 7C misses per pattern
fig = figure();  hold on; box on; set(gca,'fontsize',12);
misses = sum(BP_misses | syn_misses,2);
n_misses = histcounts(misses,0:4) / par.B;
b = bar(0:3, n_misses,'facecolor',[0 0 0.4],'facecolor','w','edgecolor','r','barwidth',1);
xlim([-0.5 3.5])
title('total misses per sequence')
ylim([0 0.7]); yticks([0 max(ylim)])
xlabel('no. misses'); ylabel('miss probability')
ax = gca; ax.YAxis(1).Color = 'r';
saveas(fig, [filename '/Figures/7C_misses-per-sequence']); close(fig);

%% 7D synfire chain LFP & weights
T1 = -125; T2 = par.NC_stim_L;
T = T1:25:T2; 
TC_W_t = nan(length(T)-1, 1); TC_n_t = nan(length(T)-1, 1); 
TC_SPKS_t = nan(length(T)-1, 1); 
BP_rep_t = nan(length(T)-1, 1);
TC_W = nanmean(dW_TC2,3);
TC_W(isnan(TC_W)) = 0;
TC_n = sum(~isnan(dW_TC2),3);
for t = 2:length(T)
    TC_W_t(t-1) = nanmean(TC_W(BP_mp <= T(t) & BP_mp >= T(t-1)));
    TC_n_t(t-1) = mean(TC_n(BP_mp <= T(t) & BP_mp >= T(t-1)));
    BP_rep_t(t-1) = mean(BP_replay(BP_mp <= T(t) & BP_mp >= T(t-1)));
    %TC_SPKS_t(t-1) = sum(TC_SPKS(:,2) <= T(t) + par.pre_stim_L & TC_SPKS(:,2) >= T(t-1) + par.pre_stim_L);
end
fig = figure();hold on; box on; set(gca,'fontsize',12);
title('synfire chain modulation')
yyaxis left
b = bar(T(2:end)-25/2,TC_W_t,'facecolor','flat', 'edgecolor','none');
b(1).FaceColor = [0.8 0.8 0.8]; 
% b = bar(TC_SPKS_t,'facecolor','flat', 'edgecolor','none');
% b(1).FaceColor = [0.8 0.8 0.8]; b(1).XData = T;
ylim([0.5 1]); ylabel('mean SC \rightarrow BP synapse (\rho)'); 
yticks(min(ylim):0.25:max(ylim)); plot([0 0], ylim, 'k-');
yyaxis right; plot((1:size(TC_LFP,2))-par.pre_stim_L, ...
    TC_LFP,'color',[0.2 0.6 0.3],'linewidth',2);
ylabel('SC LFP');
xlim([T1 T2]); xticks(0:500:T2); xlabel('time (ms)');
ax = gca; ax.YAxis(1).Color = [0.5 0.5 0.5]; ax.YAxis(2).Color = [0.2 0.6 0.3];
saveas(fig, [filename '/Figures/7D_synfire-chain-modulation']); close(fig);
%% 7E
fig = figure(); hold on; box on; set(gca,'fontsize',12);
yyaxis left
b = bar(T(2:end)-25/2, TC_n_t,'facecolor','flat', 'edgecolor','none');
b(1).FaceColor = 'r';
ylabel('n SC bindings'); plot([0 0], ylim, 'k-');
yyaxis right
b = bar(T(2:end)-25/2, BP_rep_t,'facecolor','flat', 'edgecolor','none');
b(1).FaceColor = [0.3 0.3 0.6]; 
ylabel('BP replay prob.'); xlim([T1 T2]); xticks(0:500:T2);xlabel('time (ms)');
ax = gca; ax.YAxis(1).Color = 'r'; ax.YAxis(2).Color = [0.3 0.3 0.6];
saveas(fig, [filename '/Figures/7E_synfire-chain-modulation_2']); close(fig);

%% 7F synfire chain LFP phase hits v misses
fig = figure(); hold on; box on; set(gca,'fontsize',12);
subplot(1,2,2);
TC_phase_H = TC_phase(patterns_t(BP_replay > 0.5) + par.pre_stim_L);
polarhistogram(TC_phase_H,8,'normalization','pdf','DisplayStyle','stair','edgecolor',[0.3 0.3 0.6]);
title('BP replay prob \geq 0.5'); set(gca,'fontsize',12);
ax = gca; ax.ThetaAxisUnits = 'radians'; ax.ThetaZeroLocation = 'left';
thetaticks([0:pi/2:2*pi])

subplot(1,2,1);
TC_phase_M = TC_phase(patterns_t(TC_n > 1) + par.pre_stim_L);
polarhistogram(TC_phase_M,8,'normalization','pdf','DisplayStyle','stair','edgecolor','r');
title('SC bindings \geq 2'); set(gca,'fontsize',12);
ax = gca; ax.ThetaAxisUnits = 'radians'; ax.ThetaZeroLocation = 'left';
thetaticks([0:pi/2:2*pi])
saveas(fig, [filename '/Figures/7F-G_synfire-chain-modulation_3']); close(fig);

%% 7i-iii calculate CS for specific categories
load([filename '/Processed Data/power_phase_data.mat'])
patterns_i = 1:par.B;

sp_i = [3 6 9]; h = [];
titles_i1 = {'\geq1 miss', 'binds >1', 'other'};
titles_i2 = {'has_miss', 'double_binding', 'other_hits'};

has_a_miss = misses > 0;
double_bindings = sum(TC_n > 1,2) > 0 & ~has_a_miss;
all_else = ~double_bindings & ~has_a_miss;

sim_to_same_all = cell(1,length(titles_i1));
sim_to_diff_all = cell(1,length(titles_i1));
content_spec_all = cell(1,length(titles_i1));

for i = 1:length(titles_i1)
   fprintf(['extracting data for "' regexprep(titles_i1{i},'_',' ') '" ...\n'])
    
   if(i==1); x_i = find(has_a_miss);
   elseif(i==2); x_i = find(double_bindings);
   elseif(i==3); x_i = find(all_else);
   end
   time = nan(numel(x_i),1); fpr = false;
   sim_to_same_all{i} = nan(numel(x_i), pre_stim_L*2);
   sim_to_diff_all{i} = nan(numel(x_i), pre_stim_L*2);
   content_spec_all{i} = nan(numel(x_i), pre_stim_L*2);
   
   for b = 1:numel(x_i)
       tic
       load([filename '/Processed Data/similarity_data_P' int2str(x_i(b)) '.mat'])
       sim_to_same_all{i}(b,:) = sim_to_same_cond.(titles_i2{i});
       sim_to_diff_all{i}(b,:) = sim_to_diff_cond.(titles_i2{i});
       content_spec_all{i}(b,:) = content_spec_cond.(titles_i2{i});
        time(b) = toc;
        sec_rem = (nanmean(time)) * (numel(x_i)-b);
        M = floor(sec_rem/60);
        S = floor(sec_rem - M*60);
        per_done = b/numel(x_i)*100;
        if(~fpr); fprintf('\b\b\b\b%3.0f%% time remaining %2.0fm %2.0fs\n', per_done, M, S); fpr = true;
        else; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f%% time remaining %2.0fm %2.0fs\n', per_done, M, S);
        end
   end
end
fig = figure(); set(fig,'position',[561    59   573   890])
lines_i = {':','--','-'}; T = -pre_stim_L+1:pre_stim_L;
titles_i2 = {'similar sequences', 'different sequences', 'content specificity'};
for i = 1:3
    subplot(3,1,i); hold on; box on; set(gca,'fontsize',12,'YAxisLocation','right');
    xlim([-250 2000]); xticks(-500:500:2000); plot(xlim,[0 0],'k');
    ylim([-0.03 0.08]); yticks([min(ylim) 0 max(ylim)]); plot([0 0], ylim, 'k');
    title(regexprep(titles_i2{i},'_',' '));
end

h1 = []; h2 = []; h3 = [];

for i = 1:length(titles_i1)
    win = 500;
    T1 = find(T >= -1000,1,'first');
    T2 = find(T >= -250,1,'first');
    SS = sim_to_same_all{i};
    SD = sim_to_diff_all{i};
    CS = content_spec_all{i};
    
    pre_s = mean(SS(:,T1:T2),2);
    SS = nanmean(SS - pre_s,1);
    pre_s = mean(SD(:,T1:T2),2);
    SD = nanmean(SD - pre_s,1);
    CS = filter(ones(win,1),win,nanmean(CS,1));
    SS = filter(ones(win,1),win,SS);
    SD = filter(ones(win,1),win,SD);
    
    subplot(3,1,1); h1(i) = plot(T, SS, 'color', 'r','linestyle',lines_i{i},'linewidth',2);
    
    subplot(3,1,2); h2(i) = plot(T, SD, 'color', 'k','linestyle',lines_i{i},'linewidth',2);
        
    subplot(3,1,3); h3(i) = plot(T, CS, 'color', [0.8 .4 .4],'linestyle',lines_i{i},'linewidth',2);
end
subplot(3,1,1); legend(h1,regexprep(titles_i1,'_',' '),...
    'location','northwest','orientation','vertical','fontsize',8)
subplot(3,1,2); legend(h2,regexprep(titles_i1,'_',' '),...
    'location','northwest','orientation','vertical','fontsize',8); ylabel('similarity');
subplot(3,1,3); legend(h3,regexprep(titles_i1,'_',' '),...
    'location','northwest','orientation','vertical','fontsize',8)

saveas(fig,[filename '/Figures/7i-iii_content_specificity.jpg']); close(fig);


end

