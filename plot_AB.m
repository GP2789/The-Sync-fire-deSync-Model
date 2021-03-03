function [  ] = plot_AB( filename, filenames_support )
%% load data
for i = 1:length(filenames_support)
    load([filenames_support{i} '/Data/parameters.mat'])
    load([filenames_support{i} '\Processed Data\pattern_info.mat'])
    load([filenames_support{i} '/pattern_data.mat'])
    TS_s{i+1} = int2str(par.BP_OE_TS);
    for t = 1:length(AB_lags)
        t_i = AB_lag_t(:,1) == AB_lags(t); 
        dW_t_s{i}(:,t) = nansum(dW_BP(t_i,:),1) / sum(t_i);
    end
end
load([filename '\Processed Data\pattern_info.mat'])
load([filename '/pattern_data.mat'])
load([filename '/Data/parameters.mat'])
load([filename '/Data/patterns.mat'])
TS_s{1} = int2str(par.BP_OE_TS);
load([filename '/Processed Data/similarity_data_all.mat']);
%% extract data by lag
AB_lags_2 = [150 300 600 1200]; AB_lags_2_t = [50 AB_lags_2];
X1 = -250; X2 = 1500; T1 = find(T>=X1,1,'first'); T2 = find(T>=X2,1,'first');
for t = 1:length(AB_lags)
    t_i = AB_lag_t(:,1) == AB_lags(t);
    dW_t(:,t) = nansum(dW_BP(t_i,:),1) / sum(t_i);
    W_cap_t(:,t) = mean(W_cap(t_i));
    
    [PC_CI_L_t(:,t), PC_CI_U_t(:,t), PC_AL_t(:,t)] = ...
        bootstrap( BS_power_change(t_i,2), 10000, []);
    [NS_CI_L_t(:,t), NS_CI_U_t(:,t), NS_AL_t(:,t)] = ...
        bootstrap( non_stationarity(t_i,2), 10000, []);
    
    t_i = AB_lag_t(:,1) == AB_lags(t) & is_repetition(:,2);
    if(sum(t_i)>=5); dW_R(:,t) = nansum(dW_BP(t_i,:),1) / sum(t_i);
    else; dW_R(:,t) = NaN;
    end
    W_cap_R(:,t) = mean(W_cap(t_i));
end
for t = 1:length(AB_lags)
    for i = 1:2
        if(isnan(dW_R(i,t)))
            if(t==1); dW_R(i,t) = dW_R(i,find(~isnan(dW_R(i,t+1:end)),1,'first')+1);
            elseif(t==length(AB_lags))
                dW_R(i,t) = dW_R(i,end-1);
            else
                dW_R(i,t) = mean(dW_R(i,find(~isnan(dW_R(i,t+1:end)),1,'first')+1 : ...
                    find(~isnan(dW_R(i,1:t)),1,'last')));
            end
        end
    end
    if(isnan(W_cap_R(t)))
        if(t==1); W_cap_R(t) = W_cap_R(find(~isnan(W_cap_R(t+1:end)),1,'first')+1);
        elseif(t==length(AB_lags))
            W_cap_R(t) = W_cap_R(end-1);
        else
            W_cap_R(t) = mean(W_cap_R(find(~isnan(W_cap_R(t+1:end)),1,'first')+1 : ...
                find(~isnan(W_cap_R(1:t)),1,'last')));
        end
    end
end
CS_t = zeros(length(AB_lags_2), size(sim_to_same_all,2));
for t = 1:length(AB_lags_2)
    t_i = blink_t(:,2) <= AB_lags_2(t);
    if(t>1)
        t_i = t_i & AB_lag_t(:,1) > AB_lags_2(t-1);
    end
    CS_t(t,:) = nanmean(content_spec_all(t_i,:),1);
    to_same_t(t,:) = nanmean(sim_to_same_all(t_i,:),1);
    to_diff_t(t,:) = nanmean(sim_to_diff_all(t_i,:),1);
end
%% plot BINDING AT ENCODING FOR T1
font_size = 15;
fig = figure(); set(fig, 'position',[0 0 1400 800])
subplot(3,2,1); hold on; box on; set(gca,'fontsize',font_size); 
title('binding at encoding for T1'); h1 = [];
yyaxis right; 
h1(1) = plot(AB_lags, dW_t(1,:),'k-','linewidth',2);
h1(2) = plot(AB_lags, dW_R(1,:),'k--','linewidth',2);
ylabel('weight change'); ylim([0 1]); yticks(min(ylim):0.5:1);
yyaxis left;
h1(3) = fill([AB_lags fliplr(AB_lags)], [W_cap_t zeros(size(W_cap_t))],[0.7 0.7 0.7],'edgecolor','none');  
ylabel('BP capacity (%)'); ylim([60 100]); yticks(60:20:100);
xticklabels({}); xlim([0 1250]); xticks(100:100:1200);
for t = 1:length(AB_lags_2_t); plot([1 1]*AB_lags_2_t(t), ylim,'k-'); end
legend(h1, {'T1 (+T2)', 'T1 (+T1)', 'BP cap.'}, ...
    'fontsize', 12, 'location', 'south','Orientation','horizontal')
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
%% plot BINDING AT ENCODING FOR T2
subplot(3,2,3); hold on; box on; set(gca,'fontsize',font_size); xlim([100 1200]);
title('binding at encoding for T2'); lines = {'k--','k:'}; h2 = [];
yyaxis right; 
h2(1) = plot(AB_lags, dW_t(2,:),'k-','linewidth',2);
for i = 1:length(filenames_support)
    h2(end+1) = plot(AB_lags, dW_t_s{i}(2,:), lines{i},'linewidth',2);
end
ylabel('weight change'); ylim([0 1]); yticks(min(ylim):0.5:1);
yyaxis left;
h2(end+1) = fill([AB_lags fliplr(AB_lags)], [W_cap_t zeros(size(W_cap_t))],[0.7 0.7 0.7],'edgecolor','none'); 
ylabel('BP capacity (%)'); ylim([60 100]); yticks(60:20:100);
xticklabels({}); xlim([0 1250]); xticks(100:100:1200);
for t = 1:length(AB_lags_2_t); plot([1 1]*AB_lags_2_t(t), ylim,'k-'); end
leg_lab = [cellfun(@(x,y)['\tau = ' x],TS_s,'uniformoutput',false), 'BP cap.'];
legend(h2([2 1 3 4]), leg_lab([2 1 3 4]), 'fontsize', 12, 'location', 'south','Orientation','horizontal')
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
%% plot STATIONARITY & POWER CHANGES
subplot(3,2,5); hold on; box on; set(gca,'fontsize',font_size); 
title('alpha power & stationarity at recall')
yyaxis left; 
fill([AB_lags fliplr(AB_lags)],[PC_CI_L_t fliplr(PC_CI_U_t)],'b','edgecolor','none'); alpha(0.25);
plot(AB_lags, PC_AL_t,'b-','linewidth',2); 
ylabel('power change (%)'); yticks(ceil(min(ylim):10:max(ylim)));
yyaxis right;
fill([AB_lags fliplr(AB_lags)],[NS_CI_L_t fliplr(NS_CI_U_t)],'r','edgecolor','none'); alpha(0.25);
plot(AB_lags, NS_AL_t,'r-','linewidth',2); ylabel('stationarity');
ylim([0.5 0.7]); yticks(0.5:0.1:0.7);
xlim([0 1250]); xticks([0 150 300 600 1200]); xlabel('lags (ms)'); 
for t = 1:length(AB_lags_2_t); plot([1 1]*AB_lags_2_t(t), ylim,'k-'); end
ax = gca; ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = 'r';
%% plot CONTENT SPECIFICITY BY LAG
CS_sp = [3 4 7 8];
for t = 1:4
    subplot(2,4,CS_sp(t)); hold on; box on; set(gca,'fontsize',font_size,'YAxisLocation','right'); h = [];
    title('similarity to similar patterns')
    h(1) = fill([T(T1:T2) fliplr(T(T1:T2))], [nanmean(content_spec_all(:,T1:T2),1) zeros(1, X2-X1+1)],...
        [.7 .7 .7], 'edgecolor','none'); 
    h(2) = plot(T, to_same_t(t,:),'color','r','linewidth',2); 
    h(3) = plot(T, to_diff_t(t,:),'color','k','linewidth',2); 
    h(4) = plot(T, CS_t(t,:),'color',[1 0.6 0.6],'linewidth',2); 
    xlim([-250 1500]); plot(xlim, [0 0], 'k-'); xticks(-500:500:1500); 
    ylim([-0.1 0.2]); plot([0 0], ylim, 'k'); yticks([min(ylim) 0 max(ylim)])
    if(mod(t,2)==0); ylabel('similarity'); end
    if(t>2); xlabel('time after T1 (ms)'); else; xticklabels({}); end
    try lag_1 = [int2str(AB_lags_2(t-1)) ' <']; catch; lag_1 = []; end
    title([lag_1 ' lag \leq ' int2str(AB_lags_2(t))])
    if(t==1)
        legend(h, {'specificity - all lags', 'RSA similar patterns', 'RSA different patterns', 'specificity'},...
            'fontsize',10,'location','north','Orientation','vertical');
    end
end
if(exist([filename '/Figures'],'dir')~=7); mkdir([filename '/Figures']); end
saveas(fig, [filename '/Figures/attentional_blink.jpg']); close(fig);

end

