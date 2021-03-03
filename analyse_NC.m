function [] = analyse_NC( filename, trials, simulate )
% set_par_N = {'NC_n_stim', 'NC_pre_stim_L', 'NC_stim_L'};
% set_par_x = {1,1000,500};
% var_par_N = {'NC_stim_tau', 'patterns_t'};
% %var_par_x = {[0.5 1 2 3 4], [1.5 5 10 20 35 60 90], round((0:45:360)/360*(1000/8))};
% var_par_x = {1:5:90, round(((0:pi/8:2*pi)/(2*pi))*(1000/8))};
% %% set new default parameters
% % set any additional parameters from inputs
% for n = 1:length(set_par_N)
%     if(ischar(set_par_x{n})) 
%         par.(set_par_N{n}) = {set_par_x{n}};
%     else; par.(set_par_N{n}) = set_par_x{n};
%     end
% end
% % always set to just NC @ learning phase
% par.sim_order = {'learning'}; % simulate learning / recall 
% par.sim_order_n = {'DL'}; % during learning / after learning
% par.nG = {'NC'}; % include which regions
% par.sim_type = 'learning';
% par.sim_length = par.NC_pre_stim_L*2 + par.NC_stim_L;
% par = set_parameters(par);
% %% calculate all possible parameters combinations
% var_pars = repmat(allcomb(var_par_x{:}), [trials 1]);
% %time = nan(size(var_pars,1),1);
% par.trials = size(var_pars,1);
% par.B = 1;

[ par, var_par_N, var_par_x, var_pars ] = set_pars( trials );

filename = [filename '_' int2str(length(var_par_N)) 'pars_' int2str(par.trials) 'T'];
if(exist(filename,'dir')~=7); mkdir(filename); end
if(exist([filename '/Data'],'dir')~=7); mkdir([filename '/Data']); end
if(exist([filename '/Figures'],'dir')~=7); mkdir([filename '/Figures']); end
if(simulate)
    save([filename '/Data/parameters.mat'],'par')
    save([filename '/Data/varied_parameters.mat'],'var_pars','var_par_N','var_par_x')
    %% loop through all combinations as trials, simulate network
    parfor t = 1:size(var_pars,1)
        tic;
        [ par, var_par_N, ~, var_pars ] = set_pars( trials );
        fprintf('\nsimulating trial %.0f/%.0f ... \n', t, par.trials)
        for n = 1:length(var_par_N)
            par.(var_par_N{n}) = var_pars(t, n);
        end
        if(any(contains(var_par_N, 'patterns_t'))) % if varying stim presentation time
            par.NC_stims_t = var_pars(t, contains(var_par_N, 'patterns_t'));
        else; par.NC_stims_t = 0; % else set to stim onset
        end
        [ I, weight_matrix, weight_matrix_STDP, tau_syn, wm_max, delay, par ] = create_network( par );
        [ sim_data, ~, ~ ] = simulate_network( filename, 'off', ...
            par, I, weight_matrix, weight_matrix_STDP, tau_syn, wm_max, delay);
        if(exist([filename '/Data/T' int2str(t)],'dir')~=7); mkdir([filename '/Data/T' int2str(t)]); end
        %save([filename '/Data/T' int2str(t) '/DL.mat'],'sim_data')
        parsave([filename '/Data/T' int2str(t) '/DL.mat'], sim_data);
        %time(t) = toc(x);
        %min_rem = nanmean(time)/60 * (par.trials - t);
        %H = floor(min_rem/60);
        %M = floor(min_rem - H*60); clear min_rem;
        %fprintf('estimated time remaining: %.0fh %.0fm\n', H, M);
        M = floor(toc/60);
        S = toc - M*60;
        fprintf('\nsimulation of T (%.0f/%.0f) took %.0fM %.0fS ...\n',...
            t,par.trials, M,S)
    end
    %% load & preproceses
    save([filename '/Data/parameters.mat'],'par')
    preprocess_data(filename);
end
load([filename '/Data/varied_parameters.mat'])
load([filename '/Data/parameters.mat'])
load([filename '/Processed Data/freq_data_P1.mat'])
%% calculate psp area and phase
stim_T_phase = [-250 600];% par.NC_stim_L;
stim_T_pow = [50 800];% par.NC_stim_L;

psp_area = var_par_x{1};
pi_i = -pi:pi/16:pi;
phase_comp = nan(numel(psp_area), numel(pi_i));
power_comp = nan(numel(psp_area), numel(pi_i));
count = zeros(numel(psp_area), numel(pi_i));
%% extract all results by phase & psp areas
non_stationarity = zeros(size(var_pars,1), 1);
power_diff = zeros(size(var_pars,1), 1);
pi_cycle = zeros(size(var_pars,1), 1);
phase_at_stim = zeros(size(var_pars,1), 1);
for i = 1:size(var_pars,1)
     % how many cycles to go back when checking phase 
     %(checking at stimulus onset is tricky due to onset of non-stationarity)
    dPh = PHASE.DL.NC(i, 600:par.NC_pre_stim_L-100);
    dPh = dPh(2:end) - dPh(1:end-1);
    pi_cycle(i) = round((2*pi)/mean(dPh(dPh > 0)));
    
    psp_area_t = var_pars(i,1);
    a_i = find(abs(psp_area - psp_area_t) == min(abs(psp_area - psp_area_t)));
    p_i = par.NC_pre_stim_L + ceil(var_pars(i,2)) + ceil(var_pars(i,1));
    stim_T = round(stim_T_phase * (var_pars(i,1)/max(psp_area)));
    non_stationarity(i) = mean(PHASE_NS.DL.NC(i, (stim_T(1):stim_T(2)) + p_i));
    stim_T = round(stim_T_pow * (var_pars(i,1)/max(psp_area)));
    power_diff(i) = mean(POWER.DL.NC(i, (stim_T(1):stim_T(2)) + p_i));
    
    p_i = p_i - pi_cycle(i)*2;
    phase_at_stim(i) = PHASE.DL.NC(i, p_i);
    
    p_i = find(abs(pi_i - phase_at_stim(i)) == min(abs(pi_i - phase_at_stim(i))));
    phase_comp(a_i, p_i) = nanmean([phase_comp(a_i, p_i) non_stationarity(i)]);
    power_comp(a_i, p_i) = nanmean([power_comp(a_i, p_i) power_diff(i)]);
    count(a_i, p_i) = count(a_i, p_i) + 1;
end
x_i = find(psp_area >= 20, 1, 'first');
mean_phase = nanmean(phase_comp(1:x_i,:), 1); mean_phase = smoothdata(mean_phase);
mean_power = nanmean(power_comp(1:x_i,:), 1); mean_power = smoothdata(mean_power);

most_sync = find(max(power_diff) == power_diff,1,'first');
most_desync = find(min(power_diff) == power_diff,1,'first');
most_stat = find(max(non_stationarity) == non_stationarity,1,'first');
most_nonstat = find(min(non_stationarity) == non_stationarity,1,'first');

trials = [most_stat most_nonstat; most_desync most_sync];
long_tit = {'stationary','non-stationary'; 'de-synched' 'synched'};

a_sync = find(var_pars(most_sync,1) == var_par_x{1});
p_sync = find(abs(pi_i - phase_at_stim(most_sync)) == min(abs(pi_i - phase_at_stim(most_sync))));
a_desync = find(var_pars(most_desync,1) == var_par_x{1});
p_desync = find(abs(pi_i - phase_at_stim(most_desync)) == min(abs(pi_i - phase_at_stim(most_desync))));
a_stat = find(var_pars(most_stat,1) == var_par_x{1});
p_stat = find(abs(pi_i - phase_at_stim(most_stat)) == min(abs(pi_i - phase_at_stim(most_stat))));
a_nonstat = find(var_pars(most_nonstat,1) == var_par_x{1});
p_nonstat = find(abs(pi_i - phase_at_stim(most_nonstat)) == min(abs(pi_i - phase_at_stim(most_nonstat))));

psp_diff = mode(psp_area(2:end) - psp_area(1:end-1))/2;
pi_diff = mode(pi_i(2:end) - pi_i(1:end-1))/2;
font_size = 17;
%% make figure
fig = figure(); set(fig, 'position', [0 0 1800 1000])
%% phase figure
ax(1) = subplot(7,6,[2 3 8 9 14 15 20 21 26 27]); box on; hold on;
imagesc(pi_i,psp_area,  phase_comp); 
axis xy;  set(gca, 'yaxislocation','right')
set(gca, 'fontsize', font_size); 
xticks([-pi -pi/2 0 pi/2 pi]); xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
xlim([min(pi_i)-(pi/length(pi_i)), max(pi_i)+(pi/length(pi_i))]); 
ylim([min(psp_area)-2.5 max(psp_area)+2.5]); yticks([10 20 30 60 70 80]);
c1 = colorbar('location','northoutside'); colormap(ax(1),'parula');
caxis([0.4 1]); c1.Label.String = 'stationarity';
c1.Label.FontSize = font_size + 6;
plot([-2*pi 2*pi], [1 1]*psp_area(x_i)+psp_diff, 'k-', 'linewidth', 2);

% red boxes
Y1 = psp_area(a_stat)-psp_diff; Y2 = psp_area(a_stat+1)-psp_diff; 
X1 = pi_i(p_stat)-pi_diff; X2 = pi_i(p_stat+1)-pi_diff;
plot([X1,X2,X2,X1,X1],[Y1,Y1,Y2,Y2,Y1],'r-','linewidth',2);
Y1 = psp_area(a_nonstat)-psp_diff; Y2 = psp_area(a_nonstat+1)-psp_diff; 
X1 = pi_i(p_nonstat)-pi_diff; X2 = pi_i(p_nonstat+1)-pi_diff;
plot([X1,X2,X2,X1,X1],[Y1,Y1,Y2,Y2,Y1],'r-','linewidth',2);
%% power figure
ax(2) = subplot(7,6,[4 5 10 11 16 17 22 23 28 29]); box on; hold on;
imagesc(pi_i,psp_area,  power_comp); 
axis xy; yticklabels({}); ylabel('F(\alpha) \tau (ms)'); 
set(gca, 'fontsize', font_size); 
xticks([-pi -pi/2 0 pi/2 pi]); xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
c2 = colorbar('location','northoutside'); colormap(ax(2),'parula');
caxis([-60 10]); c2.Label.String = 'power difference to baseline (%)';
c2.Label.FontSize = font_size + 6;
xlim([min(pi_i)-(pi/length(pi_i)), max(pi_i)+(pi/length(pi_i))]); 
ylim([min(psp_area)-2.5 max(psp_area)+2.5]);
plot([-2*pi 2*pi], [1 1]*psp_area(x_i)+psp_diff, 'k-', 'linewidth', 2);

% red boxes
Y1 = psp_area(a_sync)-psp_diff; Y2 = psp_area(a_sync+1)-psp_diff; 
X1 = pi_i(p_sync)-pi_diff; X2 = pi_i(p_sync+1)-pi_diff;
plot([X1,X2,X2,X1,X1],[Y1,Y1,Y2,Y2,Y1],'r-','linewidth',2);
Y1 = psp_area(a_desync-1)-psp_diff; Y2 = psp_area(a_desync)-psp_diff; 
X1 = pi_i(p_desync)-pi_diff; X2 = pi_i(p_desync+1)-pi_diff;
plot([X1,X2,X2,X1,X1],[Y1,Y1,Y2,Y2,Y1],'r-','linewidth',2);
%% upper bottom mean x axis (phase) plots
y_lims_ph = [floor(min([mean_phase ])*100)/100 ...
    ceil(max([mean_phase ])*100)/100];
subplot(7,6,[32 33]); hold on; box on; 
set(gca, 'fontsize', font_size, 'XAxisLocation','bottom'); 
yyaxis left; plot(pi_i, mean_phase, 'r','linewidth',2); 
yticks(y_lims_ph); ylim(y_lims_ph); yticklabels({}); %ylabel('stat.'); 
yyaxis right; plot(-pi:0.1:pi, cos(-pi:0.1:pi),'k','linewidth',2); yticks([-1 1]);
xlim([min(pi_i), max(pi_i)]); yticklabels({}); %xlabel('phase');
xticks([-pi -pi/2 0 pi/2 pi]); xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ax(3) = gca; ax(3).YAxis(1).Color = 'r'; ax(3).YAxis(2).Color = 'k';

y_lims_pow = [floor(min([mean_power ])) 0];% ...
    %ceil(max([mean_power ]))];
subplot(7,6,[34 35]); hold on; box on; 
set(gca, 'fontsize', font_size, 'XAxisLocation','bottom'); 
yyaxis right; plot(pi_i, mean_power, 'b','linewidth',2); 
ylim(y_lims_pow ); yticks(y_lims_pow ); yticklabels({}); %ylabel('pow.%');
yyaxis left; plot(-pi:0.1:pi, cos(-pi:0.1:pi),'k','linewidth',2); 
yticklabels({}); yticklabels({});% ylabel('cos(\pi)');
xlim([min(pi_i), max(pi_i)]); %xlabel('phase');
xticks([-pi -pi/2 0 pi/2 pi]); xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
ax(4) = gca; ax(4).YAxis(1).Color = 'k'; ax(4).YAxis(2).Color = 'b';

%% lower bottom mean y axis (TS & time)
subplot(7,6,[38 39]); hold on; box on; 
set(gca, 'fontsize', font_size); 
yyaxis left; plot(var_par_x{1},mean(phase_comp,2),'r-','linewidth',2); 
ylim([0.4 1]); yticks([0.4 1]); yticklabels({}); %ylabel('stat.')
yyaxis right; plot(var_par_x{1},mean(power_comp,2),'b-','linewidth',2); 
ylim([-60 10]); yticks([-60 0]); yticklabels({});
%plot([psp_area(1) psp_area(end)], [0 0],'k-')
ax(5) = gca; ax(5).YAxis(1).Color = 'r'; ax(5).YAxis(2).Color = 'b';
xlim([min(var_par_x{1}) max(var_par_x{1})]); xlabel('F(\alpha) \tau (ms)'); 

subplot(7,6,[40 41]); hold on; box on; 
set(gca, 'fontsize', font_size); 
stim_T = -100:1000;
yyaxis left; plot(stim_T,mean(POWER.DL.NC(:,par.NC_pre_stim_L+stim_T),1),'b-','linewidth',2);
ylim([-60 10]); yticklabels({}); yticklabels({}); %ylabel('pow.%')
%plot([stim_T(1) stim_T(end)], [0 0],'k-')
yyaxis right; plot(stim_T,mean(PHASE_NS.DL.NC(:,par.NC_pre_stim_L+stim_T),1),'r-','linewidth',2); 
ylim([0.4 1]); yticks([0.4 1]); yticklabels({}); %ylabel('stat.')
plot([0 0],[0 1],'k-','linewidth',1);
xlim([stim_T(1) stim_T(end)]); xlabel('time from stimulus onset (ms)')
ax(6) = gca; ax(6).YAxis(1).Color = 'b'; ax(6).YAxis(2).Color = 'r';

%% individual cases plots
sp_power = [1 6; 19 24];
sp_phase = [7 12; 25 30];
stim_T = par.NC_pre_stim_L - 149 : par.NC_pre_stim_L + par.NC_stim_L;
for i = 1:2
    for j = 1:2
        t = trials(i,j);
        subplot(7, 6, sp_power(j,i)); hold on; box off;
        title(long_tit{i,j})
        yyaxis left; plot(stim_T, LFP_all.DL.NC(t,stim_T),'k-'); 
        plot(stim_T, zeros(size(stim_T)), 'k-'); 
        ylim([-.5 .5]); yticks([-.5 0 .5]); yticklabels({'.5','0','.5'})
        plot(par.NC_pre_stim_L+var_pars(t,2)+var_pars(t,1)*[1 1], [-.5 .5],'k-','linewidth',3);
        yyaxis right; plot(stim_T, POWER.DL.NC(t,stim_T),'b-','linewidth',2); ylim([-100 100]);
        xlim([min(stim_T) max(stim_T)]); xticklabels({});
        ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'b'; 
        ax.Color = 'none'; ax.FontSize = font_size;   axis off;
        
        subplot(7, 6, sp_phase(j,i)); hold on; box off;
        yyaxis left; plot(stim_T, PHASE.DL.NC(t,stim_T),'k'); ylim([-pi pi]);
        %plot(stim_T, zeros(size(stim_T)), 'k:'); 
        plot(par.NC_pre_stim_L+var_pars(t,2)+var_pars(t,1)*[1 1], [-pi pi],'k-','linewidth',3);
        xlim([min(stim_T) max(stim_T)]); xticklabels({});
        yticks([-pi pi]); yticklabels({'-\pi', '\pi'});
        yyaxis right; plot(stim_T, PHASE_NS.DL.NC(t,stim_T),'r-','linewidth',2); 
        ylim([0 1]); yticks([0 1]);
        ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'r'; 
        ax.Color = 'none'; ax.FontSize = font_size; axis off;
    end
end

saveas(fig, [filename '/Figures/phase-power_par-var.jpg']); close(fig);
end

function [] = parsave(filename, sim_data); save(filename, 'sim_data'); end
function [ par, var_par_N, var_par_x, var_pars ] = set_pars( trials )
    set_par_N = {'NC_n_stim', 'NC_pre_stim_L', 'NC_stim_L'};
    set_par_x = {1,1500,500};
    var_par_N = {'NC_stim_tau', 'patterns_t'};
    %var_par_x = {[0.5 1 2 3 4], [1.5 5 10 20 35 60 90], round((0:45:360)/360*(1000/8))};
    var_par_x = {1:5:90, round(((0:pi/8:2*pi)/(2*pi))*(1000/8))};
    %% set new default parameters
    % set any additional parameters from inputs
    for n = 1:length(set_par_N)
        if(ischar(set_par_x{n})) 
            par.(set_par_N{n}) = {set_par_x{n}};
        else; par.(set_par_N{n}) = set_par_x{n};
        end
    end
    % always set to just NC @ learning phase
    par.sim_order = {'learning'}; % simulate learning / recall 
    par.sim_order_n = {'DL'}; % during learning / after learning
    par.nG = {'NC'}; % include which regions
    par.sim_type = 'learning';
    par.sim_length = par.NC_pre_stim_L*2 + par.NC_stim_L;
    par = set_parameters(par);
    %% calculate all possible parameters combinations
    var_pars = repmat(allcomb(var_par_x{:}), [trials 1]);
    %time = nan(size(var_pars,1),1);
    par.trials = size(var_pars,1);
    par.B = 1;
    par.n_ID = cell(par.network_size,1);
    par.n_ID(1:par.n_NC_E) = {'NC_E'};
    par.n_ID(par.n_NC_E+1:par.n_NC_I+par.n_NC_E) = {'NC_I'};
end