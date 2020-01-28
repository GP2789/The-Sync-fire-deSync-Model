function [  ] = make_videos( filename )
% This function creates frames for the videos shown in the supplementary
% section of the paper. These frames can be used by video editing software
% to create a video. 

%% load variables
load([filename '/Data/parameters' ]);
load([filename '/Processed Data/calcium_data.mat']);
load([filename '/Processed Data/freq_data.mat']);
par.sim_length = par.pre_stim_L * 2 + par.NC_stim_L;
if(exist([filename '/Graphs/Videos/Encoding'], 'dir')~=7); mkdir([filename '/Graphs/Videos/Encoding']); end
if(exist([filename '/Graphs/Videos/Recall'], 'dir')~=7); mkdir([filename '/Graphs/Videos/Recall']); end

A_v = 0.15; A_s = 0.5;

%% loop through patterns
load([filename '/Processed Data/spike_data_P1.mat']);
load([filename '/Processed Data/intrinsic_data_P1.mat']);
for t = 1:par.trials
    for n = 1:length(par.sim_order_n) % loop through simulation stages
        Vm = -120:1:60; x = Vm(1:end-1); % voltage range for histcounts
        Vm_O = -70:1:-50; x_O = Vm_O(1:end-1); % voltage range for histcounts

        %% get ID of synfire chain units in voltage struct
        TC_N = fieldnames(V_all.([par.sim_order_n{n}]).TC);
        TC_L1_i = contains(TC_N, 'L1'); TC_L2_i = contains(TC_N, 'L2');
        TC_E_i = contains(TC_N, 'E'); TC_I_i = contains(TC_N, 'I');
        TC_P_i = contains(TC_N, 'P'); TC_O_i = contains(TC_N, 'O');

        TC_L1_E = TC_N(TC_L1_i & TC_E_i); TC_L1_P = TC_N(TC_L1_i & TC_P_i);
        TC_L1_I = TC_N(TC_L1_i & TC_I_i); TC_L1_O = TC_N(TC_L1_i & TC_O_i);
        TC_L2_E = TC_N(TC_L2_i & TC_E_i); TC_L2_P = TC_N(TC_L2_i & TC_P_i);
        TC_L2_I = TC_N(TC_L2_i & TC_I_i); TC_L2_O = TC_N(TC_L2_i & TC_O_i);

        % get TCL1 syn input
        TC_L1_E_syn = zeros(1, par.sim_length, length(par.uTS));
        for i = 1:length(TC_L1_E); TC_L1_E_syn = TC_L1_E_syn + SYN_all.([par.sim_order_n{n}]).TC.([TC_L1_E{i}]); end
        TC_L1_E_syn = squeeze(TC_L1_E_syn / length(TC_L1_E));

         % get TCL2 syn input
        TC_L2_E_syn = zeros(1, par.sim_length, length(par.uTS));
        for i = 1:length(TC_L2_E); TC_L2_E_syn = TC_L2_E_syn + SYN_all.([par.sim_order_n{n}]).TC.([TC_L2_E{i}]); end
        TC_L2_E_syn = squeeze(TC_L2_E_syn / length(TC_L2_E));

        %%  get ID of binding pool units in voltage struct
        BP_N = fieldnames(V_all.([par.sim_order_n{n}]).BP);
        BP_E = BP_N(contains(BP_N, 'E')); BP_I = BP_N(contains(BP_N, 'I')); BP_O = BP_N(contains(BP_N, 'O'));
        BP_E_Vm = []; BP_I_Vm = []; BP_O_Vm = []; 
        for i = 1:length(BP_E); BP_E_Vm = [BP_E_Vm; V_all.([par.sim_order_n{n}]).BP.([BP_E{i}])]; end
        for i = 1:length(BP_I); BP_I_Vm = [BP_I_Vm; V_all.([par.sim_order_n{n}]).BP.([BP_I{i}])]; end
        for i = 1:length(BP_O); BP_O_Vm = [BP_O_Vm; V_all.([par.sim_order_n{n}]).BP.([BP_O{i}])]; end

        % get BP syn input
        BP_E_syn = zeros(1, par.sim_length, length(par.uTS));
        for i = 1:length(BP_E); BP_E_syn = BP_E_syn + SYN_all.([par.sim_order_n{n}]).BP.([BP_E{i}]); end
        BP_E_syn = squeeze(BP_E_syn / length(BP_E));

         %%  get NC VM
        NC_N = fieldnames(V_all.([par.sim_order_n{n}]).NC);
        NC_E = NC_N(contains(NC_N, 'E')); NC_I = NC_N(contains(NC_N, 'I')); NC_E_Vm = []; NC_I_Vm = [];
        for i = 1:length(NC_E); NC_E_Vm = [NC_E_Vm; V_all.([par.sim_order_n{n}]).NC.([NC_E{i}])]; end
        for i = 1:length(NC_I); NC_I_Vm = [NC_I_Vm; V_all.([par.sim_order_n{n}]).NC.([NC_I{i}])]; end

        % get NC syn input
        NC_E_syn = zeros(1, par.sim_length, length(par.uTS));
        for i = 1:length(NC_E); NC_E_syn = NC_E_syn + SYN_all.([par.sim_order_n{n}]).NC.([NC_E{i}]); end
        NC_E_syn = squeeze(NC_E_syn / length(NC_E));

        W = gausswin(25); W = W / sum(W); 
        W_O = gausswin(5); W_O = W_O / sum(W_O); 

        ms_start = par.TC_pre_stim_L-150; ms_end = par.TC_pre_stim_L + 3250;

        %% EXTRACT TC SPIKES
        TC_sp_E = []; TC_sp_I = []; TC_sp_P = []; TC_sp_O = []; 
        for L = 1:2
            for i = 1:length(n_G_ID{1}{L}{1})
                TC_sp_E = [TC_sp_E; SPIKES{t}{n}{1}{L}{1}{i}];
                TC_sp_P = [TC_sp_P; SPIKES{t}{n}{1}{L}{2}{i}];
                TC_sp_I = [TC_sp_I; SPIKES{t}{n}{1}{L}{3}{i}];
                TC_sp_O = [TC_sp_O; SPIKES{t}{n}{1}{L}{4}{i}];
            end
        end

        %
        f(1) = figure(1); set(f(1), 'Position', [0 0 1300 930]); % Intrinsic
        set(f(1),'defaultAxesColorOrder',[[1 0 0]; [0 0 0]]);
        Vm_BP_E = -365:60; x_BP_E = Vm_BP_E(1:end-1); % voltage range for histcounts
            Vm_BP_I = -220:60; x_BP_I = Vm_BP_I(1:end-1);
            Vm_NC = -180:60; x_NC = Vm_NC(1:end-1); 

        if(n == 2) % RETREIVAL SYNAPTIC INPUT & VOLTAGES
            for ms = ms_start : ms_end
                %% plot NC phase
                ax = subplot(10, 4, 2:4 ); box on; 
                plot(T(1:ms), PHASE.AL.NC.B1.T1(1:ms), 'k', 'linewidth', 2); hold on; 
                plot([ms ms], [-pi pi], 'k-'); 
                xlim([ms_start ms_end]); ylim([-pi pi]); 
                set(ax, 'yaxislocation', 'right'); 
                xlabel('time (ms)'); title('NC phase'); hold off;
                ax = ancestor(ax, 'axes');
                ax.YTick = [-pi pi]; ax.YTickLabel = {'-\pi','\pi'};
                yrule = ax.YAxis; yrule.FontSize = 18;

                %% plot TC L1 voltage histograms
                subplot(10, 4, [33 37]); box on; % TC L1 E
                for L = 1:2; fill([x fliplr(x)], [filter(W, 1, histcounts(V_all.DL.TC.([TC_L1_E{L}])(:,ms), Vm) / 5) zeros(size(x))], [0 0 1*(L/2)], 'edgecolor', [0 0 1*(L/2)]); hold on; alpha(A_v); end 
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
                ylabel('SC L1 E','fontweight','bold'); yticklabels('');  xlabel('V_m'); hold off;

                ax = subplot(10, 4, [34:36 38:40]); 
                for i = 1:length(par.uTS)
                    if(sum(TC_L1_E_syn(:,i)) < 0); col = 'r'; 
                    elseif(sum(TC_L1_E_syn(:,i)) > 0); col = 'b'; 
                    else; col = [];
                    end
                    if(isempty(col) ~= 1)
                        fill([T(ms_start:ms) fliplr(T(ms_start:ms))], ...
                            [TC_L1_E_syn(ms_start:ms, i)' zeros(1, ms+1-ms_start)], col, 'edgecolor', col); alpha(A_s); hold on
                    end
                end
                y_lims = [max([-max(max(TC_L1_E_syn))*3 min(min(TC_L1_E_syn))]) max(max(TC_L1_E_syn))]*1.2;
                plot([ms ms], y_lims, 'k-'); plot([ms_start ms], [0 0], 'k-'); set(ax, 'yaxislocation', 'right')
                xlim([ms_start ms_end]); ylim(y_lims); title('SC L1 E synaptic inputs'); xlabel('time (ms)'); ylabel('current');
                hold off;

                %% plot TC L2 voltage histograms
                subplot(10, 4, [25 29]); box on; % TC L2 E
                for L = 1:8; fill([x fliplr(x)], [filter(W, 1, histcounts(V_all.DL.TC.([TC_L2_E{L}])(:,ms), Vm) / 5) zeros(size(x))], [0 0 1*(L/8)], 'edgecolor', [0 0 1*(L/8)]); hold on; alpha(A_v); end
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
                ylabel('SC L2 E','fontweight','bold'); yticklabels('');   hold off

                ax = subplot(10, 4, [26:28 30:32]);
                for i = 1:length(par.uTS)
                    if(sum(TC_L2_E_syn(:,i)) < 0); col = 'r';
                    elseif(sum(TC_L2_E_syn(:,i)) > 0); col = 'b';
                    else; col = [];
                    end
                    if(isempty(col) ~= 1)
                        fill([T(ms_start:ms) fliplr(T(ms_start:ms))], ...
                            [TC_L2_E_syn(ms_start:ms, i)' zeros(1, ms+1-ms_start)], col, 'edgecolor', col); alpha(A_s); hold on
                    end
                end
                y_lims = [max([-max(max(TC_L2_E_syn))*3 min(min(TC_L2_E_syn))]) max(max(TC_L2_E_syn))]*1.2;
                plot([ms ms], y_lims, 'k-'); plot([ms_start ms], [0 0], 'k-'); set(ax, 'yaxislocation', 'right')
                xlim([ms_start ms_end]); xticklabels(''); 
                ylim(y_lims); title('SC L2 E synaptic inputs'); ylabel('current');
                hold off;

                %%  NC E
                subplot(10, 4, [9 13]); box on;
                fill([x_NC fliplr(x_NC)], [filter(W, 1, histcounts(NC_E_Vm(:,ms), Vm_NC) / par.n_NC_E) zeros(size(x_NC))], [0 0 1], 'edgecolor', [0 0 1]); alpha(A_v);
                xlim([min(Vm_NC) max(Vm_NC)]); ylim([0 max(W)]); 
                ylabel('NC E','fontweight','bold'); yticklabels(''); 

                ax = subplot(10, 4, [10:12 14:16]); 
                for i = 1:length(par.uTS)
                    A_s = 0.5;
                    if(sum(NC_E_syn(:,i)) < 0); col = 'r';
                    elseif(sum(NC_E_syn(:,i)) > 0); col = 'b';
                    else; col = [];
                    end
                    if(isempty(col) ~= 1)
                        fill([T(ms_start:ms) fliplr(T(ms_start:ms))], ...
                            [NC_E_syn(ms_start:ms, i)' zeros(1, ms+1-ms_start)], col, 'edgecolor', col); alpha(A_s); hold on
                    end
                end
                y_lims = [max([-max(max(NC_E_syn))*3 min(min(NC_E_syn))]) max(max(NC_E_syn))]*1.2;
                plot([ms ms], y_lims, 'k-'); plot([ms_start ms], [0 0], 'k-'); set(ax, 'yaxislocation', 'right')
                xlim([ms_start ms_end]); xticklabels(''); 
                ylim(y_lims); title('NC E synaptic inputs');  ylabel('current');
                hold off;

                %% BP E 
                subplot(10, 4, [17 21]); box on; % BP E
                fill([x fliplr(x)], [filter(W, 1, histcounts(BP_E_Vm(:,ms), Vm) / size(BP_E_Vm,1)) zeros(size(x))], [0 0 1], 'edgecolor', [0 0 1]); alpha(A_v);
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
                ylabel('BP E','fontweight','bold'); yticklabels(''); 

                ax = subplot(10, 4, [18:20 22:24]); 
                for i = 1:length(par.uTS)
                    if(sum(BP_E_syn(:,i)) < 0); col = 'b';
                    elseif(sum(BP_E_syn(:,i)) > 0); col = 'b';
                    else; col = [];
                    end
                    if(isempty(col) ~= 1)
                        fill([T(ms_start:ms) fliplr(T(ms_start:ms))], ...
                            [BP_E_syn(ms_start:ms, i)' zeros(1, ms+1-ms_start)], col, 'edgecolor', col); alpha(A_s); hold on
                    end
                end
                y_lims = [0 max(max(BP_E_syn))]*1.1;
                plot([ms ms], y_lims, 'k-'); plot([ms_start ms], [0 0], 'k-'); set(ax, 'yaxislocation', 'right')
                xlim([ms_start ms_end]); xticklabels(''); 
                ylim(y_lims); title('BP E synaptic inputs'); ylabel('current');
                hold off;

                saveas(f(1), [filename '/Graphs/Videos/Recall/' int2str(ms) '.jpg'])
            end
            close(f(1)); 
        elseif(n == 1) % ENCODING CALCIUM & VOLTAGES
            ms_start = par.TC_pre_stim_L - 150; ms_end = par.NC_pre_stim_L + par.NC_stim_L + 1000;
            Vm_BP_E = -365:60; x_BP_E = Vm_BP_E(1:end-1); % voltage range for histcounts
            Vm_BP_I = -220:60; x_BP_I = Vm_BP_I(1:end-1);
            Vm_NC = -180:60; x_NC = Vm_NC(1:end-1); % voltage range for histcounts

            % BP -> NC CALCIUM & WEIGHT CHANGES
            calc_t = calcium_all.B1.DL{2, 3}{1, 1}{1, 1}; calc_BP_NC = [];
            p_w_t = p_w_all.B1.DL{2, 3}{1, 1}{1, 1}; p_w_BP_NC = [];
            for i = 1:size(calc_t, 1); calc_BP_NC = [calc_BP_NC; calc_t{i}(2,:)]; p_w_BP_NC = [p_w_BP_NC; p_w_t{i}(1,:)]; end
            calc_BP_NC = max(calc_BP_NC); p_w_BP_NC = mean(p_w_BP_NC);

            % BP <-> BP CALCIUM & WEIGHT CHANGES
            calc_t = calcium_all.B1.DL{2, 2}{1, 1}{1, 1}; calc_BP_BP = [];
            p_w_t = p_w_all.B1.DL{2, 2}{1, 1}{1, 1}; p_w_BP_BP = [];
            for i = 1:size(calc_t, 1); calc_BP_BP = [calc_BP_BP; calc_t{i, i}(2,:)]; p_w_BP_BP = [p_w_BP_BP; p_w_t{i, i}(1,:)]; end
            calc_BP_BP = max(calc_BP_BP); p_w_BP_BP = mean(p_w_BP_BP);

            % TC L1 -> NC CALCIUM & WEIGHT CHANGES
            calc_t = calcium_all.B1.DL{1, 2}{1, 1}{1, 1}; calc_TCL1_BP = [];
            p_w_t = p_w_all.B1.DL{1, 2}{1, 1}{1, 1}; p_w_TCL1_BP = [];
            for i = 1:size(calc_t, 1)
                for i2 = 1:size(calc_t,2)
                    calc_TCL1_BP = [calc_TCL1_BP; calc_t{i, i2}(2,:)]; p_w_TCL1_BP = [p_w_TCL1_BP; p_w_t{i, i2}(1,:)]; 
                end
            end
            calc_TCL1_BP = max(calc_TCL1_BP); p_w_TCL1_BP = mean(p_w_TCL1_BP);

            % TC L1 -> NC CALCIUM & WEIGHT CHANGES
            calc_t = calcium_all.B1.DL{1, 2}{2, 1}{1, 1}; calc_TCL2_BP = [];
            p_w_t = p_w_all.B1.DL{1, 2}{2, 1}{1, 1}; p_w_TCL2_BP = [];
            for i = 1:size(calc_t, 1)
                for i2 = 1:size(calc_t,2)
                    calc_TCL2_BP = [calc_TCL2_BP; calc_t{i, i2}(2,:)]; p_w_TCL2_BP = [p_w_TCL2_BP; p_w_t{i, i2}(1,:)]; 
                end
            end
            calc_TCL2_BP = max(calc_TCL2_BP); p_w_TCL2_BP = mean(p_w_TCL2_BP);

            %% LOOP THROUGH MS
            for ms = ms_start : ms_end
                subplot(10, 5, [11 16]); box on; % NC E
                fill([x_NC fliplr(x_NC)], [filter(W, 1, histcounts(NC_E_Vm(:,ms), Vm_NC) / par.n_NC_E) zeros(size(x_NC))], [0 0 1], 'edgecolor', [0 0 1]); alpha(A_v);
                xlim([min(Vm_NC) max(Vm_NC)]); ylim([0 max(W)]); 
                ylabel('NC E','fontweight','bold'); yticklabels(''); 

                subplot(10, 5, [12 17]); box on; % NC I
                fill([x fliplr(x)], [filter(W, 1, histcounts(NC_I_Vm(:,ms), Vm) / par.n_NC_I) zeros(size(x))], [1 0 0], 'edgecolor', [1 0 0]); alpha(A_v);
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
                ylabel('NC I','fontweight','bold'); yticklabels(''); 

                subplot(10, 5, [21 26]); box on; % BP E
                fill([x_BP_E fliplr(x_BP_E)], [filter(W, 1, histcounts(BP_E_Vm(:,ms), Vm_BP_E) / size(BP_E_Vm,1)) zeros(size(x_BP_E))], [0 0 1], 'edgecolor', [0 0 1]); alpha(A_v);
                xlim([min(Vm_BP_E) max(Vm_BP_E)]); ylim([0 max(W)]); 
                ylabel('BP E','fontweight','bold'); yticklabels(''); 
                
                subplot(10, 5, [22 27]); box on; % BP O
                fill([x fliplr(x)], [filter(W, 1, histcounts(BP_O_Vm(:,ms), Vm) / size(BP_O_Vm,1)) zeros(size(x))], [1 0 0], 'edgecolor', [1 0 0]); alpha(A_v);
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
                ylabel('BP O','fontweight','bold'); yticklabels('');  

                subplot(10, 5, [41 46]); box on; % TC L1 E
                for L = 1:2; fill([x fliplr(x)], [filter(W, 1, histcounts(V_all.DL.TC.([TC_L1_E{L}])(:,ms), Vm) / 5) zeros(size(x))], [0 0 1*(L/2)], 'edgecolor', [0 0 1*(L/2)]); hold on; alpha(A_v); end 
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
                ylabel('SC L1 E','fontweight','bold'); yticklabels('');  xlabel('V_m'); hold off;

                subplot(10, 5, [42 47]); box on; % TC L1 O
                for L = 1:2; fill([x fliplr(x)], [filter(W, 1, histcounts(V_all.DL.TC.([TC_L1_O{L}])(:,ms), Vm)) zeros(size(x))], [1*(L/2) 0 0], 'edgecolor', [1*(L/2) 0 0]); hold on; alpha(A_v); end 
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
               ylabel('SC L1 O','fontweight','bold'); yticklabels('');  xlabel('V_m'); hold off;

                subplot(10, 5, [31 36]); box on; % TC L2 E
                for L = 1:8; fill([x fliplr(x)], [filter(W, 1, histcounts(V_all.DL.TC.([TC_L2_E{L}])(:,ms), Vm) / 5) zeros(size(x))], [0 0 1*(L/8)], 'edgecolor', [0 0 1*(L/8)]); hold on; alpha(A_v); end 
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
                ylabel('SC L2 E','fontweight','bold'); yticklabels('');  hold off;

                subplot(10, 5, [32 37]); box on; 
                for L = 1:8; fill([x fliplr(x)], [filter(W, 1, histcounts(V_all.DL.TC.([TC_L2_O{L}])(:,ms), Vm)) zeros(size(x))], [1*(L/8) 0 0], 'edgecolor', [1*(L/8) 0 0]); hold on; alpha(A_v); end 
                xlim([min(Vm) max(Vm)]); ylim([0 max(W)]); 
                ylabel('SC L2 O','fontweight','bold'); yticklabels('');  hold off;

                ax = subplot(10, 5, 3:5 ); box on; 
                plot(T(1:ms), PHASE.DL.NC.B1.T1(1:ms), 'k', 'linewidth', 2); hold on; 
                plot([ms ms], [-pi pi], 'k-'); 
                xlim([ms_start ms_end]); ylim([-pi pi]); 
                set(ax, 'yaxislocation', 'right'); 
                xlabel('time (ms)'); title('NC phase'); hold off;
                ax = ancestor(ax, 'axes');
                ax.YTick = [-pi pi]; ax.YTickLabel = {'-\pi','\pi'};
                yrule = ax.YAxis; yrule.FontSize = 18;

                subplot(10, 5, [13:15 18:20]); box on; % BP -> NC CALC & WEIGHTS
                yyaxis left; plot(T_f(1:round(ms/10)), calc_BP_NC(1:round(ms/10)), 'r', 'linewidth', 2); hold on;
                plot([1 T_f(round(ms/10))], [1 1]*par.BP_NC_T_p, 'b--'); 
                plot([1 T_f(round(ms/10))], [1 1]*par.BP_NC_T_d, 'r--'); hold off;
                ylim([0 max(calc_BP_NC) * 1.2]); ylabel('calcium');
                yyaxis right; plot(T_f(1:round(ms/10)), p_w_BP_NC(1:round(ms/10)), 'k', 'linewidth', 2); hold on;
                ylim([0 max(p_w_BP_NC) * 1.2]); xlim([ms_start ms_end]); ylabel('weights (\rho)');
                plot([ms ms], [0 max(p_w_BP_NC) * 1.2], 'k-'); hold off;
                title('BP \Rightarrow NC'); xticklabels({});

                subplot(10, 5, [23:25 28:30]); box on; % BP <-> BP CALC & WEIGHTS
                yyaxis left; plot(T_f(1:round(ms/10)), calc_BP_BP(1:round(ms/10)), 'r', 'linewidth', 2); hold on;
                plot([1 T_f(round(ms/10))], [1 1]*par.BP_BP_T_p, 'b--'); 
                plot([1 T_f(round(ms/10))], [1 1]*par.BP_BP_T_d, 'r--'); hold off;
                ylim([0 par.BP_BP_T_p * 1.2]); ylabel('calcium');
                yyaxis right; plot(T_f(1:round(ms/10)), p_w_BP_BP(1:round(ms/10)), 'k', 'linewidth', 2); hold on;
                ylim([0 min(max(p_w_BP_BP) * 1.2, 1)]); xlim([ms_start ms_end]); ylabel('weights (\rho)');
                plot([ms ms], [0 max(p_w_BP_BP) * 1.2], 'k-'); hold off;
                title('BP \Leftrightarrow BP'); xticklabels({});

                subplot(10, 5, [33:35 38:40]); box on; % TC L2 -> BP CALC & WEIGHTS
                yyaxis left; plot(T_f(1:round(ms/10)), calc_TCL2_BP(1:round(ms/10)), 'r', 'linewidth', 2); hold on;
                plot([1 T_f(round(ms/10))], [1 1]*par.TC_BP_T_p, 'b--'); 
                plot([1 T_f(round(ms/10))], [1 1]*par.TC_BP_T_d, 'r--'); hold off;
                ylim([0 max(calc_TCL2_BP) * 1.2]); ylabel('calcium');
                yyaxis right; plot(T_f(1:round(ms/10)), p_w_TCL2_BP(1:round(ms/10)), 'k', 'linewidth', 2); hold on;
                ylim([0 max(p_w_TCL2_BP) * 1.2]); xlim([ms_start ms_end]); ylabel('weights (\rho)');
                plot([ms ms], [0 max(p_w_TCL2_BP) * 1.2], 'k-'); hold off;
                title('SC L2 \Rightarrow BP'); xticklabels({});

                subplot(10, 5, [43:45 48:50]); box on; % TC L1 -> BP CALC & WEIGHTS
                yyaxis left; plot(T_f(1:round(ms/10)), calc_TCL1_BP(1:round(ms/10)), 'r', 'linewidth', 2); hold on;
                plot([1 T_f(round(ms/10))], [1 1]*par.TC_BP_T_p, 'b--'); 
                plot([1 T_f(round(ms/10))], [1 1]*par.TC_BP_T_d, 'r--'); hold off;
                ylim([0 max(calc_TCL1_BP) * 1.2]); ylabel('calcium');
                yyaxis right; plot(T_f(1:round(ms/10)), p_w_TCL1_BP(1:round(ms/10)), 'k', 'linewidth', 2); hold on;
                ylim([0 max(p_w_TCL1_BP) * 1.2]); xlim([ms_start ms_end]); ylabel('weights (\rho)');
                plot([ms ms], [0 max(p_w_TCL1_BP) * 1.2], 'k-'); hold off;
                title('SC L1 \Rightarrow BP'); xlabel('time (ms)')

                saveas(f(1), [filename '/Graphs/Videos/Encoding/' int2str(ms) '.jpg'])
            end
            close(f(1));
        end
    end
end

end

