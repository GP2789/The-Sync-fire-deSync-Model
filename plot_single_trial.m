function [ ] = plot_single_trial( filename )
% This will make Figure 4 of the paper, analysing a single trial of a
% single pattern. 

load([filename '/Data/parameters' ]);
load([filename '/Processed Data/calcium_data.mat']);
load([filename '/Processed Data/freq_data.mat']);
par.sim_length = par.pre_stim_L * 2 + par.NC_stim_L;
x_lim_W = [par.pre_stim_L - 750 par.sim_length-par.pre_stim_L+750];
T_x = ((min(x_lim_W)+1:max(x_lim_W))-par.pre_stim_L)/1000;

for b = 1:par.B
    if(par.B > 1); P_n = ['/P' int2str(b)]; else; P_n = []; end
    load([filename '/Processed Data/spike_data_P' int2str(b) '.mat'])
    for t = 1:par.trials
        fig = figure; set(fig, 'Position', [0 0 700 700]);
        for n = 1:length(par.sim_order_n) % loop through simulation stages
            %% scatter spikes and shade by region
            c = 0;
            for g=1 : length(par.nG)
                % set figure parameters
                if(g == 1) % TC
                    sp_x = [13 15 17 19]; 
                    rast_lim = [0 sum(par.n_TC) + 0.5];
                elseif(g == 2) % BP
                    sp_x = [7 9 11]; 
                    rast_lim = [sum(par.n_TC) + 0.5 sum(par.n_TC) + par.n_BP + 0.5];
                elseif(g == 3) % NC
                    sp_x = [3 5]; 
                    rast_lim = [par.network_size - par.n_NC + 0.5 par.network_size + 0.5];
                end
                if(n==1) % encoding
                    subplot(15, 2, sp_x)
                elseif(n==2) % retrieval
                    ax = subplot(15, 2, sp_x + 1); 
                    set(ax, 'yaxislocation', 'right'); ylabel('neuron ID');
                end
                ytickangle(90);
                hold all; box on;
                
                for L = 1:par.(['n_SG_' par.nG{g}]) % loop through layers
                    for g2 = 1:length(par.([par.nG{g} '_N'])) % loop through neuron type
                        if(contains(par.([par.nG{g} '_N']){g2},'E')==1 || contains(par.([par.nG{g} '_N']){g2},'P')==1); B_c=1; else; B_c = 0; end
                        if(contains(par.([par.nG{g} '_N']){g2},'I')==1 || contains(par.([par.nG{g} '_N']){g2},'O')==1); R_c=1; else; R_c = 0; end
                        if(contains(par.([par.nG{g} '_N']){g2},'O')==1 || contains(par.([par.nG{g} '_N']){g2},'P')==1); G_c=0; else; G_c = 0; end
                        patch([0 par.sim_length par.sim_length 0], ...
                            [c+0.5 c+0.5 c+par.(['n_' par.nG{g} '_' par.([par.nG{g} '_N']){g2}])(L)+0.5 c+...
                            par.(['n_' par.nG{g} '_' par.([par.nG{g} '_N']){g2}])(L)+0.5],...
                            [1 1 1], 'edgecolor','k'); alpha(0.1);
                        
                        for i=1:length(n_G_ID{g}{L}{g2}) % loop through clusters
                            if(isempty(SPIKES{t}{n}{g}{L}{g2}{i})~=1)
                                scatter(SPIKES{t}{n}{g}{L}{g2}{i}(:,2),SPIKES{t}{n}{g}{L}{g2}{i}(:,1),'Marker', '.', ...
                                    'MarkerEdgeColor', [R_c 0 B_c]); 
                            end
                            ylim(rast_lim); ylabel('neuron ID');
                            xlim(x_lim_W); xticklabels({}); 
                        end
                        c = c + par.(['n_' par.nG{g} '_' par.([par.nG{g} '_N']){g2}])(L);
                    end
                end
            end
            %% shade scatter plot with STDP at encoding
            if(n == 1)
                BP_c = []; NC_c = []; TC_L1_c = []; TC_L2_c = [];
                for i = 1:length(n_G_ID{2}{1}{1})
                    BP_c = [BP_c; calcium_all.B1.DL{2,2}{1,1}{1,1}{i,i}(2,:)];
                    NC_c = [NC_c; calcium_all.B1.DL{2,3}{1,1}{1,1}{i,1}(2,:)];

                    for j = 1:length(n_G_ID{1}{1}{1})
                        TC_L1_c = [TC_L1_c; calcium_all.B1.DL{1,2}{1,1}{1,1}{j,i}(2,:)];
                    end
                    for j = 1:length(n_G_ID{1}{2}{1})
                        TC_L2_c = [TC_L2_c; calcium_all.B1.DL{1,2}{2,1}{1,1}{j,i}(2,:)];
                    end
                end
                BP_c = max(BP_c); NC_c = max(NC_c); TC_L1_c = max(TC_L1_c); TC_L2_c = max(TC_L2_c);

                % BP CALCIUM OVER LTD THRESHOLD
                BP_t = T_f(BP_c > par.BP_BP_T_d);
                BP_t_d1 = [0 (BP_t(2:end) - BP_t(1:end-1))];
                BP_t_d2 = [(BP_t(2:end) - BP_t(1:end-1)) 0];
                BP_t = sort([BP_t(1) BP_t(BP_t_d1 > 10) BP_t(BP_t_d2 > 10) BP_t(end)]);

                % NC CALCIUM OVER LTP THRESHOLD
                NC_t = T_f(NC_c > par.BP_NC_T_p);
                NC_t_d1 = [0 (NC_t(2:end) - NC_t(1:end-1))];
                NC_t_d2 = [(NC_t(2:end) - NC_t(1:end-1)) 0];
                NC_t = sort([NC_t(1) NC_t(NC_t_d1 > 10) NC_t(NC_t_d2 > 10) NC_t(end)]);

                % TC L1 CALCIUM OVER LTP THRESHOLD
                TC_L1_t = T_f(TC_L1_c > par.TC_BP_T_p);
                TC_L1_t_d1 = [0 (TC_L1_t(2:end) - TC_L1_t(1:end-1))];
                TC_L1_t_d2 = [(TC_L1_t(2:end) - TC_L1_t(1:end-1)) 0];
                TC_L1_t = sort([TC_L1_t(1) TC_L1_t(TC_L1_t_d1 > 10) TC_L1_t(TC_L1_t_d2 > 10) TC_L1_t(end)]);

                % TC L2 CALCIUM OVER LTP THRESHOLD
                TC_L2_t = T_f(TC_L2_c > par.TC_BP_T_p);
                TC_L2_t_d1 = [0 (TC_L2_t(2:end) - TC_L2_t(1:end-1))];
                TC_L2_t_d2 = [(TC_L2_t(2:end) - TC_L2_t(1:end-1)) 0];
                TC_L2_t = sort([TC_L2_t(1) TC_L2_t(TC_L2_t_d1 > 10) TC_L2_t(TC_L2_t_d2 > 10) TC_L2_t(end)]);

                for i = 1:length(n_G_ID{2}{1}{1})
                    % BP CALCIUM > THRESHOLD
                    subplot(15,2,[7 9 11]);
                    x = [BP_t(1+(i-1)*2) BP_t(2+(i-1)*2)];
                    y1 = sum(par.n_TC)+1; y2 = y1 + par.n_BP_E - 1;
                    fill([x fliplr(x)], [y1 y1 y2 y2], [1 0.4 0], 'edgecolor', 'none'); alpha(0.3);
                    % NC CALCIUM > THRESHOLD
                    subplot(15,2,[3 5]);
                    x = [NC_t(1+(i-1)*2) NC_t(2+(i-1)*2)];
                    y1 = sum(par.n_TC)+par.n_BP+1; y2 = y1 + par.n_NC_E - 1;
                    fill([x fliplr(x)], [y1 y1 y2 y2], [1 .4 0], 'edgecolor', 'none'); alpha(0.3);
                    % TC L1 CALCIUM > THRESHOLD
                    subplot(15,2,[13 15 17 19]);
                    x = [TC_L1_t(1+(i-1)*2) TC_L1_t(2+(i-1)*2)];
                    y1 = 1; y2 = y1 + par.n_TC_E(1) - 1;
                    fill([x fliplr(x)], [y1 y1 y2 y2], [1 .4 0], 'edgecolor', 'none'); alpha(0.3);
                    % TC L2 CALCIUM > THRESHOLD
                    x = [TC_L2_t(1+(i-1)*2) TC_L2_t(2+(i-1)*2)];
                    y1 = par.n_TC(1)+1; y2 = y1 + par.n_TC_E(2) - 1;
                    fill([x fliplr(x)], [y1 y1 y2 y2], [1 .4 0], 'edgecolor', 'none'); alpha(0.3);
                end
            end
            %% show stimulus event pattern
            if(n==1) % encoding
                subplot(15, 2, 1)
            elseif(n==2) % retrieval
                subplot(15, 2, 2)
            end
            if(isempty(find(cellfun(@isempty,strfind(par.nG,'NC'))==0)) ~= 1)
                for i = 1:length(par.NC_stims_t_all{b}{t}{1})
                   hold all; box on;
                   x_t = par.NC_stims_t_all{b}{t}{1}(i)/1000;
                   plot(ones(1,2)*x_t, [0 1],'linewidth',5,'color','k'); 
                end
            end
            ylim([0 1]); yticks([]); xlim([min(T_x) max(T_x)]); 
            if(n==2); xticklabels([]); axis off; end
            
            %% plot phase patterns
            if(n == 1) % encoding
                subplot(15, 2, 21); hold on; axis off;
                LFP_DL = LFP_all.DL.NC.(['B' int2str(b)]);
                plot( (T-par.pre_stim_L)/1000, LFP_DL, 'k', 'linewidth', 1); xlim([min(T_x) max(T_x)]);

                subplot(15, 2, 23); hold on; axis off;
                P_DL = PHASE.DL.NC.(['B' int2str(b)]).(['T' int2str(t)]); 
                plot( (T-par.pre_stim_L)/1000, P_DL, 'k', 'linewidth', 1); xlim([min(T_x) max(T_x)]);
            else % recall
                subplot(15, 2, 22); hold on; axis off;
                LFP_DL = LFP_all.AL.NC.(['B' int2str(b)]);
                plot( (T-par.pre_stim_L)/1000, LFP_DL, 'k', 'linewidth', 1); xlim([min(T_x) max(T_x)]);
                
                x_f = 200; % filter similarity data
                subplot(15, 2, 24); hold on; axis off;
                P_AL = PHASE.AL.NC.(['B' int2str(b)]).(['T' int2str(t)]); 
                plot( (T-par.pre_stim_L)/1000, P_AL, 'k', 'linewidth', 1); xlim([min(T_x) max(T_x)]);
                                
                stim_L = par.pre_stim_L+min(par.NC_stims_t_all{b}{t}{1}):par.pre_stim_L+max(par.NC_stims_t_all{b}{t}{1})+x_f;
                P_DL = PHASE.DL.NC.(['B' int2str(b)]).(['T' int2str(t)])(stim_L);
                x1 = (min(stim_L)-par.pre_stim_L)/1000; x2 = x1 + length(stim_L)/1000; 
                y1 = min(P_DL); y2 = max(P_DL);
                x = [x1, x2, x2, x1, x1]; y = [y1, y1, y2, y2, y1]; 
                subplot(15, 2, 23); plot(x,y,'r:','linewidth',1);
                xy_d = phase_compare(P_AL, P_DL); xy_d = filter(ones(x_f,1)/x_f, 1, xy_d);
                
                subplot(15, 2, 26); hold on; axis off; plot(x,y,'r:','linewidth',1);
                T_t = (T-par.pre_stim_L)/1000; 
                x_t1 = find(T_t > x1, 1, 'first'); x_t2 = find(T_t >= x2, 1, 'first');
                plot( T_t(x_t1:x_t2), P_DL, 'k', 'linewidth', 1); xlim([min(T_x) max(T_x)]);
                
                ax = subplot(15, 2, [28 30]); hold on; box on; ylabel('similarity'); xlabel('time (s)');
                T_xi = T_x*1000+par.pre_stim_L;
                xy_t = ((length(P_DL)/2:par.sim_length-length(P_DL)/2)-par.pre_stim_L)/1000;
                %xy_t = linspace(min(T_x), max(T_x), length(xy_d));
                plot(xy_t, xy_d, 'k', 'linewidth', 2); plot(zeros(1,length(-4:0.1:4)), -4:0.1:4, 'k-');
                ylim([floor(min(xy_d)*10)/10 ceil(max(xy_d)*10)/10]); yticks(0:0.25:0.5);
                xlim([min(T_x) max(T_x)]); 
            end
        end
        saveas(fig,[filename '/Graphs/single_trial' P_n '_T' int2str(t) '.jpg']); close(fig);
    end
end



end

