function [  ] = content_specificity( filename )
% This will analysis similarity across and between patterns to reproduce
% Figure 5 of the paper.

load([filename '/Data/parameters' ]);

if(par.B > 1) 
    %% PHASE ANALYSIS OVER/BETWEEN PATTERNS
    if(exist([filename '/Processed Data/similarity_data.mat'],'file') ~= 2)
        % LOAD FILES & INITIALISE
        load([filename '/Processed Data/freq_data.mat']);
        par.sim_length = par.pre_stim_L * 2 + par.NC_stim_L;
        POW.DL = zeros(size(POW_all.DL.NC.B1)); POW.AL = zeros(size(POW_all.AL.NC.B1)); 
        sim_to_diff_all = []; sim_to_same_all = [];
        h2 = waitbar(0, 'Phase Compare', 'Units', 'normalized', 'Position', [0.5 0.55 0.2 0.1]);
        
        for b = 1:par.B % loop over patterns
            POW.DL = POW.DL + POW_all.DL.NC.(['B' int2str(b)]);
            POW.AL = POW.AL + POW_all.AL.NC.(['B' int2str(b)]);

            stim_L =  par.pre_stim_L: par.pre_stim_L + max(par.NC_stims_t_all{b}{1}{1}) + 500;

            X_L = round(length(stim_L)/2):round(par.sim_length-length(stim_L)/2);
            if(par.B > 1); B_n = ['/B' int2str(b)]; else; B_n = []; end
            for L = 1:par.n_SG_NC
                if(par.n_SG_NC <= 1); nL = ''; else; nL = ['_L' int2str(L)]; end
                sim_to_same = zeros(par.trials, par.sim_length);
                sim_to_diff = zeros(par.trials, par.sim_length);

                % CALCULATE RSA SIMILARITY BETWEEN PATTERNS
                for t2 = 1:par.trials % loop through all trials
                   phase_DL = PHASE.DL.(['NC' nL]).(['B' int2str(b)]).(['T' int2str(t2)])(stim_L);
                   % TO SAME PATTERNS (TRIAL VARIATION)
                   to_same_t = zeros(par.trials, par.sim_length);
                   for t3 = 1:par.trials % loop through all other trials
                       phase_AL = PHASE.AL.(['NC' nL]).(['B' int2str(b)]).(['T' int2str(t3)]);
                       to_same_t(t3, X_L) = phase_compare( phase_AL, phase_DL );
                   end
                   sim_to_same(t2, :) = mean(to_same_t, 1);
                   
                   % TO DIFFERENT PATTERNS MEAN (PATTERN VARIATION)
                   to_diff_p = zeros(par.B-1, par.sim_length); 
                   if(par.B > 1) 
                       i = 1;
                       for b2 = 1:par.B % loop through all other patterns
                           if(b2 ~= b)
                               to_diff_t = zeros(par.trials, par.sim_length);
                               for t3 = 1:par.trials % loop through all trials
                                   phase_AL_B = PHASE.AL.(['NC' nL]).(['B' int2str(b2)]).(['T' int2str(t3)]);
                                   to_diff_t(t3, X_L) = phase_compare( phase_AL_B, phase_DL );
                               end
                               to_diff_p(i, :) = mean(to_diff_t, 1);
                               i = i + 1;
                           end
                       end
                   end
                   sim_to_diff(t2, :) = mean(to_diff_p, 1);
                end
                % RECORD VARIABLES
                if(b == 1); content_spec = zeros(par.B, size(sim_to_same, 2)); end
                content_spec(b, :) = mean(sim_to_same,1) - mean(sim_to_diff, 1); % content specificity
                sim_to_same_all = [sim_to_same_all; sim_to_same]; % concatenate similarity to same pattern
                sim_to_diff_all = [sim_to_diff_all; sim_to_diff]; % concatenate similarity to other patterns
            end
            waitbar(b / par.B, h2);
        end
        close(h2);
        % SAVE VARIABLES
        save([filename '/Processed Data/similarity_data.mat'], 'freq', 'sim_to_same_all', 'sim_to_diff_all', 'content_spec', 'POW')
    else; load([filename '/Processed Data/similarity_data.mat'])
    end
    
    %% BOOTSTRAP AND FILTER DATA
    [to_same_L, to_same_U, to_same_M] = bootstrap( sim_to_same_all, 1000, []); 
    [to_diff_L, to_diff_U, to_diff_M] = bootstrap( sim_to_diff_all, 1000, []);
    to_same_L = filter(ones(1,50)/50, 1, to_same_L); to_same_M = filter(ones(1,50)/50, 1, to_same_M); to_same_U = filter(ones(1,50)/50, 1, to_same_U);
    to_diff_L = filter(ones(1,50)/50, 1, to_diff_L); to_diff_M = filter(ones(1,50)/50, 1, to_diff_M); to_diff_U = filter(ones(1,50)/50, 1, to_diff_U);
    
    stim_L2 =  par.pre_stim_L-1000: par.pre_stim_L + par.NC_stim_L + 1000;
    POW.DL = POW.DL(:,stim_L2) / par.B; POW.AL = POW.AL(:,stim_L2) / par.B;
    T2 = ((min(stim_L2):max(stim_L2))-par.pre_stim_L)/1000; 
    pow_lim = [1 15]; 
    pow_line = pow_lim(1):0.01:pow_lim(2);
    
    %% PLOT DATA
    fig = figure(1); set(fig, 'Position', [0 0 600 600]);
    for n = 1 : 2
        %% POWER DIFFERENCES
        ax(n) = subplot(9, 10, [11:15 21:25 31:35 41:45] + (n-1)*5); 
        x = find(T2 >= 0, 1, 'first'); 
        y = [find(freq.NC <= pow_lim(1), 1, 'first') find(freq.NC >= pow_lim(2), 1, 'first')];
        bs_pow = POW.([par.sim_order_n{n}]); bs_pow = bs_pow - mean(bs_pow(:, 1:x), 2);
        imagesc(T2, freq.NC(y(1):y(2)), bs_pow(y(1):y(2), :)); set(gca,'YDir','normal'); 
        box on; ylim(pow_lim);
        hold on; plot(zeros(length(pow_line), 1), pow_line, 'k-');
        ylabel('frequency (Hz)'); box on; 
        if(n==1); cb1 = colorbar; xlabel('time (s)'); else; xticklabels({}); set(ax(n), 'yaxislocation','right'); end
        xlim([min(T2) max(T2)]);
        
        %% PHASE DIFFERENCE
        if(n==2)
            ax(6) = subplot(9, 10, [56:60 66:70]); hold on; box on;
            % PLOT SIMILARITY TO OTHER PATTERNS
            fill([T2 fliplr(T2)], [to_diff_L(stim_L2) fliplr(to_diff_U(stim_L2))], 'k', 'edgecolor', 'none'); alpha(0.3);
            h=[]; h(1) = plot(T2, to_diff_M(stim_L2), 'k', 'linewidth', 2);
            % PLOT SIMILARITY TO OTHER PATTERNS
            fill([T2 fliplr(T2)], [to_same_L(stim_L2) fliplr(to_same_U(stim_L2))], 'r', 'edgecolor', 'none'); alpha(0.3);
            h(2) = plot(T2, to_same_M(stim_L2), 'color', [0.8 0.2 0.2],'linewidth', 2);
            
            plot(T2, zeros(1, length(T2)), 'k:'); plot(zeros(length(0:0.01:1), 1), 0:0.01:1, 'k-'); 
            xlim([min(T2) max(T2)]); ylabel('similarity'); xticklabels({});
            ylim([min(min([to_same_L(stim_L2) to_diff_L(stim_L2)])) max(max([to_same_U(stim_L2) to_diff_U(stim_L2)]))])
            legend(h, 'to others', 'to same', 'location','northwest'); 
            set(ax(6), 'yaxislocation','right'); %title('similarity between patterns');
        end
        
        %% STIMULUS PATTERNS  
        if(par.B > 5)
            B_f = [1 0 round(par.B/2) 0 par.B]; B_x = 5;
        else
            B_x = par.B; B_f = 1:par.B;
        end
        for b = 1 : B_x
            subplot(9, 10, b + (n-1)*B_x); hold on;
            if(mod(b,2)==1 || par.B < 5)
                for i = 1:length(par.NC_stims_t_all{B_f(b)}{1}{1})
                   plot(ones(1,2)*((par.NC_stims_t_all{B_f(b)}{1}{1}(i)+par.pre_stim_L)-par.pre_stim_L)/1000, [0 1],'linewidth',3,'color','k'); 
                end
                xlim([min(T2)+0.5 max(T2)-0.5]); yticklabels([]); xticklabels([]);
                if(n == 1); box on; yticks(''); xticks(''); title(['P' int2str(B_f(b))]); else; axis off; end
            else
                scatter([0 1 2], [0 0 0], 'ko','filled'); axis off; xlim([-1 3]);
            end
        end
        
        %% ENCODING / RETRIEVAL COMPARISON
        if(n==2)
            ax(8) = subplot(9, 10, [76:80 86:90]); box on;
            plot(T2, filter(ones(500,1)/500, 1, mean(content_spec(:,stim_L2))), 'color' ,[0.7 0.3 0.3],'linewidth', 2); hold on;
            xlabel('time (s)'); 
            plot(T2, zeros(1, length(T2)), 'k:'); xlim([min(T2) max(T2)]);
            plot(ones(length(-1:0.01:1), 1)*0, -1:0.01:1, 'k-'); 
            ylim([min(mean(content_spec(:,stim_L2))) max(mean(content_spec(:,stim_L2)))]*1.25);
            set(ax(8), 'yaxislocation','right'); ylabel('content specificity');
        end
    end
    saveas(fig,[filename '/Graphs/content_specificity.jpg']); close(fig);
    
else; disp('ANALYSIS REQUIRES SIMULATION OF MULTIPLE PATTERNS')
end


end

