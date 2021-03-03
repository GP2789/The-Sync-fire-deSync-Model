function [ sim_stats, sim_time, time ] = simulate_network( filename, timer, ...
    par, I, weight_matrix, weight_matrix_STDP, tau_syn, wm_max, delay)
%SIMULATE_NETWORK Summary of this function goes here
%   Detailed explanation goes here
tic;
%% DECLARATIONS & INITIALISATIONS
%global par; global V_m; global I;
%global weight_matrix; %global ADP;
%global tau_syn; global calcium; global weight_matrix_STDP; global p_w;
%global delay; global wm_max; global p_STDP; 

p_w = (weight_matrix./wm_max).*weight_matrix_STDP; % assign & normalise weights 0-1(max)
p_w(isnan(p_w)==1) = 0;

% IMPORTANT : only saves synapses that have STDP enabled, as all others
% will remain constant. This vastly cuts down on data storage. 
p_STDP = find(weight_matrix_STDP==1); % ^^
% storage variables are n synapses with STDP by time
calcium = zeros(length(p_STDP), par.sim_length); 

if(strcmp('on', timer)==1)
    sim_time.total_time_s = cputime;
    t = zeros(1, par.sim_length/par.d_T); x = zeros(1, par.sim_length); 
    sim_time.volts_update = t; sim_time.adnew_spikes = x; sim_time.STDP_weights = x;
    vlt.ch = t; vlt.vt = t; vlt.gt = t; vlt.sp = t;
else
    sim_time = 0;
end
spikes_t = [];
%% SET NEURON PROPERTIES BASED ON PARAMETERS
%h1 = waitbar(0, 'Simulating Network', 'Units', 'normalized', 'Position', [0.5 0.1 0.2 0.1]);
last_spikes = zeros(par.network_size, 1);
V_m = zeros(par.network_size,1);
I.L = zeros(par.network_size, par.sim_length);
if(strcmp(par.model_type,'HH')==1)
    PYR = find(cellfun(@isempty,strfind(par.n_ID,'OLM'))~=0); % find pyramidal IDs
    OLM = find(cellfun(@isempty,strfind(par.n_ID,'OLM'))==0); % find OLM IDs
    G.m = zeros(par.network_size,1); G.n = zeros(par.network_size,1); G.h = zeros(par.network_size,1);
    [G_t] = HH_UG( V_m(PYR), 'P', 'x_SS', [] ); G.m(PYR) = G_t.m; G.n(PYR) = G_t.n; G.h(PYR) = G_t.h;
    if(isempty(OLM)~=1)
        G.p = zeros(par.network_size,1); G.h_f = zeros(par.network_size,1); G.h_s = zeros(par.network_size,1); 
        [G_t] = HH_UG( V_m(OLM), 'O-LM', 'x_SS', [] ); 
        G.m(OLM) = G_t.m; G.n(OLM) = G_t.n; G.h(OLM) = G_t.h; G.p(OLM) = G_t.p; G.h_f(OLM) = G_t.h_f; G.h_s(OLM) = G_t.h_s;
    end
    V_m(PYR) = par.PYR_E_L +1; V_m(OLM) = 50;
    E_Na = zeros(par.network_size,1); E_K = zeros(par.network_size,1); E_L = zeros(par.network_size,1);
    g_Na = zeros(par.network_size,1); g_K = zeros(par.network_size,1); g_L = zeros(par.network_size,1);
    E_Na(PYR) = par.PYR_E_Na; E_K(PYR) = par.PYR_E_K; E_L(PYR) = par.PYR_E_L;
    g_Na(PYR) = par.PYR_g_bar_Na; g_K(PYR) = par.PYR_g_bar_K; g_L(PYR) = par.PYR_g_bar_L;
    if(isempty(OLM)~=1)
        E_h = zeros(par.network_size,1); g_h = zeros(par.network_size,1); g_p = zeros(par.network_size,1);
        E_Na(OLM) = par.OLM_E_Na; E_K(OLM) = par.OLM_E_K; E_L(OLM) = par.OLM_E_L;
        g_Na(OLM) = par.OLM_g_bar_Na; g_K(OLM) = par.OLM_g_bar_K; g_L(OLM) = par.OLM_g_bar_L;
        E_h(OLM) = par.OLM_E_h; g_h(OLM) = par.OLM_g_bar_h; g_p(OLM) = par.OLM_g_bar_p;
    end
    I.Na = zeros(par.network_size, par.sim_length);
    I.K = zeros(par.network_size, par.sim_length);
    if(isempty(OLM)~=1)
       I.p = zeros(par.network_size, par.sim_length);
       I.h = zeros(par.network_size, par.sim_length);
    end
elseif(strcmp(par.model_type,'IAF')==1)
    V_m = par.E; RF = zeros(par.network_size,1); 
end

% TRACKING OUTPUTS THROUGH TIME
V_m_t = zeros(par.network_size, par.sim_length);
p_w_t = zeros(length(p_STDP), par.sim_length);
MUA_t = zeros(par.network_size, par.sim_length);
ADP_t = zeros(par.network_size, par.sim_length);
SYN_t = zeros(par.network_size, par.sim_length, length(par.uTS));
MUA = zeros(par.network_size, 1/par.d_T);

%% MAIN SIMULATION
for t=1:par.sim_length/par.d_T
    %% VOLTAGE CHANGE
    if(strcmp('on', timer)==1); sim_time.volts_update(t) = cputime; vlt.ch(t) = cputime; end
    if(strcmp(par.model_type,'HH')==1)
        I_Na = (g_Na .* G.m.^3 .* G.h .* (V_m - E_Na)) * par.d_T; % Sodium current
        I_K = (g_K .* G.n.^4 .* (V_m - E_K)) * par.d_T; % Potassium current
        I_L = (g_L .* (V_m - E_L)) * par.d_T; % Leak current
        I_sum = I_Na + I_K + I_L;
        if(isempty(OLM)~=1)
            I_p = (g_p .* G.p .* (V_m - E_Na)) * par.d_T;
            I_h = (g_h .* (0.65 * G.h_f + 0.35 * G.h_s) .* (V_m - E_h)) * par.d_T;
            I_sum = I_sum + I_p + I_h;
        end
        % ADD SYNAPTIC INPUT EVERY 1 MS
        if(t*par.d_T - floor(t*par.d_T) == 0)
            I_SYN = sum(I.SYN(:,round(t*par.d_T),:),3); 
        else; I_SYN = zeros(par.network_size, 1); 
        end
        
        dv_dt = (-I_sum + I_SYN + I.ADP(:,t)) / par.c_m; 
        V_m_prev = V_m; 
    elseif(strcmp(par.model_type,'IAF')==1)
        I_L = (V_m - par.E)/par.tau_m; I_sum = I_L;
        dv_dt = ((-I_sum + sum(I.SYN(:,t,:),3) + I.ADP(:,t)) / par.c_m) .* double(RF==0); 
        RF(RF>0) = RF(RF>0)-1; % update refractory periods
    end
    
    V_m = V_m + dv_dt; %clear('dv_dt');
    if(strcmp('on', timer)==1); vlt.ch(t) = cputime - vlt.ch(t); vlt.vt(t) = cputime; end 
    
    %% UPDATE GATES
    if(strcmp('on', timer)==1); vlt.vt(t) = cputime - vlt.vt(t); vlt.gt(t) = cputime; end
    if(strcmp(par.model_type,'HH')==1)
        [dx_dt] = HH_UG( V_m, 'P', 'dx_dt', G);
        G.m(PYR) = G.m(PYR) + dx_dt.m(PYR) * par.d_T;
        G.n(PYR) = G.n(PYR) + dx_dt.n(PYR) * par.d_T; 
        G.h(PYR) = G.h(PYR) + dx_dt.h(PYR) * par.d_T;
        if(isempty(OLM)~=1)
            [dx_dt] = HH_UG( V_m, 'O-LM', 'dx_dt', G);
            G.m(OLM) = G.m(OLM) + dx_dt.m(OLM) * par.d_T;
            G.n(OLM) = G.n(OLM) + dx_dt.n(OLM) * par.d_T; 
            G.h(OLM) = G.h(OLM) + dx_dt.h(OLM) * par.d_T;
            G.p(OLM) = G.p(OLM) + dx_dt.p(OLM) * par.d_T;
            G.h_f(OLM) = G.h_f(OLM) + dx_dt.h_f(OLM) * par.d_T;
            G.h_s(OLM) = G.h_s(OLM) + dx_dt.h_s(OLM) * par.d_T;
        end
        G.m(G.m<0)=0; G.n(G.n<0)=0; G.h(G.h<0)=0; G.m(G.m>1)=1; G.n(G.n>1)=1; G.h(G.h>1)=1; % constrain probabilities to 0><1
        if(isempty(OLM)~=1); G.p(G.p<0)=0; G.h_f(G.h_f<0)=0; G.h_s(G.h_s<0)=0; G.p(G.p>1)=1; G.h_f(G.h_f>1)=1; G.h_s(G.h_s>1)=1; end
    end
    if(strcmp('on', timer)==1); vlt.gt(t) = cputime - vlt.gt(t); vlt.sp(t) = cputime; end

    if(strcmp(par.model_type,'HH')==1)
        % find spike events and add to MUA
        x = ([-10 < V_m -10 < V_m_prev] - [0 1]); %clear('V_m_prev'); % find all V that are above V_th
        x = (x(:,1) - x(:,2) .* mean(x,2)) == 0; % spike event at first V above V_th
    elseif(strcmp(par.model_type,'IAF')==1)
        x = zeros(par.network_size, 1); s = find(V_m >= par.V_th);
        if(isempty(s)~=1)
            x(s) = 1; V_m(s) = par.E; RF(s) = par.RF;
        end
    end
    MUA(:, round((t*par.d_T - floor(t*par.d_T))*(1/par.d_T))+1) = x;
    if(strcmp('on', timer)==1); vlt.sp(t) = cputime - vlt.sp(t); sim_time.volts_update(t) = cputime-sim_time.volts_update(t); end
    
    %% SUM SPIKES OVER 1 MS PERIODS
    spikes = find(MUA(:, round((t*par.d_T - floor(t*par.d_T))*(1/par.d_T))+1)>0); % get spikes
    spikes_t = [spikes_t; spikes];
    
    if(t*par.d_T - floor(t*par.d_T) == 0)
        x = round(t*par.d_T);
        %% ADD SPIKES TO FUTURE SYNAPTIC INPUT
        if(strcmp('on', timer)==1); sim_time.adnew_spikes(x) = cputime; end
        if(length(spikes_t)>=1)
            wm = weight_matrix(spikes_t,:); % get connected neurons
            TS = tau_syn(spikes_t,:,:); % get tau_S of connected neurons
            t_d = delay(spikes_t,:) + x; t_u = unique(t_d); t_u = t_u(t_u > x); 
            for k = 1:length(t_u)
                w = wm; w(t_d~=t_u(k))=0;
                s = permute(sum(w.*TS,1), [2 1 3]);
                p = I.PSP(1,1:min(par.PSP_L_mx, par.sim_length - t_u(k)+1),:);
                [c, i] = find(s~=0); i_u = unique(i);
                I.SYN(c, t_u(k) : min(par.sim_length, t_u(k) + par.PSP_L_mx-1), i_u) = ...
                    I.SYN(c, t_u(k) : min(par.sim_length, t_u(k) + par.PSP_L_mx-1), i_u) + ...
                    s(c, :, i_u) .* p(:, :, i_u);
            end
            clear('wm'); clear('TS'); clear('t_d'); clear('t_u'); 
            clear('s'); clear('p'); clear('c'); clear('i');
        end
        if(strcmp('on', timer)==1); sim_time.adnew_spikes(x) = cputime-sim_time.adnew_spikes(x); end
    
        %% STDP LEARNING
        if(strcmp('on', timer)==1); sim_time.STDP_weights(x) = cputime; end
        [ weight_matrix, p_w, calcium ] = STDP( spikes_t, x, last_spikes, weight_matrix, wm_max, weight_matrix_STDP, ...
            p_w, p_STDP, par, calcium ); 
        last_spikes( spikes_t ) = x; spikes_t = []; clear('spikes');
        if(strcmp('on', timer)==1); sim_time.STDP_weights(x) = cputime-sim_time.STDP_weights(x); sim_time.make_simdata(x) = cputime; end

        %% DATA SAMPLING
        V_m_t(:, x) = V_m; 
        I.L(:, x) = I_L;
        if(strcmp(par.model_type,'HH')==1)
            I.Na(:, x) = I_Na;
            I.K(:, x) = I_K;
            if(isempty(OLM)~=1); I.p(:, x) = I_p; I.h(:, x) = I_h; end
        end
        SYN_t(:, x, :) = I.SYN(:, x, :);
        ADP_t(:, x) = I.ADP(:, t);
        p_w_t(:, x) = p_w(p_STDP);
        MUA_t(:, x) = sum(MUA, 2);
        MUA = zeros(par.network_size, 1/par.d_T);
        if(toc>10000); stop = 0; sim_stats = []; break; else; stop = 0; end
        if(strcmp('on', timer)==1); sim_time.make_simdata(x) = cputime - sim_time.make_simdata(x); end
    end
    
%     if(rem(t, 10/par.d_T) == 0)
%         fprintf('\b\b\b\b\b%3.0f%%\n', t/(par.sim_length/par.d_T)*100);
%         %waitbar(t/(par.sim_length/par.d_T), h1); 
%     end
end
%close(h1)
if(stop ~=1)
%% SAVE SIMULATION DATA
if(strcmp('on', timer)==1); sim_time.make_simdata(t+1) = cputime; end
I.SYN = SYN_t; clear('SYN_t'); 
I.ADP = ADP_t; clear('ADP_t'); 
sim_stats.calcium = calcium; clear('calcium');
sim_stats.p_w = p_w_t; clear('p_w_t');

[s, d] = find(MUA_t>0); clear('MUA'); clear('MUA_t'); 
sim_stats.spike_detector = [s d]; clear('s'); clear('d');
sim_stats.weight_matrix = weight_matrix; clear('weight_matrix');
sim_stats.WM = p_STDP; clear('p_STDP');
sim_stats.TS = tau_syn; clear('tau_syn');
sim_stats.WM_p = p_w; clear('p_w'); 
c = 0;

% sim_stats.V_m = V_m_t;
% sim_stats.I.L = I.L;
% sim_stats.I.SYN = I.SYN;
% if(strcmp(par.model_type,'HH')==1); sim_stats.I.Na = I.Na; sim_stats.I.K = I.K; 
% if(isempty(OLM)~=1); sim_stats.I.h = I.h; sim_stats.I.p = I.p; end; end
clear('V_m'); clear('V_m_t'); clear('I');

%% TIMERS
if(strcmp('on', timer)==1) 
    sim_time.make_simdata(t+1) = cputime - sim_time.make_simdata(t+1); 
    f(1) = figure; set(f(1), 'Position', [0 0 1200 1200]); 
    T = 0:par.d_T:par.sim_length-par.d_T; X = 0:1:par.sim_length-1; 
    mx = max([sim_time.volts_update sim_time.adnew_spikes sim_time.STDP_weights]);
   
    subplot(3,1,1);hold on;
    plot(T, sim_time.volts_update, 'color', [0 1 0]); 
    plot(T, vlt.ch, 'color', [0 0.8 0]); plot(T, vlt.vt, 'color', [0 0.6 0]);
    plot(T, vlt.sp, 'color', [0 0.2 0]); plot(T, vlt.gt, 'color', [0 0.4 0]);
    x = [vlt.ch; vlt.vt; vlt.gt; vlt.sp];
    ylim([0 mx*1.2]); title('Volt Updates'); 
    
    subplot(3,1,2); hold on;
    plot(X, sim_time.adnew_spikes, 'color', [0 0 1]); 
    ylim([0 mx*1.2]); title('Add New Spikes'); 
    
    subplot(3,1,3); hold on;
    plot(X, sim_time.STDP_weights, 'color', [1 0 0]);
    ylim([0 mx*1.2]); title('STDP Updates'); saveas(f, [filename '/timers.jpg']);
    
    sim_time.make_simdata = sum(sim_time.make_simdata);
    sim_time.volts_update = sum(sim_time.volts_update); sim_time.adnew_spikes = sum(sim_time.adnew_spikes);
    sim_time.STDP_weights = sum(sim_time.STDP_weights);
    
    x = struct2cell(sim_time); 
    sim_time.total_time_s = cputime - sim_time.total_time_s;
    sim_time.remaining_fs = sim_time.total_time_s - sum([x{2:end}]);
    vlt.ch = sum(vlt.ch); vlt.vt = sum(vlt.vt);  vlt.sp = sum(vlt.sp); vlt.gt = sum(vlt.gt);
    sim_time
end

end
time = toc;
%fprintf('\b\b\b\b\bcompleted in %1.0fm %1.0fs\n', floor(time/60), time - floor(time/60)*60)
end

