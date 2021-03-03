function [ weight_matrix, p_w, calcium ] = STDP( spikes, t, last_spikes, weight_matrix, wm_max, weight_matrix_STDP, ...
    p_w, p_STDP, par, calcium )
%STDP Summary of this function goes here
%   Detailed explanation goes here
%% GLOBALS
%global weight_matrix; global wm_max; global weight_matrix_STDP;
%global p_w; global p_STDP; global par; global calcium; 

%% CALCIUM ADDITION
if(isempty(spikes) ~= 1)
    for s = 1:length(spikes)
        % add pre-syn calcium * local inhibitory current
        [~, i] = find(weight_matrix_STDP(spikes(s),:)>0);
        i = reshape(i,[length(i) 1]);
        x = spikes(s) + (i-1)*par.network_size; % find linear index for pre-syn
        if(isempty(x) ~= 1)
            [~, b] = ismember(x, p_STDP); t_D = t + par.C_D;
            if(t_D <= par.sim_length)
                C_pre = par.C_pre(p_STDP(b)) * (1 - exp(-(t - last_spikes(spikes(s))) * 0.1));
                calcium(b, t_D) = calcium(b, t_D) + C_pre; 
            end
        end

        % add post-syn calcium * local inhibitory current
        [c, ~] = find(weight_matrix_STDP(:,spikes(s))>0);
        c = reshape(c,[length(c) 1]);
        x = c + (spikes(s)-1)*size(weight_matrix,1); % find linear index for post-syn
        [~, b] = ismember(x, p_STDP);
        if(isempty(b)~=1)
            C_post = par.C_post(p_STDP(b)) * (1 - exp(-(t - last_spikes(spikes(s))) * 0.1));
            calcium(b, t) = calcium(b, t) + C_post;
        end
    end
end
 %% LTP
x_p = find(calcium(:, t) > par.T_p(p_STDP));
if(isempty(x_p)~=1) % LTP
    p_w( p_STDP(x_p) ) = p_w( p_STDP(x_p) ) ...
        + ( (par.G_p(p_STDP(x_p)) .* ( 1-p_w(p_STDP(x_p)) ) .*  heaviside( calcium(x_p, t) - par.T_p(x_p) ) ) ...
        / par.T_WC ); 
end
%% LTD
x_d = find(calcium(:, t) > par.T_d(p_STDP) & calcium(:, t) < par.T_p(p_STDP));
if(isempty(x_d)~=1) % LTD
    p_w( p_STDP(x_d) ) = p_w( p_STDP(x_d) ) ...
        - ( (par.G_d(p_STDP(x_d)) .* p_w(p_STDP(x_d)) .* heaviside( calcium(x_d, t) - par.T_d(x_d) )) ...
        / par.T_WC); 
end

%% WEIGHT DECAY
p_w(p_STDP) = p_w(p_STDP) + ...
    (((-p_w(p_STDP) .* (1-p_w(p_STDP)) .* (par.attractor - p_w(p_STDP)))) ...
    / par.T_WC);

weight_matrix(p_STDP) = p_w(p_STDP) .* wm_max(p_STDP);

%% CALCIUM DECAY
if(t < par.sim_length)
    calcium(:, t+1) = calcium(:, t+1) + calcium(:, t) - (calcium(:, t)./(par.T_Ca(p_STDP)));  
end

end

