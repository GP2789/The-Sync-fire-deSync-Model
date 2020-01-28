function [ gate ] = HH_UG( V, N_type, C_type, gate_t )
%UPDATE_GATE Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(C_type, 'tau_x') == 1)  % TAU
    gate.m = 1 ./ (alpha_m(V, N_type) + beta_m(V, N_type));
    gate.n = 1 ./ (alpha_n(V, N_type) + beta_n(V, N_type));
    gate.h = 1 ./ (alpha_h(V, N_type) + beta_h(V, N_type));
    if(strcmp(N_type,'O-LM')==1) 
        gate.p = 1 ./ (alpha_p(V) + beta_p(V));
        gate.h_f = 0.51 ./ (exp((V-1.7)/10) + exp(-(V+340)/52)) + 1;
        gate.h_s = 5.6 ./ (exp((V-1.7)/14) + exp(-(V+260)/43)) + 1;
    end
elseif(strcmp(C_type, 'x_SS') == 1) % STEADY STATE
    gate.m = alpha_m(V, N_type) ./ (alpha_m(V, N_type) + beta_m(V, N_type));
    gate.n = alpha_n(V, N_type) ./ (alpha_n(V, N_type) + beta_n(V, N_type));
    gate.h = alpha_h(V, N_type) ./ (alpha_h(V, N_type) + beta_h(V, N_type));
    if(strcmp(N_type,'O-LM')==1)
        gate.p = alpha_p(V) ./ (alpha_p(V) + beta_p(V));
        gate.h_f = 1 ./ (1 + exp((V+79.2)/9.78));
        gate.h_s = 1 ./ (1 + exp((V+2.83)/15.9)).^58; 
    end
elseif(strcmp(C_type, 'dx_dt') == 1) % GATE UPDATE
    gate.m = (alpha_m(V, N_type) .* (1-gate_t.m) - beta_m(V, N_type) .* gate_t.m); % dm_dt
    gate.n = (alpha_n(V, N_type) .* (1-gate_t.n) - beta_n(V, N_type) .* gate_t.n); % dn_dt
    gate.h = (alpha_h(V, N_type) .* (1-gate_t.h) - beta_h(V, N_type) .* gate_t.h); % dh_dt
    if(strcmp(N_type,'O-LM')==1)
        gate.p = (alpha_p(V) .* (1-gate_t.p) - beta_p(V) .* gate_t.p); % dp_dt
        gate.h_f = ((1 ./ (1 + exp((V+79.2)/9.78))) - gate_t.h_f) ./ ...
            (0.51 ./ (exp((V-1.7)/10) + exp(-(V+340)/52)) + 1);
        gate.h_s = ((1 ./ (1 + exp((V+2.83)/15.9)).^58) - gate_t.h_s) ./ ...
            (5.6 ./ (exp((V-1.7)/14) + exp(-(V+260)/43)) + 1);
    end
end

% no NaN or inf
gate.m(isnan(gate.m)==1) = 0; gate.m(isinf(gate.m)==1) = 0;
gate.n(isnan(gate.n)==1) = 0; gate.n(isinf(gate.n)==1) = 0;
gate.h(isnan(gate.h)==1) = 0; gate.h(isinf(gate.h)==1) = 0;
if(strcmp(N_type,'O-LM')==1)
    gate.p(isnan(gate.p)==1) = 0; gate.p(isinf(gate.p)==1) = 0;
    gate.h_f(isnan(gate.h_f)==1) = 0; gate.h_f(isinf(gate.h_f)==1) = 0;
    gate.h_s(isnan(gate.h_s)==1) = 0; gate.h_s(isinf(gate.h_s)==1) = 0;
end

end

%% m GATE
function [ m ] = alpha_m( V, N_type ) 
    if(strcmp(N_type,'O-LM')==1); m = (-0.1*(V + 23)) ./ (exp(-0.1*(V+23)) - 1); 
    elseif(strcmp(N_type,'P')==1); m = (0.32*(V + 54)) ./ (1 - exp(-(V+54)/4)); end
end
function [ m ] = beta_m( V, N_type ) 
    if(strcmp(N_type,'O-LM')==1); m = 4 * exp(-(V+48)/18);
    elseif(strcmp(N_type,'P')==1); m = 0.28*(V+27) ./ (exp((V+27)/5) - 1); end
end
%% n GATE
function [ n ] = alpha_n( V, N_type ) 
    if(strcmp(N_type,'O-LM')==1); n = -0.01*(V + 27) ./ (exp(-0.1*(V+27)) - 1); 
    elseif(strcmp(N_type,'P')==1); n = 0.032*(V + 52) ./ (1 - exp(-(V+52)/5)); end
end
function [ n ] = beta_n( V, N_type ) 
    if(strcmp(N_type,'O-LM')==1); n = 0.125 .* exp(-(V+37)/80); 
    elseif(strcmp(N_type,'P')==1); n = 0.5 .* exp(-(57+V)/40); end
end
%% h GATE
function [ h ] = alpha_h( V, N_type ) 
    if(strcmp(N_type,'O-LM')==1); h = 0.07 .* exp(-(V+37)/20); 
    elseif(strcmp(N_type,'P')==1); h = 0.128 .* exp(-(50+V)/18); end
end
function [ h ] = beta_h( V, N_type ) 
    if(strcmp(N_type,'O-LM')==1); h = 1 ./ (exp(-0.1*(V+7)) + 1); 
    elseif(strcmp(N_type,'P')==1); h = 4 ./ (1 + exp(-(V+27)/5)); end
end
%% p GATE
function [ p ] = alpha_p( V ); p = 1 ./ (0.15*(1 + exp(-(V+38)/6.5))); end
function [ p ] = beta_p( V ); p = exp(-(V+38)/6.5) ./ (0.15*(1 + exp(-(V+38)/6.5))); end



