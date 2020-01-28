function [  ] = create_network( )
% Changes global weight matrix with connectivity topology defined by parameters.
make.totaltime_s = cputime;
make.declare_par = cputime;

%% declarations
global par; global I;
global delay; global ADP;
global weight_matrix; global weight_matrix_STDP; 
global tau_syn; global wm_max; 

% variable declarations
T = 0:1:par.sim_length;
par.n_ID = cell(par.network_size, 1);

% unique Input stream for each tau_syn
I.SYN = zeros(par.network_size, par.sim_length, length(par.uTS));
I.PSP = zeros(1, par.sim_length, length(par.uTS));
I.ADP = zeros(par.network_size, par.sim_length/par.d_T);
ADP.max = zeros(par.network_size,1);
par.PSP = cell(length(par.uTS),1); 
for i=1:length(par.uTS) % calculate EPSP for all tau_s in network
    g_syn = (exp(1)*T./par.uTS(i)).*exp(-T./par.uTS(i)).*heaviside(T); 
    g_syn = g_syn(g_syn>0.01);
    I.PSP(1,1:length(g_syn),i) = g_syn;
    par.PSP{i} = g_syn;
end
par.PSP_L_mx = length(g_syn);

weight_matrix = zeros(par.network_size, par.network_size);
wm_max = zeros(par.network_size, par.network_size);
weight_matrix_STDP = zeros(par.network_size, par.network_size);
tau_syn = zeros(par.network_size, par.network_size, length(par.uTS));
delay = zeros(par.network_size, par.network_size);
c = 0;

%% create groups of neurons
% make neuron IDs
par.layer_map = 0;
for g=1:length(par.nG)
    for L = 1:par.(['n_SG_' par.nG{g}])
        if(par.(['n_SG_' par.nG{g}])<=1); nL = ''; else; nL = ['_L' int2str(L)]; end
        for g2=1:length(par.([par.nG{g} '_N']))
            par.n_ID(c+1:c+par.(['n_' par.nG{g} '_' par.([par.nG{g} '_N']){g2}])(L)) = {[par.nG{g} '_' par.([par.nG{g} '_N']){g2} nL]};
            c = c + par.(['n_' par.nG{g} '_' par.([par.nG{g} '_N']){g2}])(L);
            par.layer_map = [par.layer_map; c];
        end
    end
end

%% initialise STDP parameters
par.C_pre = zeros(par.network_size, par.network_size); 
par.C_post = zeros(par.network_size, par.network_size);
par.T_Ca = zeros(par.network_size, par.network_size);
par.T_d = zeros(par.network_size, par.network_size); 
par.T_p = zeros(par.network_size, par.network_size); 
par.G_p = zeros(par.network_size, par.network_size);
par.G_d = zeros(par.network_size, par.network_size); 

make.declare_par = cputime - make.declare_par;

%% CONNECT WITHIN REGION NETWORK TOPOLOGY
% loop through regions
for g=1:length(par.nG)
    G1_E = find(cellfun(@isempty,strfind(par.n_ID,[par.nG{g} '_E']))==0);
    OLM = 0;
    %% add background noise to excitatory neurons in regions
    backg_noise_t = cputime;
    if(par.([par.nG{g} '_BG']) > 0) 
       spikes = spike_train(length(G1_E), par.([par.nG{g} '_BG']), par.sim_length, 1); % make spike train
       
       p = find(par.([par.nG{g} '_BG_TS'])==par.uTS);
       g_mx = find(par.PSP{p}==max(par.PSP{p}) == 1); g_L = length(par.PSP{p});
       [i, x, c] = find(spikes>0);
       x = unique(x);
       for k=1:length(x)
           I.SYN(G1_E, max(1,x(k)-g_mx+1):min(par.sim_length, x(k)+g_L-g_mx), p) = ...
               I.SYN(G1_E, max(1,x(k)-g_mx+1):min(par.sim_length, x(k)+g_L-g_mx), p) + ...
               par.PSP{p}(max(1,g_mx-x(k)+1):min(g_L, par.sim_length - x(k) + g_mx)) .*spikes(:,x(k)) * par.([par.nG{g} '_BG_W']);
       end
    end
    make.backg_noise(g) = cputime - backg_noise_t;

    %% loop through layers in each region
    make.con_neurons(g) = cputime;
    for i = 1:par.(['n_SG_' par.nG{g}])
        if(strcmp(par.nG{g}, 'TC') == 1)
            if(i == par.(['n_SG_' par.nG{g}]))
                L = 2;
            else; L = 1;
            end
        else; L = i;
        end
        
        if(par.(['n_SG_' par.nG{g}])<=1); nL = ''; else; nL = ['_L' int2str(i)]; end
        %% CONNECT LAYERS IN TC
        if(i > 1 && strcmp(par.nG{g}, 'TC') == 1) 
            % CONNECT ALL OF L-1 TO G1 OF L (FEED-FORWARD CONNECTION)
            L1_P = find(cellfun(@isempty,strfind(par.n_ID,[par.nG{g} '_P' ['_L' int2str(i-1)]]))==0);
            L2_E = min(find(cellfun(@isempty,strfind(par.n_ID,[par.nG{g} '_E' ['_L' int2str(i)]]))==0));
            connect_neurons(L1_P, L2_E, par.([par.nG{g} '_FF_R']), par.([par.nG{g} '_FF_TS']), par.([par.nG{g} '_FF_D']), 'off',...
                par.([par.nG{g} '_FF_W']), [], 0, [], 'normal');
        end
         %% CONNECT E <-> I for all groups 
        G1_I = find(strcmp(par.n_ID,[par.nG{g} '_I' nL])==1);
        G1_E = find(strcmp(par.n_ID,[par.nG{g} '_E' nL])==1);
        
        if(strcmp(par.nG{g}, 'NC') == 1); L = 1; end

        % STDP FOR E<->E
        if(strcmp(par.nG{g},'BP')==1)
            E_E_STDP.T_Ca = par.BP_BP_T_Ca; E_E_STDP.G_p = par.BP_BP_G_p; E_E_STDP.G_d = par.BP_BP_G_d;
            E_E_STDP.T_p = par.BP_BP_T_p; E_E_STDP.T_d = par.BP_BP_T_d;
            E_E_STDP.C_pre = par.BP_BP_C_pre; E_E_STDP.C_post = par.BP_BP_C_post;
        else; E_E_STDP = []; 
        end

        % CONNECT I -> I within each group
        connect_neurons(G1_I, G1_I, par.([par.nG{g} '_II_R'])(L), par.([par.nG{g} '_II_TS'])(L), par.([par.nG{g} '_II_D'])(L), 'off',...
            par.([par.nG{g} '_II_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_II_top']), par.([par.nG{g} '_II_cos'])(L));

        %% TOPOLOGY STRATEGY DEPENDS IF THERE ARE P (E_SS) OR O (OLM) NODES OR NOT
        % this is a bit confusing but it works and allows one to change
        % parameters without changing the code. Check weight matrix to verify.
        
        if(strcmp(par.([par.nG{g} '_E_SS']), 'no')==1)
            % CONNECT E -> I within each group
            connect_neurons(G1_E, G1_I, par.([par.nG{g} '_EI_R'])(L), par.([par.nG{g} '_EI_TS'])(L), par.([par.nG{g} '_EI_D'])(L), 'off',...
                par.([par.nG{g} '_EI_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_EI_top']), par.([par.nG{g} '_EI_cos'])(L));

            % CONNECT I -> E within each group
            connect_neurons(G1_I, G1_E,  par.([par.nG{g} '_IE_R'])(L), par.([par.nG{g} '_IE_TS'])(L), par.([par.nG{g} '_IE_D'])(L), 'off',...
                par.([par.nG{g} '_IE_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_IE_top']), par.([par.nG{g} '_IE_cos'])(L));

            % CONNECT E -> E within each group 
            connect_neurons(G1_E, G1_E, par.([par.nG{g} '_EE_R'])(L), par.([par.nG{g} '_EE_TS'])(L), par.([par.nG{g} '_EE_D'])(L), 'off',...
                par.([par.nG{g} '_EE_W'])(L), par.([par.nG{g} '_EE_W'])(L), par.([par.nG{g} '_W_sigma']), E_E_STDP, par.([par.nG{g} '_EE_top']), par.([par.nG{g} '_EE_cos'])(L));

        elseif(strcmp(par.([par.nG{g} '_E_SS']), 'yes')==1)
            G1_E_SS = find(strcmp(par.n_ID,[par.nG{g} '_P' nL])==1);
            % CONNECT E_SS -> E with wide NN gaussian to propogate wave
            connect_neurons(G1_E_SS, G1_E, par.([par.nG{g} '_PE_R'])(L), par.([par.nG{g} '_PE_TS'])(L), par.([par.nG{g} '_PE_D'])(L), 'off',...
                par.([par.nG{g} '_PE_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_PE_top']), par.([par.nG{g} '_PE_cos'])(L));
            K_n = round(par.(['n_' par.nG{g} '_E'])/par.(['n_' par.nG{g} '_P']));
            for K = 1:par.(['n_' par.nG{g} '_P'])(i)
                G1_EK = G1_E((K-1)*K_n+1:K*K_n); G1_E_SSK = G1_E_SS(K); G1_IK = G1_I(K);
                G2_EK = G1_E; G2_EK((K-1)*K_n+1:K*K_n) = [];
                %G2_IK = G1_I(G1_I~=G1_IK);
                % CONNECT E -> I within each group
                connect_neurons(G1_EK, G1_IK, par.([par.nG{g} '_EI_R'])(L), par.([par.nG{g} '_EI_TS'])(L), par.([par.nG{g} '_EI_D'])(L), 'off',...
                    par.([par.nG{g} '_EI_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_EI_top']), par.([par.nG{g} '_EI_cos'])(L));
                % CONNECT I -> E within each group
                connect_neurons(G1_IK, G2_EK,  par.([par.nG{g} '_IE_R'])(L), par.([par.nG{g} '_IE_TS'])(L), par.([par.nG{g} '_IE_D'])(L), 'off',...
                    par.([par.nG{g} '_IE_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_IE_top']), par.([par.nG{g} '_IE_cos'])(L));
                % CONNECT E -> E within each group
                connect_neurons(G1_EK, G1_EK, par.([par.nG{g} '_EE_R'])(L), par.([par.nG{g} '_EE_TS'])(L), par.([par.nG{g} '_EE_D'])(L), 'off',...
                    par.([par.nG{g} '_EE_W'])(L), par.([par.nG{g} '_EE_W'])(L), par.([par.nG{g} '_W_sigma']), E_E_STDP, par.([par.nG{g} '_EE_top']), par.([par.nG{g} '_EE_cos'])(L));
                % CONNECT E -> E_SS only to select group
                connect_neurons(G1_EK, G1_E_SSK, par.([par.nG{g} '_EP_R'])(L), par.([par.nG{g} '_EP_TS'])(L), par.([par.nG{g} '_EP_D'])(L), 'off',...
                    par.([par.nG{g} '_EP_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_EP_top']), par.([par.nG{g} '_EP_cos'])(L));
%                     connect_neurons(G1_E_SSK, G1_EK, par.([par.nG{g} '_PE_R'])(L), par.([par.nG{g} '_PE_TS'])(L), par.([par.nG{g} '_PE_D'])(L), 'off',...
%                         par.([par.nG{g} '_PE_W'])(L), [], 0, [], par.([par.nG{g} '_PE_top']), par.([par.nG{g} '_PE_cos'])(L));
                if(strcmp(par.([par.nG{g} '_OLM']), 'yes')==1)
                    G1_I_SS = find(strcmp(par.n_ID,[par.nG{g} '_O' nL])==1);
                    G1_I_SSK = G1_I_SS(K); OLM = 1;
                    % CONNECT E <-> O-LM 
                    connect_neurons(G1_EK, G1_I_SSK, par.([par.nG{g} '_EO_R'])(L), par.([par.nG{g} '_EO_TS'])(L), par.([par.nG{g} '_EO_D'])(L), 'off',...
                        par.([par.nG{g} '_EO_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_EO_top']), par.([par.nG{g} '_EO_cos'])(L));
                    connect_neurons(G1_I_SSK, G1_EK, par.([par.nG{g} '_OE_R'])(L), par.([par.nG{g} '_OE_TS'])(L), par.([par.nG{g} '_OE_D'])(L), 'off',...
                        par.([par.nG{g} '_OE_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_OE_top']), par.([par.nG{g} '_OE_cos'])(L));
                    % CONNECT E_SS <-> O-LM 
                    connect_neurons(G1_E_SSK, G1_I_SSK, par.([par.nG{g} '_PO_R'])(L), par.([par.nG{g} '_PO_TS'])(L), par.([par.nG{g} '_PO_D'])(L), 'off',...
                        par.([par.nG{g} '_PO_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_PO_top']), par.([par.nG{g} '_PO_cos'])(L));
                    connect_neurons(G1_I_SSK, G1_E_SSK, par.([par.nG{g} '_OP_R'])(L), par.([par.nG{g} '_OP_TS'])(L), par.([par.nG{g} '_OP_D'])(L), 'off',...
                        par.([par.nG{g} '_OP_W'])(L), [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_OP_top']), par.([par.nG{g} '_OP_cos'])(L));
                    if(i > 1) % FOR FEEDBACK INHIBITION P2 -> P1 ( -> O1 -> E1 )
                        if(K == par.(['n_' par.nG{g} '_P'])(i))
                            %G1_I_SS_2 = find(strcmp(par.n_ID,[par.nG{g} '_O_L' int2str(i-1)])==1);
                            for i2 = 1:i-1
                                G1_I_SS_2 = find(strcmp(par.n_ID,[par.nG{g} '_O_L' int2str(i2)])==1); w_i = par.(['n_SG_' par.nG{g}])-i2;
                                connect_neurons(G1_E_SSK, G1_I_SS_2, par.([par.nG{g} '_FB_R']), par.([par.nG{g} '_FB_TS']), par.([par.nG{g} '_FB_D']), 'off',...
                                    par.([par.nG{g} '_FB_W'])/w_i, [], par.([par.nG{g} '_W_sigma']), [], par.([par.nG{g} '_FB_top']), par.([par.nG{g} '_FB_cos']));
                            end
                        end
                    end
                end
            end
            % CONNECT I <-> E_SS within each group
            connect_neurons(G1_I, G1_E_SS, par.([par.nG{g} '_IP_R'])(L), par.([par.nG{g} '_IP_TS'])(L), par.([par.nG{g} '_IP_D'])(L), 'off',...
                par.([par.nG{g} '_IP_W'])(L), [], 0, [], par.([par.nG{g} '_IP_top']), par.([par.nG{g} '_IP_cos'])(L)); 
            connect_neurons(G1_E_SS, G1_I, par.([par.nG{g} '_PI_R'])(L), par.([par.nG{g} '_PI_TS'])(L), par.([par.nG{g} '_PI_D'])(L), 'off',...
                par.([par.nG{g} '_PI_W'])(L), [], 0, [], par.([par.nG{g} '_PI_top']), par.([par.nG{g} '_PI_cos'])(L)); 
            % CONNECT E_SS -> E_SS within each group
            connect_neurons(G1_E_SS, G1_E_SS, par.([par.nG{g} '_PP_R'])(L), par.([par.nG{g} '_PP_TS'])(L), par.([par.nG{g} '_PP_D'])(L), 'off',...
                par.([par.nG{g} '_PP_W'])(L), [], 0, [], par.([par.nG{g} '_PP_top']), par.([par.nG{g} '_PP_cos'])(L)); 
        end

        if(strcmp(par.([par.nG{g} '_OLM']), 'yes')==1)
            G1_I_SS = find(strcmp(par.n_ID,[par.nG{g} '_O' nL])==1);
            if(OLM ~= 1)
            % CONNECT E <-> O-LM 
                connect_neurons(G1_E, G1_I_SS, par.([par.nG{g} '_EO_R'])(L), par.([par.nG{g} '_EO_TS'])(L), par.([par.nG{g} '_EO_D'])(L), 'off',...
                    par.([par.nG{g} '_EO_W'])(L), [], 0, [], par.([par.nG{g} '_EO_top']), par.([par.nG{g} '_EO_cos'])(L));
                connect_neurons(G1_I_SS, G1_E, par.([par.nG{g} '_OE_R'])(L), par.([par.nG{g} '_OE_TS'])(L), par.([par.nG{g} '_OE_D'])(L), 'off',...
                    par.([par.nG{g} '_OE_W'])(L), [], 0, [], par.([par.nG{g} '_OE_top']), par.([par.nG{g} '_OE_cos'])(L));
                if(strcmp(par.([par.nG{g} '_E_SS']), 'yes')==1)
                    % CONNECT E_SS <-> O-LM 
                    connect_neurons(G1_E_SS, G1_I_SS, par.([par.nG{g} '_PO_R'])(L), par.([par.nG{g} '_PO_TS'])(L), par.([par.nG{g} '_PO_D'])(L), 'off',...
                        par.([par.nG{g} '_PO_W'])(L), [], 0, [], par.([par.nG{g} '_PO_top']), par.([par.nG{g} '_PO_cos'])(L));
                    connect_neurons(G1_I_SS, G1_E_SS, par.([par.nG{g} '_OP_R'])(L), par.([par.nG{g} '_OP_TS'])(L), par.([par.nG{g} '_OP_D'])(L), 'off',...
                        par.([par.nG{g} '_OP_W'])(L), [], 0, [], par.([par.nG{g} '_OP_top']), par.([par.nG{g} '_OP_cos'])(L));
                end
            end            
            % CONNECT OLM <-> I_FS within each group
            connect_neurons(G1_I_SS, G1_I, par.([par.nG{g} '_OI_R'])(L), par.([par.nG{g} '_OI_TS'])(L), par.([par.nG{g} '_OI_D'])(L), 'off',...
                par.([par.nG{g} '_OI_W'])(L), [], 0, [], par.([par.nG{g} '_OI_top']), par.([par.nG{g} '_OI_cos'])(L));
            connect_neurons(G1_I, G1_I_SS, par.([par.nG{g} '_IO_R'])(L), par.([par.nG{g} '_IO_TS'])(L), par.([par.nG{g} '_IO_D'])(L), 'off',...
                par.([par.nG{g} '_IO_W'])(L), [], 0, [], par.([par.nG{g} '_IO_top']), par.([par.nG{g} '_IO_cos'])(L));
            % CONNECT OLM -> OLM within each group
            connect_neurons(G1_I_SS, G1_I_SS, par.([par.nG{g} '_OO_R'])(L), par.([par.nG{g} '_OO_TS'])(L), par.([par.nG{g} '_OO_D'])(L), 'off',...
                par.([par.nG{g} '_OO_W'])(L), [], 0, [], par.([par.nG{g} '_OO_top']), par.([par.nG{g} '_OO_cos'])(L));
        end
    end
    make.con_neurons(g) = cputime - make.con_neurons(g);
end

%% CREATE & CONNECT NEO-CORTICAL STIMULUS INPUT
make.input_noise = cputime;
if(sum(strcmp('NC', par.nG)) > 0)
    for L = 1 : par.n_SG_NC
        if(strcmp(par.sim_type, 'recall') ~= 1)
            if(isfield(par, 'NC_stims_t') ~= 1)
                for k = 1:10000
                    stims = sort(round(rand(1, par.NC_n_stim(L)) * par.NC_stim_L));
                    stim_d = stims(2:end) - stims(1:end-1);
                    if(sum(stim_d<150)==0); break; end
                end
                par.NC_stims_t(L, :) = stims;
            end

            stim_gamma = ones(size(par.NC_stims_t(L,:)))*par.NC_stim_tau;
            if(par.n_SG_NC <= 1); nL = ''; else; nL = ['_L' int2str(L)]; end
            G1_E = find(strcmp(par.n_ID,['NC_E' nL])==1);
            par.NC_stims_N(L,:) = G1_E(ceil(rand(1, length(par.NC_stims_t)) * length(G1_E)));
            add_input(par.NC_stim_TS, stim_gamma, ...
                par.NC_stim_BG, par.NC_stim_W, par.NC_stims_t(L,:), ...
                par.NC_pre_stim_L, G1_E, par.NC_stim_N_cos, par.NC_stims_N(L,:))
        end
    end
end

%% TC INPUT TO KICK-START TIME CELLS
G1_E = find(strcmp(par.n_ID,'TC_E_L1')==1); G1_E = G1_E(1:par.TC_stim_N);
L = par.TC_pre_stim_L/par.d_T ...
    : (par.TC_pre_stim_L + par.TC_stim_L)/par.d_T;
I.ADP(G1_E, L) = I.ADP(G1_E, L) + par.TC_stim_W;

make.input_noise = cputime - make.input_noise;

cn = cputime;

 %% CONNECT INTER-REGIONS TOPOLOGY
NC_E = find(cellfun(@isempty,strfind(par.n_ID,'NC_E'))==0);
BP_E = find(cellfun(@isempty,strfind(par.n_ID,'BP_E'))==0);

% NC <-> BP
BP_NC_STDP.T_Ca = par.BP_NC_T_Ca; BP_NC_STDP.G_p = par.BP_NC_G_p; BP_NC_STDP.G_d = par.BP_NC_G_d;
BP_NC_STDP.T_p = par.BP_NC_T_p; BP_NC_STDP.T_d = par.BP_NC_T_d;
BP_NC_STDP.C_pre = par.BP_NC_C_pre; BP_NC_STDP.C_post = par.BP_NC_C_post;
connect_neurons(NC_E, BP_E, par.NC_BP_E_R, par.NC_BP_E_TS, par.NC_BP_E_D, 'off',...
    par.NC_BP_E_W_ini, par.NC_BP_E_W_max, par.NC_BP_E_W_sigma,  [], par.NC_BP_E_top, []); 
connect_neurons(BP_E, NC_E, par.BP_NC_R, par.BP_NC_TS, par.BP_NC_D, 'off',...
    par.BP_NC_W_ini, par.BP_NC_W_max, par.BP_NC_W_sigma,  BP_NC_STDP, 'normal', []); 

% TC <-> BP
TC_BP_STDP.G_p = par.TC_BP_G_p; TC_BP_STDP.G_d = par.TC_BP_G_d;
TC_BP_STDP.T_p = par.TC_BP_T_p; TC_BP_STDP.T_d = par.TC_BP_T_d;
TC_BP_STDP.C_pre = par.TC_BP_C_pre; TC_BP_STDP.C_post = par.TC_BP_C_post; 
BP_TC_STDP.T_Ca = par.BP_TC_T_Ca; BP_TC_STDP.G_p = par.BP_TC_G_p; BP_TC_STDP.G_d = par.BP_TC_G_d;
BP_TC_STDP.T_p = par.BP_TC_T_p; BP_TC_STDP.T_d = par.BP_TC_T_d;
BP_TC_STDP.C_pre = par.BP_TC_C_pre; BP_TC_STDP.C_post = par.BP_TC_C_post;
if(sum(strcmp('BP', par.nG)) > 0)
    for L = 1:par.n_SG_TC
        TC_E = find(cellfun(@isempty,strfind(par.n_ID,['TC_E_L' int2str(L)]))==0);

        TC_BP_STDP.T_Ca = par.TC_BP_T_Ca(L); 
        connect_neurons(TC_E, BP_E, par.(['TC_L' int2str(L) '_BP_R']),  par.(['TC_L' int2str(L) '_BP_TS']), par.(['TC_L' int2str(L) '_BP_D']) , 'off',...
            par.(['TC_L' int2str(L) '_BP_W_ini']), par.(['TC_L' int2str(L) '_BP_W_max']), par.(['TC_L' int2str(L) '_BP_W_sigma']), TC_BP_STDP, 'normal', []); 

        connect_neurons(BP_E, TC_E, par.BP_TC_R, par.BP_TC_TS, par.BP_TC_D, 'off',...
            par.BP_TC_W_ini, par.BP_TC_W_max, par.BP_TC_W_sigma, BP_TC_STDP, 'normal', []);
    end
end
make.con_neurons = cputime - cn + sum(make.con_neurons);    

%% NO SELF CONNECTIONS IN WEIGHT MATRIX
for i=1:size(weight_matrix,1); weight_matrix(i,i) = 0; weight_matrix_STDP(i,i) = 0; delay(i,i) = 0; tau_syn(i,i) = 0; end 

%% function timers
make.totaltime_s = cputime - make.totaltime_s;
make.backg_noise = sum(make.backg_noise);
make.remaining_s = make.totaltime_s - sum([make.backg_noise make.con_neurons make.input_noise]);
%make

end

%% CONNECT NEURON GROUPS FUNCTION
function [] = connect_neurons(N1, N2, R, TS, D, recipricol, W_ini, W_max, W_sigma, STDP, top_type, top_width)
    global weight_matrix; global delay; global tau_syn; global weight_matrix_STDP; global wm_max; global par;
    if(strcmp(recipricol,'on')==1) % recipricol connections?
        I=2; 
        if(length(R)~=2 || length(W_ini) ~=2 || length(D) ~=2 || length(TS) ~=2)
           D = [D D]; R = [R R]; W_ini = [W_ini W_ini]; TS = [TS TS]; 
        end
    else; I=1; 
    end 
    
    for i=1:I
        G1 = eval(['N' int2str(i)]); G2 = eval(['N' int2str(3-i)]); % assign group 1 & 2
        syn = zeros(length(G1),length(G2)); % ALl to All connectivity
        for n = 1:length(G2)
           x = transpose(randperm(length(G1), round(R(i)/100*length(G1))));
           syn(x,n) = 1;
        end
        % Create weights from active synapses
        if(W_ini(i) == 0); W = zeros(size(syn)); 
        else
            W = syn .* normrnd(ones(size(syn)) * W_ini(i), W_sigma(i)); % normal distribution of weights
            if(W_ini(i) > 0); W(W<0)=0; elseif(W_ini(i) < 0); W(W>0)=0; end % ensure no sign reversal
        end
        % Define topology
        for n = 1:length(G1)
            x = n * (length(G2)/length(G1));
            if(strcmp(top_type,'normal')==1)% even distribution
                top = ones(length(G2),1);
            elseif(contains(top_type,'C')==1) % closed loop i.e. circular
                top = normpdf(G2 - min(G2)+1,x,top_width); top = top*(1/max(top));
            elseif(contains(top_type,'O')==1) % open loop i.e. linear
                top_M = normpdf(G2 - min(G2)+1, x, top_width); % N CONNECTIONS
                top_U = normpdf(G2 - min(G2)+1, x + length(G2), top_width); % N+ CONNECTIONS
                top_L = normpdf(G2 - min(G2)+1, x - length(G2), top_width); % N- CONNECTIONS
                top = top_M + top_L + top_U;
                top = top*(1/max(top));
            end
            % enable farthest neighbour as inverse of default nearest neighbour
            if(contains(top_type,'FN')==1); top = 1 - top; top = top*(1/max(top)); end
            top(sqrt(top.^2)<0.0001)=0;  % topology curves to 4 d.p.
            syn(n,:) = syn(n,:) .* transpose(top>0); % adjust possible synapses via topology map
            top = W(n, :) .* transpose(top); top(sqrt(top.^2)<0.0001)=0; W(n, :) = top; % weights * topology to 4 d.p.
        end
        weight_matrix(G1,G2) = W; % add weight initials
        [~, ts] = ismember(TS(i), par.uTS);
        tau_syn(G1, G2, ts) = syn; % add synaptic time constants
        delay(G1,G2) = syn .* D(i); % add spike delay

        if(isempty(STDP) ~= 1) % add STDP parameters
            weight_matrix_STDP(G1,G2) = syn;
            if(W_ini(i) == W_max(i)) % if initial is already max
                wm_max(G1, G2) = W .* syn;
            else % define new random max weights
                W = syn .* normrnd(ones(size(syn)) * W_max(i), W_sigma(i)); % normal distribution of weights
                if(W_max(i) > 0); W(W<0)=0; elseif(W_max(i) < 0); W(W>0)=0; end % ensure no sign reversal
                wm_max(G1, G2) = W .* syn; % add weight maximums
            end
            par.T_Ca(G1, G2) = syn * STDP.T_Ca; % calcium decay rate
            par.G_p(G1, G2) = syn * STDP.G_p; par.G_d(G1, G2) = syn * STDP.G_d; % LTP/LTD rates
            par.T_p(G1, G2) = syn * STDP.T_p; par.T_d(G1, G2) = syn * STDP.T_d; % LTP/LTD thresholds
            par.C_pre(G1, G2) = syn * STDP.C_pre; par.C_post(G1, G2) = syn * STDP.C_post; % pre/post calcium increases
        end
    end
end

%% ADD SPIKES FUNCTION
function [] = add_input(spike_ts, input_ts, rate, weight, stims, pre_stim, G_E, cos, N )
    global par; global I;
    psp = find(spike_ts == par.uTS);
    g_mx = find(par.PSP{psp}==max(par.PSP{psp}) == 1); g_L = length(par.PSP{psp});
        
    T_stim = 0:1:par.NC_pre_stim_L*2 + par.NC_stim_L - 1;
    stim = zeros(length(G_E), length(T_stim));
    for n=1:length(stims)
        stim_L = (exp(1)*T_stim./input_ts(n)) .* exp(-T_stim./input_ts(n)) .* heaviside(T_stim);
        stim_L = stim_L(stim_L>0.01);
        % choose random neuron to focus gaussian distribution around
        top_M = normpdf(G_E, N(n), cos); % N CONNECTIONS
        top_U = normpdf(G_E, N(n) + length(G_E), cos); % N+ CONNECTIONS
        top_L = normpdf(G_E, N(n) - length(G_E), cos); % N- CONNECTIONS
        top = top_M + top_L + top_U;
        top = top*(1/max(top));
        stim(:, (stims(n) + pre_stim):...
            (stims(n) + pre_stim) + length(stim_L)-1) = repmat(stim_L,[length(G_E),1]) .* top;
    end
       
    spikes = spike_train(length(G_E), rate, par.sim_length, 1); % make spike train
    spikes = spikes .* stim;
    [i, x, c] = find(spikes>0);
    x = unique(x); 
    % add spikes centred post-synaptic current to future input of neurons
   for k=1:length(x)
       I.SYN(G_E, max(1,x(k)-g_mx+1):min(par.sim_length, x(k)+g_L-g_mx), psp) = ...
           I.SYN(G_E, max(1,x(k)-g_mx+1):min(par.sim_length, x(k)+g_L-g_mx), psp) + ...
           par.PSP{psp}(max(1,g_mx-x(k)+1):min(g_L, par.sim_length - x(k) + g_mx)) ...
           .*spikes(:,x(k)) * weight;
   end
end

