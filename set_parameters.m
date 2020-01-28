function [] = set_parameters( )
%% HOW TO READ THIS SCRIPT
% REGION = NC: neo-cortex; 
%          BP: binding pool; 
%          TC: time-cells/synfire chains;
% NEURON TYPE = E/FS E: excitatory / E fast spiking; 
%               I/FS I: inhibitory / I fast spiking; 
%               O/OLM: off node / O-LM cell; 
%               P/E2/E_SS: propagation node / excitatory type2 / E slow spiking;
% SIZE GROUP =  n_SG_x: number of layers in region x;
%               n_x: number of neurons in each layer of region x; 
%               L1/L2: layer 1/layer 2 - multi layered region will have array for layer parameter i.e. [x1 x2];
% REGION PROPERTIES = BG: background noise (spikes/second);
%                     BG_TS: synaptic time constant of BG spikes;
%                     BG_W: weight of BG spikes;
%                     x_W_sigma: sample weights from around a normal distribution of sigma;
%                     x_percent: percent of neurons for each neuron type (I/OLM/E2);
%                     x_OLM/x_SS: 'yes/no' to enable/disable O/P nodes in region;
% CONNECTION PROP. = (general note) R_NN_x: Region _ Neuron type 1 to Neuron type 2 _ parameter; 
%                    R_NN_TS: synaptic time constant;
%                    R_NN_R: random % connectivity (0-100);
%                    R_NN_W: weight of connections;
%                    R_NN_top: topology type - NN_: nearest neighbour;
%                                            - FN_: farthest neighbour;
%                                            - _O: open i.e. line of nodes;
%                                            - _C: closed i.e. circle of nodes;
%                    R_NN_cos: sigma of normpdf centred over neuron (used for FN/NN conns only);
%                    R_NN_D: delay of spike event to reach connected neurons (ms);
%                    _FF/FB_: feedforward (L1-L2) / feedback (L2-L1);
%                    _NM: neuro-modulator;
% STDP PROPERTIES = (general note) R1_R2_x: Region 1 _ Region 2 _ parameter; 
%                   R1_R2_T_Ca: calcium time constant;
%                   R1_R2_C_pre: presynaptic calcium amplitude increase;
%                   R1_R2_C_post: postsynaptic calcium amplitude increase;
%                   R1_R2_T_p: threshold for LTP;
%                   R1_R2_T_d: threshold for LTD;
%                   R1_R2_G_p: rate at which LTP occurs;
%                   R1_R2_G_d: rate at which LTD occurs;
%                   C_D: calcium delay from neuron 1 to neuron 2;

%% NOTES
% These scripts contain parameters for simulation of O-LM cells &
% heterosynaptic STDP that are not used in the submitted simulations. 

%% SET PARAMETER SPACE FUNCTION
     global par;
     if(isempty(par)==1)
        par = struct;
     end

    %% CHOOSE PARAMETERS
    parameter_settings = ...
        ... %%%% NEURONS PARAMETERS %%%%
        {'model_type', 'HH', 'V_m_initial', -60, ...
        'attractor', 0.5, 'T_WC', 50, 'C_D', 13, ...
        'TC_L1_BP_NM', 0.011, 'TC_L2_BP_NM', 0.028, 'BP_NC_NM', 2,...
        'T_pre_F', 20, 'T_post_F', 20, ...
        ... % NEURON WITHIN GROUPS PARAMETERS
        ... %%%%%%%%% NEO-CORTEX %%%%%%%%
        'n_SG_NC', 1, 'n_NC', 35, 'NC_OLM', 'no', 'NC_E_SS', 'no',...
        'NC_BG', 1500, 'NC_BG_TS', 1.5, 'NC_BG_W', 0.05, 'NC_W_sigma', 0, ...
        'NC_I_percent', 30, 'NC_OLM_percent', 15,...
        ... % FS E ->
        'NC_EE_TS', 1.5, 'NC_EE_R', 100, 'NC_EE_W', 1, 'NC_EE_top', 'NN_O', 'NC_EE_cos', 1, 'NC_EE_D', 2,...
        'NC_EI_TS', 5, 'NC_EI_R', 100, 'NC_EI_W', 1, 'NC_EI_top', 'NN_O', 'NC_EI_cos', 1, 'NC_EI_D', 2,... 
        'NC_EO_TS', 5, 'NC_EO_R', 100, 'NC_EO_W', 0.1, 'NC_EO_top', 'normal', 'NC_EO_cos', 1, 'NC_EO_D', 2,... 
        ... % FS I ->
        'NC_IE_TS', 12.5, 'NC_IE_R', 100, 'NC_IE_W', -1.25, 'NC_IE_top', 'FN_O', 'NC_IE_cos', 1, 'NC_IE_D', 2,...
        'NC_II_TS', 5, 'NC_II_R', 100, 'NC_II_W', -1.25, 'NC_II_top', 'FN_O', 'NC_II_cos', 1, 'NC_II_D', 2,...
        'NC_IO_TS', 5, 'NC_IO_R', 0, 'NC_IO_W', 0, 'NC_IO_top', 'FN_O', 'NC_IO_cos', 1, 'NC_IO_D', 2,...
        ... % OLM ->
        'NC_OE_TS', 25, 'NC_OE_R', 100, 'NC_OE_W', -0.15, 'NC_OE_top', 'normal', 'NC_OE_cos', 3, 'NC_OE_D', 2,...
        'NC_OI_TS', 5, 'NC_OI_R', 0, 'NC_OI_W', 0, 'NC_OI_top', 'FN_O', 'NC_OI_cos', 1, 'NC_OI_D', 2,...
        'NC_OO_TS', 5, 'NC_OO_R', 0, 'NC_OO_W', 0, 'NC_OO_top', 'FN_O', 'NC_OO_cos', 1, 'NC_OO_D', 2,...
        ... % NC -> NC
        'NC_NC_TS', 1.5, 'NC_NC_R', 0, 'NC_NC_W', 0, 'NC_NC_D', 2, 'NC_NC_LR', 0, ...
        ... %%%%%%% BINDING POOL %%%%%%%%
        'n_SG_BP', 1, 'n_BP', 70, 'BP_OLM', 'yes', 'BP_E_SS', 'no'...
        'BP_I_percent', 30, 'BP_OLM_percent', 10, ...
        'BP_BG', 125, 'BP_BG_TS', 1.5, 'BP_BG_W', 0, 'BP_W_sigma', 0, ...
        ... % STDP
        'BP_NC_T_Ca', 20, 'BP_NC_C_pre', 0.8, 'BP_NC_C_post', 0.7, 'BP_NC_T_p', 1.3, 'BP_NC_T_d', 1, 'BP_NC_G_p', 1.5, 'BP_NC_G_d', 0.05,...
        'BP_BP_T_Ca', 20, 'BP_BP_C_pre', 1, 'BP_BP_C_post', 1, 'BP_BP_T_p', 5, 'BP_BP_T_d', 1, 'BP_BP_G_p', 0.05, 'BP_BP_G_d', 1,...
        'BP_TC_T_Ca', 20, 'BP_TC_C_pre', 0.8, 'BP_TC_C_post', 0.8, 'BP_TC_T_p', 1.3, 'BP_TC_T_d', 1, 'BP_TC_G_p', 1.5, 'BP_TC_G_d', 0.05, ...
        ... % FS E ->
        'BP_EE_TS', 5, 'BP_EE_R', 100, 'BP_EE_W', 5, 'BP_EE_top', 'NN_O', 'BP_EE_cos', 1, 'BP_EE_D', 1,... 
        'BP_EI_TS', 1.5, 'BP_EI_R', 100, 'BP_EI_W', 8, 'BP_EI_top', 'NN_O', 'BP_EI_cos', 0.5, 'BP_EI_D', 1,... 
        'BP_EO_TS', 25, 'BP_EO_R', 100, 'BP_EO_W', 0.008, 'BP_EO_top', 'normal', 'BP_EO_cos', 1, 'BP_EO_D', 1,... 
        ... % I ->
        'BP_IE_TS', 1.5, 'BP_IE_R', 100, 'BP_IE_W', -8, 'BP_IE_top', 'FN_O', 'BP_IE_cos', 3, 'BP_IE_D', 1,...
        'BP_II_TS', 1.5, 'BP_II_R', 100, 'BP_II_W', -8, 'BP_II_top', 'FN_O', 'BP_II_cos', 1.5, 'BP_II_D', 1,...
        'BP_IO_TS', 5, 'BP_IO_R', 0, 'BP_IO_W', 0, 'BP_IO_top', 'FN_O', 'BP_IO_cos', 5, 'BP_IO_D', 1,...
        ... % P ->
        'BP_OE_TS', 5, 'BP_OE_R', 100, 'BP_OE_W', -5, 'BP_OE_top', 'normal', 'BP_OE_cos', 2, 'BP_OE_D', 1,...
        'BP_OI_TS', 5, 'BP_OI_R', 0, 'BP_OI_W', 0, 'BP_OI_top', 'FN_O', 'BP_OI_cos', 5, 'BP_OI_D', 1,...
        'BP_OO_TS', 5, 'BP_OO_R', 0, 'BP_OO_W', 0, 'BP_OO_top', 'FN_O', 'BP_OO_cos', 5, 'BP_OO_D', 1,...
        ... %%%%%%%% TIME CELLS %%%%%%%%%
        'n_SG_TC', 2, 'n_TC', [16 64], 'TC_OLM', 'yes', 'TC_E_SS', 'yes',...
        'TC_I_percent', 12.5, 'TC_OLM_percent', 12.5, 'TC_P_percent', 12.5, ...
        'TC_BG', 0, 'TC_BG_TS', 1.5, 'TC_BG_W', 1, 'TC_W_sigma', 0, ...
        ... % STDP
        'TC_BP_T_Ca', [20 20], 'TC_BP_C_pre', 0.8, 'TC_BP_C_post', 0.8, 'TC_BP_T_p', 1.3, 'TC_BP_T_d', 1, 'TC_BP_G_p', 1.5, 'TC_BP_G_d', 0.05, ...
        ... % FS E ->
        'TC_EO_TS', [1.5 1.5], 'TC_EO_R', [0 0], 'TC_EO_W', [0 0], 'TC_EO_top', 'normal', 'TC_EO_cos', [1 1], 'TC_EO_D', [2 2], ...
        'TC_EI_TS', [1.5 1.5], 'TC_EI_R', [100 100], 'TC_EI_W', [1 1], 'TC_EI_top', 'normal', 'TC_EI_cos', [1 1], 'TC_EI_D', [2 2], ...
        'TC_EP_TS', [30 30], 'TC_EP_R', [100 100], 'TC_EP_W', [0.015 0.01], 'TC_EP_top', 'normal', 'TC_EP_cos', [1 1], 'TC_EP_D', [2 2],...
        'TC_EE_TS', [10 10], 'TC_EE_R', [100 100], 'TC_EE_W', [2 3], 'TC_EE_top', 'NN_O', 'TC_EE_cos', [0.5 0.5], 'TC_EE_D', [2 2],...
        ... % FS I ->
        'TC_IO_TS', [1.5 1.5], 'TC_IO_R', [0 0], 'TC_IO_W', [0 0], 'TC_IO_top', 'FN_C', 'TC_IO_cos', [10 10], 'TC_IO_D', [2 2], ...
        'TC_II_TS', [1.5 1.5], 'TC_II_R', [100 100], 'TC_II_W', [-3 -3], 'TC_II_top', 'FN_C', 'TC_II_cos', [1 1], 'TC_II_D', [2 2], ...
        'TC_IP_TS', [1.5 1.5], 'TC_IP_R', [0 0], 'TC_IP_W', [0 0], 'TC_IP_top', 'FN_C', 'TC_IP_cos', [10 10], 'TC_IP_D', [2 2], ...
        'TC_IE_TS', [1.5 1.5], 'TC_IE_R', [100 100], 'TC_IE_W', [-4 -4], 'TC_IE_top', 'normal', 'TC_IE_cos', [10 10], 'TC_IE_D', [2 2],...
        ... % SS E ->
        'TC_PO_TS', [30 30], 'TC_PO_R', [100 100], 'TC_PO_W', [0.07 0.5], 'TC_PO_top', 'normal', 'TC_PO_cos', [1 1], 'TC_PO_D', [2 2], ...
        'TC_PI_TS', [10 10], 'TC_PI_R', [0 0], 'TC_PI_W', [0 0], 'TC_PI_top', 'NN_C', 'TC_PI_cos', [10 10], 'TC_PI_D', [2 2],... 
        'TC_PP_TS', [10 10], 'TC_PP_R', [0 0], 'TC_PP_W', [0 0], 'TC_PP_top', 'NN_C', 'TC_PP_cos', [1 1], 'TC_PP_D', [2 2],... 
        'TC_PE_TS', [25 25], 'TC_PE_R', [100 100], 'TC_PE_W', [0.4 0.75], 'TC_PE_top', 'NN_C', 'TC_PE_cos', [2 2], 'TC_PE_D', [2 2],... 
        ... % O ->
        'TC_OO_TS', [30 30], 'TC_OO_R', [0 0], 'TC_OO_W', [0 0], 'TC_OO_top', 'FN_C', 'TC_OO_cos', [10 10], 'TC_OO_D', [2 2], ...
        'TC_OI_TS', [30 30], 'TC_OI_R', [0 0], 'TC_OI_W', [0 0], 'TC_OI_top', 'NN_C', 'TC_OI_cos', [1 1], 'TC_OI_D', [2 2], ...
        'TC_OP_TS', [1.5 1.5], 'TC_OP_R', [0 0], 'TC_OP_W', [-3 -3], 'TC_OP_top', 'normal', 'TC_OP_cos', [10 10], 'TC_OP_D', [2 2],...
        'TC_OE_TS', [20 20], 'TC_OE_R', [100 100], 'TC_OE_W', [-3 -3], 'TC_OE_top', 'NN_C', 'TC_OE_cos', [10 10], 'TC_OE_D', [2 2],...
        ... % LAYERS
        'TC_TC_TS', [1.5 1.5], 'TC_TC_R', [100 100], 'TC_TC_W', [1 1], 'TC_TC_W_sigma', [0 0], 'TC_TC_D', [1000 125], 'TC_TC_LR', [0 0], ...
        'TC_FF_TS', [30], 'TC_FF_R', 100, 'TC_FF_W', [0.15], 'TC_FF_D', 2, ...
        'TC_FB_TS', 30, 'TC_FB_R', 100, 'TC_FB_W', [0.17], 'TC_FB_top', 'normal', 'TC_FB_cos', 10, 'TC_FB_D', 2,...
        ... %%% NEURON BETWEEN GROUPS PARAMETERS %%%
        'NC_BP_E_R', 60, 'NC_BP_E_W_max', 0.0125, 'NC_BP_E_W_ini', 0.0125, 'NC_BP_E_top', 'normal', 'NC_BP_E_W_sigma', 0.0005, 'NC_BP_E_D', 2, 'NC_BP_E_TS', 15, ...
        'NC_BP_I_R', 0, 'NC_BP_I_W_max', 0.0, 'NC_BP_I_W_ini', 0.0, 'NC_BP_I_W_sigma', 0, 'NC_BP_I_D', 2, 'NC_BP_I_TS', 5, ...
        'BP_NC_R', 100, 'BP_NC_W_max', 0, 'BP_NC_W_ini', 0, 'BP_NC_W_sigma', 0,'BP_NC_D', 2, 'BP_NC_TS', 60, ...
        'BP_TC_R', 100, 'BP_TC_W_max', 0, 'BP_TC_W_ini', 0, 'BP_TC_W_sigma', 0.01, 'BP_TC_D', 1, 'BP_TC_TS', 1.5, ...
        'TC_L1_BP_R', 100, 'TC_L1_BP_W_max', 0, 'TC_L1_BP_W_ini', 0, 'TC_L1_BP_W_sigma', 0.01, 'TC_L1_BP_D', 1, 'TC_L1_BP_TS', 20, ...
        'TC_L2_BP_R', 100, 'TC_L2_BP_W_max', 0, 'TC_L2_BP_W_ini', 0, 'TC_L2_BP_W_sigma', 0.01, 'TC_L2_BP_D', 1, 'TC_L2_BP_TS', 1.5, ...
        ... %%% STIMULUS PRESENTATION PARAMETERS %%%
        'TC_pre_stim_L', 2050, 'TC_stim_L', 10, ... % pre stim & stim length of TC burst
        'TC_stim_W', 0.15, ... % weight of DC burst into TC cell(s) 
        'TC_stim_N', 1, ... % number of cells DC burst fed into
        'NC_n_stim', 3, ... % number of events
        'NC_stim_N_cos', 1.5, ... % width of gaussian kernel around randomly selected neuron
        'NC_pre_stim_L', 2500, ... % length of pre-stim period before events (ms)
        'NC_stim_L', 1500, ... % time over which events are randomly interspersed (ms)
        'NC_stim_type', 'gamma', 'NC_stim_tau', 60, ... % spikes multiplied by gamma curve
        'NC_stim_BG', 1000, 'NC_stim_TS', 1.5, 'NC_stim_W', 3}; % event spike input parameters
   
    % SET PARAMETERS IF NOT ALREADY SUPPLIED
    for i=1:length(parameter_settings)/2
        if(any(strcmp(parameter_settings{(i-1)*2+1},fieldnames(par)))==0)
            par.([parameter_settings{(i-1)*2+1}]) = parameter_settings{(i-1)*2+2};
        end
    end
    
    %% NEURON SIMULATION PARAMETERS
    if(strcmp(par.model_type, 'IAF')==1)
        par.V_th = -55; par.E = -70; par.c_m = 0.65; par.tau_m = 50; par.d_T = 1; par.RF = 6;
    elseif(strcmp(par.model_type, 'HH')==1)
        par.PYR_E_L = -67; par.PYR_g_bar_L = 0.1; % Leak 
        par.PYR_E_Na = 50; par.PYR_g_bar_Na = 100; % Sodium
        par.PYR_E_K = -100; par.PYR_g_bar_K = 80; % Potassium
        
        par.OLM_E_L = -65; par.OLM_g_bar_L = 0.5; % Leak 
        par.OLM_E_Na = 55; par.OLM_g_bar_Na = 52; % Sodium
        par.OLM_E_K = -90; par.OLM_g_bar_K = 11; % Potassium
        par.OLM_E_h = -20; par.OLM_g_bar_h = 1.46; % h current
        par.OLM_g_bar_p = 0.5; % Alternate Potassium
        
        par.c_m = 1; par.d_T = 0.01;
    end
    
    %% ORDER OF THE SIMULATION
    par.sim_order = {'learning','recall'}; % simulate learning / recall 
    par.sim_order_n = {'DL','AL'}; % during learning / after learning
    par.nG = {'TC','BP','NC'}; % include which regions
    
    %% DEFINE POPULATION SIZES BASED ON PARAMETERS
    par.network_size = 0;
    for i=1:length(par.nG)
        for L = 1:par.(['n_SG_' par.nG{i}])
            par.(['n_' par.nG{i} '_I'])(L) = round(par.(['n_' par.nG{i}])(L)*par.([par.nG{i} '_I_percent'])/100); 
            n = par.(['n_' par.nG{i} '_I'])(L);
            if(L==1); par.([par.nG{i} '_N']) = {'E'}; end
            if(strcmp(par.([par.nG{i} '_E_SS']), 'yes')==1)
                par.(['n_' par.nG{i} '_P'])(L) = round(par.(['n_' par.nG{i}])(L)*par.([par.nG{i} '_P_percent'])/100); 
                n = n + par.(['n_' par.nG{i} '_P'])(L);
                if(L==1); par.([par.nG{i} '_N']){end+1} = 'P'; end
            end
            if(L==1); par.([par.nG{i} '_N']){end+1} = 'I'; end
            if(strcmp(par.([par.nG{i} '_OLM']), 'yes')==1)
                par.(['n_' par.nG{i} '_O'])(L) = round(par.(['n_' par.nG{i}])(L)*par.([par.nG{i} '_OLM_percent'])/100); 
                n = n + par.(['n_' par.nG{i} '_O'])(L);
                if(L==1); par.([par.nG{i} '_N']){end+1} = 'O'; end
            end
            par.(['n_' par.nG{i} '_E'])(L) = par.(['n_' par.nG{i}])(L) - n;

            par.network_size = par.network_size + par.(['n_' par.nG{i}])(L);
        end
    end
    
    %% find unique synaptic time constants in simulation
    fN = fieldnames(par); nTS = strfind(fN,'_TS');
    iTS = find(~cellfun(@isempty,nTS)); uTS = [];
    for i=1:length(iTS)
        uTS = [uTS; transpose(par.([fN{iTS(i)}]))];
    end
    uTS = unique(uTS); par.uTS = uTS(uTS~=0);
end



