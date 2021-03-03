function [ new_filename ] = concatenate_files( filenames, new_filename )

P = 0; T = 0;
for i = 1:length(filenames)
    load([filenames{i} '/Data/parameters.mat'])
    if(par.B ~= 1); P = P + par.B; end
    if(par.trials ~= 1); T = T + par.trials; end
end
if(P == 0); P = 1; end; if(T == 0); T = 1; end
new_filename = [new_filename '_' int2str(P) 'P_' int2str(T) 'T'];

mkdir(new_filename);
mkdir([new_filename '\Data']);
pattern_c = 0;
for i = 1:length(filenames)
    %% concatenate parameters
    fprintf(['copying files from ' filenames{i} ' ... \n'])
    D_par = dir([filenames{i} '\Data']);
    D_pat = D_par(contains({D_par(:).name},'patterns'));
    D_var_par = D_par(contains({D_par(:).name},'varied_parameters'));
    D_par = D_par(strcmp({D_par(:).name},'parameters.mat'));
    load([D_par(1).folder '\' D_par(1).name])
    if(~isempty(D_pat))
        load([D_pat(1).folder '\' D_pat(1).name]);
        if(i == 1)
            rand_phase = par.rand_phase;
            NC_stims_t_all = par.NC_stims_t_all;
            NC_stims_N_all = par.NC_stims_N_all;
            patterns = patterns_t;
        else
            rand_phase = cat(3, rand_phase, par.rand_phase);
            NC_stims_t_all = cat(2, NC_stims_t_all, par.NC_stims_t_all);
            NC_stims_N_all = cat(2, NC_stims_N_all, par.NC_stims_N_all);
            patterns = [patterns; patterns_t];
        end
    end
    if(~isempty(D_var_par)) 
        load([D_var_par(1).folder '\' D_var_par(1).name]); 
        if(i == 1); var_pars_all = []; end
        var_pars_all = [var_pars_all; var_pars];
    end
    
    %% copyfiles
    D_data = dir([filenames{i} '\Data']);
    D_data = D_data([D_data(:).isdir] & ~ismember({D_data(:).name},{'.','..'}));
    for j = 1:length(D_data)
        filename_t = D_data(j).name;
        pattern = str2double(filename_t(2:end)) + pattern_c;
        if(P == 1); N = ['T' int2str(pattern) '\']; 
        else; N = ['B' int2str(pattern) '\']; 
        end
        if(exist([new_filename '\Data\' N],'dir')~=7)
            mkdir([new_filename '\Data\' N]);
        end
        D_sim = dir([D_data(j).folder '\' D_data(j).name]);
        if(P == 1); D_sim = D_sim(~[D_sim(:).isdir] & ~ismember({D_sim(:).name},{'.','..'}));
        else; D_sim = D_sim([D_sim(:).isdir] & ~ismember({D_sim(:).name},{'.','..'}));
        end
        
        for k = 1:length(D_sim)
            filename_t = D_sim(k).name;
            if(P ~= 1)
                trial = str2double(filename_t(2:end));
                if(exist([new_filename '\Data\B' int2str(pattern) '\T' int2str(trial)],'dir')~=7)
                    mkdir([new_filename '\Data\B' int2str(pattern) '\T' int2str(trial)])
                end
                 D_trial = dir([D_sim(k).folder '\' D_sim(k).name]);
                D_trial = D_trial(~ismember({D_trial(:).name},{'.','..'}));
                for t = 1:length(D_trial)
                    copyfile([D_trial(t).folder '\' D_trial(t).name], ...
                        [new_filename '\Data\B' int2str(pattern) '\T' int2str(trial) '\' D_trial(t).name])
                    fprintf('\b\b\b\b\b%3.0f%%\n', ((t + (k-1)*2 + (j-1)*par.trials*2) / (par.B*par.trials*2)) * 100)
                end
            else
                copyfile([D_sim(k).folder '\' D_sim(k).name], ...
                    [new_filename '\Data\' N D_sim(k).name])
            end
            fprintf('\b\b\b\b\b%3.0f%%\n', ((k+(j-1)*length(D_sim)) / (par.B*par.trials*length(D_sim))) * 100)
           
        end
    end
    pattern_c = pattern_c + length(D_data);
end
if(~isempty(D_pat))
    par.rand_phase = rand_phase;
    par.NC_stims_t_all = NC_stims_t_all;
    par.NC_stims_N_all = NC_stims_N_all;
    patterns_t = patterns;
    save([new_filename '\Data\patterns.mat'], 'patterns_t')
end
if(~isempty(D_var_par)) 
    var_pars = var_pars_all;
    save([new_filename '\Data\varied_parameters.mat'], 'var_par_N', 'var_par_x', 'var_pars')
end
if(P ~= 1); par.B = pattern_c;
else; par.trials = pattern_c;
end

save([new_filename '\Data\parameters.mat'], 'par')


end

