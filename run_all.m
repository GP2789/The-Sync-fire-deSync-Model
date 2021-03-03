function [ ] = run_all( filename )
% This script will run the simulations as described in the paper. If you
% would like to change parameters, i.e. the number of events per pattern,
% or any biophysical properties, then read set_parameters.m. 
%% NB please ensure FieldTrip is not on your MatLab pathway to use these scripts.

%% FIGURE 6 
% simulate paradigm for attentional blink
[ filename_def_AB ] = recall_experiment( [filename '/AB_def'], 2000, 1, ... % name / n patterns / n trials
    {'NC_n_stim','NC_stim_L','BP_OE_TS'}, {2,1200,15} ); % set parameters
[ filename_short_AB ] = recall_experiment( [filename '/AB_short'], 2000, 1, ... % name / n patterns / n trials
    {'NC_n_stim','NC_stim_L','BP_OE_TS'}, {2,1200,5} ); % set parameters
[ filename_long_AB ] = recall_experiment( [filename '/AB_long'], 2000, 1, ... % name / n patterns / n trials
    {'NC_n_stim','NC_stim_L','BP_OE_TS'}, {2,1200,30} ); % set parameters

% preprocess simulated data & extract information about patterns presented
preprocess_data( filename_def_AB );
extract_pattern_info( filename_def_AB );

preprocess_data( filename_long_AB );
preprocess_data( filename_short_AB );

extract_pattern_info( filename_long_AB );
extract_pattern_info( filename_short_AB );

% calculate content specific reinstatement for patterns presented
content_specificity( filename_def_AB );

% plot attentional blink Figure 6
plot_AB( filename_def_AB, {filename_short_AB, filename_long_AB} )

%% FIGURE 7
% simulate data for content specificity
[ filename_CS ] = recall_experiment( [filename '/multi'], 2000, 1, ... % name / n patterns / n trials
    {'NC_n_stim','NC_stim_L','BP_OE_TS'}, {3,1500,15} ); % set parameters
preprocess_data( filename_CS );
extract_pattern_info( filename_CS );
content_specificity( filename_CS );
plot_content_spec( filename_CS );

%% FIGURE 4
% analyse NC stimulus vs alpha resting state
analyse_NC([filename '/var_NC_input'], 20, true);

%% FIGURE 5
% SIMULATE AND PLOT A SINGLE PATTERN 
[ filename_single ] = recall_experiment( [filename '/single'], 1, 1, {}, {} );
preprocess_data( filename_single );
plot_single_trial( filename_single );

end

