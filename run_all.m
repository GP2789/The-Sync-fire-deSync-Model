function [ ] = run_all( filename, n_trials )
% This script will run the simulations as described in the paper. If you
% would like to change parameters, i.e. the number of events per pattern,
% or any biophysical properties, then read set_parameters.m. 
% Will reproduce figures 4 & 5 from the paper. 

%% SIMULATE AND PLOT A SINGLE PATTERN (FIGURE 4)
filename_1 = recall_experiment( filename, 'yes', 1, 'off' );
preprocess_data( filename_1 );
plot_single_trial( filename_1 );
%make_videos( filename_1 ); % for videos from supplementary section

%% SIMULATE AND PLOT SIMILARITY BETWEEN MANY PATTERNS (FIGURE 5)
filename_2 = recall_experiment( filename, 'no', n_trials, 'off' );
preprocess_data( filename_2 );
content_specificity( filename_2 );

end

