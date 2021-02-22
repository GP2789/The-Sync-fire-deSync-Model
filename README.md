IMPORTANT: UPDATE FOR NEUROPSCYHOLOGIA SUBMISSION ARRIVING ON THE WEEK OF 1ST MARCH


# The-Sync-fire-deSync-Model

Deciphering episodic content from cortical alpha oscillations.

This code will run the simulations used for the Sync-fire deSync model paper submission. 
Simply run the file, run_all.m, specifying a directory to store data in and the number of 
trials to simulate the paradigm with (10 were simulated for the paper submission).

If you want to change parameters, see set_parameters.m.
If you want to produce the video files supporting the submission, comment out the line
% make_videos(filename); in run_all.m. This will produce all the necessary frames from
a single pattern simulation, though will require video edititing software to peice together.
