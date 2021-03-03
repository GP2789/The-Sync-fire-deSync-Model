function [ LFP ] = create_LFP(spikes, hann_on, sample_on, pband, sim_length, d_T )

    spikes(:,1) = spikes(:,1) - min(spikes(:,1)) + 1;
    T = 0:d_T:sim_length;
    sim_L = round(sim_length / d_T);
   
    spikes(:,2) = round(spikes(:,2) / d_T);
    MUA = zeros(max(spikes(:,1)), sim_L);
    for i=1:size(spikes,1)
       MUA(spikes(i,1),spikes(i,2)) = MUA(spikes(i,1),spikes(i,2)) + 1; 
    end
    
    MUA_PSP = zeros(size(MUA,1),size(MUA,2));
    tau_syn = 1.5;
    g_syn = (exp(1)*T./tau_syn).*exp(-T./tau_syn).*heaviside(T); 
    g_syn = g_syn(g_syn>0.01);
    g_mx = find(g_syn==max(g_syn)); 
    g_L = length(g_syn);
    
    
    for m1=1:size(MUA,1)
        for m2 = 1:size(MUA,2)
           if(MUA(m1,m2)==1) 
               MUA_PSP(m1, max(1,m2-g_mx+1):min(sim_L, m2+g_L-g_mx)) = ...
                   MUA_PSP(m1, max(1,m2-g_mx+1):min(sim_L, m2+g_L-g_mx)) + ...
                   g_syn(max(1,g_mx-m2+1):min(g_L, sim_L - m2 + g_mx)) * MUA(m1, m2);
           end
        end
    end
    
    LFP = sum(MUA_PSP,1);
    
    if(hann_on == 1)
        filtord = round(30 / d_T);
        b = hann(filtord)./filtord;
        LFP = filtfilt(b,1,LFP);
    end
    
    if(sample_on == 1)
        Fs = 1000 / d_T;  % Sampling Frequency

        Fstop1 = 0.5;           % First Stopband Frequency
        Fstop2 = 101;           % Second Stopband Frequency
        Astop1 = 60;          % First Stopband Attenuation (dB)

        Apass  = 1;           % Passband Ripple (dB)
        Astop2 = 80;          % Second Stopband Attenuation (dB)
        match  = 'stopband';  % Band to match exactly
        Fpass1 = pband(1);           % First Passband Frequency
        Fpass2 = pband(2);           % Second Passband Frequency
        % Construct an FDESIGN object and call its BUTTER method.
        h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                              Astop2, Fs);
        Hd = design(h, 'butter', 'MatchExactly', match);
        LFP=filtfilt(Hd.sosMatrix,Hd.ScaleValues,LFP);
    else; LFP = LFP - mean(LFP);
    end
%     powhil=abs(hilbert(LFP));
%     
%     bl=[0.25 (pre_stim_length/1000)-0.25];
%     bl1=find(bl(1)==time);
%     bl2=find(bl(2)==time);
%     powhil_rc=(powhil-mean(powhil(bl1:bl2),2))./mean(powhil(bl1:bl2),2);
%     
    
    
%     LFP=filter(hanning(30)./15.5,1,LFP);% Smooth lfp a bit to get 'EEG like' waveforms
%     LFP=flipud(LFP'); LFP=filter(hanning(30)./15.5,1,LFP'); % hanning filter
%     LFP=flipud(LFP'); LFP=LFP'-mean(LFP',2);
end

