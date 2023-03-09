%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Frontal Alpha and Beta Power %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define mainfolder directory and exclude invisible folders
mainfolder = dir('E:\Education\University\Lausanne\2. VIBES\ANALYSIS\EEG\DATA\Organised\temp\final');
mainfolder = mainfolder(~ismember({mainfolder.name}, {'.','..', '.DS_Store'}));

eeglab nogui % startup the eeglab without the user interface (background)
nSubjets = numel(mainfolder); % define number of subjects
fs = 500; % define sampling rate (500 samples/s = 500 Hz in our case)
epoch = fs*2; % Epoch duration is twice the sampling rate. 
allAbsAlphas = zeros(nSubjets,3); % Initialise matrices
allRelAlphas = zeros(nSubjets,3); % Initialise matrices
allAlphasAssymetry = zeros(nSubjets,3); % Initialise matrices

%% Import data and add channel location
for i = 1:numel(mainfolder)
    subfolder = dir(strcat(mainfolder(i).folder,'\',mainfolder(i).name,'\CLEAN\old')); % Define subfolder
    subfolder = subfolder(~ismember({subfolder.name}, {'.','..', '.DS_Store'})); % exclude invisible folders
    for j = 1:numel(subfolder)
        filename = subfolder(j).name; % Define filename
        filepath = subfolder(j).folder; % Define filepath 
        EEG = pop_loadset(filename, filepath); % Load EEG data set into EEGLAB format
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 ); 
        EEG = pop_chanedit(EEG, 'load',{'E:\Education\University\Lausanne\2. VIBES\Code\EEG_Analysis\Pipeline\Helper_func\chanlocnew.ced','filetype','autodetect'}); % add chanloc file
        chanlocs = EEG.chanlocs; % Channel location extracted from the structure of the EEGLAB format 
        
        %% extract double precision data of frontal elecs
        EEG.data = double(EEG.data);  %convert to double-precision
        data = EEG.data; % ([1:4 30:32],:) Grab 7 frontal electrodes if
        % you want to work on frontal spectral power

        %% pWelch method
        allpwelch = zeros(7,251); % Initialise matrix to compute fast fourier transform in (pWelch)
        window = 1000; 
        for welchi = 1:size(data,1) % For each channel compute separately 
            welchdata = data(welchi,:);
            [pxx,f] = pwelch(welchdata,window,300,500,fs);
            allpwelch(welchi,:) = pxx';
        end
        %% Code to plot the compute fft
%             plot(f,10*log10(pxx))
%             xlim([0.5 128]);
%             xlabel('Frequency (Hz)')
%             ylabel('PSD (dB/Hz)')

 
        %% Absolute and relative (alpha) power
        % defined freq range (Basar 2013)
        % delta = 1-4; theta: 4-7 Hz; alpha: 8-13 Hz; beta, 18-25, and gamma: 30-70 Hz
        % Alpha Power
        alphabounds = ([18 25]-1); % Alpha boundaries in hz
        freqidx_a = dsearchn(f, alphabounds');  % convert to indices because hz is a linearly spaced vector
        alphapower = mean(allpwelch(:,freqidx_a(1):freqidx_a(2)),2);  % extract average absolute power
      
        % Average Alpha power over electrodes
        AvgAbsAlphapower = mean(alphapower,1);
        allAbsAlphas(i,j) = AvgAbsAlphapower;
       
        % Compute relative power
        allbounds = [0 250]-1; % All boundaries in hz
        freqidx_all = dsearchn(f, allbounds');  % convert to indices because hz is a linearly spaced vector
        allpower = mean(allpwelch(:,freqidx_all(1):freqidx_all(2)),2);  % extract average absolute power
        
        %avgpower = mean(chanpowr,2); % power of entire frequency range
        AvgRelAlphaPower = mean((alphapower./allpower),1);
        allRelAlphas(i,j) = AvgRelAlphaPower;
        
        %% Alpha assymetry
%         F3 and F4 log transformed 
        F3 = alphapower(3,:); F3 = log(F3);
        F4 = alphapower(30,:); F4 = log(F4);
        alphaAssymetryF3F4 = F4 - F3;
        allAlphasAssymetry(i,j) = alphaAssymetryF3F4;

        F7 = alphapower(4,:); F7 = log(F7);
        F8 = alphapower(31,:); F8 = log(F8);
        alphaAssymetryF7F8 = F8 - F7;
        allAlphasAssymetry(i,j) = alphaAssymetryF7F8;
    
        %% Plotting check
%         figure(j);
%         plot(hz,10*log10(allpwelch),'linew',2)
%         xlabel('Frequency (Hz)'), ylabel('PSD (dB/Hz)')
%         set(gca,'xlim',[0 30],'ylim',[0 30])
    end
end