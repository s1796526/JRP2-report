%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Pre-processing pipeline  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hello
% eegppf is a code package developed by Sareen et al. 2020 and publically available at:
% https://data.mendeley.com/datasets/fshy54ypyh/1 
% (1) Sareen, E., Singh, L., Varkey, B., Achary, K., Gupta, A. (2020). 
%     EEG dataset of individuals with intellectual and developmental disorder
%     and healthy controls under rest and music stimuli. Data in Brief, 105488, ISSN 2352-3409,
%     DOI:https://doi.org/10.1016/j.dib.2020.105488.
% (2) Sareen, E., Gupta, A., Verma, R., Achary, G. K., Varkey, B (2019), 
%     Studying functional brain networks from dry electrode EEG set during music and resting states
%     in neurodevelopment disorder, bioRxiv 759738 [Preprint]. 
%     Available from: https://www.biorxiv.org/content/10.1101/759738v1

% This code is written by Jonathan as a reference code to channel all files
% automatically through the pipeline, grabbing the individual condition
% files and saving the file in a new directory 
% questions about this code to be sent to:
% jonathan.benhaiem (at) gmail.com


mainfolder = dir('D:\xxx\'); % Define main directory 
mainfolder = mainfolder(~ismember({mainfolder.name}, {'.', '..','.DS_Store'})); % define which files to exclude, '.' is a default invisible folder in Windows systems whilr .DS_Store is Apple's counterpart. 

%% For loop grabbing, the participant folder, then the condition file in the participant folder
for i = 1: numel(mainfolder)
    subfolder = dir(strcat(mainfolder(i).folder,'\',mainfolder(i).name,'\RAW')); % concatinate the file to select in automatically as the folders progress
    subfolder = subfolder(~ismember({subfolder.name}, {'.','..', '.DS_Store'}));
    for j = 1 : numel(subfolder)
        subsub_folder = dir(strcat(subfolder(j).folder,'\',subfolder(j).name)); % Condition folder
        subsub_folder = subsub_folder(~ismember({subsub_folder.name}, {'.','..', '.DS_Store'}));

        % Arguments for EEGPPF and clean_data command
        File_path = subsub_folder.folder; 
        File_name = subsub_folder.name;
        chan_loc = 'E:\\VIBES\\Code\\EEG_Analysis\\Pipeline\\Helper_func\\chanlocnew.ced'; % Channel location file, adapted from BrainVision
        Fs = 500; % Sampling rate of EEG device (Hertz)
        lpf = 0.5; % Low-pass filter from 0.5 to 128 hz
        hpf = 128;
        linenoise = [50 100];  % Define line noise
        electrode_select = {'FP1','Fz','F3','F7','FT9','FC5','FC1','C3','T7','TP9','CP5','CP1','Pz','P3','P7','O1','Oz','O2','P4','P8','TP10', 'CP6', 'CP2','Cz','C4','T8','FT10','FC6','FC2','F4','F8','FP2'}; % Select all electrodes to be pre-processed
        rem_num = 3; % define how many components to be removed in the after the ICA. 
        cleandata = eegppf(File_name,File_path,chan_loc,Fs,lpf,hpf,linenoise,electrode_select,rem_num); % Clean data function eegppf (Sareen et al, 2020) 
        disp('Your data is clean, and will now be saved') 

        % saving cleandata in a .set file after converting it from .mat
        storepath = strcat(mainfolder(i).folder,'\',mainfolder(i).name,'\CLEAN\');
        storename = File_name;
        EEG = pop_importdata('dataformat','matlab','nbchan',32,'data',cleandata,'srate',500,'pnts',0,'xmin',0); % EEGLAB import of pre-processed file to be able to save in .set format
        EEG = pop_chanedit(EEG, 'load',{'E:\\VIBES\\Code\\EEG_Analysis\\Pipeline\\Helper_func\\chanlocnew.ced','filetype','autodetect'}); % Add chanlocation file again to save inside .set format
        pop_saveset( EEG, 'filename',storename,'filepath',storepath); % save the file into file directory. 
        disp('data saved')
    end
end