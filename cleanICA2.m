%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 2nd ICA to check quality of components  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I developed this script using mostly eeglab functions, to compute a second ICA
% all 32 components of the pre-processed data are determined and labeled, 
% to identify which files have more than a certain amount of components. 
% This not a pre-processing step but simply a code to verify the quality of the data. 
%
% questions about this specific code to be sent to:
% jonathan.benhaiem (at) gmail.com
%% dir mainfolder
mainfolder = dir('E:\VIBES\EEG\DATA\Organised'); % adapt to personal directory
mainfolder = mainfolder(~ismember({mainfolder.name}, {'.', '..', '.DS_Store'})); % exclude invisible folders

eeglab nogui
%% dir subfolder & subsubfolder
for i = 1: numel(mainfolder)
    subfolder = dir(strcat(mainfolder(i).folder,'\',mainfolder(i).name,'\CLEAN\'));
    subfolder = subfolder(~ismember({subfolder.name}, {'old','cut', '.','..', '.DS_Store'}));
    for j = 1: numel(subfolder)
        subsub_folder = dir(strcat(subfolder(j).folder,'\',subfolder(j).name));
        subsub_folder = subsub_folder(~ismember({subsub_folder.name}, {'.','..', '.DS_Store'}));
        remEpochs = [0 0 0 0 0 0 0 0]; % vector of how many epochs are good to keep
        %% load EEG set into matlab
        for filei = 1:numel(subsub_folder)
            filename = subsub_folder(filei).name; 
            filepath = subsub_folder(filei).folder;
            EEG = pop_loadset(filename, filepath); % Load EEG data set into eeglab format
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = pop_chanedit(EEG, 'load',{'E:\\VIBES\\Code\\EEG_Analysis\\Pipeline\\Helper_func\\chanlocnew.ced','filetype','autodetect'}); % add chanlocation file
            %% Extract data and extract chanlocs
            EEG.data = double(EEG.data); % make double precision data
            data = EEG.data;  % Data extracted from EEG struct
            
            %% run ICLabel on each .set file
            EEG = eeg_checkset( EEG );
            EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on'); % ICA decomposition (type runica)
            EEG = eeg_checkset( EEG );
            EEG = pop_iclabel(EEG, 'default'); % Labelling of components (ICLabel)
            EEG = eeg_checkset( EEG );
            EEG = pop_icflag(EEG, [NaN NaN;0.85 1;0.85 1;0.85 1;0.85 1;0.85 1;0.85 1]); % Setting thresholds for flagging components
            EEG = eeg_checkset( EEG );
            vector = EEG.reject.gcompreject; % vectorise all flagged components
            ncomp = sum( vector( : ) == 1 ); % Sum of all flagged components          
            
            %% Only save epochs with enough good components
            if ncomp < 10
                storepath = strcat(mainfolder(i).folder,'\',mainfolder(i).name,'\SUPERCLEAN\',subfolder(j).name);
                storename = filename;
                dest = fullfile(storepath,storename);
                pop_saveset( EEG, 'filename',storename,'filepath',storepath);
                disp('data saved')
                remEpochs(filei) = 1;
            end
        end
    end
end
