%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FUNCTIONAL CONNECTIVITY %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code is written by Jonathan as a reference code to compute
% functional connectivity in pre-processed eeg files. Both the
% Phase-locking value (PLV) and phase-lag index (PLI) are desribed here.
% For the code of the PLV and PLI functions, please contact 
% ekansh.sareen (@) epfl.ch
%
% questions about this specific code to be sent to:
% jonathan.benhaiem (at) gmail.com

%% Define mainfolder
mainfolder = dir('D:\VIBES\EEG\DATA\Organised\temp\final'); % Define main folder with directory
mainfolder = mainfolder(~ismember({mainfolder.name}, {'.','..', '.DS_Store'})); % exclude invisible folders
%% start EEGLab without GUI
eeglab nogui 
%% Initialise matrices
% constants 
windowLength = 15000; % 1 epoch for the PLV is 15000 samples or 30 seconds. 
totalLength = 120001 - 1; % define length of file
trials = totalLength / windowLength; % time divided by epoch length 
srate = 500; % Sampling rate in Hertz
numElec = 32; % 32 electrodes
cond = 3; % Three conditions in exp. design (eyes closed, open and anticipation)
subj = numel(mainfolder); % Working with 15 participant --> 15 folders inside mainfolder

% tensors
plv_tensor = zeros(numElec, numElec, cond, subj);
avg_plv_nonthresh = zeros(numElec, numElec, cond, subj);
pli_tensor = zeros(numElec, numElec, cond, subj);
avg_pli_nonthresh = zeros(numElec, numElec, cond, subj);

%% Loop to load 
for subi = 1:numel(mainfolder)
    subfolder = dir(strcat(mainfolder(subi).folder,'\',mainfolder(subi).name,'\CLEAN\old')); % Define subfolder
    subfolder = subfolder(~ismember({subfolder.name}, {'.','..', '.DS_Store'})); % exclude invisible folders
    %% Loop to load each condition
    for filei = 1:numel(subfolder)
        filename = subfolder(filei).name; 
        filepath = subfolder(filei).folder;
        EEG = pop_loadset(filename, filepath); % Load EEG data set into eeglab format
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = pop_chanedit(EEG, 'load',{'D:\\VIBES\\Code\\EEG_Analysis\\Pipeline\\Helper_func\\chanlocnew.ced','filetype','autodetect'}); % add chanloc file

        %% Extract data and create empty matrix for windows (3D)
        EEG.data = double(EEG.data); % Convert to double precision data
        data = EEG.data;  % Data extracted from EEG struct
        chanlocs = EEG.chanlocs; % Channel location extracted from the struct
        
        %% Epoch data (30 seconds/epoch (fs = 500 Hz)  --> 15000 samples/epoch)
    	allWindows = zeros(32,windowLength,trials);
        for j = 1 :size(allWindows,3)
           allWindows (:, :, j) = data(:, windowLength*(j - 1) + 1 : (windowLength*j));  % 
        end 
        
        %% PLV
        plv = zeros(numElec, numElec, trials);
        for epochi = 1:trials
            epoch = allWindows(:,:,epochi);
            plv(:,:, epochi) = PLV_1(epoch);
        end
        avgPLV = mean(plv,3);
        avg_plv_nonthresh(:,:,filei, subi) = avgPLV;

        %% Phase Lag Index
        pli = zeros(numElec, numElec, trials);
        for epochi = 1:trials
            epoch = allWindows(:,:,epochi);
            pli(:,:, epochi) = PLI_1(epoch);
        end
        avgPLI = mean(pli,3);
        avg_pli_nonthresh(:,:,filei,subi) = avgPLI;
        
    end
end

%% Saving PLV matrices
% PLV unthresholded
save('D:\VIBES\EEG\Results\Graphs\PLV\plv_tensor_nonthresh.mat', 'avg_plv_nonthresh', '-mat');

%% Average non thresholded PLV over subjects
nonthresh_avg_plv_tensor = mean(avg_plv_nonthresh, 4);
save('D:\VIBES\EEG\Results\Graphs\PLV\plv_tensor_nonThresh.mat', 'nonthresh_avg_plv_tensor', '-mat');

%% extract eyes closed and save  
nonthresh_avg_plv_closed = nonthresh_avg_plv_tensor(:,:,1); 
save('D:\VIBES\EEG\Results\Graphs\PLV\nonthresh_avg_plv_closed.mat', 'nonthresh_avg_plv_closed', '-mat');

%% extract eyes open and save
nonthresh_avg_plv_open = nonthresh_avg_plv_tensor(:,:,2); 
save('D:\VIBES\EEG\Results\Graphs\PLV\nonthresh_avg_plv_open.mat', 'nonthresh_avg_plv_open', '-mat');

%% extract anticipation and save
nonthresh_avg_plv_anti = nonthresh_avg_plv_tensor(:,:,3); 
save('D:\VIBES\EEG\Results\Graphs\PLV\nonthresh_avg_plv_anti.mat', 'nonthresh_avg_plv_anti', '-mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving PLI matrices
% PLI unthresholded
save('D:\VIBES\EEG\Results\Graphs\PLI\pli_tensor_nonthresh.mat', 'avg_pli_nonthresh', '-mat');

%% NON-thresholded averaged over subjects
nonthresh_avg_pli_tensor = mean(avg_pli_nonthresh, 4);
save('D:\VIBES\EEG\Results\Graphs\PLI\nonthresh_avg_pli_tensor.mat', 'nonthresh_avg_pli_tensor', '-mat');
%% Extract eyes closed
nonthresh_avg_pli_closed = nonthresh_avg_pli_tensor(:,:,1); 
save('D:\VIBES\EEG\Results\Graphs\PLI\nonthresh_avg_pli_closed.mat', 'nonthresh_avg_pli_closed', '-mat');
%% Extract eyes open 
nonthresh_avg_pli_open = nonthresh_avg_pli_tensor(:,:,2); 
save('D:\VIBES\EEG\Results\Graphs\PLI\nonthresh_avg_pli_open.mat', 'nonthresh_avg_pli_open', '-mat');
%% Extract anticipation 
nonthresh_avg_pli_anti = nonthresh_avg_pli_tensor(:,:,3); 
save('D:\VIBES\EEG\Results\Graphs\PLI\nonthresh_avg_pli_anti.mat', 'nonthresh_avg_pli_anti', '-mat');


%% Step in between is to apply a threshold using EEGNET toolbox 
% I did it by hand at 30% proportional threshold. 
% Then the newly acquired thresholded graphs are imported here and the
% network measures can be computed.

%% Network Measures PLI
% for all functions, they can be found in 
% the Brain Connectivity Toolbox (Rubinov, Sporns, 2010)

% Convert weight to path length to compute the clustering coefficient
w_norm_closed = weight_conversion(avg_pli_closed, 'normalize');
w_norm_open = weight_conversion(avg_pli_open, 'normalize');
w_norm_anti = weight_conversion(avg_pli_anti, 'normalize');

% clustering coefficient (wu = weighted undirected). 
clustering_coef(:,1) = clustering_coef_wu(w_norm_closed);
clustering_coef(:,2) = clustering_coef_wu(w_norm_open);
clustering_coef(:,3) = clustering_coef_wu(w_norm_anti);

% efficiency is the average of inverse shortest path length 
efficiency(:,1) = efficiency_wei(avg_pli_closed);
efficiency(:,2) = efficiency_wei(avg_pli_open);
efficiency(:,3) = efficiency_wei(avg_pli_anti);

% degrees (use for node size in brainNetViewer)
degrees(:,1) = degrees_und(avg_pli_closed);
degrees(:,2) = degrees_und(avg_pli_open);
degrees(:,3) = degrees_und(avg_pli_anti);

%% PLotting adjencency matrices 
% Chanlabels
eeglab nogui;
chanloc = readtable('chanlocnew.txt'); 
labels = chanloc.labels;
XlabelNames = [labels];
YlabelNames = [labels];

% Load files Unthresholded (which I included inside the example_file folder)
% Eyes closed
load('nonthresh_avg_plv_closed.mat');
load('nonthresh_avg_pli_closed.mat');
% load('nonthresh_AvgCorrClosed.mat');
% Eyes open
load('nonthresh_avg_plv_open.mat');
load('nonthresh_avg_pli_open.mat');
% Anticipation
load("nonthresh_avg_plv_anti.mat");
load("nonthresh_avg_pli_anti.mat");

% Thresholded (at 30%)
% PLV
% eyes closed
prop_avg_plv_closed = load('prop_avg_plv_closed.mat'); 
prop_avg_plv_closed = prop_avg_plv_closed.SavedMatrix;
prop_avg_plv_closed = squeeze(prop_avg_plv_closed(1,:,:));
% eyes open
prop_avg_plv_open = load('prop_avg_plv_open.mat'); 
prop_avg_plv_open = prop_avg_plv_open.SavedMatrix;
prop_avg_plv_open = squeeze(prop_avg_plv_open(1,:,:));
% anticipation
prop_avg_plv_anti = load('prop_avg_plv_anti.mat'); 
prop_avg_plv_anti = prop_avg_plv_anti.SavedMatrix;
prop_avg_plv_anti = squeeze(prop_avg_plv_anti(1,:,:));

% PLI
% eyes closed
prop_avg_pli_closed = load('prop_avg_pli_closed.mat'); 
prop_avg_pli_closed = prop_avg_pli_closed.SavedMatrix;
prop_avg_pli_closed = squeeze(prop_avg_pli_closed(1,:,:));

% eyes open
prop_avg_pli_open = load('prop_avg_pli_open_2.mat'); 
prop_avg_pli_open = prop_avg_pli_open.SavedMatrix;
prop_avg_pli_open = squeeze(prop_avg_pli_open(1,:,:));

% antipation
prop_avg_pli_anti = load('prop_avg_pli_anti.mat'); 
prop_avg_pli_anti = prop_avg_pli_anti.SavedMatrix;
prop_avg_pli_anti = squeeze(prop_avg_pli_anti(1,:,:));


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Thesis Plots %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non thresholded images PLI
nonthresh_anti = imagesc(nonthresh_avg_pli_anti); colorbar; title('Anticipation'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32); set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'nonthresh_anti.png');

nonthresh_open = imagesc(nonthresh_avg_pli_open); colorbar; title('Eyes open'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32);set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'nonthresh_open.png');

nonthresh_closed = imagesc(nonthresh_avg_pli_closed); colorbar; title('Eyes closed'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32); set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'nonthresh_closed.png');

% 30percent threshold images PLI
prop_anti = imagesc(prop_avg_pli_anti); colorbar; title('Anticipation'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32);set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'prop_anti.png');

prop_open = imagesc(prop_avg_pli_open); colorbar; title('Eyes open'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32); set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'prop_open.png');

prop_closed = imagesc(prop_avg_pli_closed); colorbar; title('Eyes closed'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32);set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'prop_closed.png');

% PLV
% non thresh 
nonthresh_anti_plv = imagesc(nonthresh_avg_plv_anti); colorbar; title('Anticipation'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32); set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'nonthresh_anti_plv.png');

nonthresh_open_plv = imagesc(nonthresh_avg_plv_open); colorbar; title('Eyes open'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32);set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'nonthresh_open_plv.png');

nonthresh_closed_plv = imagesc(nonthresh_avg_plv_closed); colorbar; title('Eyes closed'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32); set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'nonthresh_closed_plv.png');

% 30percent threshold images PLV
prop_anti_plv = imagesc(prop_avg_plv_anti); colorbar; title('Anticipation'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32);set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'prop_anti_plv.png');

prop_open_plv = imagesc(prop_avg_plv_open); colorbar; title('Eyes open'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32); set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'prop_open_plv.png');

prop_closed_plv = imagesc(prop_avg_plv_closed); colorbar; title('Eyes closed'); set(gca,'YTickLabel',YlabelNames,  'Ytick', 1:32);set(gca,'XTickLabel',YlabelNames,  'Xtick', 1:32); set(gca, 'FontSize', 9);
saveas(gcf, 'prop_closed_plv.png');



