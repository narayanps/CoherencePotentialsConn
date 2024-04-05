%%%%%%%%AUTHOR : NARAYAN P SUBRAMANIYAM%%%%%%%%%%%%
%%%%%%%%%%%%%2023 \ SAPIEN LABS%%%%%%%%%%%%%%%%%%%%%


clear all
close all
addpath(genpath(strcat(pwd,'/external/')))
addpath(strcat(pwd, '/func'))
addpath(genpath(strcat(pwd, '/data/PatternTN_AllTasks')))
conds = {'EyesClosedPassive', 'EyesOpenVisualTask', 'EyesOpenPatternCompletionTask', 'EyesOpenWorkingMemoryTask'};

sub=16160;
flags = [1 0]; % 1 is negative peaks, 0 positive peaks
for ii = 1:length(conds)
    
    %%%Read edf file
    [hdr,data] = read_eeg_edf(sprintf('%d_%s.edf', sub, conds{ii}));
    
    %set filtering parameters
    hpf=0.1;
    lpf=40;
    filt_eeg = eegfilt(data(3:16,:),hdr.samplingrate,hpf,lpf,length(data));
    
    %remove the mean
    filt_eeg = filt_eeg - mean(filt_eeg,2);
    
    sd = 15; % SD computed as mean of all channels in eyes closed condition
    fact = 3.0; % sd factor - set to 3 for faster computation in demo  --lesser events to detect as threshold increases
    offset=1;%in seconds %ignore 1 sec in the begining and end of the file
    Fs=hdr.samplingrate;
    all_chan_labels = cellstr(hdr.channels);
    eeg_chan_labels = all_chan_labels(3:16,:);
    num_channels = length(eeg_chan_labels);
    %detect positive and negetive deflections in every channel
    for jj=1:length(flags)
        for j=1:14
            x=filt_eeg(j,:);
            chan_features{jj,j} = get_features(x, sd, flags(jj),fact,offset, Fs, eeg_chan_labels{j});
        end
    end
    for j=1:length(eeg_chan_labels)
        if ~isempty(chan_features{1,j})
            chan_features{1,j}.event_id = strcat(chan_features{1,j}.event_id, '_neg');
        end
        if ~isempty(chan_features{2,j})
            chan_features{2,j}.event_id = strcat(chan_features{2,j}.event_id, '_pos');
        end
    end
    
    
    %align the peaks of CPs and compute correlation
    thr=0.3; % set clustering threshold here to be used later
    window_size=1;%in seconds
    [C, event_id, peak_index,l,ls] = align_and_correlate_v2(filt_eeg, chan_features, window_size, Fs, flags);
    for j=1:length(eeg_chan_labels)
        if ~isempty(chan_features{1,j})
            neg_eve = length(chan_features{1,j}.event_id);
        else
            neg_eve = 0;
        end
        
        if ~isempty(chan_features{2,j})
            pos_eve = length(chan_features{2,j}.event_id);
        else
            pos_eve = 0;
        end
        
        total_eve_per_chan(j) = neg_eve + pos_eve;%sum(cell2mat(indx));
    end
    %hierarchical clustering based on corr matrix
    if ~isempty(C)
        peak_index_s = peak_index/Fs;
        corr = cell2mat(C(:,3));
        clear C
        [ Z, I, cophenetic_CC] = make_hierarchical_cluster(corr,'average',0);
        for cc=1:length(thr)
            T_dist = cluster(Z,'cutoff',thr(cc), 'Criterion', 'distance');
            no_of_clust = length(unique(T_dist));
            for i=1:1:no_of_clust
                clust{i,1}.events = find(T_dist==i);
            end
            for i=1:1:length(clust)
                clust{i,1}.peak_index = peak_index(clust{i,1}.events);
                clust{i,1}.event_id = event_id(clust{i,1}.events);
            end
            %compute CP-based connectivity based on CP features
            % connectivity matrix size nchan X nchan X 1 X conds
            % e.g. connectivity matrix for EC based on inter-event-interval
            % EC = squeeze(A_def.A_iei(:,:,1,1))
            % EO = squeeze(A_def.A_iei(:,:,1,2))
            % PC = squeeze(A_def.A_iei(:,:,1,3))
            % WM = squeeze(A_def.A_iei(:,:,1,4))
            min_ne=2;
            Fs=hdr.samplingrate;
            [A_def.A_iei(:,:,cc,ii), A_def.A_iei_min(:,:,cc,ii), A_def.A_iei_max(:,:,cc,ii), A_def.A_min_cf(:,:,cc,ii),A_def.A_all_cf(:,:,cc,ii),...
                A_def.A_all_cf_norm(:,:,cc,ii),A_def.A_min_cf_norm(:,:,cc,ii),...
                A_def.A_min_cf_norm_clust(:,:,cc,ii),A_def.A_perc_cf(:,:,cc,ii),a{cc,ii}] ...
                = get_clust_conn_matrix_updated(clust, eeg_chan_labels,Fs, min_ne, total_eve_per_chan);
            cond_clust{ii,cc}.clust = clust;
            clear clust
            clear corr
            
        end
    else
        A_def.A_iei(:,:,1:length(thr),ii)=nan(14,14);
        A_def.A_iei_min(:,:,1:length(thr),ii)=nan(14,14);
        A_def.A_iei_max(:,:,1:length(thr),ii)=nan(14,14);
        
        A_def.A_min_cf(:,:,1:length(thr),ii)=zeros(14,14);
        A_def.A_all_cf(:,:,1:length(thr),ii)=zeros(14,14);
        A_def.A_min_cf_norm(:,:,1:length(thr),ii)=zeros(14,14);
        A_def.A_min_cf_norm_clust(:,:,1:length(thr),ii)=zeros(14,14);
        A_def.A_all_cf_norm(:,:,1:length(thr),ii)=zeros(14,14);
        A_def.A_perc_cf(:,:,1:length(thr),ii)=zeros(14,14);
    end
    clear chan_features
end
topoplot_connect(a,'emotiv.locs') %need EEGLAB