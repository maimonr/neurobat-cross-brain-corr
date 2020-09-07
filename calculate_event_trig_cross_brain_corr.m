function [cross_brain_corr,shuffled_corr,shuffled_corr_p,cross_brain_cohr,cross_brain_corr_index] = calculate_event_trig_cross_brain_corr(activation,varargin)

pnames = {'time_idx','n_shuffle_reps','shuffle_type','dT','time_bin_T','corr_index_flag','pre_calc_shuffled_corr'};
dflts  = {true(1,size(activation,4)),4,'randomPhase',NaN,NaN,false,false,[]};
[time_idx,n_shuffle_reps,shuffleType,dT,time_bin_T,corr_index_flag,pre_calc_shuffled_corr] = internal.stats.parseArgs(pnames,dflts,varargin{:});

nBat = size(activation,1);
nTrial = size(activation,2);
nChannel = size(activation,3);
n_time_wins = size(activation,4);
n_activation_bins = size(activation,5);

channel_pairs = cell(1,2);
[channel_pairs{1},channel_pairs{2}] = meshgrid(1:nChannel,1:nChannel);
channel_pairs = cellfun(@(x) reshape(x,1,[]),channel_pairs,'un',0);
channel_pairs = vertcat(channel_pairs{:})';

bat_pairs = nchoosek(1:nBat,2);
nBat_pair = size(bat_pairs,1);
nChannel_pair = size(channel_pairs,1);

cross_brain_corr = nan(nTrial,nBat_pair,nChannel_pair,n_activation_bins);
cohrFlag = false;
if ~isnan(dT) && ~isnan(time_bin_T)
    cohrFlag = true;
    df = 1/time_bin_T; %Determine the frequency resolution.
    fNQ = round(1/ dT / 2); %Determine the Nyquist frequency.
    f = (0:df:fNQ); %Construct frequency axis
    nFreq = length(f);
    cross_brain_cohr = nan(nFreq,nBat_pair,nChannel_pair,n_activation_bins);
else
    cross_brain_cohr = NaN;
end

if corr_index_flag
    cross_brain_corr_index = nan(nTrial,nBat_pair,nChannel_pair,n_activation_bins,n_time_wins);
else
    cross_brain_corr_index = NaN;
end

if strcmp(shuffleType,'random_shuffle_corr_index')
    nT = n_time_wins;
else
    nT = 1;
end

shuffled_corr_p = nan(nTrial,nBat_pair,nChannel_pair,n_activation_bins,nT);

if strcmp(shuffleType,'preCalc')
    shuffled_corr = pre_calc_shuffled_corr;
elseif strcmp(shuffleType,'none')
    shuffled_corr = [];
else
    shuffled_corr = nan(nTrial,nBat_pair,nChannel_pair,n_activation_bins,n_shuffle_reps,nT);
end

for bat_pair_k = 1:nBat_pair
    bat_idxs = bat_pairs(bat_pair_k,:);
    for f_k = 1:n_activation_bins
        bat_pair_activation = zeros(nChannel_pair,n_time_wins,nTrial,2);
        for channel_k = 1:nChannel_pair
            channel_idxs = channel_pairs(channel_k,:);
            for bat_k = 1:2
                for trial_k = 1:nTrial
                    bat_pair_activation(channel_k,:,trial_k,bat_k) = squeeze(activation(bat_idxs(bat_k),trial_k,channel_idxs(bat_k),:,f_k));
                end
            end
        end
        
        used_channel_idx = find(all(~all(isnan(bat_pair_activation),[2 3]),4))';
        
        if isempty(used_channel_idx)
            continue
        end
        
        bat_pair_activation_full = bat_pair_activation(used_channel_idx,:,:,:);
        bat_pair_activation = bat_pair_activation(used_channel_idx,time_idx,:,:);
        n_used_channels = length(used_channel_idx);
        
        current_cross_brain_corr = nan(nTrial,n_used_channels);
        
        if cohrFlag
           current_cross_brain_cohr = nan(nFreq,n_used_channels);
        end
        
        if corr_index_flag
            current_cross_brain_corr_index = nan(nTrial,n_used_channels,n_time_wins);
        end
        
        current_shuffled_corr = nan(nTrial,n_used_channels,n_shuffle_reps,nT);
        
        if strcmp(shuffleType,'preCalc')
            current_shuffled_corr = reshape(pre_calc_shuffled_corr(:,bat_pair_k,used_channel_idx,f_k,:),[nTrial n_used_channels size(pre_calc_shuffled_corr,5)]);
        end
        
        for used_channel_k = 1:n_used_channels
            current_activation = reshape(bat_pair_activation(used_channel_k,:,:,:),[sum(time_idx),nTrial,2]);
            current_activation_full = reshape(bat_pair_activation_full(used_channel_k,:,:,:),[n_time_wins,nTrial,2]);
            corr_mat = corr(current_activation(:,:,1),current_activation(:,:,2));
            current_cross_brain_corr(:,used_channel_k) = diag(corr_mat);
            
            if cohrFlag
                x = current_activation(:,:,1)';
                y = current_activation(:,:,2)';
                cohr = trial_based_coherence(x,y,dT,time_bin_T);
                current_cross_brain_cohr(:,used_channel_k) = cohr;
            end
            
            if corr_index_flag
                for trial_k = 1:nTrial
                    trial_activation = squeeze(current_activation(:,trial_k,:))';
                    current_cross_brain_corr_index(trial_k,used_channel_k,:) = get_corr_index(trial_activation);
                end
            end
            
            switch shuffleType
                case 'random_shuffle_corr_index'
                    current_shuffled_corr(:,used_channel_k,:,:) = calculate_time_shuffle_corr(current_activation_full,time_idx,n_shuffle_reps,false,'corrIndex');
                case 'randomPhase'            
                    current_shuffled_corr(:,used_channel_k,:,:) = calculate_random_phase_corr(current_activation,n_shuffle_reps,true);
                case 'randomShuffle'
                    current_shuffled_corr(:,used_channel_k,:,:) = calculate_time_shuffle_corr(current_activation_full,time_idx,n_shuffle_reps,false,'corrCoef');
                case 'none'
                    
            end
        end
        
        cross_brain_corr(:,bat_pair_k,used_channel_idx,f_k) = current_cross_brain_corr;
        
        if cohrFlag
           cross_brain_cohr(:,bat_pair_k,used_channel_idx,f_k) = current_cross_brain_cohr; 
        end
        
        if corr_index_flag
            cross_brain_corr_index(:,bat_pair_k,used_channel_idx,f_k,:) = current_cross_brain_corr_index;
        end
        
        if any(strcmp(shuffleType,{'randomPhase','randomShuffle'}))
            shuffled_corr(:,bat_pair_k,used_channel_idx,f_k,:,nT) = current_shuffled_corr;
            shuffled_corr_p(:,bat_pair_k,used_channel_idx,f_k) = sum(current_cross_brain_corr > current_shuffled_corr,3)/n_shuffle_reps;
        elseif strcmp(shuffleType,'random_shuffle_corr_index')
            shuffled_corr(:,bat_pair_k,used_channel_idx,f_k,:,nT) = current_shuffled_corr;
            current_cross_brain_corr_index = reshape(current_cross_brain_corr_index,[nTrial n_used_channels 1 n_time_wins]);
            shuffled_corr_p(:,bat_pair_k,used_channel_idx,f_k,:) = sum(current_cross_brain_corr_index > current_shuffled_corr,3)/n_shuffle_reps;
        elseif strcmp(shuffleType,'preCalc')
            n_shuffle_trials = sum(~isnan(current_shuffled_corr),3);
            shuffled_corr_p(:,bat_pair_k,used_channel_idx,f_k) = sum(current_cross_brain_corr > current_shuffled_corr,3)./n_shuffle_trials;
        end
    end
end

end

function shuffle_corr = calculate_random_phase_corr(all_time_series,n_shuffle_reps,use_trial_avg_flag)

nTrial = size(all_time_series,2);
shuffle_corr = nan(nTrial,n_shuffle_reps);

if use_trial_avg_flag
    trial_avg_time_series = squeeze(nanmean(all_time_series,2))';
end

for trial_k = 1:nTrial
    actual_time_series = squeeze(all_time_series(:,trial_k,:))';
    if ~any(isnan(actual_time_series),'all')
        if use_trial_avg_flag
            shuffle_corr(trial_k,:) = calculate_random_phase_surrogate_corr(actual_time_series,n_shuffle_reps,trial_avg_time_series);
        else
            shuffle_corr(trial_k,:) = calculate_random_phase_surrogate_corr(actual_time_series,n_shuffle_reps);
        end
    end
end

end

function shuffle_corr = calculate_time_shuffle_corr(all_time_series,time_idx,n_shuffle_reps,use_trial_avg_flag,corrType)

nTrial = size(all_time_series,2);
n_time_wins = size(all_time_series,1);

switch corrType
    case 'corrCoef'
        shuffle_corr = nan(nTrial,n_shuffle_reps);
    case 'corrIndex'
        shuffle_corr = nan(nTrial,n_shuffle_reps,n_time_wins);
end

if use_trial_avg_flag
    trial_avg_time_series = reshape(nanmean(all_time_series,2),[n_time_wins 2]);
end

for trial_k = 1:nTrial
    actual_time_series = squeeze(all_time_series(:,trial_k,:));
    if ~any(isnan(actual_time_series),'all')
        
        if use_trial_avg_flag
            actual_time_series = actual_time_series - trial_avg_time_series;
        end
        
        all_shuffle_time_series = zeros(n_time_wins,n_shuffle_reps,2);
        
        for shuffle_rep_k = 1:n_shuffle_reps
            all_shuffle_time_series(:,shuffle_rep_k,1) = circshift(actual_time_series(:,1),randi(n_time_wins));
            all_shuffle_time_series(:,shuffle_rep_k,2) = circshift(actual_time_series(:,2),randi(n_time_wins));
        end
        
        if use_trial_avg_flag
            all_shuffle_time_series = all_shuffle_time_series(time_idx,:,:) + permute(repmat(trial_avg_time_series(time_idx,:),1,1,n_shuffle_reps),[1 3 2]);
        else
            all_shuffle_time_series = all_shuffle_time_series(time_idx,:,:);
        end
        
        switch corrType
            
            case 'corrCoef'
                
                R = corr(all_shuffle_time_series(:,:,1),all_shuffle_time_series(:,:,2));
                shuffle_corr(trial_k,:) = diag(R);
                
            case 'corrIndex'
                for shuffle_rep_k = 1:n_shuffle_reps
                    activation = squeeze(all_shuffle_time_series(:,shuffle_rep_k,:))';
                    shuffle_corr(trial_k,shuffle_rep_k,:) = get_corr_index(activation);
                end
        end
    end
end

end

function corr_index = get_corr_index(activation)

n_time_wins = size(activation,2);

b = repmat(sqrt(2)/2,2,1) .* [1; -1];
bInv = -b;
bMat = repmat(b,1,n_time_wins);
bMat_Inv = repmat(bInv,1,n_time_wins);

activation = (activation - mean(activation,2))./vecnorm(activation,2,2);
activation_norm = vecnorm(activation,2,1);
theta1 = acos(dot(activation,bMat)./activation_norm);
theta2 = acos(dot(activation,bMat_Inv)./activation_norm);
corr_index = min(theta1,theta2);

end