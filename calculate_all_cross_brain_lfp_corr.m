function bat_pair_corr_info = calculate_all_cross_brain_lfp_corr(eData,include_trial_flag,varargin)

pnames = {'concatenate_sessions_flag','select_bat_nums','dataType','combine_single_units','callType'};
dflts  = {true,[],'lfp',false,'call'};
[concatenate_sessions_flag,select_bat_nums,dataType,combine_single_units,callType] = internal.stats.parseArgs(pnames,dflts,varargin{:});

baseDir = unique(eData.baseDirs);

switch dataType
    case 'lfp'
        baseDir = fullfile(baseDir{1},'data_analysis_results','lfp_data_analysis');
        nChannel = 16;
        event_trig_str = 'ps_corr';
    case 'mua'
        baseDir = fullfile(baseDir{1},'data_analysis_results','fr_corr_analysis');
        event_trig_str = 'fr_corr_MUA';
    case 'singleUnit'
        baseDir = fullfile(baseDir{1},'data_analysis_results','fr_corr_analysis');
        if ~combine_single_units
            event_trig_str = 'fr_corr_SU';
        else
            event_trig_str = 'fr_corr_combined_SU';
        end
end


switch eData.expType
    case 'adult'
        nBats_per_exp = 3;
        nBats_total = 4;
    case 'adult_operant'
        nBats_per_exp = 4;
        nBats_total = 4;
    case 'adult_social'
        nBats_per_exp = 4;
        nBats_total = 5;
end

switch callType
    case 'playback'
        event_trig_corr_fnames = dir(fullfile(baseDir,['*playback_' event_trig_str '.mat']));
        if strcmp(eData.expType,'adult')
            nBats_per_exp = 2;
            nBats_total = 2;
        end
    case 'call'
        event_trig_corr_fnames = dir(fullfile(baseDir,['*call_trig_' event_trig_str '.mat']));
    case 'operant'
        event_trig_corr_fnames = dir(fullfile(baseDir,['*call_trig_operant_box_*' event_trig_str '.mat']));
        nBats_per_exp = 2;
        nBats_total = 4;
    case 'social'
        event_trig_corr_fnames = dir(fullfile(baseDir,['*call_trig_social_' event_trig_str '.mat']));
    case 'test'
        nChannel = 4;
        event_trig_corr_fnames = dir(fullfile(baseDir,'test_call_trig_ps.mat'));
        nBats_per_exp = 3;
        nBats_total = 3;
end

n_exp_day = length(event_trig_corr_fnames);

bat_pairs_idx = nchoosek(1:nBats_per_exp,2);
all_bat_pairs_idx = nchoosek(1:nBats_total,2);

n_bat_pairs = size(bat_pairs_idx,1);
n_all_bat_pairs = size(all_bat_pairs_idx,1);

exp_bat_pairs = cell(n_bat_pairs,2,n_exp_day);
use_exp_date = true(1,n_exp_day);
for exp_k = 1:n_exp_day
    s = load(fullfile(baseDir,event_trig_corr_fnames(exp_k).name),'expParams');
    nBats_per_exp = length(s.expParams.batNums);
    if nBats_per_exp < 2
        use_exp_date(exp_k) = false;
        continue
    end
    bat_pairs_idx = nchoosek(1:nBats_per_exp,2);
    n_bat_pairs = size(bat_pairs_idx,1);
    for pair_k = 1:n_bat_pairs
        for bat_pair_k = 1:2
            bNum = s.expParams.batNums{bat_pairs_idx(pair_k,bat_pair_k)};
            if iscell(bNum)
                exp_bat_pairs(pair_k,bat_pair_k,exp_k) = bNum;
            else
                exp_bat_pairs{pair_k,bat_pair_k,exp_k} = bNum;
            end
        end
    end
end

n_exp_day = sum(use_exp_date);
event_trig_corr_fnames = event_trig_corr_fnames(use_exp_date);
exp_bat_pairs = exp_bat_pairs(:,:,use_exp_date);

all_bat_nums = unique(exp_bat_pairs(cellfun(@isstr,exp_bat_pairs)));
all_bat_pairs = cell(size(all_bat_pairs_idx));
for pair_k = 1:n_all_bat_pairs
    for bat_pair_k = 1:2
        all_bat_pairs(pair_k,bat_pair_k) = all_bat_nums(all_bat_pairs_idx(pair_k,bat_pair_k));
    end
end

n_time_bins = size(s.expParams.time_bins,1);
time_bin_T = unique(round(diff(s.expParams.time_bins,[],2),3));

if strcmp(dataType,'lfp')
    
    used_channels_by_bat = cell(1,nBats_total);
    channel_pairs = cell(1,2);
    [channel_pairs{1},channel_pairs{2}] = meshgrid(1:nChannel,1:nChannel);
    channel_pairs = cellfun(@(x) reshape(x,1,[]),channel_pairs,'un',0);
    channel_pairs = vertcat(channel_pairs{:})';
    
    for bat_k = 1:nBats_total
        used_channels_by_bat{bat_k} = eData.activeChannels{strcmp(eData.batNums,all_bat_nums{bat_k})};
    end
    
    nT = length(s.expParams.ps_time);
    n_sample_per_time_bin = round(time_bin_T/mode(diff(s.expParams.ps_time)));
    
elseif any(strcmp(dataType,{'mua','singleUnit'}))
    nT = length(s.expParams.fr_time);
    n_sample_per_time_bin = round(time_bin_T/mode(diff(s.expParams.fr_time)));
end

dT = time_bin_T/n_sample_per_time_bin;
df = 1/time_bin_T; %Determine the frequency resolution.
fNQ = 1/ dT / 2; %Determine the Nyquist frequency.
f = (0:df:fNQ); %Construct frequency axis
nFreq = length(f);

[cross_brain_corr,cross_brain_cohr,cross_brain_corr_index,shuffled_corr_p,...
    all_included_call_nums,time] = deal(cell(1,n_exp_day));
expDates = datetime([],[],[]);

for exp_k = 1:n_exp_day
    
    exp_day_str = regexp(event_trig_corr_fnames(exp_k).name,'^\d{8}','match');
    if ~isempty(exp_day_str)
        expDates(exp_k) = datetime(exp_day_str{1},'InputFormat','yyyyMMdd');
    else
        expDates(exp_k) = NaT;
    end
    
    corrData = load(fullfile(baseDir,event_trig_corr_fnames(exp_k).name),'cross_brain_corr','cross_brain_cohr','cross_brain_corr_index','expParams','shuffled_corr_p');
    
    time{exp_k} = corrData.expParams.time_bins(:,1)';
    
    if strcmp(dataType,'lfp')
        n_f_band = size(corrData.expParams.f_bins,1);
    elseif any(strcmp(dataType,{'mua','singleUnit'}))
        n_f_band = 1;
    end
    nTrial = size(corrData.cross_brain_corr,1);
    n_bat_pairs = size(corrData.cross_brain_corr,2);
    nBats_per_exp = length(corrData.expParams.batNums);
    bat_pairs_idx = nchoosek(1:nBats_per_exp,2);
    
    nChannel_pairs = size(corrData.cross_brain_corr,3);
    [cross_brain_corr{exp_k},shuffled_corr_p{exp_k}] = deal(nan(nTrial,n_bat_pairs,n_f_band,n_time_bins));
    
    if all(isnan(corrData.cross_brain_cohr))
        get_cohr_flag = false;
    else
        get_cohr_flag = true;
        cross_brain_cohr{exp_k} = nan(nFreq,n_bat_pairs,n_f_band,n_time_bins);
    end
    
    if all(isnan(corrData.cross_brain_corr_index))
        get_index_flag = false;
    else
        get_index_flag = true;
        cross_brain_corr_index{exp_k} = nan(nTrial,n_bat_pairs,n_f_band,nT);
    end
    
    if ~strcmp(callType,'playback')
        all_included_call_nums{exp_k} = corrData.expParams.included_call_IDs;
    end
    
    used_corr_idx = false(n_bat_pairs,nTrial);
    for bat_pair_k = 1:n_bat_pairs
        bat_pair_nums = exp_bat_pairs(bat_pair_k,:,exp_k);
        used_channel_idx = true(nChannel_pairs,1);
        
        if strcmp(dataType,'lfp')
            for bat_k = 1:2
                bat_used_channel_idx = used_channels_by_bat{strcmp(bat_pair_nums{bat_k},all_bat_nums)};
                used_channel_idx(~ismember(channel_pairs(:,bat_k),bat_used_channel_idx)) = false;
            end
        end
        
        batPair = corrData.expParams.batNums(bat_pairs_idx(bat_pair_k,:));
        for trial_k = 1:nTrial
            use_trial = false;
            if strcmp(callType,'playback')
                use_trial = true;
            else
                trial_bat_nums = corrData.expParams.included_bat_nums{trial_k};
                multiple_bats_idx = cellfun(@iscell,trial_bat_nums);
                if any(multiple_bats_idx)
                    trial_bat_nums = [trial_bat_nums{multiple_bats_idx} trial_bat_nums(~multiple_bats_idx)];
                end
                switch include_trial_flag
                    case 'all'
                        use_trial = true;
                    case 'either_included'
                        use_trial = contains(batPair{1},trial_bat_nums) | contains(batPair{2},trial_bat_nums);
                    case 'only_included'
                        use_trial = contains(batPair{1},trial_bat_nums) & contains(batPair{2},trial_bat_nums);
                    case 'neither_included'
                        use_trial = ~(contains(batPair{1},trial_bat_nums) | contains(batPair{2},trial_bat_nums));
                    case 'only_one_included'
                        use_trial = xor(contains(batPair{1},trial_bat_nums),contains(batPair{2},trial_bat_nums));
                    case 'select_bats'
                        use_trial = any(contains(select_bat_nums,trial_bat_nums));
                end
            end
            
            used_corr_idx(bat_pair_k,trial_k) = use_trial && ~all(isnan(corrData.cross_brain_corr(trial_k,bat_pair_k,used_channel_idx,:,:)),'all');
        end
    end
    
    for f_k = 1:n_f_band
        for bat_pair_k = 1:n_bat_pairs
            bat_pair_nums = exp_bat_pairs(bat_pair_k,:,exp_k);
            used_channel_idx = true(nChannel_pairs,1);
            if strcmp(dataType,'lfp')
                for bat_k = 1:2
                    bat_used_channel_idx = used_channels_by_bat{strcmp(bat_pair_nums{bat_k},all_bat_nums)};
                    used_channel_idx(~ismember(channel_pairs(:,bat_k),bat_used_channel_idx)) = false;
                end
            end
            for trial_k = 1:nTrial
                if used_corr_idx(bat_pair_k,trial_k)
                    for time_bin_k = 1:n_time_bins
                        
                        corr_indices = {trial_k,bat_pair_k,used_channel_idx,f_k,time_bin_k};
                        cohr_indices = {':',bat_pair_k,used_channel_idx,f_k,time_bin_k};
                        
                        trial_corr = nanmean(corrData.cross_brain_corr(corr_indices{:}));
                        trial_corr_p = sum(corrData.shuffled_corr_p(corr_indices{:})>0.95)/sum(used_channel_idx);
                        
                        if get_cohr_flag
                            trial_cohr = nanmean(corrData.cross_brain_cohr(cohr_indices{:}),3);
                            cross_brain_cohr{exp_k}(cohr_indices{[1 2 4:end]}) = trial_cohr;
                        end
                        
                        cross_brain_corr{exp_k}(corr_indices{[1 2 4:end]}) = trial_corr;
                        shuffled_corr_p{exp_k}(corr_indices{[1 2 4:end]}) = trial_corr_p;
                    end
                    if get_index_flag
                        for time_bin_k = 1:nT
                            corr_index_indices = {trial_k,bat_pair_k,used_channel_idx,f_k,time_bin_k};
                            trial_corr_index = nanmean(corrData.cross_brain_corr_index(corr_index_indices{:}));
                            cross_brain_corr_index{exp_k}(corr_index_indices{[1 2 4:end]}) = trial_corr_index;
                        end
                    end
                end
            end
        end
    end
end

[bat_pair_corr,bat_pair_cohr,bat_pair_corr_index,bat_pair_shuffled_corr_p] = deal(cell(1,n_exp_day));

for exp_k = 1:n_exp_day
    nTrial = size(cross_brain_corr{exp_k},1);
    n_bat_pairs = size(cross_brain_corr{exp_k},2);
    
    [bat_pair_corr{exp_k},bat_pair_shuffled_corr_p{exp_k}] = deal(nan(nTrial,n_all_bat_pairs,n_f_band,n_time_bins));
    
    if get_cohr_flag
        bat_pair_cohr{exp_k} = nan(nFreq,n_all_bat_pairs,n_f_band,n_time_bins);
    end
    
    if get_index_flag
        bat_pair_corr_index{exp_k} = nan(nTrial,n_all_bat_pairs,n_f_band,nT);
    end
    
    for bat_pair_k = 1:n_bat_pairs
        idx = all(ismember(all_bat_pairs,squeeze(exp_bat_pairs(bat_pair_k,:,exp_k))),2);
        
        bat_corr_indices = {':',bat_pair_k,':',':'};
        bat_pair_corr_indices = {':',idx,':',':'};
        
        bat_pair_corr{exp_k}(bat_pair_corr_indices{:}) = cross_brain_corr{exp_k}(bat_corr_indices{:});
        bat_pair_shuffled_corr_p{exp_k}(bat_pair_corr_indices{:}) = shuffled_corr_p{exp_k}(bat_corr_indices{:});
        if get_cohr_flag
            bat_pair_cohr{exp_k}(bat_pair_corr_indices{:}) = cross_brain_cohr{exp_k}(bat_corr_indices{:});
        end
        if get_index_flag
            bat_pair_corr_index{exp_k}(bat_pair_corr_indices{:}) = cross_brain_corr_index{exp_k}(bat_corr_indices{:});
        end
    end
end

if concatenate_sessions_flag
    all_exp_dates = cell(1,n_exp_day);
    for exp_k = 1:n_exp_day
        all_exp_dates{exp_k} = repmat(expDates(exp_k),size(bat_pair_corr{exp_k},1),1);
    end
    bat_pair_corr = cat(1,bat_pair_corr{:});
    bat_pair_cohr = permute(cat(5,bat_pair_cohr{:}),[5 1:4]);
    bat_pair_corr_index = cat(1,bat_pair_corr_index{:});
    bat_pair_shuffled_corr_p = cat(1,bat_pair_shuffled_corr_p{:});
    all_included_call_nums = [all_included_call_nums{:}];
    expDates = vertcat(all_exp_dates{:});
end

time = vertcat(time{:});
time = unique(time,'rows');

assert(size(time,1) == 1)

bat_pair_corr_info = struct('bat_pair_corr',bat_pair_corr,...
    'bat_pair_shuffled_corr_p',bat_pair_shuffled_corr_p,'bat_pair_cohr',bat_pair_cohr,...
    'bat_pair_corr_index',bat_pair_corr_index,'all_included_call_nums',{all_included_call_nums},...
    'all_bat_pairs',{all_bat_pairs},'expDates',expDates,'callType',callType,...
    'includedCalls',include_trial_flag,'time',time);

end