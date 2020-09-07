function [mua_lfp_corr,used_calls,expDates] = get_lfp_mua_correlation(eData)

baseDir = unique(eData.baseDirs);
baseDir = baseDir{1};

T = readtable(fullfile(baseDir,'documents','recording_logs.csv'));
lfp_analysis_dir = fullfile(baseDir,'data_analysis_results','lfp_data_analysis');
mua_analysis_dir = fullfile(baseDir,'data_analysis_results','fr_corr_analysis');

if strcmp(eData.callType,'call')
    T = T(T.usable & strcmp(T.Session,'communication'),:);
elseif strcmp(eData.callType,'operant')
    T = T(T.usable & strcmp(T.Session,'operant'),:);
end

expDates = T.Date;
nExp = length(expDates);
nBat = length(eData.batNums);

mua_lfp_corr = cell(1,nExp);
used_calls = cell(1,nExp);

exp_k = 1;
while true
    exp_date_str = datestr(expDates(exp_k),'yyyymmdd');
    mua_data_fname = fullfile(mua_analysis_dir,[exp_date_str '_call_trig_fr_corr.mat']);
    lfp_data_fname = fullfile(lfp_analysis_dir,[exp_date_str '_call_trig_ps_corr.mat']);
    if ~(exist(lfp_data_fname,'file') && exist(mua_data_fname,'file'))
        exp_k = exp_k + 1;
        continue
    end
    
    muaData = load(mua_data_fname,'expParams');
    lfpData = load(lfp_data_fname,'expParams');
    mua_fs = round(1/mode(diff(muaData.expParams.fr_time)));
    lfp_fs = round(1/mode(diff(lfpData.expParams.ps_time)));
    break
end

for exp_k = 1:nExp
    exp_date_str = datestr(expDates(exp_k),'yyyymmdd');
    switch eData.callType
        case 'call'
            mua_data_fname = fullfile(mua_analysis_dir,[exp_date_str '_call_trig_fr_corr.mat']);
            lfp_data_fname = fullfile(lfp_analysis_dir,[exp_date_str '_call_trig_ps_corr.mat']);
        case 'operant'
            boxNum = num2str(T.Box(exp_k));
            mua_data_fname = fullfile(mua_analysis_dir,[exp_date_str '_call_trig_operant_box_' boxNum '_fr_corr.mat']);
            lfp_data_fname = fullfile(lfp_analysis_dir,[exp_date_str '_call_trig_operant_box_' boxNum '_ps_corr.mat']);
    end
    
    if ~(exist(lfp_data_fname,'file') && exist(mua_data_fname,'file'))
        continue
    end
    
    lfpData = load(lfp_data_fname,'ps','expParams','n_call_artifact_times','specParams');
    muaData = load(mua_data_fname,'expFR','expParams');
    
    n_exp_bat = size(muaData.expFR,1);
    nTT = size(muaData.expFR,3);
    n_ch_per_tt = 4;
    n_csc_samp = 16665;
    
    [used_calls{exp_k},trial_idx_MUA,trial_idx_LFP] = intersect(muaData.expParams.used_call_IDs,lfpData.expParams.used_call_IDs);
    nTrial = length(trial_idx_MUA);
    
    [~,bat_idx_MUA,bat_idx_LFP] = intersect([lfpData.expParams.batNums{:}],muaData.expParams.batNums);
    batNums = muaData.expParams.batNums(bat_idx_MUA);
    
    muaData.expFR = muaData.expFR(bat_idx_MUA,trial_idx_MUA,:,:);
    lfpData.ps = lfpData.ps(bat_idx_LFP,trial_idx_LFP,:,:,:);
    
    fr_t_range = [min(muaData.expParams.fr_time) max(muaData.expParams.fr_time)];
    ps_t_range = [min(lfpData.expParams.ps_time) max(lfpData.expParams.ps_time)];
    
    t_range = [max(fr_t_range(1),ps_t_range(1)) min(fr_t_range(2),ps_t_range(2))];
    
    [fr_time,fr_t_idx] = inRange(muaData.expParams.fr_time,t_range+[-eps 0]);
    [ps_time,ps_t_idx] = inRange(lfpData.expParams.ps_time,t_range+[-eps 0]);
    
    artifact_removed_ps = prepare_ps_data_for_corr(lfpData.ps,lfpData.n_call_artifact_times,lfpData.expParams,lfpData.specParams,n_csc_samp);
    activation = get_f_bin_lfp_power(artifact_removed_ps,lfpData.specParams.freqs,lfpData.expParams.f_bins(3,:));
    
    expFR = muaData.expFR(:,:,:,fr_t_idx);
    activation = activation(:,:,:,ps_t_idx);
    
    mua_lfp_corr{exp_k} = nan(nBat,nTrial,nTT);
    for bat_k = 1:n_exp_bat
        bat_idx = strcmp(eData.batNums,batNums(bat_k));
        for tt_k = 1:nTT
            for trial_k = 1:nTrial
                chIdx = 1+(tt_k-1)*n_ch_per_tt:tt_k*n_ch_per_tt;
                chIdx = ismember(eData.activeChannels{bat_idx},chIdx);
                X = squeeze(expFR(bat_k,trial_k,tt_k,:));
                X = interp1(fr_time,X,ps_time);
                
                Y = squeeze(activation(bat_k,trial_k,chIdx,:))';
                
                mua_lfp_corr{exp_k}(bat_idx,trial_k,tt_k) = nanmean(corr(X,Y));
            end
        end
    end
end
end