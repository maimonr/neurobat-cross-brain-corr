function p_over_baseline = get_fraction_above_baseline_by_trial(baseline_corr,bat_pair_corr,all_bat_pairs,expDates)

f_k = 2;
p_over_baseline = cell(1,length(baseline_corr));
lfp_base_dir = 'E:\ephys\adult_recording\data_analysis_results\lfp_data_analysis\';
n_all_bat_pair = size(all_bat_pairs,1);

n_baseline_rep = size(baseline_corr(1).cross_brain_corr,1);

for session_k = 1:length(baseline_corr)
    
    exp_day_str = datestr(baseline_corr(session_k).expDate,'yyyymmdd');
    call_trig_lfp_fname = fullfile(lfp_base_dir,[exp_day_str '_call_trig_ps_corr.mat']);
    
    call_trig_bat_nums = load(call_trig_lfp_fname,'expParams');
    call_trig_bat_nums = call_trig_bat_nums.expParams.batNums;
    
    assert(isequal(baseline_corr(session_k).expParams.batNums,call_trig_bat_nums));
    assert(expDates(session_k) == baseline_corr(session_k).expDate);
    
    n_used_bats = length(call_trig_bat_nums);
    
    bat_pair_idxs = nchoosek(1:n_used_bats ,2);
    used_bat_pair_idx = false(1,n_all_bat_pair);
    
    for bat_pair_k = 1:size(bat_pair_idxs,1)
       bat_idxs = bat_pair_idxs(bat_pair_k,:);
       bat_pair_nums = [call_trig_bat_nums{bat_idxs}];
       used_bat_idx = all(strcmp(all_bat_pairs,repmat(bat_pair_nums,n_all_bat_pair,1)),2);
       used_bat_pair_idx(used_bat_idx) = true;
    end
    
    session_pair_corr = squeeze(bat_pair_corr{session_k}(:,:,f_k,:));
    session_baseline_corr = nan(n_baseline_rep,n_all_bat_pair);
    session_baseline_corr(:,used_bat_pair_idx) = nanmedian(baseline_corr(session_k).cross_brain_corr(:,:,:,f_k),3);
    n_used_baseline_samples = sum(~isnan(session_baseline_corr));
    
    nTrial = size(session_pair_corr,1);
    nT = size(session_pair_corr,3);
    
    session_pair_corr = repmat(session_pair_corr,[ones(1,3) n_baseline_rep]);
    session_baseline_corr = repmat(session_baseline_corr,[ones(1,2) nTrial nT]);
    session_baseline_corr = permute(session_baseline_corr,[3 2 4 1]);
    
    missing_trial_idx = all(isnan(session_pair_corr),4);
    n_used_baseline_samples = repmat(n_used_baseline_samples,[nTrial 1 nT]);
    
    p_over_baseline{session_k} = sum(session_pair_corr>session_baseline_corr,4)./n_used_baseline_samples;
    p_over_baseline{session_k}(missing_trial_idx) = NaN;
end