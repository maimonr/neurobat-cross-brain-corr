function [artifact_removed_ps,time_idx,ps_time] = prepare_ps_data_for_corr(ps,n_call_artifact_times,expParams,specParams,n_csc_samp)

n_time_bins = size(expParams.time_bins,1);
lfp_call_time_s = abs(diff(expParams.call_t_win));
call_t_length = expParams.fs*lfp_call_time_s;
max_n_artifact = call_t_length*expParams.max_artifact_frac;
artifact_trial_idx = any(n_call_artifact_times>max_n_artifact,3);
artifact_removed_ps = ps;
n_spec_bins = size(ps,4);

for bat_k = 1:size(ps,1)
    artifact_removed_ps(bat_k,artifact_trial_idx(bat_k,:),:,:,:) = NaN;
end

csc_time = linspace(-expParams.lfp_call_offset,expParams.lfp_call_offset,n_csc_samp);
sliding_win_idx = slidingWin(n_csc_samp,specParams.winSize,specParams.overlap);
ps_time = mean(csc_time(sliding_win_idx),2);
time_idx = false(n_spec_bins,n_time_bins);
for t_k = 1:n_time_bins
    [~, time_idx(:,t_k)] = inRange(ps_time,expParams.time_bins(t_k,:));
end

end