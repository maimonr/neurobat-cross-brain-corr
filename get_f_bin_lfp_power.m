function lfp_power = get_f_bin_lfp_power(ps,freqs,f_bins)

nBat = size(ps,1);
nTrial = size(ps,2);
nChannel = size(ps,3);
n_spec_wins = size(ps,4);

nFreq = length(freqs);
f_idx = false(nFreq,2);
n_freq_bins = size(f_bins,1);

for f_k = 1:n_freq_bins
    [~,f_idx(:,f_k)] = inRange(freqs,f_bins(f_k,:));
end

lfp_power = zeros(nBat,nTrial,nChannel,n_spec_wins,n_freq_bins);
for f_k = 1:n_freq_bins
    freq_band_ps = ps(:,:,:,:,f_idx(:,f_k));
    mu = nanmean(freq_band_ps,4);
    sigma = nanstd(freq_band_ps,[],4);
    freq_band_ps_zscore = (freq_band_ps - mu)./sigma;
    lfp_power(:,:,:,:,f_k) = mean(freq_band_ps_zscore,5);
end



end