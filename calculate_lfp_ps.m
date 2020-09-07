function [ps,n_call_artifact_times,specParams] = calculate_lfp_ps(event_trig_csc,csc_var_name,fs)

common_ref_flag = false;

time_half_bandwidth_product = 4;
nSamp = size(event_trig_csc(1).(csc_var_name),1);

winSize = 400;
overlap = 350;

[slepian_sequences,concentrations] = dpss(winSize,time_half_bandwidth_product,2*time_half_bandwidth_product-1);
sliding_win_idx = slidingWin(nSamp,winSize,overlap);
nWin = size(sliding_win_idx,1);
freqs = 5:150;
nFreq = length(freqs);
baseline_t_offset = -1;
call_t_win = [-1 1];
artifact_nStd_factor = 5;

specParams = struct('freqs',freqs,'winSize',winSize,'overlap',overlap,...
    'time_half_bandwidth_product',time_half_bandwidth_product,'common_ref_flag',common_ref_flag,...
    'call_t_win',call_t_win,'artifact_nStd_factor',artifact_nStd_factor,'baseline_t_offset',baseline_t_offset);

nBat = length(event_trig_csc);
nT = size(event_trig_csc(1).(csc_var_name),1);
nTrial = size(event_trig_csc(1).(csc_var_name),2);
nChannel = 16;

ps = nan(nBat,nTrial,nChannel,nWin,nFreq);
t = linspace(-event_trig_csc(1).lfp_call_offset,event_trig_csc(1).lfp_call_offset,nT);
[~,call_t_idx] = inRange(t,call_t_win);
mu = nan(nChannel,nBat);
sigma = nan(nChannel,nBat);
baseline_t_idx = t<baseline_t_offset;
for bat_k = 1:nBat
    bat_csc = event_trig_csc(bat_k).(csc_var_name);
    for channel_k = 1:size(bat_csc,3)
        mu(channel_k,bat_k) = mean(bat_csc(baseline_t_idx,:,channel_k),[1 2]);
        sigma(channel_k,bat_k) = std(bat_csc(baseline_t_idx,:,channel_k),[],[1 2]);
    end
end

mu = abs(mu);

n_call_artifact_times = zeros(nBat,nTrial,nChannel);

for bat_k = 1:nBat
    bat_csc = event_trig_csc(bat_k).(csc_var_name);
    bat_nChannel = size(bat_csc,3);
    for trial_k = 1:nTrial
        if common_ref_flag
            common_ref = repmat(squeeze(median(bat_csc(:,trial_k,:),3)),1,bat_nChannel);
        else
            common_ref = zeros(nT,bat_nChannel);
        end
        bat_trial_ps = cell(1,nChannel);
        bat_trial_csc = squeeze(bat_csc(:,trial_k,:)) - common_ref;
        parfor channel_k = 1:bat_nChannel
            current_csc = bat_trial_csc(:,channel_k);
            n_call_artifact_times(bat_k,trial_k,channel_k) = sum(abs(current_csc(call_t_idx)) > (mu(channel_k,bat_k) + artifact_nStd_factor*sigma(channel_k,bat_k)));
            win_csc = current_csc(sliding_win_idx);
            bat_trial_ps{channel_k} = pmtm(win_csc',slepian_sequences,concentrations,freqs,fs,'adapt','DropLastTaper',false)';
        end
        bat_trial_ps = cat(3,bat_trial_ps{:});
        ps(bat_k,trial_k,1:bat_nChannel,:,:) = permute(bat_trial_ps,[3 1 2]);
    end
end