function [cross_brain_corr,shuffled_corr_p,cross_brain_cohr,cross_brain_corr_index,expParams] = get_cross_brain_corr(activation,time_idx,cut_call_fname,expParams)

n_time_bins = size(time_idx,2);
n_sample_per_time_bin = mode(sum(time_idx));
time_bin_T = unique(round(diff(expParams.time_bins,[],2),3));
dT = time_bin_T/n_sample_per_time_bin;
[cross_brain_corr,shuffled_corr_p,cross_brain_cohr] = deal(cell(1,n_time_bins));
time_axis_dim = 5;

[~,shuffled_corr_pre_calc] = calculate_event_trig_cross_brain_corr(activation,'time_idx',time_idx(:,1),...
    'shuffle_type','randomShuffle','n_shuffle_reps',1e2);

for time_bin_k = 1:n_time_bins
    [cross_brain_corr{time_bin_k},~,shuffled_corr_p{time_bin_k},cross_brain_cohr{time_bin_k}] = ...
        calculate_event_trig_cross_brain_corr(activation,'time_idx',time_idx(:,time_bin_k),...
        'shuffle_type','preCalc','dT',dT,'time_bin_T',time_bin_T,...
        'pre_calc_shuffled_corr',shuffled_corr_pre_calc);
end

cross_brain_corr = cat(time_axis_dim,cross_brain_corr{:});
cross_brain_cohr = cat(time_axis_dim,cross_brain_cohr{:});
shuffled_corr_p = cat(time_axis_dim,shuffled_corr_p{:});

[~,~,~,~,cross_brain_corr_index] = calculate_event_trig_cross_brain_corr(activation,...
    'shuffle_type','none','corr_index_flag',true,'n_shuffle_reps',1e3);

if ~isempty(cut_call_fname)
    
    s = load(cut_call_fname);
    cut_call_data = s.cut_call_data;
    
    [included_bat_nums,included_call_IDs] = get_included_bat_nums(cut_call_data,expParams);
    
    if length(included_bat_nums) ~= size(activation,2)
        disp('mismatch in # of calls')
        keyboard
    end
    
    expParams.included_call_IDs = included_call_IDs;
    expParams.included_bat_nums = included_bat_nums;
end

end

function [included_bat_nums,included_call_IDs] = get_included_bat_nums(cut_call_data,expParams)

all_call_nums = [cut_call_data.uniqueID];
all_bat_nums = {cut_call_data.batNum};
nCall = length(expParams.used_call_IDs);
callPos = vertcat(cut_call_data.corrected_callpos)';

call_offset = 1e3 * expParams.lfp_call_offset .* [-1 1];

included_bat_nums = cell(1,nCall);
included_call_IDs = cell(1,nCall);

for call_k = 1:nCall
    current_call_pos = cut_call_data(all_call_nums == expParams.used_call_IDs(call_k)).corrected_callpos;
    [~,included_call_idx] = inRange(callPos(1,:),current_call_pos + call_offset);
    
    included_bat_nums{call_k} = all_bat_nums(included_call_idx);
    included_call_IDs{call_k} = all_call_nums(included_call_idx);
end

end