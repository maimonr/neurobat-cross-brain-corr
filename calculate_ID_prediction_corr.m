function bat_id_pred_corr = calculate_ID_prediction_corr(bat_pair_corr,cdList)

n_boot_rep = [1e3 1e4];
mdlType = 'glm_fit_log';
timeWin = [-1 1];
vdStrs = {'adult','adult_operant_comm'};

t = cell(1,2);
for vd_k = 1:2
    t{vd_k} = bat_pair_corr.(vdStrs{vd_k}).lfp(1).time;
end
nT = cellfun(@length,t);
[~,sortIdx] = sort(nT);
t_idx{sortIdx(2)} = ismember(t{sortIdx(2)},t{sortIdx(1)});
t_idx{sortIdx(1)} = ismember(t{sortIdx(1)},t{sortIdx(1)});
t = t{2};

f_k = 3;

excluded_bat_num = '71360';
excluded_bat_idx = cell(1,2);
for vd_k = 1:2
    excluded_bat_idx{vd_k} = ~any(ismember(bat_pair_corr.(vdStrs{vd_k}).lfp(1).all_bat_pairs,excluded_bat_num),2);
end

cross_brain_corr = cell(2);
for vd_k = 1:2
    for call_type_k = 1:2
        cross_brain_corr{call_type_k,vd_k} = squeeze(bat_pair_corr.(vdStrs{vd_k}).lfp(call_type_k).bat_pair_corr(:,excluded_bat_idx{vd_k},f_k,t_idx{vd_k}));
    end
end

bat_id_pred_corr = cell(1,2);
for vd_k = 1:2
    all_cross_brain_corr = nan(size(cross_brain_corr{1,vd_k}));
    
    for call_type_k = 1:2
        idx = ~isnan(cross_brain_corr{call_type_k,vd_k});
        all_cross_brain_corr(idx) = cross_brain_corr{call_type_k,vd_k}(idx);
    end
    
    all_bat_pairs = bat_pair_corr.(vdStrs{vd_k}).lfp(1).all_bat_pairs(excluded_bat_idx{vd_k},:);
    callNums = bat_pair_corr.(vdStrs{vd_k}).lfp(1).all_included_call_nums;
    cData = cdList(vd_k);
    included_bat_nums = cellfun(@(x) char(mode(categorical(cellflat(cData('callID',x').batNum)))),callNums,'UniformOutput',0);
    expDates = bat_pair_corr.(vdStrs{vd_k}).lfp(1).expDates;
    bat_id_pred_corr{vd_k} = predict_bat_id_corr(all_cross_brain_corr,all_bat_pairs,included_bat_nums,expDates,t,'n_boot_rep',n_boot_rep,'mdlType',mdlType,'timeWin',timeWin);
end