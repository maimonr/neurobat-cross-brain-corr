function [p,anovaTbl,stats] = calculate_interbrain_anova(bat_pair_corr)

callRange = [-0.3 0.3];
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

[~,callIdx] = inRange(t,callRange);

avgCorr = cellfun(@(x) mean(x(:,:,callIdx),3),cross_brain_corr,'un',0);
n_bat_pairs = sum(cellfun(@(x) size(x,2),avgCorr(1,:)));
avg_corr_by_bat = cell(2,n_bat_pairs);

for call_type_k = 1:2
    k = 1;
    for vd_k = 1:2
        for bat_k = 1:size(avgCorr{call_type_k,vd_k},2)
            avg_corr_by_bat{call_type_k,k} = avgCorr{call_type_k,vd_k}(:,bat_k);
            avg_corr_by_bat{call_type_k,k}  = avg_corr_by_bat{call_type_k,k}(~isnan(avg_corr_by_bat{call_type_k,k}));
            k = k + 1;
        end
    end
end

all_avg_corr = [vertcat(avg_corr_by_bat{1,:});vertcat(avg_corr_by_bat{2,:})];

prod_perceive_factor = [ones(sum(cellfun(@length,avg_corr_by_bat(1,:))),1); 2*ones(sum(cellfun(@length,avg_corr_by_bat(2,:))),1)];
bat_id_factor = cell(1,2);
for call_type_k = 1:2
    bat_id_factor{call_type_k} = cellfun(@(bat,x) bat*ones(length(x),1),num2cell(1:n_bat_pairs),avg_corr_by_bat(call_type_k,:),'un',0); 
end
bat_id_factor = vertcat(bat_id_factor{:});
bat_id_factor = [vertcat(bat_id_factor{1,:});vertcat(bat_id_factor{2,:})];

exp_group_factor = bat_id_factor <= size(avgCorr{1},2);

nested = [0 0 0; 0 0 1; 0 0 0];

[p,anovaTbl,stats] = anovan(all_avg_corr,{prod_perceive_factor,bat_id_factor,exp_group_factor},'display','off','nested',nested);

end