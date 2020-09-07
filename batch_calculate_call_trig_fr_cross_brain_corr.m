function batch_calculate_call_trig_fr_cross_brain_corr(vd,varargin)

pnames = {'combine_single_units'};
dflts  = {false};
[combine_single_units] = internal.stats.parseArgs(pnames,dflts,varargin{:});


call_base_dir = fullfile(vd.baseDirs{1},'call_data');
analysis_dir = fullfile(vd.analysisDir{1},'fr_corr_analysis');

overwrite_corr_flag = true;

call_offset = unique(abs(vd.spikeRange));

min_t = -2.75;
win_length_s = 0.5;
max_t = 2.75;
step_length_s = 0.25;
time_bins = [min_t:step_length_s:max_t-win_length_s;(min_t+win_length_s):step_length_s:max_t]';

n_shuffle_reps = 0;
minFR = 0.1;

[all_cell_fr,~,fr_exp_params,frParams] = get_multi_single_unit_spike_trains({vd},'winSize',200e-3,'overlap',170e-3,'combine_single_units_by_tetrode',combine_single_units,'minCalls',1);

all_cell_fr = all_cell_fr{1};
t = fr_exp_params.t{1};
all_used_call_nums = fr_exp_params.all_used_call_nums{1};
all_used_cell_ks = fr_exp_params.all_used_cell_ks{1};
batNums = fr_exp_params.batNums{1};
all_used_exp_dates = fr_exp_params.expDates{1};

assert(length(unique(vd.expDay)) == size(all_cell_fr,2));

T = readtable(fullfile(vd.baseDirs{1},'documents','recording_logs.csv'));

if strcmp(vd.callType,'call')
    T = T(T.usable & strcmp(T.Session,'communication'),:);
elseif strcmp(vd.callType,'operant')
    T = T(T.usable & strcmp(T.Session,'operant'),:);
end

expDates = T.Date;
nExp = length(expDates);

t1 = tic;
for exp_k = 1:nExp
    
    expDate = expDates(exp_k);
    exp_date_str = datestr(expDate,'yyyymmdd');
    exp_date_idx = all_used_exp_dates == expDate;
    assert(sum(exp_date_idx) == 1)
    switch vd.callType
        case 'call'
            cut_call_fname = fullfile(call_base_dir,[exp_date_str '_cut_call_data.mat']);
            results_fname_str = '_call_trig_fr_corr';
            results_fname = fullfile(analysis_dir,[exp_date_str '.mat']);
            exp_bat_nums = batNums;
        case 'operant'
            boxNum = num2str(T.Box(exp_k));
            cut_call_fname = fullfile(call_base_dir,[exp_date_str '_cut_call_data_operant_box_' boxNum '.mat']);
            results_fname_str = ['_call_trig_operant_box_' boxNum '_fr_corr'];
            exp_bat_nums = cellstr(num2str([T.Bat_1(exp_k); T.Bat_2(exp_k)]));
    end
    
    switch vd.cellType
        case 'multiUnit'
            results_fname = fullfile(analysis_dir,[exp_date_str results_fname_str '_MUA' '.mat']);
        case 'singleUnit'
            if ~combine_single_units
                results_fname = fullfile(analysis_dir,[exp_date_str results_fname_str '_SU' '.mat']);
            else
                results_fname = fullfile(analysis_dir,[exp_date_str results_fname_str '_combined_SU' '.mat']);
            end
            
    end
    
    used_bat_idx = ~cellfun(@isempty,all_cell_fr(:,exp_date_idx)) & ismember(batNums',exp_bat_nums);
    used_bat_nums = batNums(used_bat_idx);
    
    if sum(used_bat_idx) < 2 || (exist(results_fname,'file') && ~overwrite_corr_flag)
        continue
    end
    
    used_call_nums = all_used_call_nums(used_bat_idx,exp_date_idx);
    all_shared_call_nums = multiIntersect(used_call_nums{:});
    
    expFR = all_cell_fr(used_bat_idx,exp_date_idx);
    used_cell_ks = all_used_cell_ks(used_bat_idx,exp_date_idx);
    
    if ~all(cellfun(@(x) isempty(setxor(x,all_shared_call_nums)),used_call_nums))
        if length(all_shared_call_nums) < frParams.minCalls
           continue
        else
            used_call_idx = cellfun(@(x) ismember(x,all_shared_call_nums),used_call_nums,'un',0);
            expFR = cellfun(@(fr,callIdx) fr(:,callIdx,:),expFR,used_call_idx,'un',0);
        end        
    end
    
    [expFR,used_cell_ks,time_idx] = prepare_fr_data_for_corr(expFR,used_cell_ks,t,time_bins,minFR);
    
    expParams.used_call_IDs = all_shared_call_nums;
    expParams.lfp_call_offset = call_offset;
    expParams.n_shuffle_reps = n_shuffle_reps;
    expParams.included_cell_ks = used_cell_ks;
    expParams.time_bins = time_bins;
    expParams.fr_time = t;
    expParams.minFR = minFR;
    expParams.batNums = used_bat_nums;
    
    expFR_jittered = expFR + abs(0.01*rand(size(expFR)));
    [cross_brain_corr,shuffled_corr_p,cross_brain_cohr,cross_brain_corr_index,expParams] = get_cross_brain_corr(expFR_jittered,time_idx,cut_call_fname,expParams);
    
    frResults = struct('expFR',expFR,'cross_brain_corr',cross_brain_corr,...
        'shuffled_corr_p',shuffled_corr_p,'cross_brain_cohr',cross_brain_cohr,...
        'cross_brain_corr_index',cross_brain_corr_index,'expParams',expParams,...
        'frParams',frParams);
    
    save(results_fname,'-v7.3','-struct','frResults');
    
    fprintf('%d sessions processed of %d total',exp_k,length(expDates));
    toc(t1);
end


end

function [activation,used_fr_idx,time_idx] = prepare_fr_data_for_corr(expFR,all_used_cell_ks,t,time_bins,minFR)

used_fr_idx = cellfun(@(x) squeeze(mean(x,[1 2]))>minFR,expFR,'un',0);
expFR = cellfun(@(fr,idx) fr(:,:,idx),expFR,used_fr_idx,'un',0);
used_fr_idx = cellfun(@(cell_ks,idx) cell_ks(idx),all_used_cell_ks,used_fr_idx,'un',0);

maxCells = max(cellfun(@(x) size(x,3),expFR));
nT = size(expFR{1},1);
nTrial = size(expFR{1},2);
nBat = length(expFR);

activation = cell(1,nBat);

for bat_k = 1:nBat
    nanPad = nan(nT,nTrial,maxCells - size(expFR{bat_k},3));
    activation{bat_k} = cat(3,expFR{bat_k},nanPad);
end

activation = permute(cat(4,activation{:}),[4 2 3 1]);

n_time_bins = size(time_bins,1);
time_idx = false(nT,n_time_bins);
for t_k = 1:n_time_bins
    [~, time_idx(:,t_k)] = inRange(t,time_bins(t_k,:));
end

% assert(length(unique(sum(time_idx))) == 1)

end