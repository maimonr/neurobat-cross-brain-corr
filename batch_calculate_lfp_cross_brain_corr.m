function batch_calculate_lfp_cross_brain_corr(eData,varargin)

pnames = {'used_exp_dates','callType'};
dflts  = {[],'call'};
[used_exp_dates,callType] = internal.stats.parseArgs(pnames,dflts,varargin{:});

baseDir = unique(eData.baseDirs);
baseDir = baseDir{1};
expType = unique(eData.expType);
expType = expType{1};
overwrite_ps_flag = false;
overwrite_corr_flag = false;

[lfp_base_dir, call_base_dir, analysis_dir, csc_var_name, nExp, expDates, fs, T] = get_exp_setup(baseDir,expType,callType,used_exp_dates);

f_bins = [5 20;20 70;70 150];

min_t = -2.5;
win_length_s = 0.5;
max_t = 2.5;
step_length_s = 0.25;
time_bins = [min_t:step_length_s:max_t-win_length_s;(min_t+win_length_s):step_length_s:max_t]';
max_artifact_frac = 0.01;
n_shuffle_reps = 1;
corr_shuffle_type = 'none';
t1 = tic;

for session_k = 1:nExp
    
    [n_csc_samp, calculate_corr_flag, overwrite_ps_flag, cut_call_fname, results_fname, event_trig_csc_fnames] = ...
        get_exp_data(session_k,expType,callType,expDates,lfp_base_dir, call_base_dir, analysis_dir, csc_var_name, T, overwrite_ps_flag);
    
    if isnan(n_csc_samp)
        continue
    end
    
    [calculate_corr_flag,calculate_ps_flag,continue_flag] = determine_overwrite_status(results_fname,overwrite_ps_flag,overwrite_corr_flag,calculate_corr_flag);
    
    if continue_flag
        continue
    end
    
    if calculate_ps_flag
        [ps,n_call_artifact_times,specParams,expParams] = calculate_event_trig_lfp_ps(expType,event_trig_csc_fnames);
    else
        [ps,n_call_artifact_times,specParams,expParams] = deal([]);
    end
    
    if calculate_corr_flag
        
        expParams.f_bins = f_bins;
        expParams.fs = fs;
        expParams.time_bins = time_bins;
        expParams.max_artifact_frac = max_artifact_frac;
        expParams.n_shuffle_reps = n_shuffle_reps;
        expParams.corr_shuffle_type = corr_shuffle_type;
        [cross_brain_corr,shuffled_corr_p,cross_brain_cohr,cross_brain_corr_index,expParams,specParams,continue_flag]...
            = get_corr(results_fname,calculate_ps_flag,ps,n_call_artifact_times,expParams,specParams,n_csc_samp,cut_call_fname);
        
        if continue_flag
            continue
        end
    else
        [cross_brain_corr,shuffled_corr_p,cross_brain_cohr,cross_brain_corr_index] = deal([]);
    end
    
    if calculate_ps_flag || calculate_corr_flag
        
        lfpResults = struct('ps',ps,'cross_brain_corr',cross_brain_corr,...
            'shuffled_corr_p',shuffled_corr_p,'cross_brain_cohr',cross_brain_cohr,...
            'cross_brain_corr_index',cross_brain_corr_index,...
            'n_call_artifact_times',n_call_artifact_times,'specParams',specParams,...
            'expParams',expParams);
        
        save(results_fname,'-v7.3','-struct','lfpResults');
    end
    
    fprintf('%d sessions processed of %d total',session_k,nExp);
    toc(t1);
    clear ps corss_brain_corr trial_shuffled_corr time_shuffled_corr n_call_artifact_times specParams expParams lfpResults event_trig_csc_fnames results_fname m
end


end

function [lfp_base_dir, call_base_dir, analysis_dir, csc_var_name, nExp, expDates, fs, T] = get_exp_setup(baseDir,expType,callType,used_exp_dates)
switch expType
    case 'test'
        
        lfp_base_dir = fullfile(baseDir,'lfp_data');
        nExp = 1;
        fs = 2083;
        
    case {'adult_operant','adult','adult_social'}
        
        lfp_base_dir = fullfile(baseDir,'lfp_data');
        call_base_dir = fullfile(baseDir,'call_data');
        analysis_dir = fullfile(baseDir,'data_analysis_results','lfp_data_analysis');
        
        if strcmp(expType,{'adult_operant','adult'})
            T = readtable(fullfile(baseDir,'documents','recording_logs.csv'));
        elseif strcmp(expType,'adult_social')
            T = get_rec_logs;
        end
        csc_var_name = 'call_trig_csc';
        
        switch callType
            case 'call'
                if strcmp(expType,{'adult_operant','adult'})
                    expIdx = T.usable & strcmp(T.Session,'communication');
                elseif strcmp(expType,'adult_social')
                    expIdx = T.usable & strcmp(T.Session,'vocal');
                end
            case 'operant'
                expIdx = T.usable & strcmp(T.Session,'operant');
            case 'playback'
                csc_var_name = 'playback_csc';
                expIdx = true(1,size(T,1));
            case 'social'
                expIdx = T.usable & strcmp(T.Session,'social');
        end
        
        if ~isempty(used_exp_dates)
            idx = ismember(T.Date,used_exp_dates);
            expIdx = idx & expIdx;
        end
        
        T = T(expIdx,:);
        
        expDates = T.Date;
        nExp = length(expDates);
        fs = 2083;
        
    case 'juvenile'
        
        analysis_dir = 'E:\ephys\juvenile_recording\data_analysis_results\lfp_data_analysis\';
        
        baseDir = 'E:\ephys\juvenile_recording\';
        all_lfp_fnames = dir(fullfile(baseDir,'bat*','neurologger_recording*','lfpformat','LFP.mat'));
        
        nExp = length(all_lfp_fnames);
        fs = 1953;
        csc_var_name = 'call_trig_csc';
end
end

function [n_csc_samp, calculate_corr_flag, overwrite_ps_flag, cut_call_fname, results_fname, event_trig_csc_fnames] =...
    get_exp_data(session_k,expType,callType,expDates,lfp_base_dir, call_base_dir, analysis_dir, csc_var_name,T, overwrite_ps_flag)

[n_csc_samp, calculate_corr_flag, cut_call_fname] = deal(NaN);

switch expType
    
    case 'test'
        
        results_fname = 'E:\ephys\adult_recording\data_analysis_results\lfp_data_analysis\test_call_trig_ps.mat';
        s = load(results_fname,'expParams');
        expParams = s.expParams;
        n_csc_samp = 1+expParams.lfp_call_offset*expParams.fs*2;
        calculate_corr_flag = true;
        overwrite_ps_flag = false;
        cut_call_fname = [];
        
    case {'adult','adult_operant','adult_social'}
        expDate = expDates(session_k);
        exp_date_str = datestr(expDate,'yyyymmdd');
        
        switch callType
            case 'playback'
                cut_call_fname = [];
                event_trig_lfp_fname = [exp_date_str '_LFP_playback.mat'];
                results_fname = fullfile(analysis_dir,[exp_date_str '_playback_ps_corr.mat']);
                event_trig_csc_fnames = dir(fullfile(lfp_base_dir,['*' event_trig_lfp_fname]));
            case 'call'
                cut_call_fname = fullfile(call_base_dir,[exp_date_str '_cut_call_data.mat']);
                event_trig_lfp_fname = [exp_date_str '_LFP_call_trig.mat'];
                results_fname = fullfile(analysis_dir,[exp_date_str '_call_trig_ps_corr.mat']);
                event_trig_csc_fnames = dir(fullfile(lfp_base_dir,['*' event_trig_lfp_fname]));
            case 'operant'
                boxNum =  num2str(T.Box(session_k));
                cut_call_fname = fullfile(call_base_dir,[exp_date_str '_cut_call_data_operant_box_' boxNum '.mat']);
                event_trig_lfp_fname = [exp_date_str '_LFP_call_trig_operant.mat'];
                results_fname = fullfile(analysis_dir,[exp_date_str '_call_trig_operant_box_' boxNum '_ps_corr.mat']);
                event_trig_csc_fnames = dir(fullfile(lfp_base_dir,['*' event_trig_lfp_fname]));
                
                batNums = cellstr(num2str([T.Bat_1(session_k); T.Bat_2(session_k)]));
                session_bat_idx = false(1,length(event_trig_csc_fnames));
                for bat_k = 1:length(event_trig_csc_fnames)
                    session_bat_idx(bat_k) = contains(event_trig_csc_fnames(bat_k).name,batNums);
                end
                event_trig_csc_fnames = event_trig_csc_fnames(session_bat_idx);
        end
        
        
        
        if isempty(event_trig_csc_fnames)
            disp('No call trig lfp data')
            return
        end
        
        calculate_corr_flag = true;
        
        if isnan(n_csc_samp)
            m = matfile(fullfile(event_trig_csc_fnames(1).folder,event_trig_csc_fnames(1).name));
            n_csc_samp = size(m.(csc_var_name),1);
        end
        
    case 'juvenile'
        
        event_trig_csc_fnames = dir(fullfile(all_lfp_fnames(session_k).folder,all_lfp_fnames(session_k).name));
        calculate_corr_flag = false;
        
        s = load(fullfile(event_trig_csc_fnames(1).folder,event_trig_csc_fnames(1).name),'call_trig_csc_struct');
        warnMsg = lastwarn;
        if strcmp(warnMsg,'Variable ''call_trig_csc_struct'' not found.')
            lastwarn('')
            return
        end
        
        if isempty(n_csc_samp)
            n_csc_samp = size(s.call_trig_csc_struct.(csc_var_name),1);
        end
        
        batNums = regexp(event_trig_csc_fnames.folder,'bat\d{5}','match');
        batNums{1} = batNums{1}(length('bat')+1:end);
        
        exp_date_str = regexp(event_trig_csc_fnames.folder,'\d{8}','match');
        
        results_fname = fullfile(analysis_dir,[exp_date_str{1} '_'  batNums{1} '_' eventType '_ps_corr.mat']);
        
end
end

function [calculate_corr_flag,calculate_ps_flag,continue_flag] = determine_overwrite_status(results_fname,overwrite_ps_flag,overwrite_corr_flag,calculate_corr_flag)

calculate_ps_flag = true;
continue_flag = false;

if exist(results_fname,'file')
    if ~overwrite_ps_flag && ~overwrite_corr_flag
        disp('cross brain correlation data already saved')
        continue_flag = true;
        return
    end
    
    m = matfile(results_fname);
    
    if overwrite_corr_flag
        if ~isempty(m.cross_brain_corr)
            calculate_corr_flag = false;
        end
    else
        calculate_corr_flag = false;
    end
    
    if ~overwrite_ps_flag
        calculate_ps_flag = false;
    end
    
end

end

function [cross_brain_corr,shuffled_corr_p,cross_brain_cohr,cross_brain_corr_index,expParams,specParams,continue_flag]...
    = get_corr(results_fname,calculate_ps_flag,ps,n_call_artifact_times,expParams,specParams,n_csc_samp,cut_call_fname)

continue_flag = false;

if ~calculate_ps_flag
    lfpData = load(results_fname);
    ps = lfpData.ps;
    specParams = lfpData.specParams;
    n_call_artifact_times = lfpData.n_call_artifact_times;
    expParams = lfpData.expParams;
    clear lfpData
end

if length(expParams.batNums) <= 1
    disp('fewer than 2 bats, can''t calculate correlation')
    
    if ~calculate_ps_flag
        continue_flag = true;
    end
    
    cross_brain_corr = [];
    shuffled_corr_p = [];
    
else
    [artifact_removed_ps,time_idx,ps_time] = prepare_ps_data_for_corr(ps,n_call_artifact_times,expParams,specParams,n_csc_samp);
    activation = get_f_bin_lfp_power(artifact_removed_ps,specParams.freqs,expParams.f_bins);
    [cross_brain_corr,shuffled_corr_p,cross_brain_cohr,cross_brain_corr_index,expParams] =...
        get_cross_brain_corr(activation,time_idx,cut_call_fname,expParams);
    expParams.ps_time = ps_time;
end
end