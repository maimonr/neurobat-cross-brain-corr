function batch_calculate_lfp_cross_brain_corr(eData)

baseDir = unique(eData.baseDirs);
baseDir = baseDir{1};
expType = unique(eData.expType);
expType = expType{1};
overwrite_ps_flag = false;
overwrite_corr_flag = true;

if strcmp(expType,'test')
    
    lfp_base_dir = fullfile(baseDir,'lfp_data');
    nExp = 1;
    fs = 2083;
    
elseif any(strcmp(expType,{'adult_operant','adult'}))
    
    lfp_base_dir = fullfile(baseDir,'lfp_data');
    call_base_dir = fullfile(baseDir,'call_data');
    analysis_dir = fullfile(baseDir,'data_analysis_results','lfp_data_analysis');
    
    T = readtable(fullfile(baseDir,'documents','recording_logs.csv'));
    csc_var_name = 'call_trig_csc';
    
    if strcmp(eData.callType,'call')
        T = T(T.usable & strcmp(T.Session,'communication'),:);
    elseif strcmp(eData.callType,'operant')
        T = T(T.usable & strcmp(T.Session,'operant'),:);
    elseif strcmp(eData.callType,'playback')
        csc_var_name = 'playback_csc';
    end
    
    expDates = T.Date;
    nExp = length(expDates);
    fs = 2083;
    
elseif strcmp(expType,'juvenile')
    
    analysis_dir = 'E:\ephys\juvenile_recording\data_analysis_results\lfp_data_analysis\';
    
    baseDir = 'E:\ephys\juvenile_recording\';
    all_lfp_fnames = dir(fullfile(baseDir,'bat*','neurologger_recording*','lfpformat','LFP.mat'));
    
    nExp = length(all_lfp_fnames);
    fs = 1953;
    csc_var_name = 'call_trig_csc';
end
f_bins = [5 20;20 70;70 150];

min_t = -2.5;
win_length_s = 0.5;
max_t = 2.5;
step_length_s = 0.25;
time_bins = [min_t:step_length_s:max_t-win_length_s;(min_t+win_length_s):step_length_s:max_t]';
max_artifact_frac = 0.01;
n_shuffle_reps = 5;
n_csc_samp = [];
t1 = tic;

for session_k = 1:nExp
    
    if strcmp(expType,'test')
        
        results_fname = 'E:\ephys\adult_recording\data_analysis_results\lfp_data_analysis\test_call_trig_ps.mat';
        s = load(results_fname,'expParams');
        expParams = s.expParams;
        n_csc_samp = 1+expParams.lfp_call_offset*expParams.fs*2;
        calculate_corr_flag = true;
        overwrite_ps_flag = false;
        cut_call_fname = [];
        
    elseif any(strcmp(expType,{'adult','adult_operant'}))
        expDate = expDates(session_k);
        exp_date_str = datestr(expDate,'yyyymmdd');
        
        if strcmp(eData.callType,'playback')
            cut_call_fname = [];
            event_trig_lfp_fname = [exp_date_str '_LFP_playback.mat'];
            results_fname = fullfile(analysis_dir,[exp_date_str '_playback_ps_corr.mat']);
            event_trig_csc_fnames = dir(fullfile(lfp_base_dir,['*' event_trig_lfp_fname]));
        else
            switch eData.callType
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
            
        end
        
        if isempty(event_trig_csc_fnames)
            disp('No call trig lfp data')
            continue
        end
        
        calculate_corr_flag = true;
        
        if isempty(n_csc_samp)
            m = matfile(fullfile(event_trig_csc_fnames(1).folder,event_trig_csc_fnames(1).name));
            n_csc_samp = size(m.(csc_var_name),1);
        end
        
    elseif strcmp(expType, 'juvenile')
        
        event_trig_csc_fnames = dir(fullfile(all_lfp_fnames(session_k).folder,all_lfp_fnames(session_k).name));
        calculate_corr_flag = false;
        
        cross_brain_corr = [];
        shuffled_corr = [];
        
        s = load(fullfile(event_trig_csc_fnames(1).folder,event_trig_csc_fnames(1).name),'call_trig_csc_struct');
        warnMsg = lastwarn;
        if strcmp(warnMsg,'Variable ''call_trig_csc_struct'' not found.')
            lastwarn('')
            continue
        end
        
        if isempty(n_csc_samp)
            n_csc_samp = size(s.call_trig_csc_struct.(csc_var_name),1);
        end
        
        batNums = regexp(event_trig_csc_fnames.folder,'bat\d{5}','match');
        batNums{1} = batNums{1}(length('bat')+1:end);
        
        exp_date_str = regexp(event_trig_csc_fnames.folder,'\d{8}','match');
        
        results_fname = fullfile(analysis_dir,[exp_date_str{1} '_'  batNums{1} '_' eventType '_ps_corr.mat']);
        
    end
    
    calculate_ps_flag = true;
    
    if exist(results_fname,'file')
        if ~overwrite_ps_flag && ~overwrite_corr_flag
            disp('cross brain correlation data already saved')
            continue
        end
        
        m = matfile(results_fname);
        
        if overwrite_corr_flag
            if ~isempty(whos(m,'cross_brain_corr')) && numel(m.cross_brain_corr) == 0
                calculate_corr_flag = false;
            end
        else
            calculate_corr_flag = false;
        end
        
        if ~overwrite_ps_flag
            calculate_ps_flag = false;
        end
        
    end
    
    if calculate_ps_flag
        [ps,n_call_artifact_times,specParams,expParams] = calculate_event_trig_lfp_ps(eData,event_trig_csc_fnames);
    end
    
    if calculate_corr_flag
        
        if ~calculate_ps_flag
            lfpData = load(results_fname);
            ps = lfpData.ps;
            specParams = lfpData.specParams;
            n_call_artifact_times = lfpData.n_call_artifact_times;
            expParams = lfpData.expParams;
            clear lfpData
        end
        
        expParams.f_bins = f_bins;
        expParams.fs = fs;
        expParams.time_bins = time_bins;
        expParams.max_artifact_frac = max_artifact_frac;
        expParams.n_shuffle_reps = n_shuffle_reps;
        
        if length(expParams.batNums) <= 1
            disp('fewer than 2 bats, can''t calculate correlation')
            
            if ~calculate_ps_flag
                continue
            end
            
            cross_brain_corr = [];
            shuffled_corr_p= [];
            
        else
            [artifact_removed_ps,time_idx,ps_time] = prepare_ps_data_for_corr(ps,n_call_artifact_times,expParams,specParams,n_csc_samp);
            activation = get_f_bin_lfp_power(artifact_removed_ps,specParams.freqs,expParams.f_bins);
            [cross_brain_corr,shuffled_corr_p,cross_brain_cohr,cross_brain_corr_index,expParams] =...
                get_cross_brain_corr(activation,time_idx,cut_call_fname,expParams);
            expParams.ps_time = ps_time;
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