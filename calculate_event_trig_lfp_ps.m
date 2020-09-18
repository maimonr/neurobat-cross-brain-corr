function [ps,n_call_artifact_times,specParams,expParams] = calculate_event_trig_lfp_ps(expType,event_trig_fnames,varargin)

if ~isempty(varargin)
    callType = varargin;
else
    callType = 'call';
end

switch expType
    
    case {'adult','adult_operant','adult_social'}
        
        batNums = arrayfun(@(x) regexp(x.name,'^\d{5}','match'),event_trig_fnames,'un',0);
        batNums = batNums(~cellfun(@isempty,batNums));
        assert(length(batNums) >= 1 && length(batNums) <= 5)
        
        for bat_k = 1:length(batNums)
            event_trig_csc(bat_k) = load(fullfile(event_trig_fnames(bat_k).folder,event_trig_fnames(bat_k).name)); %#ok<AGROW>
        end
        
        [event_trig_csc.batNum] = deal(batNums{:});
        
        fs = 2083;
        
    case 'juvenile'
        
        batNums = regexp(event_trig_fnames.folder,'bat\d{5}','match');
        
        assert(length(batNums) == 1)
        
        batNums{1} = batNums{1}(length('bat')+1:end);
        event_trig_csc = load(fullfile(event_trig_fnames.folder,event_trig_fnames.name),'call_trig_csc_struct');
        event_trig_csc = event_trig_csc.call_trig_csc_struct;
        event_trig_csc.batNums = batNums{1};
        
        fs = 1953;
end

switch callType
    case 'call'
        
        csc_var_name = 'call_trig_csc';
        
        if length(batNums) > 1
            
            used_call_IDs = multiIntersect(event_trig_csc(:).used_call_IDs);
            
            for bat_k = 1:length(batNums)
                keep_used_call_idx = ismember(event_trig_csc(bat_k).used_call_IDs,used_call_IDs);
                current_used_call_IDs = event_trig_csc(bat_k).used_call_IDs(keep_used_call_idx);
                [~,sortIdx] = sort(current_used_call_IDs);
                event_trig_csc(bat_k).call_trig_csc = event_trig_csc(bat_k).call_trig_csc(:,keep_used_call_idx,:);
                event_trig_csc(bat_k).call_trig_csc = event_trig_csc(bat_k).call_trig_csc(:,sortIdx,:); 
            end
        else
            used_call_IDs = event_trig_csc.used_call_IDs;
        end
        
        used_playback_timestamps = [];
        used_playback_durations = [];
        
    case 'playback'
        
        csc_var_name = 'playback_csc';
        
        if length(batNums) > 1
            used_playback_timestamps = multiIntersect(event_trig_csc(:).playback_TTL_timestamps);
            
            for bat_k = 1:length(batNums)
                keep_used_call_idx = ismember(event_trig_csc(bat_k).playback_TTL_timestamps,used_playback_timestamps);
                event_trig_csc(bat_k).call_trig_csc = event_trig_csc(bat_k).(csc_var_name)(:,keep_used_call_idx,:);
            end
            used_playback_durations = event_trig_csc(bat_k).playback_TTL_durations(keep_used_call_idx);
        else
            used_playback_timestamps = event_trig_csc.playback_TTL_timestamps;
            used_playback_durations = event_trig_csc.playback_TTL_durations;
        end
        
        used_call_IDs = [];
        
end

[ps,n_call_artifact_times,specParams] = calculate_lfp_ps(event_trig_csc,csc_var_name,fs);

expParams = struct('batNums',{batNums},'used_call_IDs',used_call_IDs,'fs',fs,...
    'call_t_win',specParams.call_t_win,'baseline_t_offset',specParams.baseline_t_offset,...
    'artifact_nStd_factor',specParams.artifact_nStd_factor,'lfp_call_offset',event_trig_csc(1).lfp_call_offset,...
    'used_playback_durations',used_playback_durations,'used_playback_timestamps',used_playback_timestamps);


end