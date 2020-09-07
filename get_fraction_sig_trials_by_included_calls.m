include_trial_flags =  {'either_included','neither_included','playback'};
n_include_types = length(include_trial_flags);
[bat_pair_corr,all_bat_pair_shuffled_p,all_bat_pairs,expDates] = deal(cell(1,n_include_types));
for include_k = 1:n_include_types
    [bat_pair_corr{include_k},all_bat_pair_shuffled_p{include_k},all_included_call_nums,all_bat_pairs{include_k},expDates{include_k}] = calculate_all_cross_brain_lfp_corr(include_trial_flags{include_k},false);
end

%%

all_bat_pair_corr = cellfun(@(x) cat(1,x{:}),bat_pair_corr,'un',0);

multiple_h_correction = 'none';
shuffleTypes = {'trial_shuffled','time_shuffled'};
n_shuffle_types = length(shuffleTypes);
n_freq_bands = size(all_bat_pair_corr{1},3);
n_time_bins = size(all_bat_pair_corr{1},4);
n_include_types = length(include_trial_flags);

p = cell(n_include_types,n_shuffle_types);

minTrial = 50;
alpha = 0.05;

for shuffle_k = 1:n_shuffle_types
    shuffle_type_var = shuffleTypes{shuffle_k};
    bat_pair_shuffled_p = cellfun(@(x) cat(1,x.(shuffle_type_var)),all_bat_pair_shuffled_p,'un',0);
    for include_k = 1:n_include_types
        n_bat_pairs = size(all_bat_pair_corr{include_k},2);
        p{include_k,shuffle_k} = nan(n_bat_pairs,n_freq_bands,n_time_bins);
        for bat_k = 1:n_bat_pairs
            for f_k = 1:n_freq_bands
                for time_bin_k = 1:n_time_bins
                    current_indices = {bat_k,f_k,time_bin_k};
                    nTrial = size(bat_pair_shuffled_p{include_k},1);
                    n_used_trial = sum(~isnan(bat_pair_shuffled_p{include_k}(:,bat_k,f_k,time_bin_k)));
                    
                    if n_used_trial < minTrial
                        p{include_k,shuffle_k}(current_indices{:}) = NaN;
                        continue
                    end
                    p_list = 1-bat_pair_shuffled_p{include_k}(:,bat_k,f_k,time_bin_k);
                    switch multiple_h_correction
                        case 'BH'
                            
                            used_p_idx = find(~isnan(p_list));
                            used_p = p_list(used_p_idx);
                            
                            pBH = alpha*((1:n_used_trial)/n_used_trial);
                            
                            [p_sort,p_sort_idx] = sort(used_p);
                            n_sig_trial = find(p_sort<pBH',1,'last');
                            
                            if isempty(n_sig_trial)
                                n_sig_trial = 0;
                            end
                            
                            
                        case 'bonferroni'
                            
                            n_sig_trial = sum(p_list<alpha/n_used_trial);
                        case 'none'
                            
                            n_sig_trial = sum(p_list<alpha);
                    end
                    
                    
                    p{include_k,shuffle_k}(current_indices{:}) = n_sig_trial/n_used_trial;
                    
                end
            end
        end
    end
end
%%
corr_t = mean(expParams.time_bins,2);
clf

used_bat_pairs = ~any(cellfun(@(x) strcmp(x,'71360'),all_bat_pairs{1}),2);
n_used_bat_pairs = sum(used_bat_pairs);
n_plot_columns = 3;
n_plot_rows = 5;

p_ylims = 100*[0 0.2; 0 0.54; 0.5 0.84];

significance_ylabels = {'%% Sig. \n(Trial)','%% Sig. \n(Time)','%% Above \nBaseline'};

for bat_k = 1:n_used_bat_pairs
    for shuffle_k = 1:n_shuffle_types + 1
        subplot_number = bat_k + 2*n_plot_columns + (shuffle_k-1)*n_plot_columns;
        subplot(n_plot_rows,n_plot_columns,subplot_number)
        h = gca;
        hold(h,'on')
        ylim(h,p_ylims(shuffle_k,:))
        h.FontSize = 20;
        if bat_k == 1
            ylabel(h,sprintf(significance_ylabels{shuffle_k}))
        else
            h.YAxis.Visible = 'off';
        end
        if shuffle_k == n_shuffle_types + 1
            xlabel(h,'Time (s)')
        else
            h.XAxis.Visible = 'off';
        end
        h.YTick = 10*round(p_ylims(shuffle_k,:)/10);
    end
    
    subplot(n_plot_rows,n_plot_columns,[bat_k bat_k+n_plot_columns])
    h = gca;
    hold(h,'on')
    ylim([-0.05 0.3])
    set(gca,'FontSize',20)
    if bat_k == 1
        ylabel(sprintf('Median Rho'))
    else
        h.YAxis.Visible = 'off';
    end
    h = gca;
    h.XAxis.Visible = 'off';
        
    used_bat_pairs = ~any(cellfun(@(x) strcmp(x,'71360'),all_bat_pairs{1}),2);
    used_bat_pair_idx = find(used_bat_pairs);
    bat_idx = used_bat_pair_idx(bat_k);
    title(strjoin(all_bat_pairs{1}(bat_idx,:),' + '))
   
    for include_k = 1:n_include_types
        
        if bat_k > size(all_bat_pairs{include_k},1)
            continue
        end
        
        used_bat_pairs = ~any(cellfun(@(x) strcmp(x,'71360'),all_bat_pairs{include_k}),2);
        used_bat_pair_idx = find(used_bat_pairs);
        bat_idx = used_bat_pair_idx(bat_k);
        
        n_used_bat_pairs = sum(used_bat_pairs);
            
        subplot(n_plot_rows,n_plot_columns,[bat_k bat_k+n_plot_columns])
        avgR = squeeze(nanmedian(all_bat_pair_corr{include_k}(:,bat_idx,1,:),1));
        stdR = squeeze(mad(all_bat_pair_corr{include_k}(:,bat_idx,1,:),1))./squeeze(sqrt(sum(~isnan(all_bat_pair_corr{include_k}(:,bat_idx,1,:)))));
        errorbar(corr_t,avgR,stdR,'-x','LineWidth',4)
        
        for shuffle_k = 1:n_shuffle_types
            subplot_number = bat_k + 2*n_plot_columns + (shuffle_k-1)*n_plot_columns;
            subplot(n_plot_rows,n_plot_columns,subplot_number)
            plot(corr_t,100*squeeze(p{include_k,shuffle_k}(bat_idx,shuffle_k,:,:)),'-x','LineWidth',4)
        end
        
        if include_k < n_include_types
            subplot_number = bat_k + 2*n_plot_columns + n_shuffle_types*n_plot_columns;
            subplot(n_plot_rows,n_plot_columns,subplot_number)
            include_p_over_baseline = cat(1,p_over_baseline{include_k}{:});
            avg_p_over_baseline = 100*squeeze(nanmean(include_p_over_baseline(:,bat_idx,:)))';
            std_p_over_baseline = 100*squeeze(nanstd(include_p_over_baseline(:,bat_idx,:)))';
            
            errorbar(corr_t,avg_p_over_baseline,std_p_over_baseline./sqrt(squeeze(sum(~isnan(include_p_over_baseline(:,bat_idx,:))))'),'-x','LineWidth',4)
        end

    end   
   
end

subplot_number = [1 n_plot_columns+1];
subplot(n_plot_rows,n_plot_columns,subplot_number)
l = legend(cellfun(@(x) strrep(x,'_',' bat '),include_trial_flags,'un',0),'Location','northwest');
legend box off

l.Position(2) = l.Position(2) + 0.05;
