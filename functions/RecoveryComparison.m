function tbl = RecoveryComparison(fitted,recovered,names)
% as it sounds - I give data and it does the magic
% fitted & recovered is subject * variables
% names is names of the variables

table_data = cell(length(names),16);
for i=1:length(names)
    diff = recovered(:,i)-fitted(:,i);
    table_data{i,1}=names{i};
    table_data{i,2}=round([min(diff), max(diff)],2);
    table_data{i,3}=round(nanmean(diff),2);
    table_data{i,4}=round(nanstd(diff),2);
    table_data{i,5}=round(nanmedian(diff),3);
    table_data{i,6}=round(prctile(diff,75)-prctile(diff,25),3);
    table_data{i,7}=round([min(abs(diff)), max(abs(diff))],2);
    table_data{i,8}=round(nanmean(abs(diff)),2);
    table_data{i,9}=round(nanstd(abs(diff)),2);
    table_data{i,10}=round(nanmedian(abs(diff)),3);
    table_data{i,11}=round(prctile(abs(diff),75)-prctile(abs(diff),25),3);
    [recovery_corr_spearman,p_spearman]=corr(recovered(:,i),fitted(:,i),'Type','Spearman','Rows','Complete');
    [recovery_corr_pearson,p_pearson]=corr(recovered(:,i),fitted(:,i),'Type','Pearson','Rows','Complete');
    table_data{i,12}=round(recovery_corr_pearson,3);
    table_data{i,13}=round(p_pearson,3);
    table_data{i,14}=round(recovery_corr_spearman,3);
    table_data{i,15}=(p_spearman);
    table_data{i,16}=sqrt(nanmean(diff.^2));
end
tbl=cell2table(table_data,'VariableNames',{'name', 'diff_range', 'diff_mean', 'diff_std', 'diff_median','diff_interquartile_range',...
    'abs_diff_range', 'abs_diff_mean', 'abs_diff_std', 'abs_diff_median','abs_diff_interquartile_range',...
    'corr_pearson_r','corr_pearson_p','corr_spearman_r','corr_spearman_p','rmse'});