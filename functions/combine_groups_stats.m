function [p_value, true_corr_diff, perm_corr_diff] = combine_groups_stats(data, grp, n_perm)
% data is 2 columns - 1 with values of X, 2 with values of Y
% grp is group affiliation

grp_indexes = grp2idx(grp);
g_sizes = [sum(grp_indexes==1), sum(grp_indexes==2)];
true_data = {data(grp_indexes==1,:),data(grp_indexes==2,:)};

true_correlations = nan(1,3);
true_correlations(1)=corr(true_data{1}(:,1),true_data{1}(:,2),'Type','Spearman','Rows','Complete');
true_correlations(2)=corr(true_data{2}(:,1),true_data{2}(:,2),'Type','Spearman','Rows','Complete');
true_correlations(3)=corr(data(:,1),data(:,2),'Type','Spearman','Rows','Complete');
true_corr_diff = [true_correlations(1:2)-repmat(true_correlations(3),1,2),true_correlations(2)-true_correlations(1)];

perm_correlations = nan(n_perm, 3); % grp 1, 2, all
for perm = 1:n_perm
    resampled_data_all = [];
    for g = 1:2
        resampled_idx = randsample(g_sizes(g), g_sizes(g),true);
        resampled_data = true_data{g}(resampled_idx,:);
        perm_correlations(perm, g)=corr(resampled_data(:,1),resampled_data(:,2),'Type','Spearman','Rows','Complete');
        resampled_data_all = [resampled_data_all;resampled_data];
    end
    perm_correlations(perm, 3)=corr(resampled_data_all(:,1),resampled_data_all(:,2),'Type','Spearman','Rows','Complete');
end
perm_corr_diff = [perm_correlations(:,1:2)-repmat(perm_correlations(:,3),1,2), perm_correlations(:,2)-perm_correlations(:,1)];

p_value = sum(abs(perm_corr_diff)>abs(true_corr_diff))/n_perm;

end

