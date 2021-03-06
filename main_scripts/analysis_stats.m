%% load relevant stuff
close all;clc;clear all;
CODE_DIR='C:\Users\user\Google Drive\Tapping_Project';
cd(CODE_DIR);
load('raw_data.mat')
load('analyses_data.mat')

%% methods (participants)
wanted_cols = [4,16,11,5,7,8,10,17,19];
participants_data = table2array(subject_info(:,wanted_cols));
titles = subject_info.Properties.VariableNames(wanted_cols);
[participants_values,participants_stats]=CreateStatTables(participants_data,titles,grp);

%% isochronous

[basic_values,basic_stats]=CreateStatTables([mean(NMA(:,1,:),3) mean(STD(:,1,:),3)],{'NMA', 'STD'},grp);
std_effects = EffectSizeTable(nanmean(STD(:,1,:),3),{'std'},grp);

% proportion of grps outside the 2 std range (of control)
tmp_std = mean(STD(grp==1,1,:),3);
cutoff = 2*std(tmp_std) + mean(tmp_std);
[x,CHI2,P] = crosstab(grp, mean(STD(:,1,:),3) > cutoff);
percentages = x(:,2)./g_size';
%no one is less than 2 std - crosstab(grp, mean(STD(:,1,:),3) < -2*std(tmp_std) + mean(tmp_std))

%mean asynch negative for all groups
dat = mean(NMA(:,1,:),3);
for i=1:3
    [~,p]=ttest(dat(grp==i));
end

[autoreg_values,autoreg_stats]=CreateStatTables([correls_all autoreg_coeffs],{'autocorr','reg1','reg2','reg3','reg4'},grp);
autoreg_coeffs_p_fromzero=nan(3,4);
for et=1:4
    for pp=1:3
        autoreg_coeffs_p_fromzero(pp,et)=signrank(autoreg_coeffs(grp==pp,et));
    end
end

% group correlations in the histograms (pearson)
stats = nan(3,2);
for pp=1:3
    dat = isoch_all_e{pp};
    [r,p]=corr(dat(:,1),dat(:,2),'Type','Pearson','Rows','Complete');
    stats(pp,1)=r;stats(pp,2)=p;
end

[model_values,model_stats]=CreateStatTables(wing_isoch_av,{'alpha','st','sm'},grp);
alpha_effects = EffectSizeTable(wing_isoch_av(:,1),{'alpha'},grp);

[counts,chi2,p,labels] = crosstab(hierarchical_num_coeffs,grp);

% test retest
correlations_in_groups = nan(2,4,2); %NMA,std; pp1-3, all; r p
measures = {STD,NMA};
for i=1:2
    tmp = squeeze(nanmean((measures{i}-mean(measures{i}(grp==1,:,:)))./std(measures{i}(grp==1,:,:)),2));
    for pp=1:3
        [r,p]=corr(tmp(grp==pp,1),tmp(grp==pp,2),'Type','Spearman','Rows','Complete');
        correlations_in_groups(i,pp,:) = [r p];
    end
    [r,p]=corr(tmp(:,1),tmp(:,2),'Type','Spearman','Rows','Complete');
    correlations_in_groups(i,4,:) = [r,p];
end

tmp = mean(NMA,3);
nma_zscored = nanmean(tmp-nanmean(tmp(grp==1,:)) ./ nanstd(tmp(grp==1,:)),2);
tmp = mean(STD,3);
std_zscored = nanmean(tmp-nanmean(tmp(grp==1,:)) ./ nanstd(tmp(grp==1,:)),2);
Single_Participant_Graph([std_zscored grp],leg.clean,'KW')

%% alternating blocks 

[signal_detection_values,signal_detection_stats]=CreateStatTables_Blocks(individual_signal_detection,{'d_prime', 'AUC', 'nomin'},grp);
[model_changes_values,model_changes_stats]=CreateStatTables_Blocks(wing_changes_av,{'alpha', 'beta', 'st', 'sm'},grp);
all_z = [zindividual_signal_detection,zwing_changes];
names = {'zd_prime', 'zAUC', 'znomin','zalpha','zbeta','zst','zsm'};
[z_values,z_stats]=CreateStatTables(all_z,names,grp);
changes_effects = EffectSizeTable(zwing_changes,{'zalpha','zbeta','zst','zsm'},grp);

% model comparison
[aic1,bic1] = aicbic(nansum(loglik_segs(:,3:5,:,1),[2 3]),nansum(nsegs(:,3:5,:),[2 3])*4,nansum(nsegs(:,3:5,:),[2 3])*10);
[aic2,bic2] = aicbic(nansum(loglik_segs(:,3:5,:,2),[2 3]),nansum(nsegs(:,3:5,:),[2 3])*3,nansum(nsegs(:,3:5,:),[2 3])*10);
[~,lr_ps]=lratiotest(nansum(loglik_segs(:,3:5,:,1),[2 3]),nansum(loglik_segs(:,3:5,:,2),[2 3]),nansum(nsegs(:,3:5,:),[2 3]));

%% both models

[comb_values,comb_stats]=CreateStatTables(model_params_combined,{'combined update rate'},grp);
comb_effects = EffectSizeTable(model_params_combined,{'model_params_combined'},grp);

%group level correlations of update rates
var1tmp = wing_isoch_av(:,1);var2tmp = zwing_changes(:,2);
[r,p]=corr(var1tmp,var2tmp,'Type','Spearman','Rows','Complete');
txt = ['\rho',sprintf(' = %0.2g (p = %0.2g)',r,p)];

% correlations between all error types
corr_errors = nan(3,3,2); %grp; (alpha1 alpha2; alpha1 beta2; alpha2 beta2); r\p
for pp=1:3
    [r,p]=corr(wing_isoch_av(grp==pp,1),zwing_changes(grp==pp,1),'Type','Spearman','Rows','Complete');
    corr_errors(pp,1,:)=[r p];
    [r,p]=corr(wing_isoch_av(grp==pp,1),zwing_changes(grp==pp,2),'Type','Spearman','Rows','Complete');
    corr_errors(pp,2,:)=[r p];
    [r,p]=corr(zwing_changes(grp==pp,1),zwing_changes(grp==pp,2),'Type','Spearman','Rows','Complete');
    corr_errors(pp,3,:)=[r p];
end

% aq correlations
aq_corrs = nan(3,2); %con,asd,all; r,p
var2 = aq_data_austin(:,1);var1 = model_params_combined(aq_data.Subject_Index);

[r,p]=corr(var1(aq_data.Group==1),var2(aq_data.Group==1),'Type','Spearman','Rows','Complete');
aq_corrs(1,:)=[r p];
[r,p]=corr(var1(aq_data.Group==3),var2(aq_data.Group==3),'Type','Spearman','Rows','Complete');
aq_corrs(2,:)=[r p];
[r,p]=corr(var1,var2,'Type','Spearman','Rows','Complete');
aq_corrs(3,:)=[r p];
    
% aq basic info
dataq = aq_data_austin(:,3);
for pp=[1,3]
dat = dataq(aq_data.Group==pp);
med = median(dat);
intqrtl=round(prctile(dat,75)-prctile(dat,25),3);
end

p=ranksum(dataq(aq_data.Group==1),dataq(aq_data.Group==3));
Cliffs_Delta(dataq(aq_data.Group==1),dataq(aq_data.Group==3));

n_perm = 1000;
data = [aq_data_austin(:,3), model_params_combined(aq_data.Subject_Index)];
group_idx = aq_data.Group;
[p_value, true_corr_diff, perm_corr_diff] = combine_groups_stats(data, group_idx, n_perm);

%%
%load('simulation_data.mat')
fitted_data = [wing_isoch_av, squeeze(wing_changes_av(:,3,:)) squeeze(wing_changes_av(:,4,:)) squeeze(wing_changes_av(:,5,:)) zwing_changes model_params_combined];
recovered_data = [mean_recovery_isoch, mean_recovery_changes(:,:,1), mean_recovery_changes(:,:,2), mean_recovery_changes(:,:,3) mean_zrecovery_changes mean_model_params_combined_recovery];
params = {'isoch-alpha','isoch-tk','isoch-mn','50-alpha','50-beta','50-tk','50-mn',...
    '70-alpha','70-beta','70-tk','70-mn','90-alpha','90-beta','90-tk','90-mn',...
    'z-alpha','z-beta','z-tk','z-mn','combined'};

recovery_stats = RecoveryComparison(fitted_data,recovered_data,params);
recovery_stats(:,[1,8,14:16])
for pp = 1:3
    recovery_stats = RecoveryComparison(fitted_data(grp==pp,:),recovered_data(grp==pp,:),params);
    recovery_stats(:,[1,8,14:16])
end

raw_data = correls_all;
simulated_data = nanmean(correls_all_sim,2);
for pp=1:3
[r,p]=corr(simulated_data(grp==pp),raw_data(grp==pp),'Type','Spearman','Rows','Complete')
end

raw_data = correls_all;
simulated_data = correls_all_sim;
exc_subjects = any(isnan(wing_isoch_av),2);
[~,tbl] = kruskalwallis(raw_data(~exc_subjects),grp(~exc_subjects),'off');
real_effect =tbl{2,5};
effect_dist = nan(1,round(n_rep/2));
for rep = 1:round(n_rep/2)
    [~,tbl] = kruskalwallis(simulated_data(~exc_subjects,rep),grp(~exc_subjects),'off');
    effect_dist(rep) = tbl{2,5};
end

figure;histogram(effect_dist);hold on
ylm = ylim;
plot([real_effect real_effect], ylm, 'k', 'LineWidth', 2)
p_val = 2*(sum((effect_dist)>(real_effect))/round(n_rep/2));
txt = sprintf('p value %0.2g',p_val);
text(0.8,0.9,txt,'Units','normalized','FontSize',font_size.txt);