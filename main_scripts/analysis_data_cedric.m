clear

load('/datadisk/Dropbox/fromPapers/Vishne_2021/tapping_data_4bGLS.mat')

n_sub = size(tapping_data.e, 1); 
n_cond = size(tapping_data.e, 2); 
n_trials = size(tapping_data.e, 3); 

for i_sub=1:n_sub
    for i_cond=1:n_cond
        for i_trial=1:n_trials
            
            if isempty(tapping_data.e{i_sub, i_cond, i_trial})
                continue
            end
            
            tapping_data.e{i_sub, i_cond, i_trial} = ...
                tapping_data.e{i_sub, i_cond, i_trial} * 1000; 
            
            tapping_data.r{i_sub, i_cond, i_trial} = ...
                tapping_data.r{i_sub, i_cond, i_trial} * 1000; 
            
            iois = diff(tapping_data.s{i_sub, i_cond, i_trial}); 
            
            tapping_data.s{i_sub, i_cond, i_trial} = ...
                repmat(iois(1), ...
                       length(tapping_data.e{i_sub, i_cond, i_trial}),...
                       1) * 1000; 
            
        end
    end
end


%% autocorrelation

lag = 1;
start_offset=1;
correls_all = nan(n_sub, n_cond);
all_e = cell(3,1);

for i_cond=1:n_cond
    
    vec_all_e = [];
    
    for i_sub=1:n_sub
        
        subj_e = [];
        
        for i_trial = 1:2
            
            e = tapping_data.e{i_sub, i_cond, i_trial};
            
            if isempty(e)
                continue
            end
            
            e = e - nanmean(e(start_offset:end));
            
            vec_all_e = [vec_all_e; ...
                e(start_offset:(end-lag)), e((start_offset+lag):end)...
                ];
            
            subj_e = [subj_e; ...
                e(start_offset:(end-lag)), e((start_offset+lag):end)];
            
        end
        
        correls_all(i_sub, i_cond) = corr(subj_e(:,1), subj_e(:,2), ...
                                          'Rows','complete');
                            
    end
    
    all_e{i_cond} = vec_all_e;
    
end


grp1 = correls_all(:, 1); 
grp2 = correls_all(:, 2); 
grp3 = correls_all(:, 3); 

figure('color', 'w')
plot(1 + randn(length(grp1), 1) * 0.1, grp1, 'bo'); 
hold on
plot(2 + randn(length(grp2), 1) * 0.1, grp2, 'ro'); 
plot(3 + randn(length(grp3), 1) * 0.1, grp3, 'ko'); 
xlim([0, 4])




%% autoregression - single subject

hierarchical_num_coeffs = nan(n_sub, n_cond);
% Successively fit linear models trying to predict the current asynchrony
% from last n asynchronies. Start with -1, then fit -1, -2, then -1, -2,
% -3, etc. 
% Always test if the coefficient for the largest fitted lag is significant. 
% Stop once the coeficient for the largest fitted lag is not signifiacnt anymore. 

max_coef = 5; 

for i_sub=1:n_sub
    
    for i_cond=1:n_cond

        y = [];
        X = [];
        
        for i_trial=1:n_trials
            
            % prepare asynchronies (subtract mean)
            e = tapping_data.e{i_sub, i_cond, i_trial};
            if isempty(e)
                continue
            end
            e = e - nanmean(e);
            
            % append the trials into a signle long vector of the outcome
            % variable (current tap asynchrony)
            y = [y; e((max_coef+1):end)];
            Xtmp = [];
            for i = 1:max_coef
                % Concatenate the predictor vectors (i.e. lagged asynchronies) for
                % the two trials. 
                % Stack lags as columns of the design matrix. 
                Xtmp = [Xtmp, e((max_coef+1-i):end-i)];
            end
            X = [X; Xtmp];
        end
        
        for npred=2:max_coef
            % Fit models with increasing amount of lags. We will have 1
            % parameter per lag + intercept. 
            model = [eye(npred), zeros(npred,1)];

            mdl_reg = fitlm(X(:, 1:npred), y, model);
            
            % For each test if the largest fitted lag coef is significantly 
            % different from 0. If yes, continue fitting, if no, stop fitting. 
            [p, F] = coefTest(mdl_reg, [zeros(1,npred-1), 1]);

            if p > 0.1
                % test 
                hierarchical_num_coeffs(i_sub, i_cond) = npred-1;
                break
            end
        end
        
    end
end

% this vector contains the number of singnificant coefficients (i.e. lags)
% per subject
hierarchical_num_coeffs(isnan(hierarchical_num_coeffs)) = max_coef;




%% autoregression - group level


% 
% hierarchical_num_coeffs_group = nan(1,4);
% % go over each group (or all subjects when pp=4) - essentially concatenate
% % data across subjects in each group and fit linear model (i.e. this is not
% % a hieararchical model...we're just pooling all data across subjects)
% for pp=1:4
% 
%     if pp<4
%         sample_size = g_size(pp);
%         sub_inds = find(grp==pp)';
%     else
%         sample_size = n_sub;
%         sub_inds = 1:n_sub;
%     end
% 
%     yall = cell(sample_size,1); 
%     Xall = cell(sample_size,1);
%     s_ind=1;
% 
%     for i_sub=sub_inds
% 
%         % y vector contains asychronies across subjects
%         y=[];
% 
%         % X matrix contains predictor asynchronies (each column is a lag)
%         % [    |              |       ...     
%         % [    |              |       ...     
%         % [trial1_lag1    trial1_lag2 ...       
%         % [    |              |       ...        
%         % [    |              |       ...    
%         % [    |              |       ...     
%         % [    |              |       ...     
%         % [trial2_lag1    trial2_lag2 ...       
%         % [    |              |       ...        
%         % [    |              |       ...    
%         %     ...            ...
%         X=[];
%         for i_trial = 1:2
%             % prepare asynchronies (subtract mean)
%             e = tapping_data.e{i_sub,1,i_trial};
%             mean_e = nanmean(e);
%             e = e-mean_e;
%             % outcome variable is the asycnrhony of each tap (starting from
%             % lags+1 of course)...we'll concatenate the two trials into a
%             % single vector...
%             y = [y; e((max_coef+1):end)];
%             Xtmp = [];
%             for i = 1:max_coef
%                 % predictor matrix: each column is a shifted version of the
%                 % data -> by -1 for lag1, by -2 for lag2 etc...again each 
%                 % column contains data concatenated across trials
%                 Xtmp = [Xtmp, e((max_coef+1-i):end-i)];
%             end
%             % append the predictor matrix for this trial 
%             X = [X; Xtmp];
%         end
%         % put all subjects together into a cell 
%         yall{s_ind} = y; 
%         Xall{s_ind} = X;
%         s_ind = s_ind+1;
%     end
%     % concatenate the outcome variable across subjects into a one long
%     % vector
%     ycombined = cell2mat(yall);
%     % now we prepare the predictor matrix by combining across subjects
%     Xcombined = zeros(length(ycombined), max_coef*sample_size);
%     ind = 0;
%     for i_sub=1:sample_size
%         % we stack columns for each subject, so that we have a "diagonal"
%         % design matrix -> plot it with: 
%         % > imagesc(Xcombined)
%         % > colormap(gray)
%         Xcombined((ind+1):(ind+size(Xall{i_sub},1)), (i_sub-1)*max_coef+(1:max_coef)) = Xall{i_sub};
%         ind = ind + size(Xall{i_sub},1);
%     end
% 
%     % fit the models 
%     for npred=1:max_coef
%         % create model matrix: there will be one parameter per lag for each
%         % participant, + one intercept 
%         model = [eye(npred*sample_size), zeros(npred*sample_size,1)];
%         % we will take out the lags of interest (still keeping all subjets
%         % next to each other)
%         Xidx = repmat(boolean([ones(npred,1);zeros(max_coef-npred,1)]),sample_size,1);
%         mdl_reg = fitlm(Xcombined(:,Xidx), ycombined, model);
%         % prepare the linear hypothesis vector: if we set multiple
%         % parameters to 1 here, we'll essentially get ANOVA out - i.e. as
%         % if we tested main effect (model comparison with vs. without those
%         % regressors included)
%         % Here we'll test coefficient corresponding to the tested lag (and
%         % we'll set the hypothesis to 1 for this coefficient for all
%         % subjects)
%         hypothesis = repmat([zeros(1,npred-1), 1], 1, sample_size); 
%         [p, F] = coefTest(mdl_reg, hypothesis);
%         if p > 0.1
%             hierarchical_num_coeffs_group(pp) = npred-1;
%             break
%         end
%     end
% 
% end % end of group 
% 
% % this vector contains the number of singnificant coefficients (i.e. lags)
% % per group 
% hierarchical_num_coeffs_group(isnan(hierarchical_num_coeffs_group))=max_coef;
% 
% 
% 
% % here we just get the coefficients (not pvals) for each lag, separately for
% % each participant 
% num_coeffs = 4;
% autoreg_coeffs = nan(n_sub, num_coeffs);
% for i_sub = 1:n_sub
%     y=[];
%     X=[];
%     for i_trial = 1:2
%         e = tapping_data.e{i_sub,1,i_trial};
%         mean_e = nanmean(e);
%         e = e-mean_e;
%         y = [y; e(5:end)];
%         X = [X; e(4:end-1), e(3:end-2), e(2:end-3), e(1:end-4)];
%     end
%     [blk,~,resid] = regress(y,X);
%     autoreg_coeffs(i_sub,:) = blk;
% end







%% computational model


% 
% % [subject, paramter, trial]
% wing_isoch = nan(n_sub,3,2);
% 
% for i_sub = 1:n_sub
%     for i_trial = 1:2
% 
%         % take out the asynchronies for this subject/trial
%         e = tapping_data.e{i_sub,1,i_trial};
%         % take out the inter-tap intervals for this subject/trial            
%         r = tapping_data.r{i_sub,1,i_trial};
% 
%         % get the means (ignore NaNs)
%         me = nanmean(e);
%         mr = nanmean(r);
% 
%         % if more than 40% metronome sounds were missed, reject the trial
%         if ERR(i_sub,1,i_trial) > block_exc
%            continue
%         end
% 
%         % fit the model 
%         [...
%             wing_isoch(i_sub,1,i_trial), ...
%             wing_isoch(i_sub,2,i_trial), ....
%             wing_isoch(i_sub,3,i_trial)...
%         ]= ...
%             model_fit_exp1(r, e, me, mr);
%     end
% end
% 
% % average the extimated parameters across trials
% wing_isoch_av = nanmean(wing_isoch, 3);
% 
%     
% % alpha_grp1 = wing_isoch_av(grp==1, 1); 
% % alpha_grp2 = wing_isoch_av(grp==2, 1); 
% % alpha_grp3 = wing_isoch_av(grp==3, 1); 
% % 
% % figure('color', 'w')
% % plot(1, alpha_grp1, 'bo'); 
% % hold on
% % plot(2, alpha_grp2, 'ro'); 
% % plot(3, alpha_grp3, 'ko'); 
% % xlim([0, 4])
% 
% 







