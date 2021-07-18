%% load stuff
clear all; clc; close all;
CODE_DIR='C:\Users\user\Google Drive\Tapping_Project';
cd(CODE_DIR);
load('raw_data.mat')
SAVE_NAME = 'analyses_data';

% saves in different file
run_simulations = false;

% these are just to save time. if I don't rerun them it copies from the
% previous file and saves it with the new stuff
run_hierarchical_regression = false;
run_computational_modelling = false;

%% visuals

font_size = struct();
font_size.tit   = 7;
font_size.lab   = 7;
font_size.tick  = 6;
font_size.txt   = 5;
font_size.ltr   = 7;

symbols = struct();
symbols.c = {[0 114 178]/255,[213 94 0]/255,[0 158 115]/255}; % colors from Wong, 2011 (https://doi.org/10.1038/nmeth.1618)
symbols.s = {'o','^','s'};

fig_num = 0; % in case I copy a figure with this in the title to another file (otherwise it would have been in all figures)
ltrs = char("abcdefghijklmnopqrstuvwxyz");

leg = struct();
leg.size = {'CON (N=47)','DYS (N=32)','ASD (N=30)'};
leg.clean = {'CON','DYS','ASD'};

%% basic info

grp     = subject_info.Group;
n_subs  = length(grp);
tmp = tabulate(grp);
g_size  = tmp(:,2);

step_sizes = [0, 10:20:90]';
tempos = 500*ones(6,2);
tempos(:,1) = tempos(:,1)-0.5*step_sizes;
tempos(:,2) = tempos(:,2)+0.5*step_sizes;
block_exc           = 0.4;
blocks_alternating  = 3:5;

%% general

STD=nan(n_subs,6,2);NMA=nan(n_subs,6,2);ERR = nan(n_subs,6,2);
for blk = 1:6
    for part = 1:2
        for sub = 1:n_subs
            e = tapping_data.e{sub,blk,part};
            r = tapping_data.r{sub,blk,part};
            s = tapping_data.s{sub,blk,part};
            s2=s;
            for kk=1:length(s)-1
                if round(s(kk+1)/5)*5~=round(s(kk)/5)*5
                    s2(kk+1)=s(kk);
                end
            end
            asynch=e+s-s2;
            STD(sub,blk,part) = nanstd(asynch);
            NMA(sub,blk,part) = nanmean(asynch);
            ERR(sub,blk,part) = sum(isnan(asynch))/length(asynch);
        end
    end
end

%% isochronous

% autocorrelation
lag = 1;start_offset=1;
correls_all=nan(n_subs,1);
isoch_all_e=cell(3,1);
for pp=1:3
    vec_all_e=[];
    pwho= find(grp==pp);
    for sub=pwho'
        subj_e=[];
        for part = 1:2
            e = tapping_data.e{sub,1,part};
            e = e-nanmean(e(start_offset:end));
            vec_all_e=[vec_all_e;e(start_offset:(end-lag)),e((start_offset+lag):end)];
            subj_e=[subj_e;e(start_offset:(end-lag)),e((start_offset+lag):end)];
        end
        correls_all(sub)=corr(subj_e(:,1),subj_e(:,2),'Rows','complete');
    end
    isoch_all_e{pp} = vec_all_e;
end

% autoregression
if run_hierarchical_regression
    % single subject
    hierarchical_num_coeffs =nan(n_subs,1);max_coef = 5;
    for sub = 1:n_subs
        y=[];X=[];
        for part = 1:2
            e = tapping_data.e{sub,1,part};
            nma = nanmean(e);e = e-nma;
            y=[y;e((max_coef+1):end)];
            Xtmp = [];
            for i = 1:max_coef
                Xtmp = [Xtmp, e((max_coef+1-i):end-i)];
            end
            X=[X;Xtmp];
        end
        for npred=2:max_coef
            model = [eye(npred),zeros(npred,1)];
            mdl_reg = fitlm(X(:,1:npred),y,model);
            [p, F] = coefTest(mdl_reg, [zeros(1,npred-1), 1]);
            if p > 0.1
                hierarchical_num_coeffs(sub)=npred-1;
                break
            end
        end
    end
    hierarchical_num_coeffs(isnan(hierarchical_num_coeffs))=max_coef;

    % group level
    hierarchical_num_coeffs_group = nan(1,4);
    for pp = 1:4
        if pp<4
            sample_size = g_size(pp);
            sub_inds = find(grp==pp)';
        else
            sample_size = n_subs;
            sub_inds = 1:n_subs;
        end
        yall = cell(sample_size,1); Xall = cell(sample_size,1);s_ind=1;
        for sub = sub_inds
            y=[];X=[];
            for part = 1:2
                e = tapping_data.e{sub,1,part};
                nma = nanmean(e);e = e-nma;
                y=[y;e((max_coef+1):end)];
                Xtmp = [];
                for i = 1:max_coef
                    Xtmp = [Xtmp, e((max_coef+1-i):end-i)];
                end
                X=[X;Xtmp];
            end
            yall{s_ind} = y; Xall{s_ind} = X;
            s_ind = s_ind+1;
        end
        ycombined = cell2mat(yall);
        Xcombined = zeros(length(ycombined),max_coef * sample_size);
        ind = 0;
        for sub = 1:sample_size
            Xcombined((ind + 1):(ind+size(Xall{sub},1)), (sub-1)*max_coef+(1:max_coef)) = Xall{sub};
            ind = ind +size(Xall{sub},1);
        end
        for npred=1:max_coef
            model = [eye(npred*sample_size),zeros(npred*sample_size,1)];
            Xidx = repmat(boolean([ones(npred,1);zeros(max_coef-npred,1)]),sample_size,1);
            mdl_reg = fitlm(Xcombined(:,Xidx),ycombined,model);
            [p, F] = coefTest(mdl_reg, repmat([zeros(1,npred-1), 1],1,sample_size));
            if p > 0.1
                hierarchical_num_coeffs_group(pp) = npred-1;
                break
            end
        end
    end
    hierarchical_num_coeffs_group(isnan(hierarchical_num_coeffs_group))=max_coef;
else
    load('analyses_data','hierarchical_num_coeffs','hierarchical_num_coeffs_group')
end

num_coeffs = 4;
autoreg_coeffs=nan(n_subs,num_coeffs);
for sub = 1:n_subs
    y=[];X=[];
    for part = 1:2
        e = tapping_data.e{sub,1,part};
        nma = nanmean(e);e = e-nma;
        y=[y;e(5:end)];
        X=[X;e(4:end-1),e(3:end-2),e(2:end-3),e(1:end-4)];
    end
    [blk,~,resid] = regress(y,X);
    autoreg_coeffs(sub,:)=blk;
end

% computational model
if run_computational_modelling
    wing_isoch = nan(n_subs,3,2);
    for sub = 1:n_subs
        for part = 1:2
            e = tapping_data.e{sub,1,part};
            r = tapping_data.r{sub,1,part};
            ep=e;me=nanmean(ep);
            rp=r;mr=nanmean(rp);
            if ERR(sub,1,part) > block_exc
               continue
            end
            [wing_isoch(sub,1,part), wing_isoch(sub,2,part), wing_isoch(sub,3,part)]=gal_bGLS_phase_model_single_and_multiperson(rp,ep,me,mr);
        end
    end
    wing_isoch_av = nanmean(wing_isoch,3);
else
    load('analyses_data','wing_isoch','wing_isoch_av')
end

%% changes

% single participant segments
rng_t = -2:7;
e_segments = cell(n_subs,5,2,2); % sub, block, part, acc\dec
r_segments = cell(n_subs,5,2,2);
d_segments = cell(n_subs,5,2,2);
for blk=2:6
    for sub = 1:n_subs
        for part = 1:2
            e = tapping_data.e{sub,blk,part};
            s = tapping_data.s{sub,blk,part};
            r = tapping_data.r{sub,blk,part};
            [~,inds]=min(abs(s-tempos(blk,:)),[],2);
            ind_trans = find(inds(2:end)-inds(1:end-1)~=0)+1;
            if ind_trans(1)<-(rng_t(1)-1)
                ind_trans=ind_trans(2:end);
            end
            lens=ind_trans(2:end)-ind_trans(1:end-1);
            which_transition = inds(ind_trans(1:end-1)); %1 acc,2 dec
            if length(inds)-ind_trans(end)> rng_t(end)-1
                lens=[lens;length(inds)-ind_trans(end)+1];
                which_transition=[which_transition;inds(ind_trans(end))];
            end
            ind_trans=ind_trans(1:length(lens));
            tmpe = e(ind_trans+rng_t);tmpr = r(ind_trans+rng_t);tmps = s(ind_trans+rng_t);
            segs_without_nans = ~any(isnan(tmpe),2);
            all_trans_e=tmpe(segs_without_nans,:);
            all_trans_d=all_trans_e+tmps(segs_without_nans,:);
            all_trans_r=tmpr(segs_without_nans,:);
            for d=1:2
                e_segments{sub,blk-1,part,d}=all_trans_e(which_transition(segs_without_nans)==d,:);
                d_segments{sub,blk-1,part,d}=all_trans_d(which_transition(segs_without_nans)==d,:);
                r_segments{sub,blk-1,part,d}=all_trans_r(which_transition(segs_without_nans)==d,:);
            end
        end
    end
end

% averages (fig 4)
min_trans_fig4=2;
n_transitions = cellfun(@(x) size(x,1), d_segments);
n_transitions = squeeze(sum(n_transitions,3));
d_segments_both_parts = squeeze(cellfun(@(x,y) [x;y], d_segments(:,:,1,:), d_segments(:,:,2,:),'UniformOutput',false));
d_segments_both_parts(n_transitions<min_trans_fig4)= {nan(1,length(rng_t))};
d_segments_mean = cellfun(@(x) nanmean(x,1), d_segments_both_parts,'UniformOutput',false);
tmp_reshaped = cell(n_subs,1,5,2);
tmp_reshaped(:,1,:,:) = d_segments_mean;
d_segments_mean = cell2mat(tmp_reshaped);

% signal detection
exc_signal_detection    = 0:3;
grp_d                   = cell(3,5,2); %d of the entire group- grp*blk*acc\dec
individual_d            = cell(n_subs,5,2); %same just per participant
grp_signal_detection    = nan(3,5,3); %grp*blk*type - d_prime, AUC, nomins
individual_signal_detection = nan(n_subs,5,3);
for blk=2:6
    for pp=1:3
        vec_all_e=[];
        pwho= find(grp==pp);
        d1_grp=[];d2_grp=[];
        for sub=pwho'
            d1=[];d2=[];
            for part = 1:2
                e = tapping_data.e{sub,blk,part};
                s = tapping_data.s{sub,blk,part};
                [~,inds]=min(abs(s-tempos(blk,:)),[],2);                
                ind_trans = find(inds(2:end)-inds(1:end-1)~=0)+1;
                e(ind_trans+exc_signal_detection)=nan;
                d1=[d1;e(inds==1)+s(inds==1)];
                d2=[d2;e(inds==2)+s(inds==2)];
            end
            individual_d{sub,blk-1,1}=d1;individual_d{sub,blk-1,2}=d2;
            denom = (nanvar(d2)+nanvar(d1)) / 2;
            nomin = nanmean(d2)-nanmean(d1);
            d_prime = nomin/sqrt(denom);
            [~, ~, ~, AUC]= perfcurve([ones(length(d1),1);2*ones(length(d2),1)],[d1;d2], 2);
            individual_signal_detection(sub,blk-1,1)=d_prime;
            individual_signal_detection(sub,blk-1,2)=AUC;
            individual_signal_detection(sub,blk-1,3)=nomin;
            d1_grp = [d1_grp;d1];d2_grp = [d2_grp;d2];
        end
        grp_d{pp,blk-1,1}=d1_grp;grp_d{pp,blk-1,2}=d2_grp;
        denom = (nanvar(d2_grp)+nanvar(d1_grp)) / 2;
        nomin = nanmean(d2_grp)-nanmean(d1_grp);
        d_prime = nomin/sqrt(denom);
        [~, ~, ~, AUC]= perfcurve([ones(length(d1_grp),1);2*ones(length(d2_grp),1)],[d1_grp;d2_grp], 2);
        grp_signal_detection(pp,blk-1,1)=d_prime;
        grp_signal_detection(pp,blk-1,2)=AUC;
        grp_signal_detection(pp,blk-1,3)=nomin;
    end
end

tmp=(individual_signal_detection(:,blocks_alternating,:)-mean(individual_signal_detection(grp==1,blocks_alternating,:)))./std(individual_signal_detection(grp==1,blocks_alternating,:));
zindividual_signal_detection = squeeze(mean(tmp,2));

% computational model
if run_computational_modelling
    wing_changes=nan(n_subs,5,4,2);
    loglik_segs = nan(n_subs,5,2,2); %last dim - full, partial model
    nsegs=nan(n_subs,5,2);
    for blk=1:5
        for sub = 1:n_subs
            for part = 1:2
                dat1 = [e_segments{sub,blk,part,1};e_segments{sub,blk,part,2}];
                dat2 = [r_segments{sub,blk,part,1};r_segments{sub,blk,part,2}];
                if ERR(sub,blk+1,part)>block_exc
                    continue
                end
                nreps=size(dat1,1);
                nsegs(sub,blk,part)=nreps;
                if nreps<1
                    continue
                end
                fits_per_seg = nan(4,nreps);loglik_tmp = nan(2,nreps);
                for rep = 1:nreps
                     rp = dat2(rep,:)';ep = dat1(rep,:)';
                     me = NMA(sub,blk+1,part);
                     [fits_per_seg(1,rep),fits_per_seg(2,rep),fits_per_seg(3,rep),fits_per_seg(4,rep),loglik_tmp(1,rep)]=...
                      gal_bGLS_period_model_single_and_multipeson_integ(rp,ep,me);
                     [~,~,~,loglik_tmp(2,rep)]=gal_bGLS_period_model_single_and_multipeson_diff(rp,ep,me);
                end
                wing_changes(sub,blk,:,part)= nanmean(fits_per_seg,2);
                loglik_segs(sub,blk,part,:)=nansum(loglik_tmp,2);
            end
        end
    end
    wing_changes_av=nanmean(wing_changes,4);
    tmp=(wing_changes_av(:,blocks_alternating,:)-nanmean(wing_changes_av(grp==1,blocks_alternating,:)))./nanstd(wing_changes_av(grp==1,blocks_alternating,:));
    zwing_changes = squeeze(nanmean(tmp,2));
    
    tmp=(wing_changes(:,blocks_alternating,:,:)-nanmean(wing_changes(grp==1,blocks_alternating,:,:)))./nanstd(wing_changes(grp==1,blocks_alternating,:,:));    
    zwing_changes_part = squeeze(nanmean(tmp,2));
else
    load('analyses_data','wing_changes','wing_changes_av','zwing_changes','zwing_changes_part','loglik_segs','nsegs')
end

aq50_full = aq_data.AQ50_FullScore;
aq_data_austin = [sum(table2array(aq_data(:,2+[38,11,44,17,26,47,22,40,15,34,50,13])),2),...
    sum(table2array(aq_data(:,2+[23,6,19,9,12,43,5,25])),2),...
    sum(table2array(aq_data(:,2+[39,20,45,35,7,37])),2)];

ma = nanmean(wing_isoch_av(grp==1,1));sa = nanstd(wing_isoch_av(grp==1,1));
zalpha = (wing_isoch_av(:,1) - ma)./sa;
mb = nanmean(zwing_changes(grp==1,2));sb = nanstd(zwing_changes(grp==1,2));
zbeta = (zwing_changes(:,2) - mb)./sb;
model_params_combined = (zalpha + zbeta)/2;

%% and save

save(SAVE_NAME,'font_size','symbols','ltrs','fig_num','leg',...
    'grp','n_subs','g_size','tempos',...
    'STD','NMA','ERR','rng_t',...
    'hierarchical_num_coeffs','hierarchical_num_coeffs_group','autoreg_coeffs','correls_all','isoch_all_e',...
    'grp_d','individual_d','grp_signal_detection','individual_signal_detection','zindividual_signal_detection','d_segments_mean',...
    'wing_isoch','wing_isoch_av','wing_changes','wing_changes_av','zwing_changes','zwing_changes_part','loglik_segs','nsegs',...
    'model_params_combined','aq50_full','aq_data_austin');

%% simulations

if run_simulations
rng(1) % reproducability
    
n_taps  = 100; % in block
rollout = 5; % start of simulation
n_rep   = 1000;
sim_changes_blks = 3:5; % which to simulate out of 1-5 not 2-6
exclude_taps = true; % same exclusion as raw data
add_subject_nma = true; % before exclusion

%% isochronous simulation

sim_isoch = nan(n_subs, n_rep, n_taps, 2); %e,r
stimulus = 500 * ones(n_taps+rollout, 1);
for sub=1:n_subs
    if any(isnan(wing_isoch_av(sub,:)))
        continue
    else
        alpha   = wing_isoch_av(sub,1);
        tkn     = wing_isoch_av(sub,2);
        mn      = wing_isoch_av(sub,3);
        relevant_blocks = squeeze(ERR(sub,1,:)<= block_exc);
        if sum(relevant_blocks)==2
            which_nma = ones(n_rep, 1);
            which_nma(randsample(n_rep,round(n_rep/2))) = 2;
        else
            which_nma = find(relevant_blocks)*ones(n_rep, 1);
        end
        for rep=1:n_rep
            noise =randn(2,n_taps+rollout+1);
            e_sim = [0];
            for tap = 1:(n_taps+rollout-1)
                e_sim(tap+1) = (1-alpha)*e_sim(tap) + noise(1,tap)*tkn + noise(2,tap+1)*mn - noise(2,tap)*mn;
            end
            r_sim = [e_sim(1),diff(e_sim)]+stimulus';
            e_sim = e_sim((rollout+1):end);r_sim = r_sim((rollout+1):end);
            if add_subject_nma
                e_sim = e_sim-nanmean(e_sim)+NMA(sub,1,which_nma(rep));
            end
            sim_isoch(sub,rep,:,1)=e_sim;
            sim_isoch(sub,rep,:,2)=r_sim;
        end
    end
end
if exclude_taps
    for sub=1:n_subs
        for rep=1:n_rep
            e = squeeze(sim_isoch(sub,rep,:,1));
            excluded = find(e > 200 | e < -200);
            sim_isoch(sub,rep,excluded,:) = nan;
            if ~isempty(excluded) && excluded(end)==n_taps;excluded=excluded(1:end-1);end
            sim_isoch(sub,rep,excluded+1,2) = nan;
        end
    end
end

lag = 1;start_offset=1;
correls_all_sim=nan(n_subs,round(n_rep/2));
for sub=1:n_subs
    if any(isnan(wing_isoch_av(sub,:)))
        continue
    else
        relevant_sim = squeeze(sim_isoch(sub, :, :, 1));
        ind=0;
        for rep = 1:round(n_rep/2)
            subj_e=[];
            for part = 1:2
                ind=ind+1;
                e = relevant_sim(ind,:)';
                e = e-nanmean(e(start_offset:end));
                subj_e=[subj_e;e(start_offset:(end-lag)),e((start_offset+lag):end)];
            end
            correls_all_sim(sub,rep)=corr(subj_e(:,1),subj_e(:,2),'Rows','complete');
        end
    end
end

%% isochronous parameter recovery

recovery_val_isoch = nan(3, n_subs, n_rep);
for sub = 1:n_subs
    if any(isnan(wing_isoch_av(sub,:)))
        continue
    else
        for rep = 1:n_rep
            ep = squeeze(sim_isoch(sub, rep, :, 1));
            if sum(isnan(ep))/n_taps > block_exc
                continue
            end
            rp = squeeze(sim_isoch(sub, rep, :, 2));
            me = nanmean(ep); mr = nanmean(rp);
            [alpha_sim, st_sim, sm_sim] = gal_bGLS_phase_model_single_and_multiperson(rp,ep,me,mr);
            recovery_val_isoch(:,sub,rep) = [alpha_sim, st_sim, sm_sim];
        end
    end
end

%% changes blocks simulation

sim_changes = nan(n_subs, n_rep, length(sim_changes_blks), n_taps, 3); %e,r,s
for blk = sim_changes_blks
    tempo = tempos(blk+1,:);
    for sub=1:n_subs
        if any(isnan(wing_changes_av(sub,blk,:)))
            continue
        else
            alpha = wing_changes_av(sub,blk,1);
            beta = wing_changes_av(sub,blk,2);
            tkn = wing_changes_av(sub,blk,3);
            mn = wing_changes_av(sub,blk,4);
            relevant_blocks = squeeze(ERR(sub,blk+1,:)<= block_exc);
            if sum(relevant_blocks)==2
                which_part = ones(n_rep, 1);
                which_part(randsample(n_rep,round(n_rep/2))) = 2;
            else
                which_part = find(relevant_blocks)*ones(n_rep, 1);
            end
            stimulus = nan(n_rep, n_taps+rollout);
            for rep = 1:n_rep
                stimulus_tmp = tapping_data.s{sub,blk+1,which_part(rep)};
                rollout_padding = n_taps+rollout-length(stimulus_tmp);
                stimulus(rep,:) = [repmat(stimulus_tmp(1),1,rollout_padding), stimulus_tmp'];        
                [r_sim,e_sim] = Simulate_period_correction(stimulus(rep,:)',alpha,beta,tkn,mn);
                e_sim = e_sim((rollout+1):end); r_sim = r_sim((rollout+1):end);
                sim_changes(sub,rep,sim_changes_blks==blk,:,1)=e_sim;
                if add_subject_nma
                    e_sim = e_sim-nanmean(e_sim)+NMA(sub,blk+1,which_part(rep));
                end
                sim_changes(sub,rep,sim_changes_blks==blk,:,2)=r_sim;
                sim_changes(sub,rep,sim_changes_blks==blk,:,3)=stimulus(rep,(rollout+1):end);
            end
        end
    end
end
if exclude_taps
    for sub=1:n_subs
        for rep=1:n_rep
            for blk = sim_changes_blks
                e = squeeze(sim_changes(sub,rep,sim_changes_blks==blk,:,1));
                excluded = find(e > 200 | e < -200);
                sim_changes(sub,rep,sim_changes_blks==blk,excluded,1:2) = nan;
                if ~isempty(excluded) && excluded(end)==n_taps;excluded=excluded(1:end-1);end
                sim_changes(sub,rep,sim_changes_blks==blk,excluded+1,2) = nan;
            end
        end
    end
end

% for plotting figure 4 simulated
d_segments_mean_simulated = nan(n_subs,length(rng_t),length(sim_changes_blks),2);

for blk = sim_changes_blks
    for sub=1:n_subs
        if any(isnan(wing_changes_av(sub,blk,:)))
            continue
        else
            for rep=1:n_rep
                e_sim = squeeze(sim_changes(sub,rep,sim_changes_blks==blk,:,1));
                stimulus = squeeze(sim_changes(sub,rep,sim_changes_blks==blk,:,3));
                me = nanmean(e_sim);
                e_sim = e_sim-me;
                changes  = find(abs(diff(stimulus))>5)+1;
                if changes(1) < -(rng_t(1)-1)
                    changes=changes(2:end);
                end
                if length(stimulus)-changes(end) < rng_t(end)
                    changes=changes(1:end-1);
                end
                segmented_block_e = e_sim(changes+rng_t);
                segmented_block_s = stimulus(changes+rng_t);
                segmented_block_s(any(isnan(segmented_block_e),2),:)=[];
                segmented_block_e(any(isnan(segmented_block_e),2),:)=[];
                direction = nan(size(segmented_block_e,1),1);
                for k = 1:size(segmented_block_e,1)
                    if stimulus(changes(k))<stimulus(changes(k)-1)
                        direction(k) = 1;
                    else
                        direction(k) = 2;
                    end
                end
                d_segments_mean_simulated(sub,:,sim_changes_blks==blk,1) = mean(segmented_block_e(direction==1,:) + segmented_block_s(direction==1,:));
                d_segments_mean_simulated(sub,:,sim_changes_blks==blk,2) = mean(segmented_block_e(direction==2,:) + segmented_block_s(direction==2,:));
            end
        end
    end
end

%% changes blocks recovery

recovery_val_changes = nan(4,n_subs, n_rep, length(sim_changes_blks));

for blk = sim_changes_blks
    for sub=1:n_subs
        if any(isnan(wing_changes_av(sub,blk,:)))
            continue
        else
            for rep=1:n_rep
                e_sim = squeeze(sim_changes(sub,rep,sim_changes_blks==blk,:,1));
                if sum(isnan(e_sim))/n_taps > block_exc
                    continue
                end
                r_sim = squeeze(sim_changes(sub,rep,sim_changes_blks==blk,:,2));
                stimulus = squeeze(sim_changes(sub,rep,sim_changes_blks==blk,:,3));
                changes  = find(abs(diff(stimulus))>5)+1;
                if changes(1) < -(rng_t(1)-1)
                    changes=changes(2:end);
                end
                if length(stimulus)-changes(end) < rng_t(end)
                    changes=changes(1:end-1);
                end                
                me = nanmean(e_sim);
                segmented_block_e = e_sim(changes+rng_t);
                segmented_block_r = r_sim(changes+rng_t);
                segmented_block_r(any(isnan(segmented_block_e),2),:)=[];
                segmented_block_e(any(isnan(segmented_block_e),2),:)=[];
                recovered_per_seg = nan(4,size(segmented_block_e,1));
                for k = 1:size(segmented_block_e,1)
                    ep = segmented_block_e(k, :);
                    rp = segmented_block_r(k, :);
                    [alpha_sim, beta_sim, st_sim, sm_sim]=gal_bGLS_period_model_single_and_multipeson_integ(rp',ep',me);              
                    recovered_per_seg(:,k) = [alpha_sim, beta_sim, st_sim, sm_sim];
                end
                recovery_val_changes(:,sub,rep,sim_changes_blks==blk) = nanmean(recovered_per_seg,2);
            end
        end
    end
end
mean_recovery_isoch = nanmean(recovery_val_isoch,3)';
mean_recovery_changes = permute(squeeze(nanmean(recovery_val_changes,3)),[2 1 3]);
std_recovery_isoch = nanstd(recovery_val_isoch,[],3)';
std_recovery_changes = permute(squeeze(nanstd(recovery_val_changes,[],3)),[2 1 3]);

tmp=(recovery_val_changes-nanmean(recovery_val_changes(:,grp==1,:,:),2))./nanstd(recovery_val_changes(:,grp==1,:,:),[],2);
zrecovery_changes = nanmean(tmp,4);
mean_zrecovery_changes = nanmean(zrecovery_changes,3)';
std_zrecovery_changes = nanstd(zrecovery_changes,[],3)';

tmp = squeeze(recovery_val_isoch(1,:,:));
zalpha_recovery = (tmp- nanmean(tmp(grp==1,:)))./nanstd(tmp(grp==1,:));
tmp = squeeze(zrecovery_changes(2,:,:));
zbeta_recovery = (tmp- nanmean(tmp(grp==1,:)))./nanstd(tmp(grp==1,:));
model_params_combined_recovery = (zalpha_recovery + zbeta_recovery)/2;
mean_model_params_combined_recovery = nanmean(model_params_combined_recovery,2);
std_model_params_combined_recovery = nanstd(model_params_combined_recovery,[],2);

%% save 

save('simulation_data','sim_isoch','recovery_val_isoch','sim_changes','recovery_val_changes','d_segments_mean_simulated',...
    'mean_recovery_isoch','mean_recovery_changes','std_recovery_isoch','std_recovery_changes','n_rep',...
    'zrecovery_changes','mean_zrecovery_changes','std_zrecovery_changes','model_params_combined_recovery',...
    'mean_model_params_combined_recovery','std_model_params_combined_recovery','correls_all_sim');

end