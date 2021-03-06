% load relevant stuff
close all;clc;clear all;
CODE_DIR='C:\Users\user\Google Drive\Tapping_Project';
cd(CODE_DIR);
load('raw_data.mat')
load('analyses_data.mat')
load('simulation_data.mat')

%% Experiment 1

% Figure 1: NMA, SD 
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Isoch.] explanation, NMA, STD',fig_num),'units','normalized','outerposition',[0 0 0.5 0.85])
axp=[0.08 0.62 0.89 0.335];
ha(1)=axes('Units','normalized','Position',axp);
ha(2)=axes('Units','normalized','Position',[0.08 0.075 0.42 0.44]);
ha(3)=axes('Units','normalized','Position',[0.56 0.075 0.42 0.44]);

axes(ha(1))
title('Task structure','FontSize',font_size.tit)
text(0.02,0.95,ltrs(1),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')

titles = {'Mean Asynchrony', 'Standard deviation'};
for i=1:2
    if i==1
        d_val_group = mean(NMA(:,1,:),3);
    else
        d_val_group = mean(STD(:,1,:),3);
    end
    axes(ha(i+1));
    Single_Participant_Graph([d_val_group grp] ,leg.clean, 'KW')
    title(titles{i},'FontSize',font_size.tit)
    text(0.05,0.95,ltrs(i+1),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
    grid on;
end
ylabel(ha(2),'ms','FontSize',font_size.lab)

% 2 - Consecutive correlations
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Isoch.] Consec. Hist Phase',fig_num),'units','normalized','outerposition',[0 0 1 0.65])
lims = [-240,300];bins=55;
nw=4;marg_w=[.05 .025];marg_h=[0.12 0.11];gap_w=[0.008 0.008 0.05 0];
w_tmp=(1-sum(marg_w)-sum(gap_w))/4.9;size_w=[1.3*w_tmp 1.3*w_tmp 1.3*w_tmp w_tmp];
axh = (1-sum(marg_h));py = 1-marg_h(1)-axh;px = marg_w(1);
for i=1:nw
    axw = size_w(i);
    ha(i) = axes('Units','normalized','Position',[px py axw axh]);
    px = px+axw+gap_w(i);
end

coloring = (1:64)'/64;coloring2 = (129:192)'/192;color_ind = [3 1 2];
for pp=1:3
    map = zeros(64,3);
    map(:,color_ind(setdiff(1:3,pp)))=repmat(coloring,1,2);map(:,color_ind(pp))=coloring2;
    ax=ha(pp);axes(ax)
    d_val_group = isoch_all_e{pp};
    coefs=regress(d_val_group(:,2),[ones(length(d_val_group(:,1)),1),d_val_group(:,1)]);
    histogram2(d_val_group(:,1),d_val_group(:,2),bins,'DisplayStyle','tile','LineStyle','none');hold all;
    text(0.05,0.95,ltrs(pp),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
    plot(lims,[lims(1)*coefs(2,1)+coefs(1,1),lims(2)*coefs(2,1)+coefs(1,1)],'k','LineWidth',2)
    xlim(lims);ylim(lims)
    caxis([0 165])
    title(leg.size{pp},'FontSize',font_size.tit-4)
    colormap(ax,map)
    set(gca, 'FontSize',font_size.tick,'YTickLabels','');
    box off
end
axes(ha(4))
Single_Participant_Graph([correls_all grp],leg.clean,'KW')
title('Single participant correlations','FontSize',font_size.tit-3)
text(0.06,0.95,ltrs(4),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
ylabel('Correlations','FontSize',font_size.lab)
ylim([0 1])
grid on
set(gca,'FontSize',font_size.tick)
ylabel(ha(1),'Asynchrony t+1 (ms)','FontSize',font_size.lab);
set(ha(1), 'FontSize',font_size.tick,'YTickLabelMode','auto');
locs=[marg_w(1) marg_h(1)-0.01 sum(size_w(1:3))+sum(gap_w(1:2)) 0.86];
[~,h]=suplabel('Asynchrony t (ms)','x', locs);set(h,'FontSize',font_size.lab+1)
[~,h]=suplabel('Relation of consecutive asynchronies' ,'t', locs);set(h,'FontSize',font_size.tit)

% 3 - computational model
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Isoch.] Linear Sensorimotor Synchronization Model (bGLS)',fig_num),'units','normalized','outerposition',[0 0 0.7 1])
ha=tight_subplot(1, 3, 0.05, [0.075 0.5], [.05 .02]);

titles = {'Phase correction (\alpha)', 'Timekeeper noise', 'Motor noise'};
ylabels = {'a.u.', 'ms','ms'};
for i=1:3
    data = wing_isoch_av(:,i);
    tit = titles{i};ylb = ylabels{i};
    axes(ha(i))
    Single_Participant_Graph([data grp],leg.clean,'KW')
    text(0.05,0.95,ltrs(i+1),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
    title(tit,'FontSize',font_size.tit)
    ylabel(ylb,'FontSize',font_size.lab)
    set(gca,'FontSize',font_size.tick)
    if i==2 || i==3
        ylim([0, 10*ceil(max(wing_isoch_av(:,2))/10)])
    end
    grid on
end

%% Experiment 2

% 4 Change dynamics
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Alter.] Change Dynamics',fig_num),'units','normalized','outerposition',[0 0 0.5 1])
ha=tight_subplot(3, 2, [0.025 0.01], 0.075, [.14 .01]);
p=nan(1,3);lims =500+[-70,70;-70,70;-70,70];
direcs = {'Accelerating','Decelerating'};
for b=4:6
    for d=1:2
        axes(ha(d+2*(6-b)));
        tempo = tempos(b,:);
        if d==1;tempo = fliplr(tempo);end
        plot([-0.5+rng_t(1),-1,0,rng_t(end)+0.5],[tempo(1),tempo(1),tempo(2),tempo(2)],'k--');hold all;
        d_val_all = d_segments_mean(:,:,b-1,d);
        box off
        for pp=1:3
            d_val_group = d_val_all(grp==pp,:)-mean(d_val_all(grp==pp,rng_t<0),2)+tempo(1);
            count_part = sum(~all(isnan(d_val_group),2));
            p(pp)=errorbar(rng_t,nanmean(d_val_group),nanstd(d_val_group)/sqrt(count_part),symbols.sl{pp},'MarkerFaceColor',symbols.c{pp},'MarkerEdgeColor','k','LineWidth',1.25);hold all;
        end
        text(-0.9*d+1.85,0.92,ltrs(d+2*(6-b)),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
        xlim([-0.5+rng_t(1),rng_t(end)+0.5]);ylim(lims(b-3,:))
        if d==1
            ylabel(sprintf('%s %d ms','Step size: ',-diff(tempo)),'FontSize',font_size.tit+2);
            set(gca, 'XTick', rng_t,'XTickLabels','','FontSize',font_size.tick);
        else
            set(gca,'XTick', rng_t,'XTickLabels','','YTickLabels', '', 'FontSize',font_size.tick);
        end
        if b==4
            set(gca, 'XTickLabels', rng_t,'FontSize',font_size.tick);
            xlabel('Beat','FontSize',font_size.lab)
        elseif b==6
            title(direcs{d},'FontSize',font_size.tit)
        end
        grid on      
    end
end
legend(ha(6),p,leg.clean,'Location','southeast')
[~,h]=suplabel('Change dynamics' ,'t', [0.14 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit+5)
[~,h]=suplabel('Delay interval (baselined, ms)' ,'y', [0.09 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit+3)

% 5 Signal detection theory
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Alter.] Signal Detection Graphs',fig_num),'units','normalized','outerposition',[0 0 1 1])
gap=[0.018 0.001];gap2=0.045;
axh = (1-sum([0.08 0.06])-(3-1)*gap(1))/3; 
axw = (1-sum([.06 .01])-3*gap(2)-gap2)/4;
py = 1-0.06-axh; 
ii = 1;
for ih = 1:3
    px = 0.06;
    for ix = 1:3
        ha(ii) = axes('Units','normalized','Position',[px py axw axh],'XTickLabel','','YTickLabel','');
        px = px+axw+gap(2);
        ii = ii+1;
    end
    ha(ii) = axes('Units','normalized','Position',[px+gap2 py axw axh],'XTickLabel','','YTickLabel','');
    ii=ii+1;
    py = py-axh-gap(1);
end
for blk=4:6
    plt=nan(3,1);
    leg_auc=cell(3,1);
    for pp=1:3
        axes(ha((6-blk)*4+pp));
        d1 = grp_d{pp,blk-1,1};d2=grp_d{pp,blk-1,2};
        plot([nanmean(d1),nanmean(d1)],[0,1],'k','LineWidth',1.5);hold on
        plot([nanmean(d2),nanmean(d2)],[0,1],'k','LineWidth',1.5);hold on
        histogram(d2,'FaceColor',symbols.c{pp},'BinEdges',240:10:750,'Normalization','probability','EdgeColor','none');hold on
        histogram(d1,'FaceColor',symbols.c{pp}/3,'BinEdges',240:10:750,'Normalization','probability','EdgeColor','none');hold on
        xlim([240,750]);ylim([0,0.075]);box off
        txt1 = [sprintf("d' = %0.2f",grp_signal_detection(pp,blk-1,1)); sprintf("Diff. means = %0.1f",grp_signal_detection(pp,blk-1,3))];
        text(0.55,0.9,txt1,'Units','normalized','FontSize',font_size.txt+2);
        text(0.05,0.92,ltrs((6-blk)*4+pp),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold');
        if blk==4
            set(gca, 'FontSize',font_size.tick,'YTickLabels','');
        else
            set(gca, 'FontSize',font_size.tick,'YTick',0:0.01:0.07,'YTickLabels','','XTickLabels','');
        end
        if blk==6
            title(leg.size{pp},'FontSize',font_size.tit)
        end
        axes(ha((7-blk)*4));
        text(0.05,0.92,ltrs((7-blk)*4),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
        [X, Y]= perfcurve([ones(length(d1),1);2*ones(length(d2),1)],[d1;d2], 2);
        plt(pp)=plot(X,Y,'Color',symbols.c{pp},'LineWidth',1.8);hold on;
        box off
        leg_auc{pp}=[leg.clean{pp},sprintf(': %0.2g',grp_signal_detection(pp,blk-1,2))];
        if blk==6
            title('ROC curve','FontSize',font_size.tit)
        end
        ylim([0 1])
        set(gca,'XTick',0:0.25:1,'YTick',0:0.25:1);
    end
    legu=legend(plt,leg_auc,'Location','southeast','FontSize',font_size.txt+2);
    title(legu,'AUC:')
    if blk==4
        set(gca,'FontSize',font_size.tick);
    else
        set(gca,'FontSize',font_size.tick,'XTickLabels','');
    end
    ylabel(ha((6-blk)*4+1),['Step size: ',num2str(diff(tempos(blk,:))),'ms'],'FontSize',font_size.tit)
    set(ha((6-blk)*4+1),'YTick',0:0.01:0.07,'YTickLabelMode','auto');
end
xlabel(ha(10),'Delay interval (ms)','FontSize',font_size.lab)
xlabel(ha(12),'False positive rate (FPR)','FontSize',font_size.lab)
ylabel(ha(8),'True positive rate (TPR)','FontSize',font_size.lab)

% 6 Model with tempo changes
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Alter.] Linear Sensorimotor Synchronization Model (bGLS)',fig_num),'units','normalized','outerposition',[0 0 1 1])
ha=tight_subplot(2, 3, [0.07 0.03], 0.075, [.04 .02]);

axes(ha(1))
text(0.05,0.95,ltrs(1),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
title('Timekeeper dynamics','FontSize',font_size.tit)
    
titles={'Phase correction (\alpha)','Period correction (\beta)','Timekeeper noise','Motor noise'};
for i=1:4
    j=i;if i<=2;j=3-i;end
    data =zwing_changes(:,i);
    axes(ha(j+2))
    Single_Participant_Graph([data grp],leg.clean,'KW')
    text(0.05,0.95,ltrs(j+1),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
    title(titles{i},'FontSize',font_size.tit)
    set(gca,'FontSize',font_size.tick)
    grid on
    if i==3 || i==4
        ylim([-2 5])
    end
    ylabel(ha(j+2),'z-score (a.u.)','FontSize',font_size.lab)
end

% 7 Combined update
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Both.] alpha beta corr',fig_num),'units','normalized','outerposition',[0 0 1 0.55])
titles = {'Experiment 1: phase correction (\alpha)','Experiment 2: period correction (\beta) z-score'};
ha=tight_subplot(1,3,0.005,[0.15,0.05],[0.05,0.28]);
ylims = round([min(wing_isoch_av(:,1))-0.1,max(wing_isoch_av(:,1))+0.2],1);
xlims = [min(zwing_changes(:,2))-0.5,max(zwing_changes(:,2))+0.5];
plts=[];
for pp=1:3
    phase_correction = wing_isoch_av(grp==pp,1);period_correction = zwing_changes(grp==pp,2);
    [r,p]=corr(phase_correction,period_correction,'Type','Spearman','Rows','Complete');
    coefs=regress(phase_correction,[ones(g_size(pp),1),period_correction]);
    txt = ['\rho',sprintf(' = %0.2g (p = %0.2g)',r,p)];
    axes(ha(pp));
    plts(pp)= scatter(period_correction, phase_correction, 40, symbols.s{pp},'filled','MarkerEdgeColor','k');hold on;
    plot([xlims(1),xlims(2)],[xlims(1)*coefs(2)+coefs(1),xlims(2)*coefs(2)+coefs(1)],'Color',symbols.c{pp},'LineWidth',1.5)
    text(0.04, 0.06 ,txt ,'Units','normalized', 'FontSize', font_size.txt+1)
    set(gca, 'FontSize',font_size.tick,'YTickLabels','');
    xlim(xlims);ylim(ylims)
    text(0.05,0.95,ltrs(pp),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
    grid on;
end
set(ha(1),'FontSize',font_size.tick,'YTickLabelMode','auto');
ylabel(ha(1),titles{1},'FontSize',font_size.lab+2);
xlabel(ha(2),titles{2},'FontSize',font_size.lab+2);
legend(ha(3),plts,leg.clean,'Location','southeast')

ha(4) = axes('Units','normalized', 'Position',[0.77 0.15 0.2167 0.8], 'XTickLabel','', 'YTickLabel','');
axes(ha(4));
Single_Participant_Graph([model_params_combined grp],leg.clean,'KW')
text(0.05,0.95,ltrs(4),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
grid on;
ylabel('Update rate (\alpha and \beta z-score)','FontSize',font_size.lab+2)
ylim([-1.7 2])

% 8 AQ correlations
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Both.] aq update rate',fig_num),'units','normalized','outerposition',[0 0 1 0.55])
titles = {'Update rate (\alpha and \beta z-score)','Communication and mindreading skills'};
ha(1) = axes('Units','normalized', 'Position',[0.04 0.15 0.22 0.8], 'XTickLabel','', 'YTickLabel','');
ha(2:4)=tight_subplot(1,3,0.005,[0.15,0.05],[0.31,0.02]);

axes(ha(1))
for pp=[1,3]
    histogram(aq_data_austin(aq_data.Group==pp,3),'FaceColor',symbols.c{pp});hold on
end
set(gca, 'FontSize',font_size.tick);
legend(leg.clean{[1,3]})
box off
text(0.05,0.95,ltrs(1),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
ylabel('Number of subjects','FontSize',font_size.lab+2)
txt = [sprintf("Wilcoxon rank sum \ntest, p = %0.2g",ranksum(aq_data_austin(aq_data.Group==3,3),aq_data_austin(aq_data.Group==1,3))); ...
    sprintf("Cliff's Delta = %0.2g",Cliffs_Delta(aq_data_austin(aq_data.Group==3,3),aq_data_austin(aq_data.Group==1,3)))];
text(0.3,0.9,txt,'Units','normalized','FontSize',font_size.txt);

aq_scale = aq_data_austin(:,3); update_rate =model_params_combined(aq_data.Subject_Index);
ylims = round([min(update_rate)-0.4,max(update_rate)+0.2],1);
xlims = [-1, 16];
plt = [];
for pp=[1,3,0]
    if pp==0
        var1 = aq_scale;
        var2 = update_rate;
        ind = 4;
    else
        var1 = aq_scale(aq_data.Group==pp);
        var2 = update_rate(aq_data.Group==pp);
        ind=(pp+3)/2;
    end
    axes(ha(ind))
    [r,p]=corr(var1,var2,'Type','Spearman','Rows','Complete');
    coefs=regress(var2,[ones(length(var1),1),var1]);
    if pp==0
        txt = ['\rho ALL',sprintf(' = %0.2g (p = %0.2g one-sided)',r,p/2)];
        plot([xlims(1),xlims(2)],[xlims(1)*coefs(2)+coefs(1),xlims(2)*coefs(2)+coefs(1)],'k-','LineWidth',1.5)
    else
        for i=[ind,4]
            axes(ha(i))
            plt(ind-1) = scatter(var1, var2, 40, symbols.s{pp},'filled','MarkerEdgeColor','k');hold on;
        end
        axes(ha(ind))
        txt = ['\rho ',leg.clean{pp},sprintf(' = %0.2g (p = %0.2g one-sided)',r,p/2)];
        plot([xlims(1),xlims(2)],[xlims(1)*coefs(2)+coefs(1),xlims(2)*coefs(2)+coefs(1)],'Color',symbols.c{pp},'LineWidth',1.5);
    end
    text(0.04, 0.06 ,txt ,'Units','normalized', 'FontSize', font_size.txt+1)
    set(gca, 'FontSize',font_size.tick,'YTickLabels','');
    xlim(xlims);ylim(ylims)
    text(0.05,0.95,ltrs(ind),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
    grid on;
end
set(ha(2),'FontSize',font_size.tick,'YTickLabelMode','auto');
ylabel(ha(2),titles{1},'FontSize',font_size.lab+2);
[~,h]=suplabel(titles{2},'x',[.05 .15 .93 .8]);
set(h,'FontSize',font_size.lab+2)
legend(ha(4),plt,leg.clean{[1,3]},'Location','northeast')

%% Supplementary

% Test-retest
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Test retest',fig_num),'units','normalized','outerposition',[0 0 0.5 0.8])
ha=tight_subplot(1, 2, 0.1, [0.1 0.11], [.09 .05]);

trt_vars = {STD,NMA};
trt_titles = {'Standard deviation','Mean asynchrony'};
lims = [-2 4];
for i=1:2
    tmp = (trt_vars{i}-mean(trt_vars{i}(grp==1,:,:)))./std(trt_vars{i}(grp==1,:,:));
    combo_measure = squeeze(nanmean(tmp,2));
    axes(ha(i)) 
    plot(lims,lims,'--','Color',[0.5 0.5 0.5]);hold on
    xlabel('First assesment (z-score, all blocks)','FontSize',font_size.lab)
    ylabel('Second assesment (z-score, all blocks)','FontSize',font_size.lab)
    ylim(lims)
    xlim(lims)
    plts = [];txt = [];
    for pp = 1:4
        if pp<4
            x_var = combo_measure(grp==pp,1);
            y_var = combo_measure(grp==pp,2);
            col = symbols.c{pp};
            scatter(x_var,y_var,symbols.s{pp},'filled','MarkerEdgeColor','k');hold on
            r = corr(x_var,y_var,'Type','Spearman','Rows','Complete');
            txt = [txt, newline, '\rho ',leg.clean{pp},sprintf(' = %0.2g',r)];
        else
            x_var = combo_measure(:,1); y_var = combo_measure(:,2);
            col = 'k';
            r = corr(x_var,y_var,'Type','Spearman','Rows','Complete');
            txt = [txt, newline, '\rho ALL',sprintf(' = %0.2g',r)];
        end
        coefs=regress(y_var,[ones(length(x_var),1),x_var]);
        plts(pp)=plot(lims,[lims(1)*coefs(2)+coefs(1),lims(2)*coefs(2)+coefs(1)],'Color',col,'LineWidth',1.5);
    end
    text(0.05,0.95,ltrs(i),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
    text(0.05,0.87,txt,'Units', 'Normalized','FontSize',font_size.txt,'HorizontalAlignment','left')
    grid on
    title(trt_titles{i})
    set(gca, 'FontSize',font_size.tick);
end 
legend(plts,[leg.size,'ALL (N=109)'],'Location','southeast');
[~,h]=suplabel('Test retest reliability','t',[.1 .09 .86 .86]);set(h,'FontSize',font_size.tit)

% All blocks SD and NMA
stat_style = 'KW';
fig_num = fig_num+1;
for i=1:2
    if i==1
        figure_title = 'Standard deviation';
        variable = mean(STD,3);
        ylimit = [10,135];
        ticks = 20:20:120;
    else
        figure_title = 'Mean asynchrony';
        variable = mean(NMA,3);
        ylimit = [-150,30];
        ticks = -140:20:20;
    end

    figure('Name',sprintf('[%d][Supp.] %s - Conditions',fig_num,figure_title),'units','normalized','outerposition',[0 0 1 1])
    ha=tight_subplot(2, 3, [0.045 0.006], [0.05 0.075], [.05 .015]);
    [~,h]=suplabel(figure_title,'t',[.08 .08 .875 .893]);
    set(h,'FontSize',font_size.tit+2)
    for blk = 1:6
        axes(ha(blk))
        Single_Participant_Graph([variable(:,blk) grp],leg.size,stat_style)
        tit = sprintf('Step size: %d ms',diff(tempos(blk,:)));
        title(tit,'FontSize',font_size.tit-2)
        text(0.05,0.95,ltrs(blk),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
        ylim(ylimit);
        if blk<4
            set(gca, 'FontSize',font_size.tick,'XTickLabels','');
        end
        if mod(blk,3)~=1
            set(gca, 'FontSize',font_size.tick,'YTickLabels','','YTick',ticks);
        else
            set(gca, 'FontSize',font_size.tick,'YTick',ticks);
            ylabel('ms','FontSize',font_size.lab)
        end
        grid on;
    end
end

% Autoregressive model
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Autoregressive model',fig_num),'units','normalized','outerposition',[0 0 0.5 0.6])
ha=tight_subplot(1, 2, 0.1, [0.15 0.15], [.08 .05]);

axes(ha(1))
[counts,~,~,labels] = crosstab(hierarchical_num_coeffs,grp);
ba=bar(cellfun(@(x) str2double(x),labels(:,1)), counts);
for pp=1:3
    ba(pp).FaceColor = symbols.c{pp};
end
title('Hierarchical regression','FontSize',font_size.tit)
ylabel('Number of subjects','FontSize',font_size.lab)
xlabel('Number of regressors in final model','FontSize',font_size.lab)
legend(leg.size)
ylim([0 36])
set(gca, 'XTick', 1:5, 'XTickLabels', {'1','2','3','4','\geq 5'},'FontSize',font_size.tick);
text(0.05,0.95,ltrs(1),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')

axes(ha(2))
significance= {'***','**','*'};
thresh = [0.001,0.01,0.05,20];
num_coeffs = size(autoreg_coeffs,2);
plot([0.5 num_coeffs+0.5],[0 0],'--'); hold on
for pp=1:3
    y_values = autoreg_coeffs(grp==pp,:);
    x_values = repmat(1:num_coeffs,g_size(pp),1)+(rand(size(y_values))-0.5)/6 +(2*pp-5)/12;
    p(pp)=scatter(x_values(:),y_values(:),40,symbols.s{pp},'filled','MarkerEdgeColor','k');
    hold on
    errorbar((1:num_coeffs)+1/4,mean(y_values),std(y_values)/sqrt(g_size(pp)),symbols.sl{pp},'MarkerFaceColor',symbols.c{pp}, 'MarkerEdgeColor','k','LineWidth',0.8,'CapSize',12);
end
for coef=1:num_coeffs
    [~,tbl] = anovan(autoreg_coeffs(:,coef),{grp},'varnames',{'group'},'Display','off');
    j=find(thresh==min(thresh(thresh-tbl{2,7} > 0)));
    if j<4
        text(coef+1/4,1,significance{j},'FontSize',font_size.txt+4,'HorizontalAlignment','center','FontWeight','bold');
    end
end
xlim([0.5 num_coeffs+0.5])
set(gca, 'XTick', 1:num_coeffs, 'XTickLabels', {'b1','b2','b3','b4'},'FontSize',font_size.tick); %will need change for more than 4 coeffs
grid on
box off
title('Four predictor model','FontSize',font_size.tit)
ylabel('Values','FontSize',font_size.lab)
xlabel('Coefficient','FontSize',font_size.lab)
legend(p,leg.size)
text(0.05,0.95,ltrs(2),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')

[~,h]=suplabel('Autoregressive Model','t',[.08 .08 .84 .88]);
set(h,'FontSize',font_size.tit)

% Step sizes 10-30
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Change Dynamics (10-30)',fig_num),'units','normalized','outerposition',[0 0 0.5 0.8])
ha=tight_subplot(2, 2, [0.025 0.01], 0.09, [.14 .01]);
p=nan(1,2);lims =500+ [-35,35;-35,35];
for b=2:3
    for d=1:2
        axes(ha(d+2*(3-b)));
        tempo = tempos(b,:);
        if d==1;tempo = fliplr(tempo);end
        plot([-0.5+rng_t(1),-1,0,rng_t(end)+0.5],[tempo(1),tempo(1),tempo(2),tempo(2)],'k--');hold all;
        box off
        d_val_all = d_segments_mean(:,:,b-1,d);
        for pp=1:3
            d_val_group = d_val_all(grp==pp,:)-mean(d_val_all(grp==pp,rng_t<0),2)+tempo(1);
            count_part = sum(~all(isnan(d_val_group),2));
            p(pp)=errorbar(rng_t,nanmean(d_val_group),nanstd(d_val_group)/sqrt(count_part),symbols.sl{pp},'MarkerFaceColor',symbols.c{pp},'MarkerEdgeColor','k','LineWidth',1.25);hold all;
        end
        text(-0.9*d+1.85,0.92,ltrs(d+2*(3-b)),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
        xlim([-0.5+rng_t(1),rng_t(end)+0.5]);ylim(lims(b-1,:))
        if d==1
            ylabel(sprintf('%s %d ms','Step size: ',-diff(tempo)),'FontSize',font_size.tit+2);
            set(gca, 'XTick', rng_t,'XTickLabels','','FontSize',font_size.tick);
        else
            set(gca,'XTick', rng_t,'XTickLabels','','YTickLabels', '', 'FontSize',font_size.tick);
        end
        if b==2
            set(gca, 'XTickLabels', rng_t,'FontSize',font_size.tick);
            xlabel('Beat','FontSize',font_size.lab)
        elseif b==3
            title(direcs{d},'FontSize',font_size.tit)
        end
        grid on      
    end
end
legend(ha(4),p,leg.clean,'Location','northeast')
[~,h]=suplabel('Change dynamics' ,'t', [0.14 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit+5)
[~,h]=suplabel('Delay interval (baselined, ms)' ,'y', [0.09 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit+3)

% Parameter recovery
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Parameter recovery',fig_num),'units','normalized','outerposition',[0 0 1 0.6])
title_panels = {'Experiment 1:','Experiment 2:','Experiment 2:','Combined update rate ';...
    'phase correction (\alpha)','phase correction (\alpha) z-score','period correction (\beta) z-score','(\alpha and \beta z-score)'};
ha=tight_subplot(1, 4, 0.01, [0.12 0.14], [.05 .03]);
limits_all = [0 1; -3 3; -3 3; -3 3];
fitted = [wing_isoch_av(:,1), zwing_changes(:,1), zwing_changes(:,2), model_params_combined];
recovered = [mean_recovery_isoch(:,1), mean_zrecovery_changes(:,1), mean_zrecovery_changes(:,2), mean_model_params_combined_recovery];
recovered_sem = [std_recovery_isoch(:,1), std_zrecovery_changes(:,1), std_zrecovery_changes(:,2), std_model_params_combined_recovery]/sqrt(n_rep);
for i=1:4
    axes(ha(i))
    limits = limits_all(i,:);
    errorbar(fitted(:,i), recovered(:,i),recovered_sem(:,i),recovered_sem(:,i),'o','MarkerFaceColor','auto');hold on
    plot(limits,limits,'k--');hold on
    set(gca, 'FontSize',font_size.tick,'YTickLabel','');
    text(0.05,0.95,ltrs(i),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
    r=corr(fitted(:,i),recovered(:,i),'Type','Spearman','Rows','Complete');
    txt = ['\rho ALL',sprintf(' = %0.2g',r)];
    text(0.7, 0.06 ,txt ,'Units','normalized', 'FontSize', font_size.txt+1)
    title(sprintf('%s\n%s',title_panels{1,i},title_panels{2,i}),'FontSize',font_size.tit)
    xlim(limits);ylim(limits);
    box off;grid on
end
set(ha(1),'YTickLabelMode','auto');
ylabel(ha(1),'Recovered value','FontSize',font_size.lab)
[~,h]=suplabel('Simulation input value' ,'x', [0.05 0.12 0.92 0.9]);set(h,'FontSize',font_size.lab)

% Simulated consecutive correlations
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Simulated consecutive correlations',fig_num),'units','normalized','outerposition',[0 0 1 0.65])
ha=tight_subplot(1, 3, 0.05, [0.12 0.16], [.05 .03]);
title_panels = {'Relation of experiment and simulated data','Simulated data','Experimental data'};

raw_data = correls_all;
simulated_data = nanmean(correls_all_sim,2);lim = [0 1];
exc_subjects = any(isnan(wing_isoch_av),2);
plts = [];
for i = 1:3
    axes(ha(i))
    if i==3
        Single_Participant_Graph([raw_data(~exc_subjects) grp(~exc_subjects)],leg.clean,'KW');
        title('Experimental data','FontSize',font_size.tit);
    elseif i==2
        Single_Participant_Graph([simulated_data grp],leg.clean,'KW');
        title('Simulated data','FontSize',font_size.tit);
    else
        for pp = 1:3
            scatter(simulated_data(grp==pp),raw_data(grp==pp),symbols.s{pp},'filled','MarkerEdgeColor','k');hold on
        end
        plot(lim,lim,'k--')
        r=corr(simulated_data,raw_data,'Type','Spearman','Rows','Complete');
        txt = ['\rho ',sprintf(' = %0.3g',r)];
        text(0.8, 0.28 ,txt ,'Units','normalized', 'FontSize', font_size.txt+2)
        ylabel('Experimental data','FontSize',font_size.lab)
        xlabel('Simulated data','FontSize',font_size.lab)
        xlim(lim);
        legend(plts, leg.clean, 'location','southeast')
    end
    set(gca, 'FontSize',font_size.tick);
    title(title_panels{i},'FontSize',font_size.tit)
    ylim(lim)
    grid on
    text(0.05,0.95,ltrs(i),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
end
[~,h]=suplabel('Correlation of consecutive asynchronies' ,'t', [0.05 0.12 0.92 0.83]);set(h,'FontSize',font_size.tit+3)

% Simulated change dynamics
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Simulated Change Dynamics',fig_num),'units','normalized','outerposition',[0 0 0.5 1])
ha=tight_subplot(3, 2, [0.025 0.01], 0.075, [.14 .01]);
p=nan(1,3);lims =500+[-70,70;-70,70;-70,70];
direcs = {'Accelerating','Decelerating'};
for b=4:6
    for d=1:2
        axes(ha(d+2*(6-b)));
        tempo = tempos(b,:);
        if d==1;tempo = fliplr(tempo);end
        plot([-0.5+rng_t(1),-1,0,rng_t(end)+0.5],[tempo(1),tempo(1),tempo(2),tempo(2)],'k--');hold all;
        d_val_all = d_segments_mean_simulated(:,:,b-3,d);
        box off
        for pp=1:3
            d_val_group = d_val_all(grp==pp,:)-mean(d_val_all(grp==pp,rng_t<0),2)+tempo(1);
            count_part = sum(~all(isnan(d_val_group),2));
            p(pp)=errorbar(rng_t,nanmean(d_val_group),nanstd(d_val_group)/sqrt(count_part),symbols.sl{pp},'MarkerFaceColor',symbols.c{pp},'MarkerEdgeColor','k','LineWidth',1.25);hold all;
        end
        text(-0.9*d+1.85,0.92,ltrs(d+2*(6-b)),'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold')
        xlim([-0.5+rng_t(1),rng_t(end)+0.5]);ylim(lims(b-3,:))
        if d==1
            ylabel(sprintf('%s %d ms','Step size: ',-diff(tempo)),'FontSize',font_size.tit+2);
            set(gca, 'XTick', rng_t,'XTickLabels','','FontSize',font_size.tick);
        else
            set(gca,'XTick', rng_t,'XTickLabels','','YTickLabels', '', 'FontSize',font_size.tick);
        end
        if b==4
            set(gca, 'XTickLabels', rng_t,'FontSize',font_size.tick);
            xlabel('Beat','FontSize',font_size.lab)
        elseif b==6
            title(direcs{d},'FontSize',font_size.tit)
        end
        grid on      
    end
end
legend(ha(6),p,leg.clean,'Location','southeast')
[~,h]=suplabel('Simulated change dynamics' ,'t', [0.14 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit+5)
[~,h]=suplabel('Delay interval (baselined, ms)' ,'y', [0.09 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit+3)

%% save all

FolderName = 'C:\Users\user\Google Drive\Tapping_Project\Figures\';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = get(FigHandle, 'Name');
    saveas(FigHandle, [FolderName, FigName, '.jpg']);
    saveas(FigHandle, [FolderName, FigName, '.svg']);
end