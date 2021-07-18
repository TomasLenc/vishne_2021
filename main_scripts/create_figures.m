% load relevant stuff
close all;clc;clear all;
CODE_DIR='C:\Users\user\Google Drive\Tapping_Project';
cd(CODE_DIR);
load('raw_data.mat')
load('analyses_data.mat')
load('simulation_data.mat')

ltr_args = {'Units', 'Normalized','FontSize',font_size.ltr,'HorizontalAlignment','center','FontWeight', 'Bold'};

%% Experiment 1

% Figure 1: NMA, SD 
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Isoch.] explanation, NMA, STD',fig_num),'units','centimeters','outerposition',[0 0 8.8 10],'Renderer','Painters')

axp=[0.08 0.62 0.89 0.335];
ha(1)=axes('Units','normalized','Position',axp);
ha(2)=axes('Units','normalized','Position',[0.08 0.075 0.42 0.44]);
ha(3)=axes('Units','normalized','Position',[0.56 0.075 0.42 0.44]);

axes(ha(1))
t=title('Task structure','FontSize',font_size.tit);
text(0.02,0.95,ltrs(1),ltr_args{:})

titles = {'Mean Asynchrony', 'Standard deviation'};
for i=1:2
    if i==1
        dat = mean(NMA(:,1,:),3);
    else
        dat = mean(STD(:,1,:),3);
    end
    axes(ha(i+1));
    Single_Participant_Graph([dat grp] ,leg.clean, symbols, 'KW')
    title(titles{i},'FontSize',font_size.tit)
    text(0.05,0.95,ltrs(i+1),ltr_args{:})
end
ylabel(ha(2),'ms','FontSize',font_size.lab)

% 2 - Consecutive correlations
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Isoch.] Consec. Hist Phase',fig_num),'units','centimeters','outerposition',[0 0 18 7.5],'Renderer','Painters')
lims = [-240,300];bins=55;
nw=4;marg_w=[.05 .025];marg_h=[0.12 0.12];gap_w=[0.008 0.008 0.045 0];
w_tmp=(1-sum(marg_w)-sum(gap_w))/4.9;size_w=[1.3*w_tmp 1.3*w_tmp 1.3*w_tmp w_tmp];
axh = (1-sum(marg_h));py = 1-marg_h(1)-axh;px = marg_w(1);
for i=1:nw
    axw = size_w(i);
    ha(i) = axes('Units','normalized','Position',[px py axw axh]);
    px = px+axw+gap_w(i);
end

for pp=1:3
    col_steps = linspace(0,1,64);
    map = (1-col_steps').*repmat(symbols.c{pp},64,1) + (col_steps').*repmat([1 1 1],64,1);
    ax=ha(pp);axes(ax)
    asynch_group = isoch_all_e{pp};
    coefs=regress(asynch_group(:,2),[ones(length(asynch_group(:,1)),1),asynch_group(:,1)]);
    histogram2(asynch_group(:,1),asynch_group(:,2),bins,'DisplayStyle','tile','LineStyle','none');hold all;
    text(0.05,0.95,ltrs(pp),ltr_args{:})
    plot(lims,[lims(1)*coefs(2,1)+coefs(1,1),lims(2)*coefs(2,1)+coefs(1,1)],'k','LineWidth',1.25)
    xlim(lims);ylim(lims)
    caxis([0 165])
    title(leg.size{pp},'FontSize',font_size.lab)
    colormap(ax,map)
    set(gca, 'FontSize',font_size.tick,'YTickLabels','','TickLength',[0.006 0.006]);
    box off;grid off
end
axes(ha(4))
Single_Participant_Graph([correls_all grp],leg.clean,symbols,'KW')
title('Single participant correlations','FontSize',font_size.lab)
text(0.06,0.95,ltrs(4),ltr_args{:})
ylabel('Correlations','FontSize',font_size.lab)
ylim([0 1])
set(gca,'FontSize',font_size.tick,'TickLength',[0.006 0.006])
ylabel(ha(1),'Asynchrony {\itt+1} (ms)','FontSize',font_size.lab);
set(ha(1), 'FontSize',font_size.tick,'YTickLabelMode','auto','TickLength',[0.006 0.006]);
locs=[marg_w(1) marg_h(1)-0.01 sum(size_w(1:3))+sum(gap_w(1:2)) 0.86];
xlabel(ha(2),'Asynchrony {\itt} (ms)','FontSize',font_size.lab)
[~,h]=suplabel('Relation of consecutive asynchronies' ,'t', locs);set(h,'FontSize',font_size.tit)

% 3 - computational model
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Isoch.] Linear Sensorimotor Synchronization Model (bGLS)',fig_num),'units','centimeters','outerposition',[0 0 18 12],'Renderer','Painters')
ha=tight_subplot(1, 3, 0.05, [0.075 0.5], [.05 .02]);

titles = {'Phase correction (\alpha)', 'Timekeeper noise', 'Motor noise'};
ylabels = {'arb. units', 'ms','ms'};
for i=1:3
    data = wing_isoch_av(:,i);
    tit = titles{i};ylb = ylabels{i};
    axes(ha(i))
    Single_Participant_Graph([data grp],leg.clean,symbols,'KW')
    text(0.05,0.95,ltrs(i+1),ltr_args{:})
    title(tit,'FontSize',font_size.tit)
    ylabel(ylb,'FontSize',font_size.lab)
    set(gca,'FontSize',font_size.tick,'TickLength',[0.006 0.006])
    if i==2 || i==3
        ylim([0, 5*ceil(max(wing_isoch_av(:,2))/5)])
    end
end

%% Experiment 2

% 4 Change dynamics
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Alter.] Change Dynamics',fig_num),'units','centimeters','outerposition',[0 0 8.8 13],'Renderer','Painters')
ha=tight_subplot(3, 2, [0.025 0.01], 0.075, [.13 .01]);
p=nan(1,3);lims =500+[-70,70];
direcs = {'Accelerating','Decelerating'};
for b=4:6
    for d=1:2
        axes(ha(d+2*(6-b)));
        tempo = tempos(b,:);
        if d==1;tempo = fliplr(tempo);end
        plot([-0.5+rng_t(1),-1,0,rng_t(end)+0.5],[tempo(1),tempo(1),tempo(2),tempo(2)],'k--');hold all;
        d_val_all = d_segments_mean(:,:,b-1,d);
        for pp=1:3
            d_val_group = d_val_all(grp==pp,:)-mean(d_val_all(grp==pp,rng_t<0),2)+tempo(1);
            count_part = sum(~all(isnan(d_val_group),2));
            p(pp)=errorbar(rng_t,nanmean(d_val_group),nanstd(d_val_group)/sqrt(count_part),symbols.s{pp},'MarkerFaceColor',symbols.c{pp},...
               'LineWidth',0.6,'LineStyle','-','Color',symbols.c{pp},'MarkerSize',2.5, 'CapSize', 4);hold all;
        end
        text(-0.9*d+1.85,0.92,ltrs(d+2*(6-b)),ltr_args{:})
        xlim([-0.5+rng_t(1),rng_t(end)+0.5]);ylim(lims)
        if d==1
            ylabel(sprintf('%s %d ms','Step size: ',-diff(tempo)),'FontSize',font_size.tit);
            set(gca, 'XTick', rng_t,'XTickLabels','','FontSize',font_size.tick,'TickLength',[0.006 0.006]);
        else
            set(gca,'XTick', rng_t,'XTickLabels','','YTickLabels', '', 'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
        end
        if b==4
            set(gca, 'XTickLabels', rng_t,'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
            xlabel('Beat','FontSize',font_size.lab)
        elseif b==6
            title(direcs{d},'FontSize',font_size.tit)
        end
        box off
    end
end
[ll,icons] = legend(ha(6),p,leg.clean,'Location','southeast');
for obj = 1:3
    icons(obj).Position = icons(obj).Position + [0.1 0 0];
end
for obj = 4:6
    icons(obj).Children.Children(1).XData = icons(obj).Children.Children(1).XData + 0.14;
    icons(obj).Children.Children(2).XData = icons(obj).Children.Children(2).XData + [0.2 0.08];
end
ll.Box = 'off';

[~,h]=suplabel('Change dynamics' ,'t', [0.14 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit)
[~,h]=suplabel('Delay interval (baselined, ms)' ,'y', [0.12 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit)

% 5 Signal detection theory
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Alter.] Signal Detection Graphs',fig_num),'units','centimeters','outerposition',[0 0 18 14],'Renderer','Painters')
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
        plot([nanmean(d1),nanmean(d1)],[0,1],'Color',ones(1,3)*0.02,'LineWidth',0.75);hold on
        plot([nanmean(d2),nanmean(d2)],[0,1],'Color',ones(1,3)*0.02,'LineWidth',0.75);hold on
        histogram(d2,'FaceColor',symbols.c{pp},'BinEdges',240:10:750,'Normalization','probability','EdgeColor','none');hold on
        histogram(d1,'FaceColor',symbols.c{pp}/3,'BinEdges',240:10:750,'Normalization','probability','EdgeColor','none');hold on
        xlim([240,750]);ylim([0,0.075]);box off
        txt1 = [sprintf("d' = %0.2f",grp_signal_detection(pp,blk-1,1)); sprintf("Diff. means = %0.1f",grp_signal_detection(pp,blk-1,3))];
        text(0.58,0.9,txt1,'Units','normalized','FontSize',font_size.txt);
        text(0.05,0.92,ltrs((6-blk)*4+pp),ltr_args{:});
        if blk==4
            set(gca, 'FontSize',font_size.tick,'YTickLabels','','TickLength',[0.006 0.006]);
        else
            set(gca, 'FontSize',font_size.tick,'YTick',0:0.01:0.07,'YTickLabels','','XTickLabels','','TickLength',[0.006 0.006]);
        end
        if blk==6
            title(leg.size{pp},'FontSize',font_size.tit)
        end
        axes(ha((7-blk)*4));
        text(0.05,0.92,ltrs((7-blk)*4),ltr_args{:})
        [X, Y]= perfcurve([ones(length(d1),1);2*ones(length(d2),1)],[d1;d2], 2);
        plt(pp)=plot(X,Y,'Color',symbols.c{pp},'LineWidth',1);hold on;
        box off
        leg_auc{pp}=[leg.clean{pp},sprintf(': %0.2g',grp_signal_detection(pp,blk-1,2))];
        if blk==6
            title('ROC curve','FontSize',font_size.tit)
        end
        ylim([0 1])
        set(gca,'XTick',0:0.25:1,'YTick',0:0.25:1,'TickLength',[0.006 0.006]);
    end
    [legu, icons]= legend(plt,leg_auc,'Location','southeast','FontSize',font_size.txt);
    legu.Box = 'off';
    for obj = 1:3
        icons(obj*2 + 2).XData = icons(obj*2 + 2).XData + [0.3 0];
    end
    text(0.59, 0.35, 'AUC:', 'FontSize', font_size.txt)
    if blk==4
        set(gca,'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
    else
        set(gca,'FontSize',font_size.tick,'XTickLabels','','TickLength',[0.006 0.006]);
    end
    ylabel(ha((6-blk)*4+1),['Step size: ',num2str(diff(tempos(blk,:))),'ms'],'FontSize',font_size.tit)
    set(ha((6-blk)*4+1),'YTick',0:0.01:0.07,'YTickLabelMode','auto','TickLength',[0.006 0.006]);
end
xlabel(ha(10),'Delay interval (ms)','FontSize',font_size.lab)
xlabel(ha(12),'False positive rate (FPR)','FontSize',font_size.lab)
ylabel(ha(8),'True positive rate (TPR)','FontSize',font_size.lab)

% 6 Model with tempo changes
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Alter.] Linear Sensorimotor Synchronization Model (bGLS)',fig_num),'units','centimeters','outerposition',[0 0 18 13],'Renderer','Painters')
ha=tight_subplot(2, 3, [0.07 0.03], 0.075, [.04 .02]);

axes(ha(1))
text(0.05,0.95,ltrs(1),ltr_args{:})
title('Timekeeper dynamics','FontSize',font_size.tit)
    
titles={'Phase correction (\alpha)','Period correction (\beta)','Timekeeper noise','Motor noise'};
for i=1:4
    j=i;if i<=2;j=3-i;end
    data =zwing_changes(:,i);
    axes(ha(j+2))
    Single_Participant_Graph([data grp],leg.clean,symbols,'KW')
    text(0.05,0.95,ltrs(j+1),ltr_args{:})
    title(titles{i},'FontSize',font_size.tit)
    set(gca,'FontSize',font_size.tick,'TickLength',[0.006 0.006])
    if i==3 || i==4
        ylim([-2 5])
    end
    ylabel(ha(j+2),'z-score (arb. units)','FontSize',font_size.lab)
end

% 7 Combined update
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Both.] alpha beta corr',fig_num),'units','centimeters','outerposition',[0 0 18 7],'Renderer','Painters')
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
    plts(pp)= scatter(period_correction, phase_correction, 17,symbols.c{pp},'filled','Marker',symbols.s{pp},'MarkerFaceAlpha',0.45);hold on;
    plot([xlims(1),xlims(2)],[xlims(1)*coefs(2)+coefs(1),xlims(2)*coefs(2)+coefs(1)],'Color',symbols.c{pp},'LineWidth',1.25)
    text(0.04, 0.06 ,txt ,'Units','normalized', 'FontSize', font_size.txt)
    set(gca, 'FontSize',font_size.tick,'YTickLabels','','TickLength',[0.006 0.006]);
    xlim(xlims);ylim(ylims)
    text(0.05,0.95,ltrs(pp),ltr_args{:})
end
set(ha(1),'FontSize',font_size.tick,'YTickLabelMode','auto','TickLength',[0.006 0.006]);
ylabel(ha(1),titles{1},'FontSize',font_size.lab);
xlabel(ha(2),titles{2},'FontSize',font_size.lab);
legu = legend(ha(3),plts,leg.clean,'Location','northeast','LineWidth',0.25);
legu.Position(2) = legu.Position(2) + 0.03;

ha(4) = axes('Units','normalized', 'Position',[0.77 0.15 0.2167 0.8], 'XTickLabel','', 'YTickLabel','');
axes(ha(4));
Single_Participant_Graph([model_params_combined grp],leg.clean,symbols,'KW')
text(0.05,0.95,ltrs(4),ltr_args{:})
ylabel('Update rate (\alpha and \beta z-score)','FontSize',font_size.lab)
ylim([-1.7 2])

% 8 AQ correlations
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Both.] aq update rate',fig_num),'units','centimeters','outerposition',[0 0 18 7],'Renderer','Painters')
titles = {'Update rate (\alpha and \beta z-score)','Communication and mindreading skills'};
ha(1) = axes('Units','normalized', 'Position',[0.04 0.15 0.22 0.8], 'XTickLabel','', 'YTickLabel','');
ha(2:4)=tight_subplot(1,3,0.005,[0.15,0.05],[0.31,0.02]);

axes(ha(1))
histogram(aq_data_austin(aq_data.Group==1,3),'FaceColor',symbols.c{1},'FaceAlpha',0.4,'EdgeAlpha',0);hold on
histogram(aq_data_austin(aq_data.Group==3,3),'EdgeColor',symbols.c{3},'FaceAlpha',0,'EdgeAlpha',0.4,'LineWidth',1.2);hold on
set(gca, 'FontSize',font_size.tick,'TickLength',[0.006 0.006]);box off
legend(leg.clean{[1,3]},'LineWidth',0.25)
text(0.05,0.95,ltrs(1),ltr_args{:})
ylabel('Number of subjects','FontSize',font_size.lab)
txt = [sprintf("Wilcoxon rank sum \ntest, p = %0.2g",ranksum(aq_data_austin(aq_data.Group==3,3),aq_data_austin(aq_data.Group==1,3))); ...
    sprintf("Cliff's Delta = %0.2g",Cliffs_Delta(aq_data_austin(aq_data.Group==3,3),aq_data_austin(aq_data.Group==1,3)))];
text(0.51,0.65,txt,'Units','normalized','FontSize',font_size.txt);

aq_scale = aq_data_austin(:,3); update_rate =model_params_combined(aq_data.Subject_Index); %Index
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
        txt = ['\rho ALL',sprintf(' = %0.2g\np = %0.2g one-sided',r,p/2)];
        plot([xlims(1),xlims(2)],[xlims(1)*coefs(2)+coefs(1),xlims(2)*coefs(2)+coefs(1)],'k-','LineWidth',1.25)
    else
        for i = [4, ind]  
            axes(ha(i))
            if pp == 1
                plt(ind-1) =scatter(var1, var2, 17,symbols.c{pp},'filled','Marker',symbols.s{pp},'MarkerFaceAlpha',0.45);hold on;
            else
                plt(ind-1) =scatter(var1, var2, 17,symbols.c{pp},'Marker',symbols.s{pp},'MarkerEdgeAlpha',0.45,'LineWidth',1.2);hold on;
            end
        end
    
        txt = ['\rho ',leg.clean{pp},sprintf(' = %0.2g\np = %0.2g one-sided',r,p/2)];
        plot([xlims(1),xlims(2)],[xlims(1)*coefs(2)+coefs(1),xlims(2)*coefs(2)+coefs(1)],'Color',symbols.c{pp},'LineWidth',1.25);
    end
    text(0.04, 0.084 ,txt ,'Units','normalized', 'FontSize', font_size.txt)
    set(gca, 'FontSize',font_size.tick,'YTickLabels','','TickLength',[0.006 0.006]);
    xlim(xlims);ylim(ylims)
    text(0.05,0.95,ltrs(ind),ltr_args{:})
end
set(ha(2),'FontSize',font_size.tick,'YTickLabelMode','auto','TickLength',[0.006 0.006]);
ylabel(ha(2),titles{1},'FontSize',font_size.lab);
[~,h]=suplabel(titles{2},'x',[.05 .18 .93 .8]);
set(h,'FontSize',font_size.lab)
legend(ha(4),plt,leg.clean{[1,3]},'Location','northeast','LineWidth',0.25)

%% Supplementary

% Test-retest
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Test retest',fig_num),'units','centimeters','outerposition',[0 0 8.8 7.4],'Renderer','Painters')
ha=tight_subplot(1, 2, 0.1, [0.15 0.15], [.09 .05]);

trt_vars = {STD,NMA};
trt_titles = {'Standard deviation','Mean asynchrony'};
limits = [-2 6; -4 4];
for i=1:2
    tmp = (trt_vars{i}-mean(trt_vars{i}(grp==1,:,:)))./std(trt_vars{i}(grp==1,:,:));
    combo_measure = squeeze(nanmean(tmp,2));
    axes(ha(i)) 
    lims = limits(i,:);
    plts = [];txt = [];
    for pp = 1:4
        if pp<4
            x_var = combo_measure(grp==pp,1);
            y_var = combo_measure(grp==pp,2);
            col = symbols.c{pp};
            scatter(x_var,y_var,17,symbols.c{pp},'filled','Marker',symbols.s{pp},'MarkerFaceAlpha',0.45);hold on
            r = corr(x_var,y_var,'Type','Spearman','Rows','Complete');
            txt = [txt, newline, '\rho ',leg.clean{pp},sprintf(' = %0.2g',r)];
        else
            x_var = combo_measure(:,1); y_var = combo_measure(:,2);
            col = 'k';
            r = corr(x_var,y_var,'Type','Spearman','Rows','Complete');
            txt = [txt, newline, '\rho ALL',sprintf(' = %0.2g',r)];
        end
        coefs=regress(y_var,[ones(length(x_var),1),x_var]);
        plts(pp)=plot(lims,[lims(1)*coefs(2)+coefs(1),lims(2)*coefs(2)+coefs(1)],'Color',col,'LineWidth',1);hold on;
    end
    plot(lims,lims,'--','Color',[0.3 0.3 0.3],'LineWidth',1);hold on
    text(0.05,0.95,ltrs(i),ltr_args{:})
    text(0.05,0.83,txt,'Units', 'Normalized','FontSize',font_size.txt,'HorizontalAlignment','left')
    ylim(lims);xlim(lims)
    title(trt_titles{i})
    xlabel('First assesment (z-score, all blocks)','FontSize',font_size.lab)
    ylabel('Second assesment (z-score, all blocks)','FontSize',font_size.lab)
    set(gca, 'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
end 
[legu, icons] = legend(plts,[leg.clean,'ALL'],'Location','southeast','box','off');
for obj = 1:4
    icons(obj*2 + 3).XData = icons(obj*2 + 3).XData + [0.35 0.08];
    icons(obj).Position = icons(obj).Position + [0.08 0 0];
end
[~,h]=suplabel('Test retest reliability','t',[.1 .09 .86 .86]);set(h,'FontSize',font_size.tit)

% All blocks SD and NMA
stat_style = 'KW';
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
    fig_num = fig_num+1;
    figure('Name',sprintf('[%d][Supp.] %s - Conditions',fig_num,figure_title),'units','centimeters','outerposition',[0 0 18 12],'Renderer','Painters')
    ha=tight_subplot(2, 3, [0.045 0.006], [0.05 0.075], [.05 .015]);
    [~,h]=suplabel(figure_title,'t',[.08 .08 .875 .893]);
    set(h,'FontSize',font_size.tit)
    for blk = 1:6
        axes(ha(blk))
        Single_Participant_Graph([variable(:,blk) grp],leg.size,symbols,stat_style)
        tit = sprintf('Step size: %d ms',diff(tempos(blk,:)));
        title(tit,'FontSize',font_size.tit-2)
        text(0.05,0.95,ltrs(blk),ltr_args{:})
        ylim(ylimit);
        if blk<4
            set(gca, 'FontSize',font_size.tick,'XTickLabels','','TickLength',[0.006 0.006]);
        end
        if mod(blk,3)~=1
            set(gca, 'FontSize',font_size.tick,'YTickLabels','','YTick',ticks,'TickLength',[0.006 0.006]);
        else
            set(gca, 'FontSize',font_size.tick,'YTick',ticks,'TickLength',[0.006 0.006]);
            ylabel('ms','FontSize',font_size.lab)
        end
    end
end

% Autoregressive model
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Autoregressive model',fig_num),'units','centimeters'	,'outerposition',[0 0 8.8 8],'Renderer','Painters')
ha=tight_subplot(1, 2, 0.07, [0.12 0.12], [.08 .05]);

axes(ha(1))
[counts,~,~,labels] = crosstab(hierarchical_num_coeffs,grp);
ba=bar(cellfun(@(x) str2double(x),labels(:,1)), counts, 'EdgeColor', 'none');
for pp=1:3
    ba(pp).FaceColor = symbols.c{pp};
end
box off
title('Hierarchical regression','FontSize',font_size.tit)
ylabel('Number of subjects','FontSize',font_size.lab)
xlabel('Number of regressors in final model','FontSize',font_size.lab)
legend(leg.size,'LineWidth',0.25)
ylim([0 36])
set(gca, 'XTick', 1:5, 'XTickLabels', {'1','2','3','4','\geq 5'},'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
text(0.05,0.95,ltrs(1),ltr_args{:})

axes(ha(2))
significance= {'***','**','*'};
thresh = [0.001,0.01,0.05,20];
num_coeffs = size(autoreg_coeffs,2);
plot([0.5 num_coeffs+0.5],[0 0],'k--','LineWidth',0.25); hold on
for pp=1:3
    y_values = autoreg_coeffs(grp==pp,:);
    x_values = repmat(1:num_coeffs,g_size(pp),1)+(rand(size(y_values))-0.5)/6 +(2*pp-5)/12;
    p(pp)=scatter(x_values(:),y_values(:),11,symbols.c{pp},'filled','Marker',symbols.s{pp},'MarkerFaceAlpha',0.45);
    hold on
    errorbar((1:num_coeffs)+1/4,mean(y_values),std(y_values)/sqrt(g_size(pp)),symbols.s{pp},'MarkerFaceColor',symbols.c{pp},...
        'LineStyle','-','Color',symbols.c{pp},'MarkerSize',2,'LineWidth',0.5,'CapSize',4);
end
for coef=1:num_coeffs
    [~,tbl] = anovan(autoreg_coeffs(:,coef),{grp},'varnames',{'group'},'Display','off');
    j=find(thresh==min(thresh(thresh-tbl{2,7} > 0)));
    if j<4
        text(coef+1/4,1,significance{j},'FontSize',9,'HorizontalAlignment','center');
    end
end
xlim([0.5 num_coeffs+0.5])
set(gca, 'XTick', 1:num_coeffs, 'XTickLabels', {'b1','b2','b3','b4'},'FontSize',font_size.tick,'TickLength',[0.006 0.006]); %will need change for more than 4 coeffs
box off
title('Four predictor model','FontSize',font_size.tit)
ylabel('Values','FontSize',font_size.lab)
xlabel('Coefficient','FontSize',font_size.lab)
legend(p,leg.size,'LineWidth',0.25)
text(0.05,0.95,ltrs(2),ltr_args{:})

[~,h]=suplabel('Autoregressive Model','t',[.08 .12 .87 .84]);
set(h,'FontSize',font_size.tit)

% Step sizes 10-30
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Change Dynamics (10-30)',fig_num),'units','centimeters','outerposition',[0 0 8.8 10],'Renderer','Painters')
ha=tight_subplot(2, 2, [0.025 0.01], 0.09, [.13 .01]);
p=nan(1,2);lims =500+ [-35,35];
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
            p(pp)=errorbar(rng_t,nanmean(d_val_group),nanstd(d_val_group)/sqrt(count_part),symbols.s{pp},'MarkerFaceColor',symbols.c{pp},...
                'LineWidth',0.6,'LineStyle','-','Color',symbols.c{pp},'MarkerSize',2.5,'CapSize',4);hold all;

        end
        text(-0.9*d+1.85,0.92,ltrs(d+2*(3-b)),ltr_args{:})
        xlim([-0.5+rng_t(1),rng_t(end)+0.5]);ylim(lims)
        if d==1
            ylabel(sprintf('%s %d ms','Step size: ',-diff(tempo)),'FontSize',font_size.tit);
            set(gca, 'XTick', rng_t,'XTickLabels','','FontSize',font_size.tick,'TickLength',[0.006 0.006]);
        else
            set(gca,'XTick', rng_t,'XTickLabels','','YTickLabels', '', 'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
        end
        if b==2
            set(gca, 'XTickLabels', rng_t,'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
            xlabel('Beat','FontSize',font_size.lab)
        elseif b==3
            title(direcs{d},'FontSize',font_size.tit)
        end
    end
end
[ll,icons] = legend(ha(4),p,leg.clean,'Location','northeast');
for obj = 1:3
    icons(obj).Position = icons(obj).Position + [0.1 0 0];
end
for obj = 4:6
    icons(obj).Children.Children(1).XData = icons(obj).Children.Children(1).XData + 0.14;
    icons(obj).Children.Children(2).XData = icons(obj).Children.Children(2).XData + [0.2 0.08];
end
ll.Box = 'off';
[~,h]=suplabel('Change dynamics' ,'t', [0.14 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit)
[~,h]=suplabel('Delay interval (baselined, ms)' ,'y', [0.12 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit)

% Parameter recovery
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Parameter recovery',fig_num),'units','centimeters','outerposition',[0 0 18 7.5],'Renderer','Painters')
title_panels = {'Experiment 1:','Experiment 2:','Experiment 2:','Combined update rate ';...
    'phase correction (\alpha)','phase correction (\alpha) z-score','period correction (\beta) z-score','(\alpha and \beta z-score)'};
ha=tight_subplot(1, 4, 0.01, [0.12 0.14], [.05 .03]);
limits_all = [0 1; -3 3; -3 3; -3 3];
fitted = [wing_isoch_av(:,1), zwing_changes(:,1), zwing_changes(:,2), model_params_combined];
recovered = [mean_recovery_isoch(:,1), mean_zrecovery_changes(:,1), mean_zrecovery_changes(:,2), mean_model_params_combined_recovery];
for i=1:4
    axes(ha(i))
    limits = limits_all(i,:);
    scatter(fitted(:,i),recovered(:,i),17,'filled','MarkerFaceAlpha',0.5);  hold on
    plot(limits,limits,'k--');hold on
    set(gca, 'FontSize',font_size.tick,'YTickLabel','','TickLength',[0.006 0.006]);
    text(0.05,0.95,ltrs(i),ltr_args{:})
    r=corr(fitted(:,i),recovered(:,i),'Type','Spearman','Rows','Complete');
    txt = ['\rho ALL',sprintf(' = %0.2g',r)];
    text(0.65, 0.06 ,txt ,'Units','normalized', 'FontSize', font_size.txt)
    title(sprintf('%s\n%s',title_panels{1,i},title_panels{2,i}),'FontSize',font_size.tit)
    xlim(limits);ylim(limits);
    box off;
end
set(ha(1),'YTickLabelMode','auto','TickLength',[0.006 0.006]);
ylabel(ha(1),'Recovered value','FontSize',font_size.lab)
[~,h]=suplabel('Simulation input value' ,'x', [0.05 0.17 0.92 0.9]);set(h,'FontSize',font_size.lab)

% Simulated consecutive correlations
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Simulated consecutive correlations',fig_num),'units','centimeters','outerposition',[0 0 18 8],'Renderer','Painters')
ha=tight_subplot(1, 3, 0.03, [0.12 0.12], [.05 .03]);
title_panels = {'Relation of experiment and simulated data','Simulated data','Experimental data'};

experiment_data = correls_all;
simulated_data = nanmean(correls_all_sim,2);lim = [0 1];
exc_subjects = any(isnan(wing_isoch_av),2);
plts = [];
for i = 1:3
    axes(ha(i))
    if i==3
        Single_Participant_Graph([experiment_data(~exc_subjects) grp(~exc_subjects)],leg.clean,symbols,'KW');
        title('Experimental data','FontSize',font_size.tit);
    elseif i==2
        Single_Participant_Graph([simulated_data grp],leg.clean,symbols,'KW');
        title('Simulated data','FontSize',font_size.tit);
    else
        for pp = 1:3
            scatter(simulated_data(grp==pp),experiment_data(grp==pp),17,symbols.c{pp},'filled','Marker',symbols.s{pp},'MarkerFaceAlpha',0.45);hold on
        end
        plot(lim,lim,'k--')
        r=corr(simulated_data,experiment_data,'Type','Spearman','Rows','Complete');
        txt = ['\rho ',sprintf(' = %0.3g',r)];
        text(0.74, 0.3 ,txt ,'Units','normalized', 'FontSize', font_size.txt)
        ylabel('Experimental data','FontSize',font_size.lab)
        xlabel('Simulated data','FontSize',font_size.lab)
        xlim(lim);
        legend(plts, leg.clean, 'location','southeast','LineWidth',0.25)
    end
    set(gca, 'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
    title(title_panels{i},'FontSize',font_size.tit)
    ylim(lim);
    text(0.05,0.95,ltrs(i),ltr_args{:})
end
[~,h]=suplabel('Correlation of consecutive asynchronies' ,'t', [0.05 0.12 0.92 0.85]);set(h,'FontSize',font_size.tit)

% Simulated change dynamics
fig_num = fig_num+1;
figure('Name',sprintf('[%d][Supp.] Simulated Change Dynamics',fig_num),'units','centimeters','outerposition',[0 0 8.8 13],'Renderer','Painters')
ha=tight_subplot(3, 2, [0.025 0.01], 0.075, [.13 .01]);
p=nan(1,3);lims =500+[-70,70];
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
            p(pp)=errorbar(rng_t,nanmean(d_val_group),nanstd(d_val_group)/sqrt(count_part),symbols.s{pp},'MarkerFaceColor',symbols.c{pp},...
               'LineWidth',0.6,'LineStyle','-','Color',symbols.c{pp},'MarkerSize',2.5,'CapSize',4);hold all;
        end
        text(-0.9*d+1.85,0.92,ltrs(d+2*(6-b)),ltr_args{:})
        xlim([-0.5+rng_t(1),rng_t(end)+0.5]);ylim(lims)
        if d==1
            ylabel(sprintf('%s %d ms','Step size: ',-diff(tempo)),'FontSize',font_size.tit);
            set(gca, 'XTick', rng_t,'XTickLabels','','FontSize',font_size.tick,'TickLength',[0.006 0.006]);
        else
            set(gca,'XTick', rng_t,'XTickLabels','','YTickLabels', '', 'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
        end
        if b==4
            set(gca, 'XTickLabels', rng_t,'FontSize',font_size.tick,'TickLength',[0.006 0.006]);
            xlabel('Beat','FontSize',font_size.lab)
        elseif b==6
            title(direcs{d},'FontSize',font_size.tit)
        end
    end
end
[ll,icons] = legend(ha(6),p,leg.clean,'Location','southeast');
for obj = 1:3
    icons(obj).Position = icons(obj).Position + [0.1 0 0];
end
for obj = 4:6
    icons(obj).Children.Children(1).XData = icons(obj).Children.Children(1).XData + 0.14;
    icons(obj).Children.Children(2).XData = icons(obj).Children.Children(2).XData + [0.2 0.08];
end
ll.Box = 'off';
[~,h]=suplabel('Simulated change dynamics' ,'t', [0.14 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit)
[~,h]=suplabel('Delay interval (baselined, ms)' ,'y', [0.12 0.075 0.85 0.9]);set(h,'FontSize',font_size.tit)

%% save all

FolderName = 'C:\Users\user\Google Drive\Tapping_Project\Figures\';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName = get(FigHandle, 'Name');
    saveas(FigHandle, [FolderName, FigName, '.jpg']);
    saveas(FigHandle, [FolderName, FigName, '.pdf'],'pdf');
    %    exportgraphics(FigHandle, [FolderName, FigName, '.pdf']);
end