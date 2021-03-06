function [tbl1,tbl2] = CreateStatTables(data,names,grp)
% as it sounds - I give data and it does the magic
% data is subject*variables
% names is names of the variables
% tbl1 - name, grp, mean, median, std, range (min,max), #in group
% (might have been nans)
% tbl2 - name, typestats(KW\AN), stat_value,df,p_all,p1-2,p1-3,p2-3
types_stats = {'ANOVA','Kruskal-Wallis'};
leg = {'CON','DYS','ASD'};

stats = cell(3*length(names),9);
for pp=1:3
    for i=1:length(names)
        ind = pp+3*(i-1);
        dat = data(grp==pp,i);
        stats{ind,1}=names{i};stats{ind,2}=leg{pp};
        stats{ind,3}=round(nanmean(dat),2);stats{ind,4}=round(nanmedian(dat),3);stats{ind,5}=round(nanstd(dat),2);
        stats{ind,6}=round(nanstd(dat)/sqrt(sum(~isnan(dat))),2);
        stats{ind,7}=round([min(dat), max(dat)],2);stats{ind,8}=sum(~isnan(dat));
        stats{ind,9}=round(prctile(dat,75)-prctile(dat,25),3);
    end
end
tbl1=cell2table(stats,'VariableNames',{'type', 'group', 'mean', 'median', 'std', 'sem','range', 'num','interquartile_range'});
stats = cell(2*length(names),8);
for i=1:length(names)
    for tp=1:2
        ind = tp+2*(i-1);
        dat = data(:,i);
        stats{ind,1}=names{i};stats{ind,2}=types_stats{tp};
        if tp==1
            [~,tbl,statistics] = anovan(dat,{grp},'varnames',{'group'},'Display','off');
            stats{ind,3}=round(tbl{2,6},3);stats{ind,4}=[tbl{2,3},tbl{3,3}];stats{ind,5}=tbl{2,7};
        elseif tp==2
            [~,tbl,statistics] = kruskalwallis(dat,grp,'off');
            stats{ind,3}=round(tbl{2,5},3);stats{ind,4}=tbl{2,3};stats{ind,5}=tbl{2,6};
        end
        mltcmp = multcompare(statistics,'Display','off');
        stats{ind,6}=mltcmp(1,6);stats{ind,7}=mltcmp(2,6);stats{ind,8}=mltcmp(3,6);
    end
end
tbl2=cell2table(stats,'VariableNames',{'type','stat_type', 'stat_value', 'df', 'p_all', 'p12','p13','p23'});
end