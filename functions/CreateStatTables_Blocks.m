function [tbl1,tbl2] = CreateStatTables_Blocks(data,names,grp)
% as it sounds - I give data and it does the magic
% data is subject*blocks(5)*other stuff (changes)
% names is names of the other stuff
% tbl1 - name, step, grp, mean, median, std, range (min,max), #in group
% (might have been nans)
% tbl2 - name, step, typestats(KW\AN), stat_value,df,p_all,p1-2,p1-3,p2-3
types_stats = {'ANOVA','Kruskal-Wallis'};
leg = {'CON','DYS','ASD'};
vals = [500,505,515,525,535,545];

stats = cell(5*3*size(data,3),9);
for blk=2:6
    for pp=1:3
        for i=1:size(data,3)
            ind = pp+3*(blk-2)+15*(i-1);
            dat = data(grp==pp,blk-1,i);
            stats{ind,1}=names{i};stats{ind,2}=2*vals(blk)-1000;stats{ind,3}=leg{pp};
            stats{ind,4}=round(nanmean(dat),3);stats{ind,5}=round(nanmedian(dat),3);stats{ind,6}=round(nanstd(dat),3);
            stats{ind,7}=round(nanstd(dat)/sqrt(sum(~isnan(dat))),3);
            stats{ind,8}=round([min(dat), max(dat)],3);stats{ind,9}=sum(~isnan(dat));
        end
    end
end
tbl1=cell2table(stats,'VariableNames',{'type', 'step', 'group', 'mean', 'median', 'std', 'sem','range', 'num'});
stats = cell(5*2*size(data,3),9);
for blk=2:6
    for i=1:size(data,3)
        for tp=1:2
            ind = tp+2*(blk-2)+10*(i-1);
            dat = data(:,blk-1,i);
            stats{ind,1}=names{i};stats{ind,2}=2*vals(blk)-1000;stats{ind,3}=types_stats{tp};
            if tp==1
                [~,tbl,statistics] = anovan(dat,{grp},'varnames',{'group'},'Display','off');
                stats{ind,4}=round(tbl{2,6},3);stats{ind,5}=[tbl{2,3},tbl{3,3}];stats{ind,6}=round(tbl{2,7},5);
            elseif tp==2
                [~,tbl,statistics] = kruskalwallis(dat,grp,'off');
                stats{ind,4}=round(tbl{2,5},3);stats{ind,5}=tbl{2,3};stats{ind,6}=round(tbl{2,6},5);
            end
            mltcmp = multcompare(statistics,'Display','off');
            stats{ind,7}=round(mltcmp(1,6),5);stats{ind,8}=round(mltcmp(2,6),5);stats{ind,9}=round(mltcmp(3,6),5);
        end
    end
end
tbl2=cell2table(stats,'VariableNames',{'type', 'step', 'stat_type', 'stat_value', 'df', 'p_all', 'p12','p13','p23'});
end