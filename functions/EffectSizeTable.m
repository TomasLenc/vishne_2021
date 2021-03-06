function tbl = EffectSizeTable(data,names,grp)
% as it sounds - I give data and it does the magic
% data is subject*variables
% names is names of the variables
% tbl2 - Cliffs delta for 1-2, 1-3, 2-3 (in the future I can add cohen's d)

stats = cell(length(names),4);
for i=1:length(names)
    dat = data(:,i);
    stats{i,1}=names{i};
    for pp=1:2
        for pp2=(pp+1):3
            stats{i, pp+pp2-1}= round(Cliffs_Delta(dat(grp==pp),dat(grp==pp2)),3);
        end
    end
end
tbl=cell2table(stats,'VariableNames',{'type','effect12','effect13','effect23'});
end