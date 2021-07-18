function fig = Single_Participant_Graph(data, labels, symbols, stats_choice)
% data: n_subjects * 2 - column 1 is the data to plot, column 2 is group
%   affiliation (1, 2, 3 etc.)
% labels: cell array of strings corresponding to the groups
% symbols: struct with two fields - c - cell array with colors for each
% group and s - cell array with symbols for each group
% stats_choice: optional argument: 'AN' - ANOVA, 'KW' - Kruskal-Wallis test,
%   'NO' - no stats (default)
% Center statistic default is mean, switched to median if stats are not
% parametric ('KW'). Dispersion is SEM, switched to 25,75% percentiles if stats are
% not parametric ('KW')

if ~exist('stats_choice','var')
    stats_choice = 'NO';
end

num_groups = size(labels,2);

center_stat = nan(1,num_groups);
for i=1:num_groups
    grp_values = find(data(:,2)==i);
    rand_add = (rand(length(grp_values),1)-0.5)/3;
    if ~strcmp(stats_choice,'KW')
        center_stat(i) = nanmean(data(grp_values,1));
        spread_stat([1,2]) = nanstd(data(grp_values,1))/sqrt(length(grp_values(~isnan(grp_values))));
    else
        center_stat(i) = nanmedian(data(grp_values,1));
        spread_stat(1) = center_stat(i)-prctile(data(grp_values,1),25);
        spread_stat(2) = prctile(data(grp_values,1),75)-center_stat(i);
    end
    scatter(data(grp_values,2)+rand_add,data(grp_values,1),17,symbols.c{i},'filled','Marker',symbols.s{i},'MarkerFaceAlpha',0.45);hold on
    errorbar(i,center_stat(i),spread_stat(1),spread_stat(2),'-.','MarkerSize',4,'Color','k', 'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',0.65,'CapSize',6);hold on
end

set(gca, 'XTick', 1:num_groups ,'XTickLabel',labels,'FontSize',5,'TickLength',[0.006 0.006]);
xlim([0,num_groups+1])

loc = [1.5,3.5,2.5]; %for 3 comparisons only!
significance_long = {'*** %0.3f','** %0.3f','* %0.3f','ns %0.3f'};
significance_short = {'*** %0.2f','** %0.2f','* %0.2f','ns %0.2f'};
thresh = [0.001,0.01,0.05];

if ~strcmp(stats_choice,'NO')
    if strcmp(stats_choice,'AN')
        [~,tbl,stats] = anovan(data(:,1),{data(:,2)},'varnames',{'group'},'Display','off');
        Fstring = sprintf('F(%d,%d) = %0.2f',tbl{2,3},tbl{3,3},tbl{2,6});
        if tbl{2,7}<0.001
            pstring = sprintf('p = %0.4f',tbl{2,7});
        elseif tbl{2,7}<0.1
            pstring = sprintf('p = %0.3f',tbl{2,7});
        else
            pstring = sprintf('p = %0.2f',tbl{2,7});
        end
    elseif strcmp(stats_choice,'KW')
        [~,tbl,stats] = kruskalwallis(data(:,1),data(:,2),'off');
        Fstring = sprintf('H(%d) = %0.2f',tbl{2,3},tbl{2,5});
        if tbl{2,6}<0.001
            pstring = sprintf('p = %0.4f',tbl{2,6});
        elseif tbl{2,6}<0.1
            pstring = sprintf('p = %0.3f',tbl{2,6});
        else
            pstring = sprintf('p = %0.2f',tbl{2,6});
        end
    end

    text(0.032,0.02,[Fstring newline pstring],'Units','normalized','FontSize',5,'VerticalAlignment','bottom');
    mltcmp = multcompare(stats,'Display','off');

    %specific comparisons
    for i=1:size(mltcmp,1)
        a = center_stat(mltcmp(i,1));
        b = center_stat(mltcmp(i,2));
        plot([loc(i) loc(i)],[a b],'Color',[1 1 1]/2,'LineWidth',0.35);
    end
    for i=1:num_groups
        j = i+0.75;
        if i==1
            j = 3.75;
        end
        plot([i-0.75,j],[center_stat(i),center_stat(i)],'LineWidth',1.25,'Color',symbols.c{i});hold on
    end
    for i=1:size(mltcmp,1)
        a = center_stat(mltcmp(i,1));
        b = center_stat(mltcmp(i,2));
        p = mltcmp(i,6);
        j=1;
        while j<4
            if p<thresh(j)
                break
            else
                j=j+1;
            end
        end
        y_loc = (a+b)/2;
       % if p>0.4
            y_loc = max(a,b);%+abs((max(data(:,1))-min(data(:,1)))/10);
      %  end
        if p<0.001
            str = sprintf('*** %0.4f',p);
        elseif p<0.1
            str = sprintf(significance_long{j},p);
        elseif p>0.99
            str = 'ns 0.99';
        else
            str = sprintf(significance_short{j},p);
        end
        tmp = strsplit(str);
        str1 = tmp{1}; str2 = tmp{2};
        if j == 4
            text(loc(i), y_loc, [str1, newline, str2],'FontSize',5,'HorizontalAlignment','center','VerticalAlignment','bottom');
        else
            text(loc(i), y_loc, str2,'FontSize',5,'HorizontalAlignment','center','VerticalAlignment','bottom');
            text(loc(i), y_loc, str1,'FontSize',9,'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
    end
else
    for i=1:num_groups
        j = i+0.75;
        if i==1
            j = 3.75;
        end
        plot([i-0.75,j],[center_stat(i),center_stat(i)],'LineWidth',1.25,'Color',symbols.c{i});hold on
    end    
end

box off
end