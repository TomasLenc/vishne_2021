function fig = Single_Participant_Graph(data, labels, stats_choice)
% data num of participants * 2 - first column is the data, second is group (1, 2, etc.)
% labels -cell array of strings corresponding to groups
% stats is an optional argument - default is no stats (also using 'NO'), 'AN' - anova, 'KW' -
% kruskal wallis
% default for center stat is mean and dispersion is sem, if 'KW' it's
% medians and 25,75% percentiles

if ~exist('stats_choice','var')
    stats_choice = 'NO';
end

SYMBS={'bo','r^','gs'};
num_groups = size(labels,2);

center_stat = nan(1,num_groups);
for i=1:num_groups
    grp_values = find(data(:,2)==i);
    rand_add = (rand(length(grp_values),1)-0.5)/5;
    if ~strcmp(stats_choice,'KW')
        center_stat(i) = nanmean(data(grp_values,1));
        spread_stat([1,2]) = nanstd(data(grp_values,1))/sqrt(length(grp_values(~isnan(grp_values))));
    else
        center_stat(i) = nanmedian(data(grp_values,1));
        spread_stat(1) = center_stat(i)-prctile(data(grp_values,1),25);
        spread_stat(2) = prctile(data(grp_values,1),75)-center_stat(i);
    end
    j = i+0.75;color = SYMBS{i}(1);
    if i==1
        j = 3.75;
    end
    scatter(data(grp_values,2)+rand_add,data(grp_values,1),40,SYMBS{i},'filled','MarkerEdgeColor','k');
    hold on
    errorbar(i,center_stat(i),spread_stat(1),spread_stat(2),'-.','MarkerSize',9,'Color','k', 'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',10);
    hold on
    plot([i-0.75,j],[center_stat(i),center_stat(i)],'LineWidth',2,'Color',color);
    hold on
end

set(gca, 'XTick', 1:num_groups ,'XTickLabel',labels,'FontSize',13);
xlim([0,num_groups+1])

loc = [1.5,3.5,2.5]; %for 3 comparisons only!
significance_verylong = {'*** %0.4f','** %0.3f','* %0.3f','ns %0.3f'};
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

    text(0.05,0.02,[Fstring newline pstring],'Units','normalized','FontSize',11,'VerticalAlignment','bottom');
    mltcmp = multcompare(stats,'Display','off');

    %specific comparisons
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
        if p>0.4
            y_loc = max(a,b)+abs((max(data(:,1))-min(data(:,1)))/25);
        end
        if p<0.001
            str = sprintf(significance_verylong{j},p);
        elseif p<0.1
            str = sprintf(significance_long{j},p);
        elseif p>0.99
            str = 'ns 0.99';
        else
            str = sprintf(significance_short{j},p);
        end
        plot([loc(i) loc(i)],[a b],'Color','k');
        text(loc(i)-0.25, y_loc, str,'FontSize',11);
    end
end

box off
end