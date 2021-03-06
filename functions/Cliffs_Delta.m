function delta = Cliffs_Delta(data1,data2)
% calculates cliff's delta (see wiki
% https://en.wikipedia.org/wiki/Effect_size#Effect_size_for_ordinal_data)

n = length(data1);
m = length(data2);

delta = 0;
for i=1:n
    for j=1:m
        if data1(i)>data2(j)
            delta = delta+1;
        elseif data1(i)<data2(j)
            delta = delta-1;
        end
    end
end

delta = delta/(n*m);
end