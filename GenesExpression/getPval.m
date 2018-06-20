function [vHist, timeAxis] = getPval(vitalStat,month,clust)

maxMonth = max(month(~vitalStat));
Oj = zeros(maxMonth,1);

for ii = 1:length(vitalStat)
    if (vitalStat(ii) == 0)
        Oj(month(ii)) = Oj(month(ii))+1;
    end
end
Nj = length(vitalStat) - cumsum(Oj);

Oj1 = zeros(maxMonth,1);
vitalStat1 = vitalStat(clust == 1);
month1 = month(clust == 1);
for ii = 1:length(vitalStat1)
    if (vitalStat1(ii) == 0)
        Oj1(month1(ii)) = Oj1(month1(ii))+1;
    end
end
Nj1 = length(vitalStat1) - cumsum(Oj1);
l

vHist = vHist/max(vHist);
timeAxis = 1:maxMonth;