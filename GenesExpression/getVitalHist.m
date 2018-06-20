function [vHist, timeAxis] = getVitalHist(vitalStat,month)

maxMonth = max(month(~vitalStat));
temp = zeros(maxMonth,1);

for ii = 1:length(vitalStat)
    if (vitalStat(ii) == 0)
        temp(month(ii)) = temp(month(ii))+1;
    end
end

vHist = length(vitalStat) - cumsum(temp);
vHist = vHist/max(vHist);
timeAxis = 1:maxMonth;