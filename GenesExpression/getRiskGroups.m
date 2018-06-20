function [x1,x2] = getRiskGroups(month,vitalStat,groups)


maxMonth = 60;

group1 = (month <= maxMonth | vitalStat) & groups == 1; %groupsPca == 1;
group2 = (month <= maxMonth | vitalStat) & groups == 2; %groupsPca == 2;
month1 = month(group1);
month1(month1 > maxMonth & vitalStat(group1)) = maxMonth;
month2 = month(group2);
month2(month2 > maxMonth & vitalStat(group2)) = maxMonth;

x1 = [ month1 , vitalStat(group1) ];
x2 = [ month2 , vitalStat(group2) ];