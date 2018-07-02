function pVal = logrankWrap(x1,x2)

strStat                 = evalc('logrank(x1,x2);');
B                       = regexp(strStat,'\d.\d*','Match');
pVal                    = str2double(B{12});