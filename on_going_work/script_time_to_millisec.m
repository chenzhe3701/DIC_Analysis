% 2019-03-23 to 2019-03-25
% 11:00, 11:30, 12:00, 22:00, 22:30, 23:00

days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31];
tYear = 2019;
tMonth = 3;
tDay = 24;
tHour = 12;
tMin = 51;
tSec = 0;

nHourZone = +8;

nDayExtra_1 = 0;
for year = 1970:tYear-1
   if (mod(year,4)==0)&&(mod(year,100)~=0) || (mod(year,400)==0)
       nDayExtra_1 = nDayExtra_1+1;
   end
end

nDayExtra_2 = 0;
if ((mod(tYear,4)==0)&&(mod(tYear,100)~=0) || (mod(tYear,400)==0)) && (tMonth > 2)
    nDayExtra_2 = 1;
end

nDays = (tYear-1970)*365 + sum(days_per_month(1:tMonth-1)) + (tDay-1) + nDayExtra_1 + nDayExtra_2;

nSec = ((nDays * 24 + tHour - nHourZone) * 60 + tMin) * 60 + tSec;

nMilliSec = nSec*1000;

sprintf('%.0f',nMilliSec)