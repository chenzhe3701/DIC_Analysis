fun = @(M) nanmean(M,1);
e_temp = colfilt(exx_t,[5,5],'sliding',fun);
myplot(e_temp);