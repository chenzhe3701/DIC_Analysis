function s = round_odd(s)
% chenzhe, 2018-09-10
% round to nearest odd number

% convert to positive
ind_neg = s<0;
s(ind_neg) = s(ind_neg) * -1;

% find parts
part_dec = mod(s,2);
part_int = s - part_dec;

% add 1
s = part_int + 1;

% convert back
s(ind_neg) = s(ind_neg) * -1;
