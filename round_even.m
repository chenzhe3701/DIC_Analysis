function s = round_even(s)
% chenzhe, 2018-09-10
% round to nearest even number

% convert to positive
ind_neg = s<0;
s(ind_neg) = s(ind_neg) * -1;

% find parts
part_dec = mod(s,2);
part_int = s - part_dec;

% if the mod part >=1, add 2 to the integer part
ind = part_dec >= 1;
s = part_int + 2*ind;

% convert back
s(ind_neg) = s(ind_neg) * -1;

