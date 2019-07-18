% chenzhe, 2018-09-17
% write a function to calculate effective strain
% main purpose is to record the method/equation
%
% the input is a n-by-3 matrix

function eEff = effective_strain_nBy3(M)
    eEff =  sqrt(2/3*(M(:,1).^2 + M(:,3).^2 + 2*M(:,2).^2));
end