% chenzhe, 2019-07-07
% write a function to calculate effective strain
% main purpose is to record the method/equation
%
% the 3 inputs are exx, exy, eyy

function eEff = effective_strain(exx,exy,eyy)
    eEff =  sqrt(2/3*(exx.^2 + eyy.^2 + 2*exy.^2));
end