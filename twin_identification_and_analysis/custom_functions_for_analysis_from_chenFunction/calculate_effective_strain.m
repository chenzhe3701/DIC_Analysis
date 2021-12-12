% chenzhe, 2018-03-25
% write a function to calculate effective strain
% main purpose is to record the method/equation

function eEff = calculate_effective_strain(exx,exy,eyy)
    eEff =  sqrt(2/3*(exx.^2 + eyy.^2 + 2*exy.^2));
end