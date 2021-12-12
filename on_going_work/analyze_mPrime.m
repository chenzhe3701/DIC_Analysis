function [mPrime_each_twin, rank_mPrime_each_twin, ssn_nb_for_each_ss, SF_nb_for_each_ss] = analyze_mPrime(ssn_allowed, ssn_nb_allowed, euler, euler_nb, stressTensor, gbNormal, sampleMaterial, str_for_twin)
% chenzhe, 2019-10-10. good to code, but have not used.  maybe save it for
% Ti analysis as well.
% 
% Need input:
% ssn allowed/considered in the grain considered
% ssn allowed/considered neighbor grain
% Example:
% To analyze the effect of neighbor on twin activation, it's good to set
% in_ssn_allowed = [19:24].
%
%

[schmidFactorG1, schmidFactorG2, mPrimeMatrix, resBurgersMatrix, mPrimeMatrixAbs, resBurgersMatrixAbs] ...
    = calculate_mPrime_and_resB(euler, euler_nb, stressTensor, gbNormal, sampleMaterial, str_for_twin);

SFs_nb = schmidFactorG2(ssn_nb_allowed);

mPrime_local = mPrimeMatrixAbs(ssn_allowed, ssn_nb_allowed);

[mPrime_each_twin, ind] = max(mPrime_local, [], 2);     % output: m'

[sorted_mPrime_each_twin,~] = sort(mPrime_each_twin,'descend');
[~,rank_mPrime_each_twin] = ismember(mPrime_each_twin, sorted_mPrime_each_twin);    % output: m'-rank

ssn_nb_for_each_ss = ssn_nb_allowed(ind);   % output: the ssn in neighbor, for each ss/ts in the grain of interest
SF_nb_for_each_ss = SFs_nb(ind);    % output: the SF of that selected ssn in neighbor

end