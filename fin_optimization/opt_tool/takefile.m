function [CoeffsF,CoeffsE,Alphas,Betas,Altitudes,Machs] = takefile (data)


% Coefficients in full configuration
CoeffsF = data.full.Coeffs;

% Coefficients in empty configuration
CoeffsE = data.empty.Coeffs;

s = data.full.State;
Alphas = s.Alphas';
Betas = s.Betas';
Altitudes = s.Altitudes';
Machs = s.Machs';

end