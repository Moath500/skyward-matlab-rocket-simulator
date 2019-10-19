function [settings] = ManageData (data, settings)

% Coefficients in full configuration
settings.CoeffsF = data.full.Coeffs;

% Coefficients in empty configuration
settings.CoeffsE = data.empty.Coeffs;

s = data.full.State;
settings.Alphas = s.Alphas';
settings.Betas = s.Betas';
settings.Altitudes = s.Altitudes';
settings.Machs = s.Machs';

end