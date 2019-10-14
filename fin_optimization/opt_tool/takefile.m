function [CoeffsF,CoeffsE,Alphas,Betas,Altitudes,Machs] = takefile (r_name)

DATA_PATH = 'data/';
filename = strcat(DATA_PATH, r_name);

% Coefficients in full configuration
filename_full = strcat(filename,'_full.mat');
CoeffsF = load(filename_full,'Coeffs');
CoeffsF = CoeffsF.Coeffs;


% Coefficients in empty configuration
filename_empty = strcat(filename,'_empty.mat');
CoeffsE = load(filename_empty,'Coeffs');
CoeffsE = CoeffsE.Coeffs;

s = load(filename_full,'State');
Alphas = s.State.Alphas';
Betas = s.State.Betas';
Altitudes = s.State.Altitudes';
Machs = s.State.Machs';

end