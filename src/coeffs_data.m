
load for006_empty.mat

iAlpha0 = find(State.Alphas == 0);
iBeta0 = find(State.Betas == 0);

M = State.Machs';
altitudes = State.Altitudes';

CD = [];
for i = 1:length(altitudes)
    CD(:, i) = Coeffs.CD(iAlpha0, :, iBeta0, i)';
end

CA = [];
for i = 1:length(altitudes)
    CA(:, i) = Coeffs.CA(iAlpha0, :, iBeta0, i)';
end