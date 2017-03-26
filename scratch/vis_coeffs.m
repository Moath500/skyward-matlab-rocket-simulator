% Author: Francesco Colombi
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: francesco.colombi@skywarder.eu
% Release date: 16/04/2016

close all
clear, clc

load for006_empty.mat

iAlpha0 = find(State.Alphas == 0);
iBeta0 = find(State.Betas == 0);

% for iMach = 1:length(State.Machs)
%     figure();
%     for iAlt = 1:length(State.Altitudes)
%         subplot(3, 1, 1), hold on, grid on;
%         title(['CL vs Alpha for Beta = 0, Mach = ', num2str(State.Machs(iMach))]);
%         plot(State.Alphas, Coeffs.CL(:, iMach, iBeta0, iAlt));
%     end
%
%     for iAlt = 1:length(State.Altitudes)
%         subplot(3, 1, 2), hold on, grid on;
%         title(['CD vs Alpha for Beta = 0, Mach = ', num2str(State.Machs(iMach))]);
%         plot(State.Alphas, Coeffs.CD(:, iMach, iBeta0, iAlt));
%     end
%
%     for iAlt = 1:length(State.Altitudes)
%         subplot(3, 1, 3), hold on, grid on;
%         title(['Cm vs Alpha for Beta = 0, Mach = ', num2str(State.Machs(iMach))]);
%         plot(State.Alphas, Coeffs.CM(:, iMach, iBeta0, iAlt));
%     end
% end

figure(), grid on, hold on;
for iMach = 1:length(State.Machs)
    for iAlt = 1
        title('X_{CP} vs Alpha for Beta = 0, different Mach');
        plot(State.Alphas, Coeffs.X_C_P(:, iMach, iBeta0, iAlt));
    end
end
legend(num2str(State.Machs(1:end)'))

figure();
plot(State.Machs, Coeffs.X_C_P(iAlpha0, :, iBeta0, 1)), grid on;
title('X_{CP} vs Mach, with alpha and beta equal to 0');
xlabel('Mach [-]')
ylabel('X_{CP} from moment center in ref. lenght, neg. AFT of moment center')
