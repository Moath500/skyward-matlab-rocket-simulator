function [all_steps] = RecallOdeFcn(fun,T,Y,varargin)
% RECALLODEFCN - This function allows to compute some parameters used
% inside the ODE integrations
%
% OUTPUTS:
%             - all_steps: structure which contains all the parameters needed 
%
% Author: Adriano Filippo Inno
% Skyward Experimental Rocketry | AFD Dept | crd@skywarder.eu
% email: adriano.filippo.inno@skywarder.eu
% Release date: 16/11/2018

NT = length(T);
fun_info = functions(fun);

for i = 1:NT
    
    [~,single_step] = fun(T(i),Y(i,:),varargin{:});
    
     all_steps.coeff.XCP(i) = single_step.coeff.XCP;
end