function [ ] = print_license()
% Function that prints at first run the LICENSE saved in the LICENSE folder
% in the same directory of this file

% Author: Ruben Di Battista
% Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
% email: ruben.dibattista@skywarder.eu
% Website: http://www.skywarder.eu
% April 2014; Last revision: 25.IV.2014
% License:  2-clause BSD

if exist('LICENSE','dir')
    if not(exist('LICENSE/run.txt','file'))
        fh = fopen('LICENSE/LICENSE','r');
        while not(feof(fh))
            l = fgetl(fh);
            fprintf(strcat(l,'\n'));
        end
        fclose(fh);
        fh = fopen('LICENSE/run.txt','wt');
        fprintf(fh,strcat('License accepted:',datestr(clock,0)));
        fclose(fh);
        
        agree = input('Write ''AGREE'' to accept the terms\n\n','s');
        while not(strcmpi(agree,'agree'))
            agree = input('Write ''AGREE'' to accept the terms\n\n','s');
        end
    end
end


end

