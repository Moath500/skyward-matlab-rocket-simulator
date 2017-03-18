function [T,Y] = MAIN()
    %MAIN SCRIPT - This function retrieves all the parameters and run the
    %simulation
    %
    % All the parameters are configured in config.m

    % Author: Ruben Di Battista
    % Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
    % email: ruben.dibattista@skywarder.eu
    % Website: http://www.skywarder.eu
    % April 2014; Last revision: 31.XII.2014
    % License: 2-clause BSD

    %Printing License
    print_license();
    tic
    
    %Retrieving Parameters
    run('config.m');
    

    %Checking if stochastic or standard simulation needed
    if settings.stoch.N > 1
        fprintf('Stochastic Simulation Fired...\n\n');
        if settings.stoch.parallel == 1
            [LP,Z]=stoch_run_p(settings);
        else
            [LP,Z]=stoch_run(settings);
        end
    else
        fprintf('Standard Simulation Fired...\n\n');
        [T,Y, Ta,Ya]=std_run(settings);
    end

    toc
end
