function [T,Y, Ta,Ya] = MAIN(settings)
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
    
    % Author: Francesco Colombi
    % Skyward Experimental Rocketry | CRD Dept | crd@skywarder.eu
    % email: francesco.colombi@skywarder.eu
    % Release date: 16/04/2016
    
    %Printing License
    print_license();
    tic
    
    %Retrieving Parameters
%     run('config.m');
    

    %Checking if stochastic or standard simulation needed
    if settings.ballistic    
        fprintf('Standard Ballistic Simulation Fired...\n\n');
            [T,Y, Ta,Ya]=std_run_ballistic(settings);
    else
        fprintf('Standard Simulation Fired...\n\n');
            [T,Y, Ta,Ya]=std_run(settings);
    end
    

    toc
end
