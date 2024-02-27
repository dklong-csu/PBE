function Cost = GOLDSIM_objectiveFunction(parameters, data, model_reduction)

%--------------------------------------------------------------------------
%   Some parameters are powers of 10 so transform those
%       Parameters allowed to vary over orders of magnitude can be
%       difficult to optimize, so it can sometimes be easier to estimate
%       the parameter to be: p=10^q
%       and estimate q instead of p
%--------------------------------------------------------------------------
%   parameter K (index=1)
%   parameter c (index=4)
% parameters = pre_parameters;
% parameters(1) = 10^pre_parameters(1);
% parameters(4) = 10^pre_parameters(4);

%--------------------------------------------------------------------------
%   First perform the simulation
%--------------------------------------------------------------------------

ic = [0.0001, 0.0003];
[sol, mySettings] = GOLDSIM_simulateGoldParticles(parameters,model_reduction,ic, data.times(end));

%--------------------------------------------------------------------------
%   L2 norm between simulation and data for each time
%   Add up the difference at every time for the cost function
%--------------------------------------------------------------------------
Cost = 0;
try
    [sim_diams, sim_PSDs] = PBElib_getPSDs(sol,mySettings,data.times);
    for ttt=1:length(data.times)
        %   simulations at that time
        d = sim_diams(:,ttt);
        psd = sim_PSDs(:,ttt);
        psd = psd .* (d.^3);    % volume weight
    
        %   Interpolate the simulation so we can compare to the data
        %   Data is volume weighted so scale by diam^3
        sim = griddedInterpolant(d,psd,'makima','linear');
        fcn_siminterp = @(diam) max(0, sim(diam));
    
        %   Normalize area under curve to be equal to 1
        sim_area = integral(fcn_siminterp,d(1),d(end));
        fcn_qsim = @(diam) fcn_siminterp(diam) / sim_area;
    
        %   Compare to the data
        data_d = data.diam;
        data_q = data.PSDs(:,ttt);
        sim_q = fcn_qsim(data_d);
    
        data_cdf = cumtrapz(data_d, data_q);
        sim_cdf = cumtrapz(data_d, sim_q);
        Cost = Cost + norm(data_q-sim_q)^2;
    end
catch
    for ttt=1:length(data.times)
        data_d = data.diam;
        data_q = data.PSDs(:,ttt);
    
        data_cdf = cumtrapz(data_d, data_q);
        Cost = Cost + norm(data_q)^2;
    end
end
end