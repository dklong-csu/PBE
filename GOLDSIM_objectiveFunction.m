function Cost = GOLDSIM_objectiveFunction(parameters, data, model_reduction)
% compare_data.diam = data_diam;
% compare_data.PSDs = data_PSDs;
% compare_data.times = data_times;

%--------------------------------------------------------------------------
%   First perform the simulation
%--------------------------------------------------------------------------
ic = [0.0001, 0.0003];
[sol, mySettings] = GOLDSIM_simulateGoldParticles(parameters,model_reduction,ic);

%--------------------------------------------------------------------------
%   L2 norm between simulation and data for each time
%   Add up the difference at every time for the cost function
%--------------------------------------------------------------------------
Cost = 0;
[sim_diams, sim_PSDs] = PBElib_getPSDs(sol,mySettings,data.times);
for ttt=1:length(data.times)
    %   simulations at that time
    d = sim_diams(:,ttt);
    psd = sim_PSDs(:,ttt);

    %   Interpolate the simulation so we can compare to the data
    sim = griddedInterpolant(d,psd,'makima','linear');
    fcn_siminterp = @(diam) (diam>=d(1)) .* (diam<=d(end)) .* sim(diam);

    %   Normalize area under curve to be equal to 1
    sim_area = integral(fcn_siminterp,d(1),d(end));
    fcn_qsim = @(diam) sim(diam) / sim_area;

    %   Compare to the data
    data_d = data.diam;
    data_q = data.PSDs(:,ttt);
    sim_q = fcn_qsim(data_d);
    Cost = Cost + norm(data_q-sim_q)^2;
end
end