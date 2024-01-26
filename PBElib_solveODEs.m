function sol = PBElib_solveODEs(mySettings, y0, T,useJac)
%SOLVEODES solves odes
%   FIXME documentation

%   FIXME arguments

F = @(t,y) PBElib_RHS(t,y,mySettings);

if useJac
    odeOpts = odeset('Stats','off',...
        'RelTol',1e-8,...
        'AbsTol',1e-14,...
        'Jacobian',@(t,y) PBElib_Jacobian(t,y,mySettings),...
        'NonNegative',1);
else
    odeOpts = odeset('Stats','off',...
        'RelTol',1e-8,...
        'AbsTol',1e-14,...
        'NonNegative',1);
end

sol = ode15s(F,[0, T],y0, odeOpts);

end