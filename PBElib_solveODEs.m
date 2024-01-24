function sol = PBElib_solveODEs(mySettings, y0, T, useMex,useJac)
%SOLVEODES solves odes
%   FIXME documentation

%   FIXME arguments


if useMex
    F = @(t,y) PBElib_RHS_mex(t,y,mySettings);
else
    F = @(t,y) PBElib_RHS(t,y,mySettings);
end

if useJac
    if useMex
        odeOpts = odeset('Stats','off',...
            'RelTol',1e-8,...
            'AbsTol',1e-14,...
            'Jacobian',@(t,y) PBElib_Jacobian_mex(t,y,mySettings),...
            'NonNegative',1);
    else
        odeOpts = odeset('Stats','off',...
        'RelTol',1e-8,...
        'AbsTol',1e-14,...
        'Jacobian',@(t,y) PBElib_Jacobian(t,y,mySettings),...
        'NonNegative',1);
    end
else
    odeOpts = odeset('Stats','off',...
        'RelTol',1e-8,...
        'AbsTol',1e-14,...
        'NonNegative',1);
end

sol = ode15s(F,[0, T],y0, odeOpts);

end