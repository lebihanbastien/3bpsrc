function [ default ] = parameters_default_init(cst)
% Initializes the parameters at default value

%-------------------------------------------------------------------------%
%Integration
%-------------------------------------------------------------------------%
default.ode45.RelTol = 1e-12;
default.ode45.AbsTol = 1e-12;
default.ode87.RelTol = 1e-12;
default.ode87.AbsTol = 1e-12;

%-------------------------------------------------------------------------%
%Differential correction
%-------------------------------------------------------------------------%
default.diff_corr.precision = 1e-10;
default.diff_corr.type = cst.corr.Z0_FIXED;   %Rq: Z0_FIXED should be used with small Az, X0_FIXED otherwise
default.diff_corr.isON = cst.FALSE;

%-------------------------------------------------------------------------%
%Precision in libration points computation
%-------------------------------------------------------------------------%
default.libp.precision = 1e-12;

%-------------------------------------------------------------------------%
%Plotting or not?
%-------------------------------------------------------------------------%
default.plot.diff_corr       = cst.FALSE; %during differential correction
default.plot.XY              = cst.TRUE;  %Plot the XY view
default.plot.YZ              = cst.FALSE; %Plot the YZ view
default.plot.XZ              = cst.FALSE; %Plot the XZ view
default.plot.TD              = cst.TRUE;  %Plot the 3D view
default.plot.halo_orbit      = cst.TRUE;  %during halo orbit computation
default.plot.lyap_orbit      = cst.TRUE;  %during lyap orbit computation
default.plot.manifold_branch = cst.TRUE;  %during manifold computation

end

