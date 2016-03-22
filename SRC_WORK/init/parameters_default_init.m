%-------------------------------------------------------------------------%
% parameters_default_init(cst)
%
% Initializes the parameters at default values.
%       
% Author: BLB
% Version: 2.0
% Year: 2016
%-------------------------------------------------------------------------%
function [ default ] = parameters_default_init(cst)
%-------------------------------------------------------------------------%
%Integration
%-------------------------------------------------------------------------%
default.ode45.RelTol = 3e-14; %minimum allowed by ode45
default.ode45.AbsTol = 1e-14; %arbitrary (already very slow!)
default.ode87.RelTol = 3e-14; %minimum allowed by ode45
default.ode87.AbsTol = 1e-14; %arbitrary (already very slow!)

%-------------------------------------------------------------------------%
%Differential correction
%-------------------------------------------------------------------------%
default.diff_corr.precision = 1e-10;
default.diff_corr.type = cst.corr.Z0_FIXED;   %Rq: Z0_FIXED should be used with small Az, X0_FIXED otherwise.
%[Deprecated]
default.diff_corr.isON = cst.FALSE;

%-------------------------------------------------------------------------%
%Precision in libration points computation
%-------------------------------------------------------------------------%
default.libp.precision = 1e-12;

%-------------------------------------------------------------------------%
%Plotting or not?
%-------------------------------------------------------------------------%
default.plot.orbit           = cst.TRUE;  %Global switch for plotting
default.plot.diff_corr       = cst.FALSE; %during differential correction
default.plot.XY              = cst.TRUE;  %Plot the XY view
default.plot.YZ              = cst.FALSE; %Plot the YZ view
default.plot.XZ              = cst.FALSE; %Plot the XZ view
default.plot.TD              = cst.TRUE;  %Plot the 3D view
default.plot.manifold_branch = cst.TRUE;  %during manifold computation
default.plot.LineSmoothing   = 'off';     %during manifold computation

default.plot.firstPrimDisp   = cst.FALSE;  %is the first primary (e.g. the Sun in the Sun-Earth system) displayed?
default.plot.allLibPoints    = cst.FALSE;  %are all libration points displayed?
default.plot.names           = cst.FALSE;  %are the names displayed?
default.plot.tdAxes          = cst.FALSE;  %are the pretty 3D axes displayed?
default.plot.bigPrimFac      = 1.0;         %the primaries appear bigPrimFac x bigger than they actually are (easier to see on screen)

%[Deprecated]
default.plot.halo_orbit      = cst.TRUE;  %during halo orbit computation (deprecated)
default.plot.lyap_orbit      = cst.TRUE;  %during lyap orbit computation (deprecated)



%-------------------------------------------------------------------------%
%Computation type
%-------------------------------------------------------------------------%
default.computation.type = cst.computation.MEX; %use MEX routines by default

end

