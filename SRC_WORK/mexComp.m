%--------------------------------------------------------------------------
% Script for mex compilation. Use with care. 
%
% With great power comes great responsibility.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%--------------------------------------------------------------------------
% Define some paths
gslpath         = '-I/usr/local/include/gsl';
incpath         = '-I/C';
libgslpath      = 'lib/libgsl.a';
libgslcblaspath = 'lib/libgslcblas.a';
                    

                                         
 %% Build cmic.cpp
mex(gslpath, 'cmic.cpp', libgslpath,... 
                         libgslcblaspath, ...
                         'Cpp/FTA.cpp',...
                         'Cpp/Oftsc.cpp',...
                         'Cpp/Ofsc.cpp',...
                         'Cpp/pmcoc.cpp',...
                         'Cpp/env.cpp');                    
 
 %% Build cumo.cpp
mex(gslpath, 'cumo.cpp', libgslpath,... 
                          libgslcblaspath, ...
                          'Cpp/FTA.cpp',...
                          'Cpp/Oftsc.cpp',...
                          'Cpp/Ofsc.cpp',...
                          'Cpp/pmcoc.cpp',...
                          'Cpp/vf.cpp',...
                          'Cpp/ode.cpp',...
                          'Cpp/nrutil.c',...
                          'Cpp/single_orbit.cpp',...
                          'Cpp/env.cpp');
                      
 %% Build cmo.cpp
mex(gslpath, 'cmo.cpp', libgslpath,... 
                          libgslcblaspath, ...
                          'Cpp/FTA.cpp',...
                          'Cpp/Oftsc.cpp',...
                          'Cpp/Ofsc.cpp',...
                          'Cpp/pmcoc.cpp',...
                          'Cpp/vf.cpp',...
                          'Cpp/ode.cpp',...
                          'Cpp/nrutil.c',...
                          'Cpp/single_orbit.cpp',...
                          'Cpp/env.cpp');
                      
                      
%% Build ode78_cr3bp.c
mex(gslpath, 'ode78_cr3bp.c', libgslpath,... 
                              libgslcblaspath,... 
                              'C/ode78.c',... 
                              'C/cr3bp_derivatives.c',... 
                              'C/custom.c',... 
                              'C/custom_ode.c');
                          
%% Build ode78_cr3bp_event.c
mex(gslpath, 'ode78_cr3bp_event.c', libgslpath,...
                                    libgslcblaspath, ...
                                    'C/ode78.c',... 
                                    'C/custom_odezero.c', ...
                                    'C/cr3bp_derivatives.c', ...
                                    'C/custom.c', ...
                                    'C/custom_ode.c');
                                
%% Build ode78_bcp.c                               
mex(gslpath, 'ode78_bcp.c', libgslpath,... 
                              libgslcblaspath,... 
                              'C/ode78.c',... 
                              'C/cr3bp_derivatives.c',... 
                              'C/custom.c',... 
                              'C/custom_ode.c');
                          
%% Build ode78_bcp_event.c                          
mex(gslpath, 'ode78_bcp_event.c', libgslpath,...
                                    libgslcblaspath, ...
                                    'C/ode78.c',... 
                                    'C/custom_odezero.c', ...
                                    'C/cr3bp_derivatives.c', ...
                                    'C/custom.c', ...
                                    'C/custom_ode.c');
                                
%% Build test.cpp
mex(gslpath, 'test.cpp', libgslpath,... 
                         libgslcblaspath);
                                