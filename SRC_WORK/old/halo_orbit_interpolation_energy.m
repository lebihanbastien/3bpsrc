%-------------------------------------------------------------------------%
% Polynomial interpolation in abacus to recover a good guess of the
% the initial conditions orbit.y0, given either an energy level orbit.C or
% an estimate of the vertical extension orbit.Az.
% The temporary structure fit contains all the elements 
% for this interpolation.
%
% WARNING 1: the interpolation alone should give a good enough guess to
% directly produce a periodic orbit. Hence the differential correction
% process performed in the routine orbit_refinement should be very fast, or
% even useless.
%
% WARNING 2: the interpolation possibilies are limited: see energy bounds
% in halo_init.
%
% RM: the data files contain the NORTHERN family.
%
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function orbit = halo_orbit_interpolation_energy(cr3bp, orbit, halo_init, params, cst)
%-------------------------------------------------------------------------%
% First guess from abacus
%-------------------------------------------------------------------------%
if(isfield(orbit, 'C')) %if an energy level is provided
    
    if(orbit.C > halo_init.Cjaclimit(1) && orbit.C < halo_init.Cjaclimit(2))
        fit.half = 2;
        fit.degree = 2*fit.half;
        [~, array_position] = min(abs(halo_init.matrix(:,8) - orbit.C));
        fit.x =  halo_init.matrix(array_position - fit.half: array_position + fit.half,8);
        for count =1:6
            fit.y(:,count) =  halo_init.matrix(array_position - fit.half: array_position + fit.half,count);
            %Fitting for every dimension of the state (6)
            [fit.p, ~, fit.mu] = polyfit(fit.x,fit.y(:,count),fit.degree);
            %Evaluation
            fit.f(count) = polyval(fit.p,orbit.C,[],fit.mu);
        end
    else
        disp('WARNING: the desired jacobi cst is out of bounds in halo_init.matrix');
    end
    
elseif(isfield(orbit, 'Az')) %if a vertical extension is provided
    
    if(orbit.Az < halo_init.Azlimit)
        fit.half = 2;
        fit.degree = 2*fit.half;
        [~, array_position] = min(abs(halo_init.matrix(:,7) - orbit.Az));
        fit.x =  halo_init.matrix(array_position - fit.half: array_position + fit.half,7);
        for count =1:6
            fit.y(:,count) =  halo_init.matrix(array_position - fit.half: array_position + fit.half,count);
            %Fitting for every dimension of the state (6)
            [fit.p, ~, fit.mu] = polyfit(fit.x,fit.y(:,count),fit.degree);
            %Evaluation
            fit.f(count) = polyval(fit.p,orbit.Az,[],fit.mu);
        end
    else
        disp('WARNING: the vertical extension Az is out of range in halo_init.matrix');
    end
    
    
end

% First guess
yv0 = fit.f(1:6);   

% Specific case of the SOUTHERN family
if(strcmp(orbit.family,cst.orbit.family.SOUTHERN))
    yv0(3) = -yv0(3);    
end

%-------------------------------------------------------------------------%
% Refinement
%-------------------------------------------------------------------------%     
orbit = orbit_refinement(cr3bp, orbit, params, yv0, cst);

end