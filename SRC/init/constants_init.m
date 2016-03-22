%-------------------------------------------------------------------------%
% Initialisation of various constants for halo orbit computation
%-------------------------------------------------------------------------%
% @return the structure cst
function cst = constants_init()
%-------------------------------------------------------------------------%
% Misc
%-------------------------------------------------------------------------%
cst.TRUE  = 1;
cst.FALSE = 0;

%-------------------------------------------------------------------------%
% Environment
%-------------------------------------------------------------------------%
cst.env.G = 6.67428e-11; %gravitational constant
cst.env.AU = 1.49597871e8;        %astronomical unit in km
cst.env.julian.y2015 = 2457023.5; %Julian date of 01/01/2015 at 00:00
cst.env.julian.y2000 = 2451544.5; %Julian date of 01/01/2015 at 00:00
cst.env.hours = 3600;             %in seconds
cst.env.days  = 86400;            %in seconds
cst.env.years = 31556926;         %in seconds (true value)

%-------------------------------------------------------------------------%
% Orbit
%-------------------------------------------------------------------------%
% Family
cst.orbit.NORTHERN = 'NORTHERN';
cst.orbit.SOUTHERN = 'SOUTHERN';

% Computation type used to built the orbit ('EMPTY' initially)
cst.orbit.EMPTY        =  'EMPTY';
cst.orbit.APPROXIMATED =  'APPROXIMATED';
cst.orbit.REAL         =  'REAL';

% Initialisation of the STM as the identity matrix
cst.orbit.STM0 = eye(6);

%-------------------------------------------------------------------------%
% Type of differential correction
%-------------------------------------------------------------------------%
cst.corr.X0_FIXED = 1;
cst.corr.Z0_FIXED = 2;

%-------------------------------------------------------------------------%
% Type of plot
%-------------------------------------------------------------------------%
cst.plot.ADIM = 0;
cst.plot.DIM  = 1;

%-------------------------------------------------------------------------%
%Manifold
%-------------------------------------------------------------------------%
% EXTERIOR == towards east, INTERIOR == towards west
cst.manifold.EXTERIOR = 1;
cst.manifold.INTERIOR = -1;

% Type
cst.manifold.STABLE = 1;
cst.manifold.UNSTABLE = -1;

% Type of termination (by event or freely)
cst.manifold.event.type.FREE      = 'FREE';                       %no termination event
cst.manifold.event.type.X_SECTION = 'X_SECTION';                  %termination on a section x = cst
cst.manifold.event.type.Y_SECTION = 'Y_SECTION';                  %termination on a section y = cst
cst.manifold.event.type.Z_SECTION = 'Z_SECTION';                  %termination on a section z = cst
cst.manifold.event.type.ANGLE_SECTION = 'ANGLE_SECTION';          %termination on a section phi = cst, with phi a given angle
cst.manifold.event.type.FLIGHT_PATH_ANGLE = 'FLIGHT_PATH_ANGLE';  %termination on a section fpa = cst

% Is the event terminal?
cst.manifold.event.isterminal.YES = 1;
cst.manifold.event.isterminal.NO = 0;

% Which zeros trigger the event?
cst.manifold.event.direction.INCREASING =  1;
cst.manifold.event.direction.DECREASING = -1;
cst.manifold.event.direction.ALL        =  0;


end


