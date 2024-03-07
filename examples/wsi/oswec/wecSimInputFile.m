%% Simulation Data
simu = simulationClass();               % Initialize Simulation Class
simu.modelName = 'OSWEC testing';
simu.simMechanicsFile = 'OSWEC.slx';    % Specify Simulink Model File
simu.mode = 'normal';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer = 'on';                   % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                     % Simulation Start Time [s]
simu.rampTime = 10;                    % Wave Ramp Time [s]
simu.endTime = 100;                     % Simulation End Time [s]        
simu.solver = 'ode4';                   % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.03;                         % Simulation Time-Step [s]
simu.cicEndTime = 30;                   % Specify CI Time [s]

% % Regular Waves 
waves = waveClass('regular');           % Initialize Wave Class and Specify Type                                 
waves.height = 0.1;                     % Wave Height [m]
waves.period = 10;                      % Wave Period [s]

%% Body Data
% Flap
body(1) = bodyClass('hydroData/oswec.h5');      % Initialize bodyClass for Flap
body(1).geometryFile = 'geometry/flap.obj';     % Geometry File
body(1).mass = 127000;                          % User-Defined mass [kg]
body(1).inertia = [1.85e6 1.85e6 1.85e6];       % Moment of Inertia [kg-m^2]
body(1).position = [0.0, 0.0, -3.9];

% Base
body(2) = bodyClass('hydroData/oswec.h5');      % Initialize bodyClass for Base
body(2).geometryFile = 'geometry/base.obj';     % Geometry File
body(2).mass = 999;                             % Placeholder mass for a fixed body
body(2).inertia = [999 999 999];                % Placeholder inertia for a fixed body
body(2).position = [0.0, 0.0, -10.9];
body(2).fixed = 'True';

% Rotational PTO
pto(1) = ptoClass('RotationalPTO');                % Initialize ptoClass for PTO1
pto(1).stiffness = 0;                              % PTO Stiffness Coeff [Nm/rad]
pto(1).damping = 0;                                % PTO Damping Coeff [Nsm/rad]
pto(1).bodies = [body(1), body(2)];
pto(1).location = [0, 0, -8.9];