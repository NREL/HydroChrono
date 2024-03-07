%% Simulation Data
simu = simulationClass();               % Initialize Simulation Class
simu.modelName = 'RM3 regular wave testing';
simu.simMechanicsFile = 'RM3.slx';      % Specify Simulink Model File
simu.mode = 'normal';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer = 'on';                   % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                     % Simulation Start Time [s]
simu.rampTime = 100;                    % Wave Ramp Time [s]
simu.endTime = 400;                     % Simulation End Time [s]
simu.solver = 'ode4';                   % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.06; 							% Simulation time-step [s]

%% Wave Information 
% % noWaveCIC, no waves with radiation CIC  
% waves = waveClass('noWaveCIC');       % Initialize Wave Class and Specify Type  

% % Regular Waves  
waves = waveClass('regular');           % Initialize Wave Class and Specify Type                                 
waves.height = 0.1;                     % Wave Height [m]
waves.period = 11;                      % Wave Period [s]

% % Regular Waves with CIC
% waves = waveClass('regularCIC');          % Initialize Wave Class and Specify Type                                 
% waves.height = 2.5;                       % Wave Height [m]
% waves.period = 8;                         % Wave Period [s]

% % Irregular Waves using PM Spectrum 
%  waves = waveClass('irregular');           % Initialize Wave Class and Specify Type
%  waves.height = 2.5;                       % Significant Wave Height [m]
%  waves.period = 8;                         % Peak Period [s]
%  waves.spectrumType = 'PM';                % Specify Wave Spectrum Type
%  waves.direction=[0];

% % Irregular Waves using JS Spectrum with Equal Energy and Seeded Phase
% waves = waveClass('irregular');           % Initialize Wave Class and Specify Type
% waves.height = 2.5;                       % Significant Wave Height [m]
% waves.period = 8;                         % Peak Period [s]
% waves.spectrumType = 'JS';                % Specify Wave Spectrum Type
% waves.bem.option = 'EqualEnergy';         % Uses 'EqualEnergy' bins (default) 
% waves.phaseSeed = 1;                      % Phase is seeded so eta is the same

% % Irregular Waves using PM Spectrum with Traditional and State Space 
% waves = waveClass('irregular');           % Initialize Wave Class and Specify Type
% waves.height = 2.5;                       % Significant Wave Height [m]
% waves.period = 8;                         % Peak Period [s]
% waves.spectrumType = 'PM';                % Specify Wave Spectrum Type
% simu.stateSpace = 1;                      % Turn on State Space
% waves.bem.option = 'Traditional';         % Uses 1000 frequnecies

% % Irregular Waves with imported spectrum
% waves = waveClass('spectrumImport');      % Create the Wave Variable and Specify Type
% waves.spectrumFile = 'spectrumData.mat';  % Name of User-Defined Spectrum File [:,2] = [f, Sf]

% % Waves with imported wave elevation time-history  
% waves = waveClass('elevationImport');          % Create the Wave Variable and Specify Type
% waves.elevationFile = 'elevationData.mat';     % Name of User-Defined Time-Series File [:,2] = [time, eta]

%% Body Data
% Float
body(1) = bodyClass('hydroData/rm3.h5');      
    % Create the body(1) Variable, Set Location of Hydrodynamic Data File 
    % and Body Number Within this File.   
body(1).geometryFile = 'geometry/float.obj';    % Location of Geomtry File
body(1).mass = 725834;                   
    % Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    % Weight.
body(1).inertia = [20907301 21306090.66 37085481.11];  % Moment of Inertia [kg*m^2]     
body(1).position = [0.0, 0.0, -0.72];

% Spar/Plate
body(2) = bodyClass('hydroData/rm3.h5'); 
body(2).geometryFile = 'geometry/plate.obj'; 
body(2).mass = 886691;                   
body(2).inertia = [94419614.57 94407091.24 28542224.82];
body(2).position = [0, 0, -21.29];

%% PTO and Constraint Parameters
% Floating (3DOF) Joint
constraint(1) = constraintClass('Constraint1'); % Initialize Constraint Class for Constraint1
constraint(1).location = [0 0 0];               % Constraint Location [m]

% PTO 1 - Translational
pto(1) = ptoClass('TranslationalPTO');                      % Initialize PTO Class for PTO1
pto(1).stiffness = 0;                           % PTO Stiffness [N/m]
pto(1).damping = 1200000;                       % PTO Damping [N/(m/s)]
pto(1).location = [0 0 0];                      % PTO Location [m]
pto(1).bodies = [body(1), body(2)];
pto(1).attachments = [[0, 0, -0.72], [0, 0, -21.29]];

% PTO 2 - Linear Spring-Damper
pto(2) = ptoClass('LinSpringDamper');             % Initialize PTO Class for PTO2
pto(2).stiffness = 0;                        % Spring Stiffness [N/m]
pto(2).damping = 50000;                           % Damping Coefficient [N/(m/s)]
pto(2).restLength = 1.0;                          % Spring Rest Length [m]
pto(2).pto_radius = 0.05;                         % Spring Visual Radius [m]
pto(2).pto_resolution = 20;                       % Spring Visual Resolution
pto(2).pto_spring_turns = 10;                     % Number of Visual Coils in the Spring
pto(2).location = [0 0 0];                        % PTO Location [m] (if needed)
pto(2).bodies = [body(1), body(2)];               % Bodies connected by PTO
pto(2).attachments = [[0.0, 0, -0.5], [0.0, 0, -2.5]]; % Attachment Points [m]
