%% Simulation Data
simu = simulationClass();                          
simu.modelName = 'TALOS';
simu.explorer = 'off';                   % Turn SimMechanics Explorer (on/off)
simu.dt = 0.01;                                
simu.endTime = 100.0;

%% Wave Information - Regular Waves
waves = waveClass('regular');           % Initialize Wave Class and Specify Type                                 
waves.height = 0.1;                     % Wave Height [m]
waves.period = 8;                      % Wave Period [s]

%% Body Data - Hull
body(1) = bodyClass('hydroData/talos.h5');
body(1).geometryFile = 'geometry/talos_hull.obj';            
body(1).mass = 1873881;                              
body(1).inertia = [237600000, 237600000, 237600000]; 
body(1).position = [0.0, 0.0, -7.96];                
body(1).orientation = [0.0, 0.0, 0.0];               

%% Body Data - Internal Reaction Mass (IRM)
body(2) = bodyClass('nonhydro');
body(2).type = 'nonhydro-sphere'
body(2).radius = 4.0;
body(2).mass = 1873881;
body(2).inertia = [106150000, 106150000, 106150000];
body(2).position = [0.0, 0.0, -5.7];

%% PTO Definitions - Lower PTOs
pto(1) = ptoClass('LinSpringDamper');
pto(1).stiffness = 500000;
pto(1).damping = 2000000;
pto(1).rest_length = 5.0;
pto(1).bodies = [body(1), body(2)]; % Adjust the body indices as appropriate
pto(1).attachments = [[5.000, 0.000, -8.660], [2.500, 0.000, -4.330]];

pto(2) = ptoClass('LinSpringDamper');
pto(2).stiffness = 500000;
pto(2).damping = 2000000;
pto(2).rest_length = 5.0;
pto(2).bodies = [body(1), body(2)];
pto(2).attachments = [[-2.500, 4.330, -8.660], [-1.250, 2.165, -4.330]];

pto(3) = ptoClass('LinSpringDamper');
pto(3).stiffness = 500000;
pto(3).damping = 2000000;
pto(3).rest_length = 5.0;
pto(3).bodies = [body(1), body(2)];
pto(3).attachments = [[-2.500, -4.330, -8.660], [-1.250, -2.165, -4.330]];

%% PTO Definitions - Upper PTOs
pto(4) = ptoClass('LinSpringDamper');
pto(4).stiffness = 2000000;
pto(4).damping = 250000;
pto(4).rest_length = 2.0;
pto(4).bodies = [body(1), body(2)];
pto(4).attachments = [[5.000, 0.000, 8.660], [2.500, 0.000, 4.330]];

pto(5) = ptoClass('LinSpringDamper');
pto(5).stiffness = 2000000;
pto(5).damping = 250000;
pto(5).rest_length = 2.0;
pto(5).bodies = [body(1), body(2)];
pto(5).attachments = [[-2.500, 4.330, 8.660], [-1.250, 2.165, 4.330]];

pto(6) = ptoClass('LinSpringDamper');
pto(6).stiffness = 2000000;
pto(6).damping = 250000;
pto(6).rest_length = 2.0;
pto(6).bodies = [body(1), body(2)];
pto(6).attachments = [[-2.500, -4.330, 8.660], [-1.250, -2.165, 4.330]];
