%% Simulation Data
simu = simulationClass();
simu.modelName = 'TALOS_testing';
simu.explorer = 'on';
simu.dt = 0.010005002501250625;
simu.endTime = 150.0;

%% Wave Information - Regular Waves
waves = waveClass('regular');
waves.height = 0.2;
waves.period = 9;

%% Body Data - 
body(1) = bodyClass('../hydroData/talos.h5');
body(1).geometryFile = '../geometry/talos_hull.obj';
body(1).mass = 1873881;
body(1).inertia = [237600000, 237600000, 237600000];
body(1).position = [0.0, 0.0, -7.96];
body(1).orientation = [0.0, 0.0, 0.0];

%% Body Data - nonhydro-sphere
body(2) = bodyClass('');
body(2).mass = 1873881;
body(2).inertia = [106150000, 106150000, 106150000];
body(2).position = [0.0, 0.0, -5.7];
body(2).type = 'nonhydro-sphere';
body(2).radius = 4.0;

%% PTO Definitions
pto(1) = ptoClass('LinSpringDamper');
pto(1).stiffness = 30000.0;
pto(1).damping = 5000.0;
pto(1).rest_length = 5.0;
pto(1).bodies = [body(1), body(2)];
pto(1).attachments = [[5.0, 0.0, -8.66], [2.5, 0.0, -4.33]];

%% PTO Definitions
pto(2) = ptoClass('LinSpringDamper');
pto(2).stiffness = 30000.0;
pto(2).damping = 5000.0;
pto(2).rest_length = 5.0;
pto(2).bodies = [body(1), body(2)];
pto(2).attachments = [[-2.5, 4.33, -8.66], [-1.25, 2.165, -4.33]];

%% PTO Definitions
pto(3) = ptoClass('LinSpringDamper');
pto(3).stiffness = 30000.0;
pto(3).damping = 5000.0;
pto(3).rest_length = 5.0;
pto(3).bodies = [body(1), body(2)];
pto(3).attachments = [[-2.5, -4.33, -8.66], [-1.25, -2.165, -4.33]];

%% PTO Definitions
pto(4) = ptoClass('LinSpringDamper');
pto(4).stiffness = 30000.0;
pto(4).damping = 5000.0;
pto(4).rest_length = 2.0;
pto(4).bodies = [body(1), body(2)];
pto(4).attachments = [[5.0, 0.0, 8.66], [2.5, 0.0, 4.33]];

%% PTO Definitions
pto(5) = ptoClass('LinSpringDamper');
pto(5).stiffness = 30000.0;
pto(5).damping = 5000.0;
pto(5).rest_length = 2.0;
pto(5).bodies = [body(1), body(2)];
pto(5).attachments = [[-2.5, 4.33, 8.66], [-1.25, 2.165, 4.33]];

%% PTO Definitions
pto(6) = ptoClass('LinSpringDamper');
pto(6).stiffness = 30000.0;
pto(6).damping = 5000.0;
pto(6).rest_length = 2.0;
pto(6).bodies = [body(1), body(2)];
pto(6).attachments = [[-2.5, -4.33, 8.66], [-1.25, -2.165, 4.33]];

