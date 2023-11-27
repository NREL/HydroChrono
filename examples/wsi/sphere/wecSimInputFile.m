%% Simulation Data
simu = simulationClass();
simu.modelName = 'Sphere testing';                              
simu.simMechanicsFile = 'sphere.slx';
simu.explorer = 'on';                   % Turn SimMechanics Explorer (on/off)
simu.endTime = 40;                        
simu.dt = 0.01;                         
simu.cicEndTime = 15;

%% Wave Information  
% Regular Waves  
waves = waveClass('noWaveCIC');    

%% Body Data
% Sphere
body(1) = bodyClass('hydroData/sphere.h5');    	
body(1).geometryFile = 'geometry/sphere.obj';    
body(1).mass = 261800;                  
body(1).inertia = [20907301 21306090.66 37085481.11];     
body(1).position = [0.0, 0.0, 0.0];

%% PTO and Constraint Parameters
% Floating (3DOF) Joint
constraint(1) = constraintClass('Constraint1'); 
constraint(1).location = [0 0 0];                    