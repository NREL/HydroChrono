import subprocess

class SimulationConfig:
    def __init__(self, model_name='TALOS', explorer='on', dt=0.01, end_time=200.0):
        self.model_name = model_name
        self.explorer = explorer
        self.dt = dt
        self.end_time = end_time

class WaveConfig:
    def __init__(self, height=0.1, period=8):
        self.height = height
        self.period = period

class BodyConfig:
    def __init__(self, hydroDataFile, geometryFile, mass, inertia, position, orientation, radius, type):
        self.hydroDataFile = hydroDataFile
        self.geometryFile = geometryFile
        self.mass = mass
        self.inertia = inertia
        self.position = position
        self.orientation = orientation
        self.radius = radius
        self.type = type

class PTOConfig:
    def __init__(self, pto_type, stiffness, damping, rest_length, bodies, attachments):
        self.pto_type = pto_type
        self.stiffness = stiffness
        self.damping = damping
        self.rest_length = rest_length
        self.bodies = bodies
        self.attachments = attachments

def generate_input_file(filename, sim_config, wave_config, body_configs, pto_configs):
    with open(filename, 'w') as f:
        f.write("%% Simulation Data\n")
        f.write(f"simu = simulationClass();\n")
        f.write(f"simu.modelName = '{sim_config.model_name}';\n")
        f.write(f"simu.explorer = '{sim_config.explorer}';\n")
        f.write(f"simu.dt = {sim_config.dt};\n")
        f.write(f"simu.endTime = {sim_config.end_time};\n\n")

        f.write("%% Wave Information - Regular Waves\n")
        f.write("waves = waveClass('regular');\n")
        f.write(f"waves.height = {wave_config.height};\n")
        f.write(f"waves.period = {wave_config.period};\n\n")

        # Write body data
        for i, body in enumerate(body_configs, start=1):
            f.write(f"%% Body Data - {body.type}\n")
            f.write(f"body({i}) = bodyClass('{body.hydroDataFile}');\n")
            if body.geometryFile:
                f.write(f"body({i}).geometryFile = '{body.geometryFile}';\n")
            f.write(f"body({i}).mass = {body.mass};\n")
            f.write(f"body({i}).inertia = {body.inertia};\n")
            f.write(f"body({i}).position = {body.position};\n")
            if body.orientation:
                f.write(f"body({i}).orientation = {body.orientation};\n")
            if body.type == 'nonhydro-sphere':
                f.write(f"body({i}).type = '{body.type}';\n")
                f.write(f"body({i}).radius = {body.radius};\n")
            f.write("\n")

        # Write PTO data
        for i, pto in enumerate(pto_configs, start=1):
            f.write(f"%% PTO Definitions\n")
            f.write(f"pto({i}) = ptoClass('{pto.pto_type}');\n")
            f.write(f"pto({i}).stiffness = {pto.stiffness};\n")
            f.write(f"pto({i}).damping = {pto.damping};\n")
            f.write(f"pto({i}).rest_length = {pto.rest_length};\n")
            f.write(f"pto({i}).bodies = [body({pto.bodies[0]}), body({pto.bodies[1]})];\n")
            f.write(f"pto({i}).attachments = {pto.attachments};\n")
            f.write("\n")

# Creating simulation & wave configurations
sim_config = SimulationConfig(model_name='TALOS_testing', explorer='on', dt=(20.0/1999), end_time=150.0)
wave_config = WaveConfig(height=0.2, period=9)

# Define bodies in the system (hydrodynamic bodies first)
body1 = BodyConfig('../hydroData/talos.h5', '../geometry/talos_hull.obj', 1873881, [237600000, 237600000, 237600000], 
                   [0.0, 0.0, -7.96], [0.0, 0.0, 0.0], None, '')
body2 = BodyConfig('', '', 1873881, [106150000, 106150000, 106150000], 
                   [0.0, 0.0, -5.7], None, 4.0, 'nonhydro-sphere')

# Define 'lower' 3 PTOs
pto1 = PTOConfig('LinSpringDamper', 3e4, .5e4, 5.0, [1, 2], [[5.000, 0.000, -8.660], [2.500, 0.000, -4.330]])
pto2 = PTOConfig('LinSpringDamper', 3e4, .5e4, 5.0, [1, 2], [[-2.500, 4.330, -8.660], [-1.250, 2.165, -4.330]])
pto3 = PTOConfig('LinSpringDamper', 3e4, .5e4, 5.0, [1, 2], [[-2.500, -4.330, -8.660], [-1.250, -2.165, -4.330]])

# # Define 'upper' 3 PTOs
pto4 = PTOConfig('LinSpringDamper', 3e4, .5e4, 2.0, [1, 2], [[5.000, 0.000, 8.660], [2.500, 0.000, 4.330]])
pto5 = PTOConfig('LinSpringDamper', 3e4, .5e4, 2.0, [1, 2], [[-2.500, 4.330, 8.660], [-1.250, 2.165, 4.330]])
pto6 = PTOConfig('LinSpringDamper', 3e4, .5e4, 2.0, [1, 2], [[-2.500, -4.330, 8.660], [-1.250, -2.165, 4.330]])

# Assembles bodies, PTOs & write the input file
body_configs = [body1, body2]
pto_configs = [pto1, pto2, pto3, pto4, pto5, pto6]
input_file_path = "./wecSimInputFile.m"
print("Generating the input file...")
generate_input_file(input_file_path, sim_config, wave_config, body_configs, pto_configs)
print(f"Input file '{input_file_path}' generated successfully.")


# Path to your hydrochrono_wsi executable
executable_path = "[build_location]/demo_hydrochrono_wsi.exe"

# Call the executable with the input file as an argument
print(f"Calling the hydrochrono executable '{executable_path}' with the input file '{input_file_path}'.")
subprocess.run([executable_path, input_file_path, './results'])
print("Executable finished running.")