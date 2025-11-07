import numpy as np
import os
import gzip

# Create data directory if it doesn't exist
os.makedirs('data', exist_ok=True)

# Function to read the stars.dat file with error handling
def read_stars_data(filename):
    """
    Read stars.dat file handling potential formatting issues
    """
    # Check if file is gzipped
    if filename.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    data_list = []
    with opener(filename, mode) as f:
        for i, line in enumerate(f):
            # Skip empty lines
            if not line.strip():
                continue
            
            # Split the line and convert to float
            values = line.strip().split()
            
            # Check if we have the expected number of columns (6)
            if len(values) == 6:
                try:
                    row = [float(v) for v in values]
                    data_list.append(row)
                except ValueError as e:
                    print(f"Warning: Could not parse line {i+1}: {line.strip()}")
                    print(f"Error: {e}")
            else:
                print(f"Warning: Line {i+1} has {len(values)} columns instead of 6: {line.strip()}")
                # If you want to pad with zeros or handle differently:
                # if len(values) < 6:
                #     values.extend([0.0] * (6 - len(values)))
                #     row = [float(v) for v in values[:6]]
                #     data_list.append(row)
    
    return np.array(data_list)

# Try to read the file
try:
    # First try the regular file
    if os.path.exists('stars_cleaned.dat'):
        data = read_stars_data('stars_cleaned.dat')
    elif os.path.exists('stars.dat.gz'):
        data = read_stars_data('stars.dat.gz')
    else:
        raise FileNotFoundError("Neither 'stars.dat' nor 'stars.dat.gz' found")
    
    print(f"Successfully loaded {len(data)} rows of data")
    print(f"Data shape: {data.shape}")
    
    # If you prefer to use np.loadtxt with error handling:
    # Alternative approach using genfromtxt which is more forgiving
    # data = np.genfromtxt('stars.dat', invalid_raise=False, filling_values=0.0)
    # # Remove any rows with NaN values
    # data = data[~np.isnan(data).any(axis=1)]
    
except Exception as e:
    print(f"Error reading file: {e}")
    # Alternative: try using genfromtxt which is more flexible
    try:
        print("Trying alternative method with genfromtxt...")
        data = np.genfromtxt('stars.dat', invalid_raise=False)
        # Remove rows with any NaN values
        data = data[~np.isnan(data).any(axis=1)]
        print(f"Successfully loaded {len(data)} rows using genfromtxt")
    except Exception as e2:
        print(f"Alternative method also failed: {e2}")
        raise

# Continue with the rest of your processing...

# The file format is: z, y, x, vz, vy, vx
# We need to reorder to: x, y, z, vx, vy, vz
z = data[:, 0]
y = data[:, 1]
x = data[:, 2]
vz = data[:, 3]
vy = data[:, 4]
vx = data[:, 5]

# Number of particles (should be 43802)
n_bodies = len(x)
print(f"Read {n_bodies} particles (expected 43802)")

# Since masses are not provided, we'll assign equal masses
# Total galaxy mass is 2.77 x 10^10 solar masses
# But only 70% is in particles (30% is in rigid sphere)
total_mass_normalized = 1.0
particle_mass_fraction = 0.7
mass_per_particle = (total_mass_normalized * particle_mass_fraction) / n_bodies
masses = np.full(n_bodies, mass_per_particle)

# Set radius for all particles
radius = 0.05

# 1. Generate nstars.dat
with open('data/nstars.dat', 'w') as f:
    f.write(f"{n_bodies} // n_bodies \n")
    f.write(f"{np.mean(masses):.15f} // avg_mass \n")

# 2. Generate positions_ascii.dat
with open('data/positions_ascii.dat', 'w') as f:
    for i in range(n_bodies):
        f.write(f"{x[i]:.15f} {y[i]:.15f} {z[i]:.15f}\n")

# 3. Generate velocities_ascii.dat
with open('data/velocities_ascii.dat', 'w') as f:
    for i in range(n_bodies):
        f.write(f"{vx[i]:.15f} {vy[i]:.15f} {vz[i]:.15f}\n")

# 4. Generate masses_ascii.dat
with open('data/masses_ascii.dat', 'w') as f:
    for i in range(n_bodies):
        f.write(f"{masses[i]:.15e}\n")

# 5. Generate radii_ascii.dat
with open('data/radii_ascii.dat', 'w') as f:
    for i in range(n_bodies):
        f.write(f"{radius:.15f}\n")

# 6. Generate input.dat
# Calculate simulation parameters
dt = 0.05
total_time = 1000.0
n_frames = 400
nsteps = int(total_time / dt)  # 20000 steps
nprint = nsteps // n_frames    # 50 steps per frame
ndump = nprint                  # Same as nprint

# Find bounding box for initial conditions
x_min, x_max = np.min(x), np.max(x)
y_min, y_max = np.min(y), np.max(y)
z_min, z_max = np.min(z), np.max(z)

print(f"\nData extent:")
print(f"  x: [{x_min:.2f}, {x_max:.2f}]")
print(f"  y: [{y_min:.2f}, {y_max:.2f}]")
print(f"  z: [{z_min:.2f}, {z_max:.2f}]")

# Add some padding to the bounding box
padding = 5.0
L_bound_x = max(abs(x_min), abs(x_max)) + padding
L_bound_y = max(abs(y_min), abs(y_max)) + padding
L_bound_z = max(abs(z_min), abs(z_max)) + padding

with open('input.dat', 'w') as f:
    f.write("ndims   n_bodies        Force_Diagnose (only first timestep)\n")
    f.write(f"3       {n_bodies}             0\n")
    
    f.write("Gravitational_Constant \n")
    f.write("1.0\n")
    
    f.write("mass_0            radius_0        radius_var (0->0.9)\n")
    f.write("1                 0.05            0.0\n")
    
    f.write("dt                nsteps          nprint      ndump\n")
    f.write(f"{dt}             {nsteps}            {nprint}          {ndump}\n")
    
    f.write("MASS_INITIALIZATION_TYPE ! 0 = read from file, 1 = normalized\n")
    f.write("0\n")
    
    f.write("INITIALIZATION_TYPE         ! 0 = read from file \n")
    f.write("0\n")
    
    f.write("nbx_init(1)       nbx_init(2)      nbx_init(3) \n")
    f.write("1                 1                1\n")
    
    f.write("L_init(1)         L_init(2)        L_init(3) \n")
    f.write("1                 1                1\n")
    
    f.write("VELOCITY_INITIALIZATION_TYPE ! 0 = read from file \n")
    f.write("0\n")
    
    f.write("vel_var       omega_init\n")
    f.write("0.0           0.0\n")
    
    f.write("BOUNDARY_CONDITION_TYPE ! 1=open, 2=reflective, 3=periodic (not implemented yet) \n")
    f.write("1\n")
    
    f.write("iReflective_BC(1)     iReflective_BC(2)     iReflective_BC(3) ! is real so keep either 0.0 or 1.0 \n")
    f.write("0.0                   0.0                   0.0\n")
    
    f.write("L_bound(1)        L_bound(2)       L_bound(3) \n")
    f.write(f"{L_bound_x:.1f}                 {L_bound_y:.1f}                {L_bound_z:.1f}\n")

print(f"\nGenerated files for {n_bodies} particles")
print(f"Mass information:")
print(f"  Total galaxy mass: 2.77e10 solar masses")
print(f"  Particle mass fraction: 70%")
print(f"  Mass per particle: {mass_per_particle:.6e} solar masses")
print(f"\nSimulation parameters:")
print(f"  dt = {dt}")
print(f"  Total time = {total_time}")
print(f"  Number of steps = {nsteps}")
print(f"  Frames to output = {n_frames}")
print(f"  Steps per frame = {nprint}")
print(f"  Boundary box: [{-L_bound_x:.2f}, {L_bound_x:.2f}] x [{-L_bound_y:.2f}, {L_bound_y:.2f}] x [{-L_bound_z:.2f}, {L_bound_z:.2f}]")
print(f"\nNote: The original units are:")
print(f"  Length: 1.33 kpc")
print(f"  Time: 4.33 Myr")
print(f"  Velocity: Length/Time = 307 km/s")
