
import importlib.util
import subprocess
import os




# === Default flags
FFLAGS = ''
FFLAGS = '-O3'
OMP_THREADS = str(os.cpu_count()-2)
OMP_SCHEDULE= "static,50"
ACC_DEVICE_NUM="0"
CUDA_VISIBLE_DEVICES="0"



# === src files
src_cpu = 'src/MT_OpenMPsim.f90'
src_gpu = 'src/GPU_OpenACCsim.f90'

# === hidden var, set in .f90 directly
hidden_nstars = 2000



os.system("clear")

def is_installed(package_name):
    """Check if a package is installed."""
    return importlib.util.find_spec(package_name) is not None

# ANSI color codes
VYAN = "\033[96m"
MAGENTA = "\033[95m"
YELLOW = "\033[93m"
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"
CHECK = "\u2714"  # ✔
CROSS = "\u2718"  # ✘


print(MAGENTA, " === Python dependencies check ===", RESET)

# Python packages checks
packages = ["numpy", "matplotlib"]

deps_met = True

for pkg in packages:
    if is_installed(pkg):
        print(f"{GREEN}{CHECK} {pkg} {RESET}")
    else:
        print(f"{RED}{CROSS} {pkg} {RESET}")
        deps_met = False



print(MAGENTA, " === Compilers check ===", RESET)

# compilers check
def check_compiler(compiler):
    try:
        subprocess.run(["which", compiler], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False



compiler = "gfortran"
if check_compiler(compiler):
    print(f"{GREEN}{CHECK} {compiler} \n     => CPU and OpenMP compilation is available{RESET}")
else:
    print(f"{RED}{CROSS} {compiler} \n     => CPU and OpenMP compilation is NOT available. THIS IS BLOCKING.{RESET}")
    deps_met = True


is_gpu_available = False
compiler = "pgfortran"
if check_compiler(compiler):
    print(f"{GREEN}{CHECK} {compiler} \n     => Legacy GPU compilation is available{RESET}")
    is_gpu_available = True
else:
    print(f"{RED}{CROSS} {compiler} \n     => Legacy GPU compilation is NOT available{RESET}")

# compiler = "nvfortran"
# if check_compiler(compiler):
#     print(f"{GREEN}{CHECK} {compiler} \n     => Untested GPU compilation is available{RESET}")
#     is_gpu_available = True
# else:
#     print(f"{RED}{CROSS} {compiler} \n     => Untested GPU compilation is NOT available{RESET}")
# print(YELLOW, """
#  => nvfortran is non blocking as it is not used.  <=
#  => You will need to change the compiler files if <=
#  =>        you do not have pgfortran              <=""")
print(YELLOW," ====> GPU dpendencies are non-blocking here.\n\n\n", RESET)




if deps_met:
    print(GREEN, " All dependencies seem to be installed!")
else:
    print(RED, "  You are missing essential dependencies  \n  Please install the missing deps or modify main.py to skip deps check.", RESET)
    quit()





ansi_art = """
\033[49m       \033[38;5;15;49m▄▄▄\033[48;5;15m  \033[38;5;15;49m▄▄▄▄\033[49m                                            \033[m
\033[49m    \033[38;5;15;49m▄\033[38;5;15;48;5;15m▄\033[49;38;5;15m▀\033[49m        \033[49;38;5;15m▀▀\033[38;5;15;48;5;15m▄\033[38;5;15;49m▄▄\033[49m                               \033[38;5;33;49m▄\033[48;5;33m  \033[38;5;33;49m▄\033[49m     \033[m
\033[49m   \033[38;5;15;48;5;15m▄\033[49;38;5;15m▀\033[49m              \033[49;38;5;15m▀▀\033[38;5;15;49m▄\033[49m                    \033[38;5;15;49m▄▄\033[48;5;15m \033[49;38;5;15m▀▀▀▀▀\033[48;5;33m      \033[49m    \033[m
\033[49m \033[38;5;15;49m▄\033[48;5;15m \033[49m                   \033[49;38;5;15m▀\033[48;5;15m \033[38;5;15;49m▄\033[49m             \033[38;5;15;49m▄▄\033[49;38;5;15m▀▀\033[49m        \033[49;38;5;33m▀\033[48;5;33m    \033[38;5;15;48;5;33m▄\033[38;5;15;49m▄\033[49m   \033[m
\033[49m \033[48;5;15m \033[49m                       \033[49;38;5;15m▀\033[38;5;15;49m▄▄\033[49m       \033[38;5;15;49m▄\033[38;5;15;48;5;15m▄\033[49;38;5;15m▀\033[49m                  \033[49;38;5;15m▀\033[48;5;15m \033[49m  \033[m
\033[49m \033[48;5;15m \033[49m                         \033[49;38;5;15m▀\033[38;5;9;48;5;202m▄\033[48;5;9m    \033[38;5;15;48;5;15m▄\033[49;38;5;15m▀\033[49m                       \033[48;5;15m \033[49m \033[m
\033[49m \033[48;5;15m \033[49m                         \033[49;38;5;9m▀\033[48;5;9m     \033[49;38;5;9m▀\033[49m                        \033[48;5;15m \033[38;5;15;49m▄\033[m
\033[49m \033[48;5;15m \033[38;5;15;49m▄\033[49m                       \033[38;5;15;49m▄\033[48;5;15m \033[49;38;5;9m▀▀\033[38;5;9;48;5;9m▄\033[49;38;5;9m▀▀\033[48;5;15m \033[38;5;15;49m▄\033[49m                       \033[48;5;15m \033[49;38;5;15m▀\033[m
\033[49m  \033[49;38;5;15m▀\033[38;5;15;49m▄\033[49m                   \033[38;5;15;49m▄\033[48;5;15m \033[49;38;5;15m▀\033[49m         \033[49;38;5;15m▀\033[48;5;15m \033[38;5;15;49m▄\033[49m                    \033[48;5;15m \033[49m \033[m
\033[49m    \033[38;5;47;48;5;15m▄\033[48;5;47m    \033[38;5;47;49m▄\033[49m         \033[38;5;15;49m▄▄\033[49;38;5;15m▀▀\033[49m               \033[49;38;5;15m▀\033[48;5;15m \033[38;5;15;49m▄\033[49m               \033[38;5;15;49m▄\033[38;5;15;48;5;15m▄\033[49m  \033[m
\033[49m    \033[48;5;47m      \033[38;5;15;49m▄▄▄▄▄\033[38;5;15;48;5;15m▄\033[49;38;5;15m▀▀▀\033[49m                      \033[49;38;5;15m▀▀\033[38;5;15;49m▄▄▄\033[49m      \033[38;5;15;49m▄▄▄\033[49;38;5;15m▀▀\033[49m   \033[m
\033[49m     \033[49;38;5;47m▀▀▀▀\033[49m                                     \033[49;38;5;15m▀▀▀▀▀▀▀\033[49m       \033[m
\033[49m                                                            \033[m

\033[96m        -====-      N-Body Simulation      -====- \033[0m
  \033[0m                 Author - Léo BECHET
                  Licence - CC-BY-NC-SA
  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \033[0m
"""
print(ansi_art)





def SetDefaultSettings():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES

    FFLAGS = '-O3'
    OMP_THREADS = str(os.cpu_count()-2)
    OMP_SCHEDULE= "static,50"
    ACC_DEVICE_NUM="0"
    CUDA_VISIBLE_DEVICES="0"


def ShowEnvVar():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES

    print(f' FFLAGS _____________ : {FFLAGS}')
    print(f' OMP_THREADS          : {OMP_THREADS}')
    print(f' OMP_SCHEDULE _______ : {OMP_SCHEDULE}')
    print(f' ACC_DEVICE_NUM       : {ACC_DEVICE_NUM}')
    print(f' CUDA_VISIBLE_DEVICES : {CUDA_VISIBLE_DEVICES}')
    Menu()

def EditEnvVar():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES

    print(f' 0) Cancel')
    print(f' 1) FFLAGS _____________ : {FFLAGS}')
    print(f' 2) OMP_THREADS          : {OMP_THREADS}')
    print(f' 3) OMP_SCHEDULE _______ : {OMP_SCHEDULE}')
    print(f' 4) ACC_DEVICE_NUM       : {ACC_DEVICE_NUM}')
    print(f' 5) CUDA_VISIBLE_DEVICES : {CUDA_VISIBLE_DEVICES}')

    PASS = True
    while PASS:
        PASS = False
        inp = int(input(f'\n Which Env Var do you wish to edit? (number)\n > '))

        match inp:
            case 1:
                FFLAGS=str(input("value > "))
            case 2:
                OMP_THREADS=str(input("value > "))
            case 3:
                OMP_SCHEDULE=str(input("value > "))
            case 4:
                ACC_DEVICE_NUM=str(input("value > "))
            case 5:
                CUDA_VISIBLE_DEVICES=str(input("value > "))

            case _:
                PASS = True
                print(" ID not recognised")
    print("\n")
    Menu()



import numpy as np

def generate_points_sphere(n_points, radius=1.0, offset=(0,0,0)):
    n_points = int(n_points)

    # Ensure even number of points for pairing
    if n_points % 2 != 0:
        print(f"Warning: n_points must be even for COM at origin. Using {n_points-1} points instead.")
        n_points = n_points - 1

    points = []
    n_pairs = n_points // 2

    # Calculate representative distance based on volume per particle
    sphere_volume = (4/3) * np.pi * radius**3
    volume_per_particle = sphere_volume / n_points
    # Approximate distance between particles (cube root of volume per particle)
    min_distance = 0.5 * (volume_per_particle ** (1/3))

    max_attempts = n_pairs * 100  # Prevent infinite loop
    attempts = 0
    pairs_added = 0

    while pairs_added < n_pairs and attempts < max_attempts:
        attempts += 1

        # Generate random point in sphere
        x, y, z = np.random.uniform(-radius, radius, 3)
        if x**2 + y**2 + z**2 <= radius**2:
            # Create the pair: original and diametrically opposite
            point1 = [x, y, z]
            point2 = [-x, -y, -z]  # Diametrically opposite

            # Check distance to all existing points for both points in the pair
            too_close = False

            for existing_point in points:
                # Check distance from point1
                dist1 = np.sqrt((point1[0] - existing_point[0])**2 +
                               (point1[1] - existing_point[1])**2 +
                               (point1[2] - existing_point[2])**2)
                # Check distance from point2
                dist2 = np.sqrt((point2[0] - existing_point[0])**2 +
                               (point2[1] - existing_point[1])**2 +
                               (point2[2] - existing_point[2])**2)

                if dist1 < min_distance or dist2 < min_distance:
                    too_close = True
                    break

            # Also check distance between the pair itself
            pair_distance = np.sqrt((point1[0] - point2[0])**2 +
                                   (point1[1] - point2[1])**2 +
                                   (point1[2] - point2[2])**2)
            if pair_distance < min_distance:
                too_close = True

            if not too_close:
                # Add both points with offset
                points.append([point1[0] + offset[0],
                              point1[1] + offset[1],
                              point1[2] + offset[2]])
                points.append([point2[0] + offset[0],
                              point2[1] + offset[1],
                              point2[2] + offset[2]])
                pairs_added += 1

    if pairs_added < n_pairs:
        print(f"Warning: Could only place {pairs_added*2} out of {n_points} points with minimum distance {min_distance:.3f}")

    # Verify COM is at origin (plus offset)
    if len(points) > 0:
        points_array = np.array(points)
        com = np.mean(points_array, axis=0)
        expected_com = np.array(offset)
        com_error = np.linalg.norm(com - expected_com)
        if com_error > 1e-10:
            print(f"Warning: COM error = {com_error:.2e}")
        else:
            print(f"COM verified at {com} (offset: {offset})")

    return points

def generate_velocity(pos_data, omega_z=1):
    x, y, z = pos_data[:, 0], pos_data[:, 1], pos_data[:, 2]

    vx = -omega_z * y
    vy = omega_z * x
    vz = np.zeros_like(z)

    velocities = np.column_stack((vx, vy, vz))

    return velocities

def calculate_mass_for_virial_equilibrium(points, velocities, G=1.0, softening=0.0):
    """
    Calculate the mass per particle needed to achieve KE = PE (virial equilibrium).

    For a system in virial equilibrium: 2*KE + PE = 0, which means KE = -PE/2
    For gravitational systems, PE is negative, so KE = |PE|/2

    Parameters:
    -----------
    points : array_like, shape (n, 3)
        Positions of particles
    velocities : array_like, shape (n, 3)
        Velocities of particles
    G : float
        Gravitational constant (default=1.0)
    softening : float
        Softening length to avoid singularities (default=0.0)

    Returns:
    --------
    mass : float
        Mass per particle to achieve KE = PE
    """
    points = np.array(points)
    velocities = np.array(velocities)
    n_particles = len(points)

    if n_particles < 2:
        raise ValueError("Need at least 2 particles")

    # Calculate total kinetic energy coefficient (without mass)
    # KE = 0.5 * m * sum(v^2) for all particles
    # Since all particles have same mass: KE = 0.5 * m * n * sum(v^2)
    v_squared_sum = np.sum(velocities**2)
    ke_coefficient = 0.5 * v_squared_sum

    # Calculate potential energy coefficient (without mass^2)
    # PE = -G * sum_i sum_j>i (m_i * m_j / r_ij)
    # Since all masses are equal: PE = -G * m^2 * sum_i sum_j>i (1 / r_ij)
    pe_coefficient = 0.0

    for i in range(n_particles):
        for j in range(i + 1, n_particles):
            # Calculate distance between particles i and j
            r_vec = points[i] - points[j]
            r_squared = np.sum(r_vec**2) + softening**2
            r = np.sqrt(r_squared)
            pe_coefficient += 1.0 / r

    pe_coefficient *= -G

    # For virial equilibrium: KE = -PE/2
    # 0.5 * m * v_squared_sum = -(-G * m^2 * pe_coefficient) / 2
    # 0.5 * m * v_squared_sum = G * m^2 * pe_coefficient / 2
    # m * v_squared_sum = G * m^2 * pe_coefficient
    # v_squared_sum = G * m * pe_coefficient
    # m = v_squared_sum / (G * pe_coefficient)

    # However, for the condition KE = PE (not virial equilibrium):
    # 0.5 * m * v_squared_sum = |G * m^2 * pe_coefficient|
    # Since pe_coefficient is already negative:
    # 0.5 * m * v_squared_sum = -G * m^2 * pe_coefficient
    # 0.5 * v_squared_sum = -G * m * pe_coefficient
    # m = -0.5 * v_squared_sum / (G * pe_coefficient)

    if pe_coefficient >= 0:
        raise ValueError("Potential energy coefficient should be negative")

    mass = -0.5 * v_squared_sum / (G * pe_coefficient)

    # Verify the calculation
    total_ke = 0.5 * mass * v_squared_sum
    total_pe = G * mass**2 * pe_coefficient

    print(f"Calculated mass per particle: {mass:.6f}")
    print(f"Total KE: {total_ke:.6f}")
    print(f"Total PE: {total_pe:.6f}")
    print(f"KE/|PE| ratio: {total_ke/abs(total_pe):.6f}")

    return mass

def GenerateData():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES, N_STARS
    
    # Create data folder if it doesn't exist
    os.makedirs('data', exist_ok=True)
    
    print(YELLOW, " Generating test data", RESET)
    
    pos_data = np.array(generate_points_sphere(hidden_nstars, 1.0))
    # Save as binary (keep transpose for backward compatibility)
    pos_data.T.tofile('data/positions.dat')
    
    # Save as ASCII - NO TRANSPOSE, each row is one particle
    np.savetxt('data/positions_ascii.dat', pos_data, fmt='%.16e', delimiter=' ')
    
    vel_data = generate_velocity(pos_data)
    
    # Save velocities
    vel_data.T.tofile('data/velocities.dat')
    
    # Save as ASCII - NO TRANSPOSE, each row is one particle  
    np.savetxt('data/velocities_ascii.dat', vel_data, fmt='%.16e', delimiter=' ') 
    mass = calculate_mass_for_virial_equilibrium(pos_data, vel_data)
    # Save number of stars to separate file
    with open('data/nstars.dat', 'w') as f:
        f.write(f"{hidden_nstars}\n{mass}")
    print("No. of stars = ", hidden_nstars)
    print("inv of stars = ", 1/hidden_nstars)

def CompileAndExec_CPU():

    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES, N_STARS

    print(YELLOW, " Generating test data", RESET)
    pos_data = np.array(generate_points_sphere(hidden_nstars, 1.0, offset=(0., 0., 0.)))
    pos_data.T.tofile('data/positions.dat')
    np.savetxt('data/positions_ascii.dat', pos_data.T, fmt='%.16e', delimiter=' ')
    vel_data = generate_velocity(pos_data)
    vel_data.T.tofile('data/velocities.dat')
    np.savetxt('data/velocities_ascii.dat', vel_data.T, fmt='%.16e', delimiter=' ')
    print(" Compiling program")
    i = os.system(f'gfortran -fopenmp {FFLAGS} -o bin/CPU.o {src_cpu}')

    if i == 0 :
        os.environ['OMP_SCHEDULE'] = OMP_SCHEDULE
        os.environ['OMP_NUMT_THREADS'] = OMP_THREADS
        os.system("bin/CPU.o")
        plot_energy()
        anim_realTime('out_data/position.dat')
    else:
        print(RED, "Error during compilation, aborting.", RESET)

    Menu()









def CompileAndExec_GPU():
    global FFLAGS, OMP_THREADS, OMP_SCHEDULE, ACC_DEVICE_NUM, CUDA_VISIBLE_DEVICES, N_STARS
    print(YELLOW, "  This program does not check CUDA capabilities on your machine and assume it is installed and configured")
    print("  Please specify your GPU architecture.")
    print("  With XXX as your GPU architecture. For QUADRO RTX 5000, use cc75 (default, enter to use)", RESET)
    GPU_arch = input(MAGENTA+"default=cc75 > "+RESET)
    if GPU_arch == "" : GPU_arch="cc75"

    print(YELLOW, " Generating test data", RESET)
    pos_data = np.array(generate_points_sphere(hidden_nstars, 1.0, offset=(0, 0.25, 0.25)))
    pos_data.T.tofile('data/positions.dat')
    vel_data = generate_velocity(pos_data)
    vel_data.T.tofile('data/velocities.dat')
    print(YELLOW, " === Compiling program ===", RESET)
    i = os.system(f'pgfortran -acc {FFLAGS} -gpu={GPU_arch} -o bin/GPU.o {src_gpu}') # Modify nvfortran here if needed

    if i == 0 :
        os.environ['ACC_DEVICE_NUM'] = ACC_DEVICE_NUM
        os.environ['CUDA_VISIBLE_DEVICES'] = CUDA_VISIBLE_DEVICES
        os.system("bin/GPU.o")
        anim_realTime('out_data/position.dat')

    else:
        print(RED, "Error during compilation, aborting.", RESET)

    Menu()













import matplotlib.pyplot as plt
import matplotlib.animation

def plot_energy():
    data = np.loadtxt('out_data/energy.dat')
    potential_energy = data[:, 0]
    kinetic_energy = data[:, 1]
    total_energy = kinetic_energy + potential_energy
    plt.figure(figsize=(10, 6))
    plt.plot(kinetic_energy, label='Kinetic Energy', color='blue')
    plt.plot(potential_energy, label='Potential Energy', color='red')
    plt.plot(total_energy, label='Total Energy', color='green', linestyle='--')
    plt.xlabel('Time Step (T)')
    plt.ylabel('Energy')
    plt.title('Kinetic, Potential, and Total Energy vs Time')
    plt.legend()
    plt.grid(True)
    plt.show()





















def load_file(filename: str, nbline: int = 0, timestep: int = -1, param: list = []):
    with open(filename, "rb") as file:
        if timestep == -1:
            record_size = np.fromfile(file, dtype=np.int32, count=1)  # Skip record size marker
            metadata = np.fromfile(file, dtype=np.int32, count=3)    # Read metadata
            record_size = np.fromfile(file, dtype=np.int32, count=1)  # Skip end record marker
            return metadata
        else:
            if param is None or len(param) < 3:
                raise ValueError("Parameter list 'param' must include [num_particles, num_timesteps, float_precision].")
            num_particles, num_timesteps, float_precision = param
            dtype = np.float32 if float_precision == 4 else np.float64
            record_bytes = 4 + (num_particles * 3 * dtype().nbytes) + 4
            file.seek(4 + (record_bytes * timestep), 0)  # 4 bytes for record size + N records
            _ = np.fromfile(file, dtype=np.int32, count=1)  # Skip record size marker
            data = np.fromfile(file, dtype=dtype, count=num_particles * 3)  # Read positions
            _ = np.fromfile(file, dtype=np.int32, count=1)  # Skip end record marker
            data = data.reshape(3, nbline)
            return data

def anim_realTime(file: str):
    figs = (720, 720)
    def update_graph(t: int):
        data = load_file(file, param[0], t, param)
        graph._offsets3d = (data[0, :], data[1, :], data[2, :])
    param = load_file(file)
    data = load_file(file, param[0], 0, param)
    dpi_v = 100
    sfig = (figs[0]/dpi_v, figs[1]/dpi_v)
    fig = plt.figure(figsize = sfig, facecolor = 'k', num = 1)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1, 1, 1])
    ax.axes.set_xlim3d(left=-3, right=3)
    ax.axes.set_ylim3d(bottom=-3, top=3)
    ax.axes.set_zlim3d(bottom=-3, top=3)
    ax.tick_params(axis='x', colors='w', which='both')
    ax.tick_params(axis='y', colors='w', which='both')
    ax.tick_params(axis='z', colors='w', which='both')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]['color'] =  (0.25,0.25,0.25,1)
    ax.yaxis._axinfo["grid"]['color'] =  (0.25,0.25,0.25,1)
    ax.zaxis._axinfo["grid"]['color'] =  (0.25,0.25,0.25,1)
    ax.set_facecolor("k")
    plt.tight_layout(pad=1)
    graph = ax.scatter(data[0, :], data[1, :], data[2, :], s=2, c='w', linewidths=0)
    ani = matplotlib.animation.FuncAnimation(fig, update_graph, param[1],
                            interval=1, blit=False)
    plt.show()









def Menu():
    print(YELLOW, "\n ======   M E N U   ======", RESET)
    print(f'  0) Quit')
    print(f'  1) Show Env Var')
    print(f'  2) Edit Env Var')
    print(GREEN, f' 3) Compile and exec CPU version', RESET)

    if is_gpu_available:
        print(GREEN,f' 4) Compile and exec GPU version',RESET)
    else:
        print(RED, f' 4) Compile and exec GPU version [UNAVAILABLE]', RESET)

    print(f'  5) Reset Env Var')
    print(f'  6) Generate Data')

    print("\n")

    PASS = True
    while PASS:
        inp = int(input(" > "))
        match inp:
            case 0:
                quit()
            case 1:
                ShowEnvVar()
            case 2:
                EditEnvVar()
                ShowEnvVar()
            case 3:
                CompileAndExec_CPU()
            case 4:
                if is_gpu_available:
                    CompileAndExec_GPU()
                else:
                    print(RED, "   GPU is not available", RESET)
            case 5:
                SetDefaultSettings()
                ShowEnvVar()

            case 6: 
                GenerateData() 
            case _:
                PASS = True
                print(" ID not recognised")


ShowEnvVar()

Menu()















