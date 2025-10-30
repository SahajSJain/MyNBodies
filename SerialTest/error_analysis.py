import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from scipy.integrate import cumulative_trapezoid as cumtrapz
import warnings
warnings.filterwarnings('ignore')

# Create output directory
os.makedirs('./Plots', exist_ok=True)

# Set matplotlib parameters for better plots
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.figsize'] = (6, 6)

def read_time_history():
    """Read the Time_History.dat file"""
    data = []
    with open('Time_History.dat', 'r') as f:
        # Skip header
        f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) > 0:
                data.append([int(parts[0]), float(parts[1])] + [float(p) for p in parts[2:]])
    
    df = pd.DataFrame(data, columns=['step', 'time', 'COM_x', 'COM_y', 'COM_z', 
                                     'COM_vx', 'COM_vy', 'COM_vz', 'KE_total', 
                                     'PE_total', 'E_total'])
    return df

def read_csv_files():
    """Read all nb.*.csv files and organize by timestep"""
    csv_files = sorted(glob.glob('NBDUMPCSV/nb.*.csv'))
    
    all_data = {}
    for file in csv_files:
        # Extract file number
        file_num = int(file.split('.')[-2])
        
        # Read CSV
        df = pd.read_csv(file)
        all_data[file_num] = df
    
    return all_data

def theoretical_trajectory(t, particle_id):
    """Calculate theoretical trajectory for particle 1 or 2"""
    if particle_id == 1:
        x = np.cos(t)
        y = np.sin(t)
        vx = -np.sin(t)
        vy = np.cos(t)
    else:  # particle 2
        x = -np.cos(t)
        y = -np.sin(t)
        vx = np.sin(t)
        vy = -np.cos(t)
    
    z = 0
    vz = 0
    
    # Acceleration (centripetal)
    r = 2  # distance between particles
    a_mag = 16 / (r**2 * 4)  # GM/r^2 where M=4, G=1
    
    if particle_id == 1:
        ax = -a_mag * np.cos(t)
        ay = -a_mag * np.sin(t)
    else:
        ax = a_mag * np.cos(t)
        ay = a_mag * np.sin(t)
    
    az = 0
    
    # Force = mass * acceleration
    mass = 4
    fx = mass * ax
    fy = mass * ay
    fz = mass * az
    
    # Kinetic energy per particle
    ke = 0.5 * mass * (vx**2 + vy**2)
    
    # Potential energy (shared between particles, so half for each)
    pe = -8 / 2  # -GM^2/r / 2
    
    # Total energy per particle
    te = ke + pe
    
    return {
        'x': x, 'y': y, 'z': z,
        'vx': vx, 'vy': vy, 'vz': vz,
        'ax': ax, 'ay': ay, 'az': az,
        'fx': fx, 'fy': fy, 'fz': fz,
        'ke': ke, 'pe': pe, 'te': te
    }

def calculate_errors(sim_data, theory_data, times):
    """Calculate mean and RMS integrated errors"""
    diff = sim_data - theory_data
    
    # Integrate absolute error over time
    integrated_error = cumtrapz(np.abs(diff), times, initial=0)
    mean_integrated_error = integrated_error[-1] / (times[-1] - times[0])
    
    # RMS integrated error
    integrated_squared_error = cumtrapz(diff**2, times, initial=0)
    rms_integrated_error = np.sqrt(integrated_squared_error[-1] / (times[-1] - times[0]))
    
    return mean_integrated_error, rms_integrated_error, diff

def plot_comparison(times, sim_data, theory_data, diff, ylabel, filename, title):
    """Create comparison plot with dual y-axes"""
    fig, ax1 = plt.subplots(figsize=(6, 6))
    
    # Left y-axis: theoretical vs simulation
    color1 = 'tab:blue'
    color2 = 'tab:red'
    ax1.set_xlabel('Time')
    ax1.set_ylabel(ylabel, color='black')
    line1 = ax1.plot(times, theory_data, color=color1, linewidth=2, label='Theoretical')
    line2 = ax1.plot(times, sim_data, color=color2, linewidth=2, linestyle='--', label='Simulation')
    ax1.tick_params(axis='y')
    ax1.grid(True, alpha=0.3)
    
    # Right y-axis: difference
    ax2 = ax1.twinx()
    color3 = 'tab:green'
    ax2.set_ylabel('Difference', color=color3)
    line3 = ax2.plot(times, diff, color=color3, linewidth=1, alpha=0.7, label='Difference')
    ax2.tick_params(axis='y', labelcolor=color3)
    
    # Combine legends
    lines = line1 + line2 + line3
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='best')
    
    plt.title(title)
    plt.tight_layout()
    plt.savefig(f'./Plots/{filename}', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    print("Reading Time_History.dat...")
    time_history = read_time_history()
    
    print("Reading CSV files...")
    csv_data = read_csv_files()
    
    # Initialize storage for time series data
    times = []
    particle1_data = {key: [] for key in ['x', 'y', 'z', 'vx', 'vy', 'vz', 
                                          'ax', 'ay', 'az', 'fx', 'fy', 'fz', 
                                          'ke', 'pe', 'te']}
    particle2_data = {key: [] for key in particle1_data.keys()}
    
    # Extract time series data
    for file_num in sorted(csv_data.keys()):
        # Get corresponding time from time_history
        if file_num < len(time_history):
            time = time_history.iloc[file_num]['time']
            times.append(time)
            
            # Get particle data
            df = csv_data[file_num]
            
            # Particle 1 (id=1)
            p1 = df[df['id'] == 1].iloc[0]
            for key in particle1_data.keys():
                particle1_data[key].append(p1[key])
            
            # Particle 2 (id=2)
            p2 = df[df['id'] == 2].iloc[0]
            for key in particle2_data.keys():
                particle2_data[key].append(p2[key])
    
    times = np.array(times)
    
    # Convert to numpy arrays
    for key in particle1_data.keys():
        particle1_data[key] = np.array(particle1_data[key])
        particle2_data[key] = np.array(particle2_data[key])
    
    # Calculate theoretical values
    theory1_data = {key: [] for key in particle1_data.keys()}
    theory2_data = {key: [] for key in particle2_data.keys()}
    
    for t in times:
        theory1 = theoretical_trajectory(t, 1)
        theory2 = theoretical_trajectory(t, 2)
        
        for key in theory1_data.keys():
            theory1_data[key].append(theory1[key])
            theory2_data[key].append(theory2[key])
    
    # Convert to numpy arrays
    for key in theory1_data.keys():
        theory1_data[key] = np.array(theory1_data[key])
        theory2_data[key] = np.array(theory2_data[key])
    
    # Open error analysis file
    with open('./Plots/error_analysis.dat', 'w') as f:
        f.write("Error Analysis for N-body Simulation vs Theoretical Predictions\n")
        f.write("="*70 + "\n\n")
        
        # Analyze and plot for each quantity
        quantities = [
            ('x', 'Position X (m)', 'position_x'),
            ('y', 'Position Y (m)', 'position_y'),
            ('vx', 'Velocity X (m/s)', 'velocity_x'),
            ('vy', 'Velocity Y (m/s)', 'velocity_y'),
            ('ke', 'Kinetic Energy', 'kinetic_energy'),
            ('pe', 'Potential Energy', 'potential_energy'),
            ('te', 'Total Energy', 'total_energy'),
            ('fx', 'Force X', 'force_x'),
            ('fy', 'Force Y', 'force_y')
        ]
        
        for key, ylabel, filename in quantities:
            # Particle 1
            mean_err1, rms_err1, diff1 = calculate_errors(
                particle1_data[key], theory1_data[key], times
            )
            plot_comparison(times, particle1_data[key], theory1_data[key], 
                          diff1, ylabel, f'{filename}_particle1.png', 
                          f'Particle 1: {ylabel}')
            
            # Particle 2
            mean_err2, rms_err2, diff2 = calculate_errors(
                particle2_data[key], theory2_data[key], times
            )
            plot_comparison(times, particle2_data[key], theory2_data[key], 
                          diff2, ylabel, f'{filename}_particle2.png', 
                          f'Particle 2: {ylabel}')
            
            # Write errors to file
            f.write(f"{ylabel}:\n")
            f.write(f"  Particle 1 - Mean Integrated Error: {mean_err1:.6e}, RMS Integrated Error: {rms_err1:.6e}\n")
            f.write(f"  Particle 2 - Mean Integrated Error: {mean_err2:.6e}, RMS Integrated Error: {rms_err2:.6e}\n")
            f.write(f"  Average   - Mean Integrated Error: {(mean_err1+mean_err2)/2:.6e}, RMS Integrated Error: {(rms_err1+rms_err2)/2:.6e}\n")
            f.write("\n")
            
            print(f"Processed {ylabel}")
        
        # Plot trajectories in xy plane
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        
        # Particle 1 trajectory
        ax1.plot(theory1_data['x'], theory1_data['y'], 'b-', linewidth=2, label='Theoretical')
        ax1.plot(particle1_data['x'], particle1_data['y'], 'r--', linewidth=2, label='Simulation')
        ax1.set_xlabel('X Position')
        ax1.set_ylabel('Y Position')
        ax1.set_title('Particle 1 Trajectory')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal')
        
        # Particle 2 trajectory
        ax2.plot(theory2_data['x'], theory2_data['y'], 'b-', linewidth=2, label='Theoretical')
        ax2.plot(particle2_data['x'], particle2_data['y'], 'r--', linewidth=2, label='Simulation')
        ax2.set_xlabel('X Position')
        ax2.set_ylabel('Y Position')
        ax2.set_title('Particle 2 Trajectory')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig('./Plots/trajectories_xy.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Calculate trajectory errors
        traj_err1 = np.sqrt((particle1_data['x'] - theory1_data['x'])**2 + 
                           (particle1_data['y'] - theory1_data['y'])**2)
        traj_err2 = np.sqrt((particle2_data['x'] - theory2_data['x'])**2 + 
                           (particle2_data['y'] - theory2_data['y'])**2)
        
        mean_traj_err1 = np.mean(traj_err1)
        mean_traj_err2 = np.mean(traj_err2)
        
        f.write("\nTrajectory Error (Euclidean distance in xy-plane):\n")
        f.write(f"  Particle 1 - Mean Error: {mean_traj_err1:.6e}\n")
        f.write(f"  Particle 2 - Mean Error: {mean_traj_err2:.6e}\n")
        f.write(f"  Average   - Mean Error: {(mean_traj_err1+mean_traj_err2)/2:.6e}\n")
    
    print("\nAnalysis complete! Results saved in ./Plots/")
    print("- error_analysis.dat contains numerical error analysis")
    print("- PNG files show visual comparisons")

if __name__ == "__main__":
    main()
