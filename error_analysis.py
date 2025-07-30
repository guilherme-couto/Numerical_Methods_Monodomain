import os
import numpy as np
import matplotlib.pyplot as plt
import csv
import re

BASE_DIR = 'outputs/error_analysis_aniso_rotation'
METHODS = ['fe', 'osadi', 'ssiadi', 'do', 'cs', 'hv', 'cg']
METHODS_FILE = 'Vm_50000.vtk'  # Solution file for each method
REFERENCE_PATH = os.path.join(BASE_DIR, 'reference', 'frames', 'Vm_500000.vtk')


def read_vtk_scalar_field(filepath):
    """Read a VTK ASCII file and return dimensions, spacing, and flattened scalar field."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    dims = tuple(map(int, lines[4].strip().split()[1:4]))
    spacing = tuple(map(float, lines[6].strip().split()[1:4]))
    data_start = lines.index("LOOKUP_TABLE default\n") + 1
    data = np.array([float(x.strip()) for x in lines[data_start:]])
    
    return dims, spacing, data


def reshape_field(data, dims):
    """Reshape 1D data to 2D (assuming z = 1)."""
    nx, ny, _ = dims
    return data.reshape((ny, nx))  # VTK orders as (y, x)


def sample_reference_to_coarse(ref_data, ref_dims, ref_spacing, coarse_dims, coarse_spacing):
    """Sample fine reference solution to match coarse grid points."""
    factor_x = int(round(coarse_spacing[0] / ref_spacing[0]))
    factor_y = int(round(coarse_spacing[1] / ref_spacing[1]))
    sampled_ref = ref_data[::factor_y, ::factor_x]
    expected_shape = (coarse_dims[1], coarse_dims[0])  # (ny, nx)
    assert sampled_ref.shape == expected_shape, f"Shape mismatch: {sampled_ref.shape} vs {expected_shape}"
    return sampled_ref


def compute_l2_error(u_numeric, u_ref):
    """Compute relative L2 norm error."""
    diff = u_numeric - u_ref
    N = diff.size
    l2_error = np.linalg.norm(diff) / np.sqrt(N)
    # l2_error = np.linalg.norm(diff) / np.linalg.norm(u_ref)
    return l2_error


def parse_simulation_info(filepath):
    """Parse simulation_infos.txt and extract timing data."""
    timings = {}
    pattern = re.compile(r'(.*)\s*=\s*([0-9.]+)\s*s')
    with open(filepath, 'r') as f:
        for line in f:
            match = pattern.match(line.strip())
            if match:
                key = match.group(1).strip()
                value = float(match.group(2))
                timings[key] = value
    return timings


def plot_total_execution_times(total_times):
    """Plot total execution times for each method."""
    methods = list(total_times.keys())
    values = [total_times[m] for m in methods]
    
    plt.figure(figsize=(8, 5))
    plt.bar([m.upper() for m in methods], values, color='orange')
    plt.ylabel('Total Execution Time (s)')
    plt.title('Execution Time Comparison Between Methods')
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    out_path = os.path.join(BASE_DIR, 'execution_time_plot.png')
    plt.savefig(out_path)
    plt.close()
    print(f"Saved execution time plot to {out_path}")


def save_summary_csv(errors, all_timings):
    """Save summary table combining errors and execution timings."""
    all_keys = set()
    for method_timings in all_timings.values():
        all_keys.update(method_timings.keys())
    all_keys = sorted(all_keys)

    summary_path = os.path.join(BASE_DIR, 'summary_table.csv')
    with open(summary_path, 'w', newline='') as f:
        writer = csv.writer(f)
        header = ['Method', 'Relative_L2_Error'] + all_keys
        writer.writerow(header)
        for method in METHODS:
            row = [method.upper(), errors.get(method, '')]
            timings = all_timings.get(method, {})
            for key in all_keys:
                row.append(timings.get(key, ''))
            writer.writerow(row)
    print(f"Saved summary table to {summary_path}")


def main():
    ref_dims, ref_spacing, ref_flat = read_vtk_scalar_field(REFERENCE_PATH)
    ref_data = reshape_field(ref_flat, ref_dims)
    
    errors = {}
    total_times = {}
    all_timings = {}
    
    for method in METHODS:
        print(f"\nProcessing method: {method.upper()}")
        method_frame_path = os.path.join(BASE_DIR, method, 'frames', METHODS_FILE)
        coarse_dims, coarse_spacing, coarse_flat = read_vtk_scalar_field(method_frame_path)
        coarse_data = reshape_field(coarse_flat, coarse_dims)
        
        sampled_ref = sample_reference_to_coarse(ref_data, ref_dims, ref_spacing, coarse_dims, coarse_spacing)
        error = compute_l2_error(coarse_data, sampled_ref)
        errors[method] = error
        print(f"{method.upper()} relative L2 error: {error:.6e}")
        
        sim_info_path = os.path.join(BASE_DIR, method, 'simulation_infos.txt')
        if os.path.exists(sim_info_path):
            timings = parse_simulation_info(sim_info_path)
            all_timings[method] = timings
            total_times[method] = timings.get('SIMULATION TOTAL EXECUTION TIME', 0.0)
        else:
            print(f"Warning: Missing {sim_info_path}")
    
    # Save error CSV
    error_csv_path = os.path.join(BASE_DIR, 'l2_errors.csv')
    with open(error_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Method', 'Relative_L2_Error'])
        for method, error in errors.items():
            writer.writerow([method.upper(), error])
    print(f"Saved error CSV to {error_csv_path}")
    
    # Save error plot
    plt.figure(figsize=(8, 5))
    method_names = [m.upper() for m in errors.keys()]
    error_values = [errors[m] for m in errors.keys()]
    plt.bar(method_names, error_values)
    plt.ylabel('Relative L2 Error')
    plt.title('Error Analysis of Numerical Methods')
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    error_plot_path = os.path.join(BASE_DIR, 'l2_error_plot.png')
    plt.savefig(error_plot_path)
    plt.close()
    print(f"Saved error plot to {error_plot_path}")
    
    # Plot execution times
    plot_total_execution_times(total_times)
    
    # Save summary CSV
    save_summary_csv(errors, all_timings)


if __name__ == '__main__':
    main()
