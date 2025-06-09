import os
import numpy as np
import matplotlib.pyplot as plt
import csv

BASE_DIR = 'outputs/error_analysis'
METHODS = ['fe', 'osadi', 'ssiadi', 'do', 'hv']
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

    # Extract matching points
    sampled_ref = ref_data[::factor_y, ::factor_x]
    expected_shape = (coarse_dims[1], coarse_dims[0])  # (ny, nx)
    assert sampled_ref.shape == expected_shape, f"Shape mismatch: {sampled_ref.shape} vs {expected_shape}"
    
    return sampled_ref

def compute_l2_error(u_numeric, u_ref):
    """Compute relative L2 norm error."""
    diff = u_numeric - u_ref
    N = diff.size
    # l2_error = np.linalg.norm(diff) / np.linalg.norm(u_ref)
    l2_error = np.linalg.norm(diff) / np.sqrt(N)  # Alternative normalization
    return l2_error

def main():
    ref_dims, ref_spacing, ref_flat = read_vtk_scalar_field(REFERENCE_PATH)
    ref_data = reshape_field(ref_flat, ref_dims)
    
    errors = {}
    
    for method in METHODS:
        method_path = os.path.join(BASE_DIR, method, 'frames', METHODS_FILE)
        coarse_dims, coarse_spacing, coarse_flat = read_vtk_scalar_field(method_path)
        coarse_data = reshape_field(coarse_flat, coarse_dims)
        
        sampled_ref = sample_reference_to_coarse(ref_data, ref_dims, ref_spacing, coarse_dims, coarse_spacing)
        
        error = compute_l2_error(coarse_data, sampled_ref)
        errors[method] = error
        print(f"{method.upper()} relative L2 error: {error:.6e}")
    
    # Save to CSV
    csv_path = os.path.join(BASE_DIR, 'l2_errors.csv')
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Method', 'Relative_L2_Error'])
        for method, error in errors.items():
            writer.writerow([method.upper(), error])
    
    # Plot
    plt.figure(figsize=(8, 5))
    methods = list(errors.keys())
    error_values = [errors[m] for m in methods]
    methods = [method.upper() for method in list(errors.keys())]
    plt.bar(methods, error_values)
    plt.ylabel('Relative L2 Error')
    plt.title('Error Analysis of Numerical Methods')
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    plot_path = os.path.join(BASE_DIR, 'l2_error_plot.png')
    plt.savefig(plot_path)
    plt.close()
    print(f"Saved CSV to {csv_path}")
    print(f"Saved plot to {plot_path}")

if __name__ == '__main__':
    main()
