import os
import subprocess
import configparser
import shutil
import matplotlib.pyplot as plt
import numpy as np
import csv

# === CONFIGURATION FOR RUNNING ===
METHODS = ['fe', 'osadi', 'ssiadi', 'do', 'cs', 'dg', 'dy', 'hv', 'cg']
INI_TEMPLATE = 'configs/stability/template.ini'
WORKING_INI = 'configs/stability/current_test.ini'
SIM_EXECUTABLE = './bin/monodomain_simulation'
CHECK_FRAME = 'lastframe.vtk'
OUTPUT_BASE_DIR = 'outputs/stability_test/fiber45'
TOTAL_TIME = 500.0  # ms
INITIAL_DT = 0.01

# === CONFIGURATION FOR ANALYSIS ===
REFERENCE_FILE = 'outputs/error_analysis_aniso_rotation/reference_45/frames/Vm_500000.vtk'

ERROR_CSV = os.path.join(OUTPUT_BASE_DIR, 'stability_errors.csv')
TIME_CSV = os.path.join(OUTPUT_BASE_DIR, 'stability_timings.csv')
ERROR_PLOT = os.path.join(OUTPUT_BASE_DIR, 'stability_error_plot.png')
TIME_PLOT = os.path.join(OUTPUT_BASE_DIR, 'stability_time_plot.png')

def generate_dt_values(total_time, start_dt):
    dt = start_dt
    values = []
    while dt <= 1.0:
        remainder = total_time % dt
        if abs(remainder) < 1e-10 or abs(remainder - dt) < 1e-10:
            values.append(dt)
        dt = round(dt + 0.0001, 4)
    return values

DT_VALUES = generate_dt_values(TOTAL_TIME, INITIAL_DT)
print(f"DT_VALUES = {DT_VALUES}\n")

def has_nan_in_vtk(filepath):
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if 'NaN' in line or 'nan' in line:
                    return True
        return False
    except Exception as e:
        return True


def prepare_ini(template_path, output_path, method, dt, output_dir):
    config = configparser.ConfigParser()
    config.read(template_path)

    config['simulation']['numerical_method'] = method
    config['simulation']['dt'] = str(dt)
    config['simulation']['output_dir'] = output_dir

    with open(output_path, 'w') as configfile:
        config.write(configfile)


def run_simulation(ini_path):
    result = subprocess.run([SIM_EXECUTABLE, ini_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.returncode


def log_method_result(method, dt, status, log_dir):
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"{method}.log")
    with open(log_path, 'a') as f:
        f.write(f"dt = {dt:.5f}: {'STABLE ✓' if status else 'UNSTABLE ✗'}\n")


def test_method_stability(method):
    print(f"\n[INFO] Testing method: {method.upper()}")
    max_dt = None
    stability_map = {}

    for dt in DT_VALUES:
        print(f"  Trying dt = {dt:.5f}... ", end='')
        method_dir = os.path.join(OUTPUT_BASE_DIR, method, f'dt_{dt:.5f}')
        frames_dir = os.path.join(method_dir, 'frames')
        vtk_file = os.path.join(frames_dir, CHECK_FRAME)

        shutil.rmtree(method_dir, ignore_errors=True)
        prepare_ini(INI_TEMPLATE, WORKING_INI, method, dt, method_dir)

        returncode = run_simulation(WORKING_INI)

        if returncode != 0:
            print("✗ Failed (return code)")
            log_method_result(method, dt, False, os.path.join(OUTPUT_BASE_DIR, 'logs'))
            stability_map[dt] = False
            break

        if not os.path.isfile(vtk_file):
            print("✗ Failed (no output)")
            log_method_result(method, dt, False, os.path.join(OUTPUT_BASE_DIR, 'logs'))
            stability_map[dt] = False
            break

        if has_nan_in_vtk(vtk_file):
            print("✗ Failed (NaN)")
            log_method_result(method, dt, False, os.path.join(OUTPUT_BASE_DIR, 'logs'))
            stability_map[dt] = False
            break

        print("✓ Stable")
        log_method_result(method, dt, True, os.path.join(OUTPUT_BASE_DIR, 'logs'))
        stability_map[dt] = True
        max_dt = dt

    return max_dt, stability_map


def plot_stability_map(all_results):
    fig, ax = plt.subplots(figsize=(10, 6))
    for method, result in all_results.items():
        dt_vals = list(result.keys())
        statuses = ['✓' if result[dt] else '✗' for dt in dt_vals]
        y = [method.upper()] * len(dt_vals)
        colors = ['green' if s == '✓' else 'red' for s in statuses]
        ax.scatter(dt_vals, y, color=colors, label=method.upper(), s=100, marker='o')

    ax.set_xlabel("dt (ms)")
    ax.set_title("Stability Regions by Method")
    ax.set_xscale("log")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    output_path = os.path.join(OUTPUT_BASE_DIR, 'stability_region.png')
    plt.savefig(output_path)
    print(f"[INFO] Stability region plot saved to: {output_path}")


# === UTILITIES FOR VTK HANDLING ===
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
    return data.reshape((ny, nx))  # VTK stores as (y, x)


def sample_reference_to_coarse(ref_data, ref_dims, ref_spacing, coarse_dims, coarse_spacing):
    """Sample fine reference solution to match coarse grid points."""
    factor_x = int(round(coarse_spacing[0] / ref_spacing[0]))
    factor_y = int(round(coarse_spacing[1] / ref_spacing[1]))
    sampled_ref = ref_data[::factor_y, ::factor_x]
    expected_shape = (coarse_dims[1], coarse_dims[0])  # (ny, nx)
    assert sampled_ref.shape == expected_shape, f"Shape mismatch: {sampled_ref.shape} vs {expected_shape}"
    return sampled_ref


# === NUMERICAL TOOLS ===
def compute_l2_error(u_numeric, u_ref):
    """Compute relative L2 norm error."""
    diff = u_numeric - u_ref
    N = diff.size
    l2_error = np.linalg.norm(diff) / np.sqrt(N)
    return l2_error


def parse_simulation_info(filepath):
    """Extract timing values from simulation_infos.txt using '=' separator."""
    timings = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if '=' in line:
                key, val = line.split('=', 1)
                key = key.strip()
                val = val.strip().split()[0]  # pega só o número antes da unidade 's'
                try:
                    timings[key] = float(val)
                except ValueError:
                    continue
    return timings


# === ANALYSIS FUNCTIONS ===
def collect_all_simulations(base_dir):
    """Return a dictionary mapping methods to their dt-simulation folders."""
    simulations = {}
    for method in os.listdir(base_dir):
        method_path = os.path.join(base_dir, method)
        if os.path.isdir(method_path):
            dt_dirs = [
                os.path.join(method_path, d)
                for d in os.listdir(method_path)
                if os.path.isdir(os.path.join(method_path, d)) and d.startswith('dt_')
            ]
            simulations[method] = dt_dirs
    return simulations


def analyze_simulations():
    print("[INFO] Reading reference solution...")
    ref_dims, ref_spacing, ref_flat = read_vtk_scalar_field(REFERENCE_FILE)
    ref_2d = reshape_field(ref_flat, ref_dims)

    simulations = collect_all_simulations(OUTPUT_BASE_DIR)
    error_results = []
    time_results = []

    print("[INFO] Analyzing simulations...")
    for method, dt_dirs in simulations.items():
        for dt_path in dt_dirs:
            dt_str = os.path.basename(dt_path).replace('dt_', '')
            try:
                dt = float(dt_str)
            except ValueError:
                continue

            vtk_path = os.path.join(dt_path, 'frames', CHECK_FRAME)
            info_path = os.path.join(dt_path, 'simulation_infos.txt')

            if not os.path.isfile(vtk_path):
                print(f"[WARN] Missing VTK for {method} at dt={dt}")
                continue

            try:
                dims, spacing, flat = read_vtk_scalar_field(vtk_path)
                data_2d = reshape_field(flat, dims)

                # Skip fields with NaNs or Infs
                if not np.isfinite(data_2d).all():
                    print(f"[WARN] Invalid values (NaN/Inf) in {vtk_path} for {method} at dt={dt}")
                    continue

                sampled_ref = sample_reference_to_coarse(ref_2d, ref_dims, ref_spacing, dims, spacing)
                error = compute_l2_error(data_2d, sampled_ref)

            except Exception as e:
                print(f"[ERROR] Could not process {vtk_path}: {e}")
                continue

            time = None
            if os.path.isfile(info_path):
                timings = parse_simulation_info(info_path)
                time = timings.get('SIMULATION TOTAL EXECUTION TIME', None)

            # Discard invalid time as well
            if time is not None and (not np.isfinite(time) or time <= 0):
                time = None

            error_results.append((method, dt, error))
            time_results.append((method, dt, time))

    return error_results, time_results


# === EXPORTING AND PLOTTING ===
def save_csv(data, csv_path, headers):
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(data)
    print(f"[INFO] Saved CSV to {csv_path}")


def plot_error_graph(errors, out_path):
    plt.figure(figsize=(10, 6))
    methods = sorted(set(row[0] for row in errors))

    for method in methods:
        points = sorted((dt, err) for m, dt, err in errors if m == method)
        dts = [dt for dt, err in points if err is not None and err > 0]
        errs = [err for dt, err in points if err is not None and err > 0]

        failed_dts = [dt for dt, err in points if err is None or err <= 0]

        plt.plot(dts, errs, marker='o', label=method.upper())

        # Marcar os pontos que falharam
        if failed_dts:
            plt.plot(failed_dts, [1e0]*len(failed_dts), 'x', color='red', label=f"{method.upper()} FAIL")

    plt.xlabel('dt (ms)')
    plt.ylabel('Relative L2 Error')
    plt.title('Stability Error vs dt')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"[INFO] Saved error plot to {out_path}")


def plot_time_graph(times, out_path):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    methods = sorted(set(row[0] for row in times))

    for method in methods:
        points = sorted((dt, t) for m, dt, t in times if m == method and t is not None and np.isfinite(t) and t > 0)
        if not points:
            continue

        dts = [dt for dt, _ in points]
        times_vals = [t for _, t in points]

        plt.plot(dts, times_vals, marker='o', label=method.upper())

    plt.xlabel('dt (ms)')
    plt.ylabel('Execution Time (s)')
    plt.title('Execution Time vs dt')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"[INFO] Saved timing plot to {out_path}")

def main():
    print("[INFO] Starting numerical stability test...\n")
    os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)
    all_results = {}
    summary = {}

    for method in METHODS:
        max_stable_dt, stability_map = test_method_stability(method)
        summary[method] = max_stable_dt
        all_results[method] = stability_map

    # Print summary
    print("\n[SUMMARY]")
    for method, dt in summary.items():
        if dt:
            print(f"{method.upper()}: maximum stable dt = {dt:.5f} ms")
        else:
            print(f"{method.upper()}: unstable even at dt = {INITIAL_DT}")

    # Plot stability regions
    plot_stability_map(all_results)
    
    print("\n[INFO] Finished running all methods!\n")
    
    print("\n[INFO] Starting stability error and timing analysis...")
    errors, times = analyze_simulations()

    save_csv(errors, ERROR_CSV, ['Method', 'dt', 'L2 Error'])
    save_csv(times, TIME_CSV, ['Method', 'dt', 'Execution Time'])

    plot_error_graph(errors, ERROR_PLOT)
    plot_time_graph(times, TIME_PLOT)
    print("[INFO] Analysis complete.")


if __name__ == '__main__':
    main()
