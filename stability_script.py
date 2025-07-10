import os
import subprocess
import configparser
import shutil
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
METHODS = ['fe', 'osadi', 'ssiadi', 'do', 'hv', 'cg']
INI_TEMPLATE = 'configs/stability/template.ini'
WORKING_INI = 'configs/stability/current_test.ini'
SIM_EXECUTABLE = './bin/monodomain_simulation'
CHECK_FRAME = 'lastframe.vtk'
OUTPUT_BASE_DIR = 'outputs/stability_test/fiber60'
TOTAL_TIME = 500.0  # ms
INITIAL_DT = 0.01

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


if __name__ == '__main__':
    main()
