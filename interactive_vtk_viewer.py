import os
import glob
import numpy as np
import plotly.express as px
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def read_vtk_ascii(filepath):
    """Reads a single VTK file and returns the Vm data as a 2D NumPy array."""
    reader = vtk.vtkStructuredPointsReader()
    reader.SetFileName(filepath)
    reader.Update()

    image_data = reader.GetOutput()
    dims = image_data.GetDimensions()
    scalars = vtk_to_numpy(image_data.GetPointData().GetScalars())
    
    Nx, Ny = dims[0], dims[1]
    return scalars.reshape((Ny, Nx))  # VTK uses Fortran order

def load_vtk_stack(frames_dir):
    """Loads all VTK frames from the specified directory and returns a 3D NumPy array."""
    vtk_files = sorted(glob.glob(os.path.join(frames_dir, "Vm_*.vtk")))
    
    if not vtk_files:
        raise FileNotFoundError(f"No VTK files found in {frames_dir}")

    print(f"[INFO] Found {len(vtk_files)} VTK frames. Loading data...")
    data_stack = np.array([read_vtk_ascii(f) for f in vtk_files])
    print(f"[INFO] Data loaded successfully. Stack shape: {data_stack.shape}")
    return data_stack

def generate_plotly_animation(data_stack, output_path):
    """Generates and saves an interactive HTML animation from the data stack."""
    n_frames = data_stack.shape[0]
    print(f"[INFO] Generating interactive HTML animation with {n_frames} frames...")

    vmin = data_stack.min()
    vmax = data_stack.max()
    print(f"[INFO] Fixed color scale: vmin = {vmin:.3f}, vmax = {vmax:.3f}")

    fig = px.imshow(
        data_stack,
        animation_frame=0,
        labels={"animation_frame": "Frame"},
        color_continuous_scale="magma",
        origin="lower",
        aspect="equal",
        zmin=vmin,
        zmax=vmax
    )

    fig.update_layout(
        title="Vm Animation",
        coloraxis_showscale=True,
        margin=dict(l=0, r=0, t=40, b=0)
    )

    fig.write_html(output_path)
    print(f"[SUCCESS] HTML animation saved to: {output_path}")

def main(output_dir):
    frames_dir = os.path.join(output_dir, "frames")
    if not os.path.exists(frames_dir):
        raise FileNotFoundError(f"The directory '{frames_dir}' does not exist.")

    data_stack = load_vtk_stack(frames_dir)

    output_html = os.path.join(output_dir, "interactive_viewer.html")
    generate_plotly_animation(data_stack, output_html)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate an interactive HTML animation from VTK frames.")
    parser.add_argument("output_dir", type=str, help="Directory containing a 'frames/' subdirectory with VTK files.")
    args = parser.parse_args()

    try:
        main(args.output_dir)
    except Exception as e:
        print(f"[ERROR] {e}")

# Example usage:
# python3 interactive_vtk_viewer.py /path/to/simulation/output