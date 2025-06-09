import os
import glob
import numpy as np
import imageio
import matplotlib.pyplot as plt

import vtk
from vtk.util.numpy_support import vtk_to_numpy


def read_vtk_ascii(filepath):
    """
    Read scalar field from a VTK file generated in ASCII format.

    Parameters:
        filepath (str): Path to the .vtk file.

    Returns:
        tuple: (2D NumPy array of scalar values, spacing (dx, dy), dimensions (Nx, Ny))
    """
    reader = vtk.vtkStructuredPointsReader()
    reader.SetFileName(filepath)
    reader.Update()

    image_data = reader.GetOutput()
    dims = image_data.GetDimensions()
    spacing = image_data.GetSpacing()

    point_data = image_data.GetPointData()
    if point_data.GetScalars() is None:
        raise ValueError(f"No scalar data found in file {filepath}")

    scalars = vtk_to_numpy(point_data.GetScalars())
    Nx, Ny = dims[0], dims[1]
    Vm = scalars.reshape((Ny, Nx))
    
    return Vm, spacing[:2], (Nx, Ny)

def create_gif_from_vtk_frames(output_dir, gif_name="simulation.gif", cmap='plasma', fps=20):
    """
    Read all VTK frames from the 'frames' subdirectory and generate a GIF.

    Parameters:
        output_dir (str): Path to the simulation output directory (containing 'frames/').
        gif_name (str): Filename of the output GIF.
        cmap (str): Colormap to use for the image.
        fps (int): Frames per second for the GIF.
    """
    frames_dir = os.path.join(output_dir, "frames")
    vtk_files = sorted(glob.glob(os.path.join(frames_dir, "Vm_*.vtk")))

    if not vtk_files:
        raise FileNotFoundError(f"No VTK files found in {frames_dir}")

    images = []
    print(f"Reading {len(vtk_files)} frames...")

    # Optional: extract global min/max for consistent color scaling
    global_min, global_max = float("inf"), float("-inf")
    data_list = []

    for file in vtk_files:
        Vm, _, _ = read_vtk_ascii(file)
        data_list.append(Vm)
        global_min = min(global_min, Vm.min())
        global_max = max(global_max, Vm.max())

    for Vm in data_list:
        fig, ax = plt.subplots(figsize=(5, 5))
        im = ax.imshow(Vm, cmap=cmap, vmin=global_min, vmax=global_max, origin='lower')
        ax.axis('off')
        fig.colorbar(im, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
        fig.tight_layout(pad=0.5)

        # Save frame to temporary image in memory
        fig.canvas.draw()
        image = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
        image = image.reshape(fig.canvas.get_width_height()[::-1] + (4,))
        images.append(image)

        plt.close(fig)

    gif_path = os.path.join(output_dir, gif_name)
    imageio.mimsave(gif_path, images, fps=fps, loop=0)
    print(f"GIF saved to {gif_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate GIF from VTK frames.")
    parser.add_argument("output_dir", type=str, help="Path to the simulation output directory.")
    parser.add_argument("--gif_name", type=str, default="simulation.gif", help="Name of the output GIF file.")
    parser.add_argument("--cmap", type=str, default="plasma", help="Colormap used for plotting.")
    parser.add_argument("--fps", type=int, default=20, help="Frames per second for the GIF.")

    args = parser.parse_args()
    create_gif_from_vtk_frames(args.output_dir, args.gif_name, args.cmap, args.fps)

# Example usage:
# python generate_vtk_gif.py /path/to/simulation/output --gif_name my_simulation.gif --cmap plasma --fps 20
# Make sure to have the 'frames' directory with VTK files in the specified output directory.
# The script will read all VTK files, generate images, and save them as a GIF.