#include "../include/save_functions.h"

// Map of save functions
static const struct
{
    const char *name;
    save_function_t function;
    const char *file_extension;
} save_function_map[] = {
    {"save_as_vtk", save_as_vtk, "vtk"},
    {"save_as_vtu", save_as_vtu, "vtu"},
    {"save_as_txt", save_as_txt, "txt"},
    {NULL, NULL, NULL}};

// Function to get the appropriate save function based on the name
save_function_t get_save_function(const char *name)
{
    if (name == NULL)
    {
        printf("Error: save function name is NULL\n");
        return NULL;
    }

    for (int i = 0; save_function_map[i].name != NULL; i++)
        if (strcmp(save_function_map[i].name, name) == 0)
            return save_function_map[i].function;
    return NULL;
}

// Function to get the appropriate file extension based on the save function name
const char *get_file_extension(const char *name)
{
    if (name == NULL)
    {
        printf("Error: save function name is NULL\n");
        return NULL;
    }

    for (int i = 0; save_function_map[i].name != NULL; i++)
        if (strcmp(save_function_map[i].name, name) == 0)
            return save_function_map[i].file_extension;
    return NULL;
}

// Implementation of save functions
inline void save_as_vtk(const char *file_path, const real *data, const int Nx, const int Ny, const real delta_x, const real delta_y)
{
    FILE *file = fopen(file_path, "w");
    if (file == NULL)
    {
        printf("Error opening file %s\n", file_path);
        exit(1);
    }

    // Write the VTK file - Header for 2D data
    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Monodomain Data 2D\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %d %d 1\n", Nx, Ny);
    fprintf(file, "ORIGIN 0 0 0\n");
    fprintf(file, "SPACING %.6g %.6g 1.0\n", delta_x, delta_y);
    fprintf(file, "POINT_DATA %d\n", Nx * Ny);
    fprintf(file, "SCALARS Vm %s 1\n", REAL_TYPE);
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int index = 0; index < Nx * Ny; index++)
        fprintf(file, "%e\n", data[index]);

    fclose(file);
}

// Writes a scalar field as a VTK Unstructured Grid (.vtu) for 2D finite difference data
inline void save_as_vtu(const char *file_path, const real *data, const int Nx, const int Ny, const real dx, const real dy)
{
    FILE *file = fopen(file_path, "w");
    if (file == NULL)
    {
        printf("Error opening file %s\n", file_path);
        exit(1);
    }

    const int num_points = Nx * Ny;
    const int num_cells = (Nx - 1) * (Ny - 1);

    // VTK header
    fprintf(file,
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            "  <UnstructuredGrid>\n"
            "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
            num_points, num_cells);

    // Write point coordinates
    fprintf(file, "      <Points>\n");
    fprintf(file, "        <DataArray type=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", REAL_TYPE_NAME);
    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            fprintf(file, " %.6f %.6f 0.0\n", i * dx, j * dy);
        }
    }
    fprintf(file, "        </DataArray>\n");
    fprintf(file, "      </Points>\n");

    // Write cell definitions (quads)
    fprintf(file, "      <Cells>\n");

    // Connectivity: list of vertex indices for each cell
    fprintf(file, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int j = 0; j < Ny - 1; ++j)
    {
        for (int i = 0; i < Nx - 1; ++i)
        {
            int p0 = j * Nx + i;
            int p1 = p0 + 1;
            int p2 = p1 + Nx;
            int p3 = p0 + Nx;
            fprintf(file, "%d %d %d %d\n", p0, p1, p2, p3);
        }
    }
    fprintf(file, "        </DataArray>\n");

    // Offsets: cumulative count of points per cell (always 4 for quads)
    fprintf(file, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for (int i = 1; i <= num_cells; ++i)
    {
        fprintf(file, "%d\n", i * 4);
    }
    fprintf(file, "        </DataArray>\n");

    // Types: VTK cell type ID (9 = quad)
    fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int i = 0; i < num_cells; ++i)
    {
        fprintf(file, "9\n");
    }
    fprintf(file, "        </DataArray>\n");
    fprintf(file, "      </Cells>\n");

    // Scalar data associated with points (Vm field)
    fprintf(file, "      <PointData Scalars=\"Vm\">\n");
    fprintf(file, "        <DataArray type=\"%s\" Name=\"Vm\" format=\"ascii\">\n", REAL_TYPE_NAME);
    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file, "%e\n", data[i]);
    }
    fprintf(file, "        </DataArray>\n");
    fprintf(file, "      </PointData>\n");

    // Footer
    fprintf(file,
            "    </Piece>\n"
            "  </UnstructuredGrid>\n"
            "</VTKFile>\n");

    fclose(file);
}

inline void save_as_txt(const char *file_path, const real *data, const int Nx, const int Ny, const real delta_x, const real delta_y)
{
    FILE *file = fopen(file_path, "w");
    if (file == NULL)
    {
        printf("Error opening file %s\n", file_path);
        exit(1);
    }

    for (int index = 0; index < Nx * Ny; index++)
        fprintf(file, "%e\n", data[index]);

    fclose(file);
}

// Writes a .pvd master file for ParaView time series
// Inputs:
//   filename     → output .pvd file name (e.g., "series.pvd")
//   num_frames   → number of time steps
//   dt           → time step size (in simulation units)
//   prefix       → filename prefix used in .vtu files (e.g., "frame_")
// Output:
//   Creates a .pvd file that references frame_0000.vtu, frame_0001.vtu, ..., with correct timesteps
inline void write_pvd_file(const char *file_path, int num_frames, real dt, const char *prefix)
{
    FILE *file = fopen(file_path, "w");
    if (file == NULL)
    {
        printf("Error opening file %s\n", file_path);
        exit(1);
    }

    fprintf(file,
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            "  <Collection>\n");

    for (int i = 0; i < num_frames; ++i)
    {
        fprintf(file,
                "    <DataSet timestep=\"%.6f\" group=\"\" part=\"0\" file=\"%s%04d.vtu\"/>\n",
                i * dt, prefix, i);
    }

    fprintf(file,
            "  </Collection>\n"
            "</VTKFile>\n");

    fclose(file);
}

inline void save_simulation_state(const char *file_path, const real *Vm, const real *sV, const int Nx, const int Ny, const int n_sv)
{
    FILE *file = fopen(file_path, "wb");
    if (file == NULL)
    {
        printf("Error opening save state file for writing: %s\n", file_path);
        exit(1);
    }

    // Write metadata
    fwrite(&Nx, sizeof(int), 1, file);
    fwrite(&Ny, sizeof(int), 1, file);
    fwrite(&n_sv, sizeof(int), 1, file);

    size_t total_points = Nx * Ny;

    // Write Vm data
    if (fwrite(Vm, sizeof(real), total_points, file) != total_points)
    {
        printf("Error writing Vm data.\n");
        exit(1);
    }

    // Write state variables data
    if (fwrite(sV, sizeof(real), total_points * n_sv, file) != total_points * n_sv)
    {
        printf("Error writing sV data.\n");
        exit(1);
    }

    fclose(file);
}


inline void restore_simulation_state(const char *file_path, real *Vm, real *sV, const int Nx, const int Ny, const int n_sv)
{
    FILE *file = fopen(file_path, "rb");
    if (file == NULL)
    {
        printf("Error opening restore file for reading: %s\n", file_path);
        exit(1);
    }

    // Read and validate metadata
    int file_Nx, file_Ny, file_n_sv;
    fread(&file_Nx, sizeof(int), 1, file);
    fread(&file_Ny, sizeof(int), 1, file);
    fread(&file_n_sv, sizeof(int), 1, file);

    if (file_Nx != Nx || file_Ny != Ny || file_n_sv != n_sv)
    {
        printf("Error: grid or state size mismatch in restored file.\n");
        exit(1);
    }

    size_t total_points = Nx * Ny;

    printf("Restoring simulation state from: %s\n", file_path);

    // Read Vm
    if (fread(Vm, sizeof(real), total_points, file) != total_points)
    {
        printf("Error reading Vm data.\n");
        exit(1);
    }

    // Read state variables
    if (fread(sV, sizeof(real), total_points * n_sv, file) != total_points * n_sv)
    {
        printf("Error reading sV data.\n");
        exit(1);
    }

    fclose(file);
}