#pragma once

#include "binary_decoder.h"


typedef struct histo {
  histo();
  histo(std::string, std::string, size_t, aladyn_float, aladyn_float);
  std::string basename;
  std::string data_to_bin;
  size_t nbin;
  aladyn_float min_value, max_value;
  aladyn_float get_bin_size();
  size_t get_column_to_bin();
  bool enabled;
  aladyn_float *data = new aladyn_float[nbin + 3];
  void write_binned_data();
  std::string get_filename_out();
} histo;

typedef struct densityplot {
  densityplot();
  densityplot(std::string, std::string, std::string, size_t, size_t, aladyn_float, aladyn_float, aladyn_float, aladyn_float);
  std::string basename;
  std::string data_to_bin_x, data_to_bin_y;
  size_t nbin_x, nbin_y;
  aladyn_float min_value_x, max_value_x;
  aladyn_float min_value_y, max_value_y;
  aladyn_float get_bin_size_x();
  aladyn_float get_bin_size_y();
  size_t get_x_column_to_bin();
  size_t get_y_column_to_bin();
  bool enabled;
  aladyn_float **data = new aladyn_float*[nbin_x + 3];
  void write_binned_data();
  std::string get_filename_out();
} densityplot;

typedef struct Parameters {
  Parameters();
  int argc;
  std::string * argv;
  std::string filebasename;
  size_t file_version;
  bool we_dont_know_file_version;
  bool sim_is_2d;
  bool we_dont_know_if_sim_is_2d;
  bool file_has_weight;
  bool we_dont_know_if_file_has_weight;
  bool file_has_charge;
  bool we_dont_know_if_file_has_charge;
  bool we_have_to_do_swap;
  bool we_dont_know_if_we_have_to_do_swap;
  bool out_ppg, out_json, out_csv, out_xyze, out_cutx, out_cuty, out_cutz, out_grid2d, out_clean_bin, out_lineoutx, out_vtk, out_vtk_nostretch;
  unsigned int ncpu_x, ncpu_y, ncpu_z, ncpu;
  unsigned int ndv, resampling_factor;
  size_t nparts;
  size_t npx, npy, npz, npx_per_cpu, npy_per_cpu, npz_per_cpu;
  size_t npx_resampled, npy_resampled, npz_resampled, npx_resampled_per_cpu, npy_resampled_per_cpu, npz_resampled_per_cpu;
  size_t header_size_bytes;
  long long int nptot;
  int endianness;
  bool fixed_aladyn_version;
  bool multifile;
  bool stretched_grid;  //nb: the tools assumes that, if a grid if found anywhere (inside the .dat file or at the end of the .bin file), then it is stretched, even if it's not
  int stretched_along_x;
  aladyn_float mass_MeV;
  int endian_file, endian_machine;
  aladyn_float tnow, xmin, xmax, ymin, ymax, zmin, zmax;
  std::vector<aladyn_float> where_to_cut_grid_along_x, where_to_cut_grid_along_y, where_to_cut_grid_along_z;
  std::vector<aladyn_float> xcoord, ycoord, zcoord;
  std::vector<float> realpar;
  std::vector<int> intpar;
  std::vector<std::string> phasespace_file_labels;
  std::vector<std::string> grid_file_labels;
  std::vector<histo> histograms;
  std::vector<densityplot> densityplots;
  bool overwrite_weight, overwrite_charge;
  aladyn_float overwrite_weight_value, overwrite_charge_value;
  int subsample;
  int span;
  bool grid_file, phasespace_file;
  void check_forced_version(const int, const char **);
  void parse_command_line();
  void debug_read_parameters();
  void check_swap();
  void check_filename(const char *);
  void read_params_from_dat_file(std::ifstream &);
  void read_params_from_bin_file(const char *);
  void man(const char[]);
  void parse_json();
} Parameters;

