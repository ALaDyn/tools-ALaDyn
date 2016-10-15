#ifndef __LEGGI_ALADYN_FORTRAN
#define __LEGGI_ALADYN_FORTRAN

//#define ENABLE_DEBUG

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#define MAJOR_RELEASE  7
#define MINOR_RELEASE  0
#define BUGFIX_RELEASE 2

typedef float aladyn_float;

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <limits>
#include <cfloat>
#include <ios>
#include <cstdarg>

#if (defined CINECA)
#include <inttypes.h>
#include <stdint.h>
#endif

#if (!defined CINECA) && (defined _MSC_VER)
#include<cstdint>
#endif

#if (!defined CINECA) && (defined __GNUC__)
/* Test for GCC > 4.6.0 */
#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ > 6))
#include<cstdint>
#else
#include <inttypes.h>
#include <stdint.h>
#endif
#endif

#if defined (_MSC_VER)
#include<wchar.h>
#define fseeko _fseeki64
#define ftello _ftelli64
#endif

#if defined (__MINGW32__)
#define fseeko fseeko64
#define ftello ftello64
#endif

#if defined (__CYGWIN__)
#define fseeko fseek
#define ftello ftell
#endif


#define MAX_NUM_OF_PARTICLES_IN_MEMORY 10000000            // in reality we store double this number -1
#define UMA_G                          1.660538921E-24     // from uma to grams
#define C                              2.99792458E+10      // cm / s
#define ME_G                           9.10938291E-28      // electron mass [g]
#define MP_G                           1.6726231E-24       // proton mass [g]
#define MP_MEV                         938.272013          // proton mass [MeV/c^2]
#define ME_MEV                         0.510998928         // electron mass [MeV/c^2]
#define MHI_UMA                        26.981538           // atomic weight of Aluminum in atomic mass units
#define MLI_UMA                        12.0107             // atomic weight of Carbon in atomic mass units
//#define CHARGE                       4.80320425e-10      // statC  - official value
#define CHARGE                         4.803262e-10        // statC  - Turchetti's value
#define MAX_NUMBER_OF_CPUS             32768
#define NUMBER_OF_PARAMS_IN_DAT_FILE   20
#define NUMBER_OF_BIN_BY_DEFAULT       120
#define FORCE_PRINTF_BUFFER_SIZE       1024
#define INVALID_COLUMN                -1
#define COLUMN_X                       0
#define COLUMN_Y                       1
#define COLUMN_Z                       2
#define COLUMN_PX                      3
#define COLUMN_PY                      4
#define COLUMN_PZ                      5
#define COLUMN_GAMMA                   6
#define COLUMN_THETA                   7
#define COLUMN_E                       8
#define COLUMN_THETAT                  9
#define COLUMN_TY                     10
#define COLUMN_TZ                     11
#define COLUMN_W                      12
#define COLUMN_CH                     13



struct Parameters
{
  Parameters();
  int argc;
  std::string * argv;
  std::string filebasename;
  int file_version;
  bool we_dont_know_file_version;
  bool sim_is_2d;
  bool we_dont_know_if_sim_is_2d;
  bool we_have_to_find_minmax;
  bool we_dont_know_if_we_have_to_find_minmax;
  bool we_have_to_do_binning;
  bool we_dont_know_if_we_have_to_do_binning;
  bool file_has_weight;
  bool we_dont_know_if_file_has_weight;
  bool file_has_charge;
  bool we_dont_know_if_file_has_charge;
  bool we_have_to_do_swap;
  bool we_dont_know_if_we_have_to_do_swap;
  bool out_vtk, out_ppg, out_params, out_csv, out_xyze, out_cutx, out_cuty, out_cutz, out_grid2d, out_clean_bin, out_lineoutx, out_vtk_nostretch;
  bool minmax_found;
  unsigned int ncpu_x, ncpu_y, ncpu_z, ncpu;
  unsigned int ndv, resampling_factor;
  size_t nparts;
  size_t npx, npy, npz, npx_per_cpu, npy_per_cpu, npz_per_cpu;
  size_t npx_resampled, npy_resampled, npz_resampled, npx_resampled_per_cpu, npy_resampled_per_cpu, npz_resampled_per_cpu;
  size_t header_size_bytes;
  size_t nptot;
  int endianness;
  bool fixed_aladyn_version;
  bool multifile;
  bool stretched_grid;  //nb: the tools assumes that, if a grid if found anywhere (inside the .dat file or at the end of the .bin file), then it is stretched, even if it's not
  int stretched_along_x;
  aladyn_float mass_MeV;
  int endian_file, endian_machine;
  int last_cpu;
  int nparams;
  std::string support_label;
  aladyn_float tnow, xmin, xmax, pxmin, pxmax, ymin, ymax, pymin, pymax, zmin, zmax, pzmin, pzmax, wmin, wmax, chmin, chmax, Emin, Emax, gammamin, gammamax, thetamin, thetamax, thetaTmin, thetaTmax, tymin, tymax, tzmin, tzmax;
  int nbin, nbin_x, nbin_y, nbin_z, nbin_px, nbin_py, nbin_pz, nbin_w, nbin_ch, nbin_E, nbin_gamma, nbin_theta, nbin_thetaT, nbin_ty, nbin_tz;
  bool xmin_b, xmax_b, pxmin_b, pxmax_b, ymin_b, ymax_b, pymin_b, pymax_b, zmin_b, zmax_b, pzmin_b, pzmax_b, wmin_b, wmax_b, chmin_b, chmax_b, Emin_b, Emax_b,
    tymin_b, tymax_b, tzmin_b, tzmax_b, gammamin_b, gammamax_b, thetamin_b, thetamax_b, thetaTmin_b, thetaTmax_b, nbin_b, nbin_x_b, nbin_y_b,
    nbin_z_b, nbin_ty_b, nbin_tz_b, nbin_px_b, nbin_py_b, nbin_pz_b, nbin_w_b, nbin_ch_b, nbin_E_b, nbin_theta_b, nbin_thetaT_b, nbin_gamma_b;
  std::vector<aladyn_float> where_to_cut_grid_along_x, where_to_cut_grid_along_y, where_to_cut_grid_along_z;
  std::vector<aladyn_float> xcoord, ycoord, zcoord;
  std::vector<float> realpar;
  std::vector<int> intpar;
  bool overwrite_weight, overwrite_charge;
  bool do_not_ask_missing;
  aladyn_float overwrite_weight_value, overwrite_charge_value;
  bool do_plot_wspec, do_plot_Espec, do_plot_thetaspec, do_plot_thetaTspec, do_plot_chspec, do_plot_Etheta, do_plot_EthetaT;
  bool do_plot_xy, do_plot_xz, do_plot_yz, do_plot_xpx, do_plot_xpy, do_plot_xpz, do_plot_ypx;
  bool do_plot_ypy, do_plot_ypz, do_plot_zpx, do_plot_zpy, do_plot_zpz, do_plot_pxpy, do_plot_pxpz, do_plot_pypz;
  bool do_plot_xw, do_plot_rcf;
  int subsample;
  int span;
  bool grid_file, phasespace_file;
  aladyn_float get_size_x();
  aladyn_float get_size_y();
  aladyn_float get_size_z();
  aladyn_float get_size_ty();
  aladyn_float get_size_tz();
  aladyn_float get_size_px();
  aladyn_float get_size_py();
  aladyn_float get_size_pz();
  aladyn_float get_size_w();
  aladyn_float get_size_ch();
  aladyn_float get_size_gamma();
  aladyn_float get_size_theta();
  aladyn_float get_size_thetaT();
  aladyn_float get_size_E();
  aladyn_float get_size_(int);
  aladyn_float get_min(int);
  aladyn_float get_max(int);
  int get_nbin(int);
  void check_forced_version(const int, const char **);
  void parse_command_line();
  void debug_read_parameters();
  void ask_file_endianness();
  void ask_file_dims();
  bool check_params();
  void check_swap();
  void check_filename(const char *);
  void read_params_from_dat_file(std::ifstream &);
  void read_params_from_bin_file(const char *);
  void man(const char[]);
  };



#endif
