#ifndef __LEGGI_ALADYN_FORTRAN
#define __LEGGI_ALADYN_FORTRAN

//#define ENABLE_DEBUG

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#define _SCL_SECURE_NO_WARNINGS

#define MAJOR_RELEASE  6
#define MINOR_RELEASE  2
#define BUGFIX_RELEASE 2

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


#define MAX_NUM_OF_PARTICLES_PER_SHOT 10000000      // in realta' il numero massimo che viene caricato in memoria e' il doppio di questo -1

#define MAX(x,y)                  ((x)>(y)?(x):(y))
#define MIN(x,y)                  ((x)<(y)?(x):(y))
#define TRUE                      1
#define FALSE                     0

#define UMA_G                     1.660538921E-24     // fattore di conversione uma -> grammi
#define C                         2.99792458E+10      // cm / s
#define ME_G                      9.10938291E-28      // electron mass [g]
#define MP_G                      1.6726231E-24       // proton mass [g]
#define MP_MEV                    938.272013          // proton mass [MeV/c^2]
#define ME_MEV                    0.510998928         // electron mass [MeV/c^2]
#define MHI_UMA                   26.981538           // atomic weight of Aluminum in atomic mass units
#define MLI_UMA                   12.0107             // atomic weight of Carbon in atomic mass units
//#define CHARGE                  4.80320425e-10      // statC  commentata perche' Turchetti la usa un po' diversa
#define CHARGE                    4.803262e-10        // statC    valore usato da Turchetti; nb: e' impreciso negli ultimi due decimali
#define FROM_TESLA_TO_GAUSS       1.0e+4
//#define DA_ERG_A_MEV            6.2415097523028e+5  // conversione Servizi
#define DA_ERG_A_MEV              6.241509744512e+5   // conversione Sinigardi
#define FROM_VOLT_TO_STATVOLT     3.335640951982e-3   // 1 statvolt = 299.792458 volts.

#define NUMERO_MASSIMO            1.0e30
#define MAX_LENGTH_FILENAME       200
#define MAX_NUMBER_OF_CPUS        32768
#define NUMERO_PARAMETRI_FILE_DAT 20
#define NUMBER_OF_BIN_BY_DEFAULT  120

#define NPARAMETRI                17
#define WEIGHT                    0
#define SWAP                      1
#define OUT_VTK                   2
#define OUT_PROPAGA               3
#define FIND_MINMAX               4
#define DO_BINNING                5
#define OUT_PARAMS                6
#define OUT_CSV                   7
#define NCOLONNE                  8 
/*  per i dump degli spazi delle fasi qui memorizziamo il numero di colonne presenti nel file binario [dovrebbe coincidere con ndv]
(valori pari a 4 o 5 significano simulazioni 2D con o senza weight, valori tipo 6 o 7 invece sono per sim 3D con o senza weight)
per i dump dei dati su griglia qui invece memorizziamo quanti sono i punti (ricampionati) lungo z (> 1 significa che la griglia e' 3D) */
#define OUT_XYZE                  9
#define OUT_CUTX                  10
#define OUT_CUTY                  11
#define OUT_CUTZ                  12
#define OUT_GRID2D                13
#define OUT_CLEAN_BINARY          14
#define OUT_LINEOUT_X             15
#define OUT_VTK_NOSTRETCH         16

#define SEI_DIMENSIONI            6                   // x, y, z, px, py, pz
#define ALTRI_PARAMETRI           8                   // gamma, theta, thetaT, E, ty, tz, w, ch

#define FORCE_PRINTF_BUFFER_SIZE  1024


struct Parametri
{
  int argc;
  std::string * argv;
  std::string filebasename;
  unsigned int ncpu_x, ncpu_y, ncpu_z, ncpu;
  unsigned int ndv, fattore_ricampionamento;
  size_t npx, npy, npz, npx_per_cpu, npy_per_cpu, npz_per_cpu;
  size_t npx_ricampionati, npy_ricampionati, npz_ricampionati, npx_ricampionati_per_cpu, npy_ricampionati_per_cpu, npz_ricampionati_per_cpu;
  size_t header_size_bytes;
  size_t nptot;
  int endianness;
  unsigned int file_version;
  bool fixed_aladyn_version;
  bool multifile;
  bool stretched_grid;  //nb: il programma al momento ragiona che se trova una griglia da qualche parte (nel dat o in fondo al bin), allora la assume stretchata anche se magari non lo e'
  int stretched_along_x;
  float massa_particella_MeV;
  int nbin, nbin_x, nbin_y, nbin_z, nbin_px, nbin_py, nbin_pz, nbin_w, nbin_ch, nbin_E, nbin_gamma, nbin_theta, nbin_thetaT, nbin_ty, nbin_tz;
  int endian_file, endian_machine;
  int last_cpu;
  int nparams;
  int p[NPARAMETRI];
  bool p_b[NPARAMETRI];
  char support_label[MAX_LENGTH_FILENAME];
  float minimi[SEI_DIMENSIONI + ALTRI_PARAMETRI], massimi[SEI_DIMENSIONI + ALTRI_PARAMETRI];
  float tnow, xmin, xmax, pxmin, pxmax, ymin, ymax, pymin, pymax, zmin, zmax, pzmin, pzmax, wmin, wmax, chmin, chmax, Emin, Emax, gammamin, gammamax, thetamin, thetamax, thetaTmin, thetaTmax, tymin, tymax, tzmin, tzmax;
  std::vector<float> posizioni_taglio_griglia_x, posizioni_taglio_griglia_y, posizioni_taglio_griglia_z;
  std::vector<float> xcoord, ycoord, zcoord, realpar;
  std::vector<int> intpar;
  bool xmin_b, xmax_b, pxmin_b, pxmax_b, ymin_b, ymax_b, pymin_b, pymax_b, zmin_b, zmax_b, pzmin_b, pzmax_b, wmin_b, wmax_b, chmin_b, chmax_b, Emin_b, Emax_b,
    tymin_b, tymax_b, tzmin_b, tzmax_b, gammamin_b, gammamax_b, thetamin_b, thetamax_b, thetaTmin_b, thetaTmax_b, nbin_b, nbin_x_b, nbin_y_b,
    nbin_z_b, nbin_ty_b, nbin_tz_b, nbin_px_b, nbin_py_b, nbin_pz_b, nbin_w_b, nbin_ch_b, nbin_E_b, nbin_theta_b, nbin_thetaT_b, nbin_gamma_b;
  bool overwrite_weight, overwrite_charge;
  bool do_not_ask_missing;
  float overwrite_weight_value, overwrite_charge_value;
  int fai_plot_wspec, fai_plot_Espec, fai_plot_thetaspec, fai_plot_thetaTspec, fai_plot_chspec, fai_plot_Etheta, fai_plot_EthetaT;
  int fai_plot_xy, fai_plot_xz, fai_plot_yz, fai_plot_xpx, fai_plot_xpy, fai_plot_xpz, fai_plot_ypx;
  int fai_plot_ypy, fai_plot_ypz, fai_plot_zpx, fai_plot_zpy, fai_plot_zpz, fai_plot_pxpy, fai_plot_pxpz, fai_plot_pypz;
  int fai_plot_xw, fai_plot_rcf;
  int subsample;
  int span;
  bool file_griglia, file_spaziofasi;
  bool file_particelle_P, file_particelle_E, file_particelle_HI, file_particelle_LI, file_particelle_generic_ion;
  bool file_campi_Ex, file_campi_Ey, file_campi_Ez, file_campi_Bx, file_campi_By, file_campi_Bz;
  bool file_densita_elettroni, file_densita_protoni, file_densita_HI, file_densita_LI, file_densita_generic_ion, file_densita_driver;
  bool file_densita_energia_griglia_elettroni, file_densita_energia_griglia_protoni, file_densita_energia_griglia_HI, file_densita_energia_griglia_LI, file_densita_energia_griglia_generic_ion;
  Parametri();
  float dimmi_dimx();
  float dimmi_dimy();
  float dimmi_dimz();
  float dimmi_dimty();
  float dimmi_dimtz();
  float dimmi_dimpx();
  float dimmi_dimpy();
  float dimmi_dimpz();
  float dimmi_dimw();
  float dimmi_dimch();
  float dimmi_dimgamma();
  float dimmi_dimtheta();
  float dimmi_dimthetaT();
  float dimmi_dimE();
  float dimmi_dim(int);
  int dimmi_nbin(int);
  void check_forced_version(const int, const char **);
  void parse_command_line();
  void leggi_file_dat(std::ifstream &);
  void debug_read_parameters();
  void chiedi_endian_file();
  void chiedi_numero_colonne();
  void chiedi_2Do3D();
  bool check_parametri();
  void check_filename(const char *);
  void organizza_minimi_massimi();
  void leggi_parametri_da_file_bin(const char *);
};



#endif
