#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdint>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

class Diag_data
{
public:
  Diag_data();

  int versione;
  std::string nomefile;
  char tipofile;

  bool multifile;
  int32_t nfl, nfl_counter;
  int32_t nstot;

  bool btmax;
  double tmax;

  int32_t mod_id, dmodel_id, LP_ord, der_ord;
  int32_t Z_i, A_i, Z1_i, A1_i, Z2_i, A2_i, iform, ibeam, str;
  double xmax, xmin, ymax, ymin;
  double lam0, w0x, w0y, chann_rad;
  double a0, lp_int, lp_pow;
  double targ_x1, targ_x2, n_over_nc, el_lp;
  double np1, np2, lx1, lx2, lx3, lx4, lx5;
  double ompe2, nmacro, np_over_nmacro, np_per_cell;
  int32_t Nx, Ny, Nz, n_cell, Nsp, Nsb;
  int32_t iter, nst, sp_step, nvar, npvar;

  double * timesteps;
  double ** Etot;
  double ** Emean;
  double ** Emax;
  double ** px;
  double ** py;
  double ** pz;
  double ** sigma_px;
  double ** sigma_py;
  double ** sigma_pz;
  double ** Jz;
  double ** mean_charge;
  double ** charge_per_cell;
  double * charge_tot;
  double * Ex2;
  double * Ey2;
  double * Ez2;
  double * Bx2;
  double * By2;
  double * Bz2;
  double * Ex2_on_solid_target;
  double * Ey2_on_solid_target;
  double * Ez2_on_solid_target;
  double * Bx2_on_solid_target;
  double * By2_on_solid_target;
  double * Bz2_on_solid_target;
  double * Ex_max;
  double * Ey_max;
  double * Ez_max;
  double * Bx_max;
  double * By_max;
  double * Bz_max;

  void allocate_arrays();
  void deallocate_arrays();
  void start(std::ifstream &);
  void print_debug_v1_v2();
  void print_debug_v3();
  void print_debug_v4();
  void read_header_v1_v2(std::ifstream &);
  void read_header_v3(std::ifstream &);
  void read_header_v4(std::ifstream &);
  void decode_diag_v1_v2(std::ifstream &);
  void decode_diag_v3(std::ifstream &);
  void decode_diag_v4(std::ifstream &);
  void decode_spec_v1_v2(std::ifstream &);
  void decode_spec_v3(std::ifstream &);
  void decode_spec_v4(std::ifstream &);
};

Diag_data::Diag_data() {
  timesteps = nullptr;
  charge_tot = nullptr;
  multifile = false;
  nstot = 0;
  nfl = 1;
  nfl_counter = 0;
  btmax = false;
  tmax = 0.0;
  Etot = nullptr;
  Emean = nullptr;
  Emax = nullptr;
  px = nullptr;
  py = nullptr;
  pz = nullptr;
  sigma_px = nullptr;
  sigma_py = nullptr;
  sigma_pz = nullptr;
  Jz = nullptr;
  mean_charge = nullptr;
  charge_per_cell = nullptr;
  Ex2 = nullptr;
  Ey2 = nullptr;
  Ez2 = nullptr;
  Bx2 = nullptr;
  By2 = nullptr;
  Bz2 = nullptr;
  Ex2_on_solid_target = nullptr;
  Ey2_on_solid_target = nullptr;
  Ez2_on_solid_target = nullptr;
  Bx2_on_solid_target = nullptr;
  By2_on_solid_target = nullptr;
  Bz2_on_solid_target = nullptr;
  Ex_max = nullptr;
  Ey_max = nullptr;
  Ez_max = nullptr;
  Bx_max = nullptr;
  By_max = nullptr;
  Bz_max = nullptr;
}

void Diag_data::start(std::ifstream &infile) {
  if (versione < 3) read_header_v1_v2(infile);
  else if (versione == 3) read_header_v3(infile);
  else read_header_v4(infile);
}

void Diag_data::allocate_arrays() {
  timesteps = new double[nst];

  Etot = new double*[Nsp];
  Emean = new double*[Nsp];
  Emax = new double*[Nsp];
  px = new double*[Nsp];
  py = new double*[Nsp];
  pz = new double*[Nsp];
  sigma_px = new double*[Nsp];
  sigma_py = new double*[Nsp];
  sigma_pz = new double*[Nsp];
  Jz = new double*[Nsp];
  mean_charge = new double*[Nsp];
  charge_per_cell = new double*[Nsp];
  for (int32_t i = 0; i < Nsp; i++) {
    Etot[i] = new double[nst];
    Emean[i] = new double[nst];
    Emax[i] = new double[nst];
    px[i] = new double[nst];
    py[i] = new double[nst];
    pz[i] = new double[nst];
    sigma_px[i] = new double[nst];
    sigma_py[i] = new double[nst];
    sigma_pz[i] = new double[nst];
    Jz[i] = new double[nst];
    mean_charge[i] = new double[nst];
    charge_per_cell[i] = new double[nst];
  }
  charge_tot = new double[nst];
  Ex2 = new double[nst];
  Ey2 = new double[nst];
  Ez2 = new double[nst];
  Bx2 = new double[nst];
  By2 = new double[nst];
  Bz2 = new double[nst];
  Ex2_on_solid_target = new double[nst];
  Ey2_on_solid_target = new double[nst];
  Ez2_on_solid_target = new double[nst];
  Bx2_on_solid_target = new double[nst];
  By2_on_solid_target = new double[nst];
  Bz2_on_solid_target = new double[nst];
  Ex_max = new double[nst];
  Ey_max = new double[nst];
  Ez_max = new double[nst];
  Bx_max = new double[nst];
  By_max = new double[nst];
  Bz_max = new double[nst];
}

void Diag_data::deallocate_arrays() {
  delete[] timesteps; timesteps = nullptr;

  for (int32_t i = 0; i < Nsp; i++) {
    delete[] Etot[i]; Etot[i] = nullptr;
    delete[] Emean[i]; Emean[i] = nullptr;
    delete[] Emax[i]; Emax[i] = nullptr;
    delete[] px[i]; px[i] = nullptr;
    delete[] py[i]; py[i] = nullptr;
    delete[] pz[i]; pz[i] = nullptr;
    delete[] sigma_px[i]; sigma_px[i] = nullptr;
    delete[] sigma_py[i]; sigma_py[i] = nullptr;
    delete[] sigma_pz[i]; sigma_pz[i] = nullptr;
    delete[] Jz[i]; Jz[i] = nullptr;
    delete[] mean_charge[i]; mean_charge[i] = nullptr;
    delete[] charge_per_cell[i]; charge_per_cell[i] = nullptr;
  }
  delete[] Etot; Etot = nullptr;
  delete[] Emean; Emean = nullptr;
  delete[] Emax; Emax = nullptr;
  delete[] px; px = nullptr;
  delete[] py; py = nullptr;
  delete[] pz; pz = nullptr;
  delete[] sigma_px; sigma_px = nullptr;
  delete[] sigma_py; sigma_py = nullptr;
  delete[] sigma_pz; sigma_pz = nullptr;
  delete[] Jz; Jz = nullptr;
  delete[] mean_charge; mean_charge = nullptr;
  delete[] charge_per_cell; charge_per_cell = nullptr;

  delete[] charge_tot; charge_tot = nullptr;
  delete[] Ex2; Ex2 = nullptr;
  delete[] Ey2; Ey2 = nullptr;
  delete[] Ez2; Ez2 = nullptr;
  delete[] Bx2; Bx2 = nullptr;
  delete[] By2; By2 = nullptr;
  delete[] Bz2; Bz2 = nullptr;
  delete[] Ex2_on_solid_target; Ex2_on_solid_target = nullptr;
  delete[] Ey2_on_solid_target; Ey2_on_solid_target = nullptr;
  delete[] Ez2_on_solid_target; Ez2_on_solid_target = nullptr;
  delete[] Bx2_on_solid_target; Bx2_on_solid_target = nullptr;
  delete[] By2_on_solid_target; By2_on_solid_target = nullptr;
  delete[] Bz2_on_solid_target; Bz2_on_solid_target = nullptr;
  delete[] Ex_max; Ex_max = nullptr;
  delete[] Ey_max; Ey_max = nullptr;
  delete[] Ez_max; Ez_max = nullptr;
  delete[] Bx_max; Bx_max = nullptr;
  delete[] By_max; By_max = nullptr;
  delete[] Bz_max; Bz_max = nullptr;
}

void Diag_data::print_debug_v1_v2() {
  std::cout << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;
  std::cout << Z_i << "\t" << A_i << "\t" << iform << "\t" << ibeam << std::endl;
  std::cout << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
  std::cout << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
  std::cout << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
  std::cout << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
  std::cout << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << std::endl;
  std::cout << ompe2 << "\t" << nmacro << "\t" << np_over_nmacro << std::endl;
  std::cout << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
  std::cout << iter << "\t" << nst << "\t" << sp_step << "\t" << nvar << "\t" << npvar << std::endl;
}

void Diag_data::print_debug_v3() {
  std::cout << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;
  std::cout << Z1_i << "\t" << A1_i << "\t" << Z2_i << "\t" << A2_i << "\t" << iform << "\t" << str << std::endl;
  std::cout << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
  std::cout << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
  std::cout << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
  std::cout << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
  std::cout << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << std::endl;
  std::cout << ompe2 << "\t" << nmacro << "\t" << np_per_cell << std::endl;
  std::cout << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
  std::cout << iter << "\t" << nst << "\t" << nvar << "\t" << npvar << std::endl;
}

void Diag_data::print_debug_v4() {
  std::cout << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;
  std::cout << Z1_i << "\t" << A1_i << "\t" << Z2_i << "\t" << A2_i << "\t" << iform << "\t" << str << std::endl;
  std::cout << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
  std::cout << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
  std::cout << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
  std::cout << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
  std::cout << lx1 << "\t" << lx2 << "\t" << lx3 << "\t" << lx4 << "\t" << lx5 << std::endl;
  std::cout << ompe2 << "\t" << nmacro << "\t" << np_per_cell << std::endl;
  std::cout << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
  std::cout << iter << "\t" << nst << "\t" << nvar << "\t" << npvar << std::endl;
}

void Diag_data::read_header_v1_v2(std::ifstream &infile) {
  std::string riga;
  std::vector<std::string> tokens;

  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  if (tokens[0] == "nfl") {
    multifile = true;
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    nfl = boost::lexical_cast<int32_t>(tokens[0]);
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    nstot = boost::lexical_cast<int32_t>(tokens[0]);
    riga.clear(), std::getline(infile, riga);
  }

  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  mod_id = boost::lexical_cast<int32_t>(tokens[0]);
  dmodel_id = boost::lexical_cast<int32_t>(tokens[1]);
  LP_ord = boost::lexical_cast<int32_t>(tokens[2]);
  der_ord = boost::lexical_cast<int32_t>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  Z_i = boost::lexical_cast<int32_t>(tokens[0]);
  A_i = boost::lexical_cast<int32_t>(tokens[1]);
  iform = boost::lexical_cast<int32_t>(tokens[2]);
  ibeam = boost::lexical_cast<int32_t>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  xmax = boost::lexical_cast<double>(tokens[0]);
  xmin = boost::lexical_cast<double>(tokens[1]);
  ymax = boost::lexical_cast<double>(tokens[2]);
  ymin = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  lam0 = boost::lexical_cast<double>(tokens[0]);
  w0x = boost::lexical_cast<double>(tokens[1]);
  w0y = boost::lexical_cast<double>(tokens[2]);
  chann_rad = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  a0 = boost::lexical_cast<double>(tokens[0]);
  lp_int = boost::lexical_cast<double>(tokens[1]);
  lp_pow = boost::lexical_cast<double>(tokens[2]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  targ_x1 = boost::lexical_cast<double>(tokens[0]);
  targ_x2 = boost::lexical_cast<double>(tokens[1]);
  n_over_nc = boost::lexical_cast<double>(tokens[2]);
  el_lp = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  np1 = boost::lexical_cast<double>(tokens[0]);
  lx1 = boost::lexical_cast<double>(tokens[1]);
  lx3 = boost::lexical_cast<double>(tokens[2]);
  np2 = boost::lexical_cast<double>(tokens[3]);
  lx5 = boost::lexical_cast<double>(tokens[4]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  ompe2 = boost::lexical_cast<double>(tokens[0]);
  nmacro = boost::lexical_cast<double>(tokens[1]);
  np_over_nmacro = boost::lexical_cast<double>(tokens[2]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  Nx = boost::lexical_cast<int32_t>(tokens[0]);
  Ny = boost::lexical_cast<int32_t>(tokens[1]);
  Nz = boost::lexical_cast<int32_t>(tokens[2]);
  n_cell = boost::lexical_cast<int32_t>(tokens[3]);
  Nsp = boost::lexical_cast<int32_t>(tokens[4]);
  Nsb = boost::lexical_cast<int32_t>(tokens[5]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  iter = boost::lexical_cast<int32_t>(tokens[0]);
  nst = boost::lexical_cast<int32_t>(tokens[1]);
  sp_step = boost::lexical_cast<int32_t>(tokens[2]);
  nvar = boost::lexical_cast<int32_t>(tokens[3]);
  npvar = boost::lexical_cast<int32_t>(tokens[4]);

  allocate_arrays();
  if (tipofile == 'd') decode_diag_v1_v2(infile);
  if (tipofile == 's') decode_spec_v1_v2(infile);
}

void Diag_data::read_header_v3(std::ifstream &infile) {
  std::string riga;
  std::vector<std::string> tokens;

  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  if (tokens[0] == "nfl") {
    multifile = true;
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    nfl = boost::lexical_cast<int32_t>(tokens[0]);
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    nstot = boost::lexical_cast<int32_t>(tokens[0]);
    riga.clear(), std::getline(infile, riga);
  }

  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  mod_id = boost::lexical_cast<int32_t>(tokens[0]);
  dmodel_id = boost::lexical_cast<int32_t>(tokens[1]);
  LP_ord = boost::lexical_cast<int32_t>(tokens[2]);
  der_ord = boost::lexical_cast<int32_t>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  Z1_i = boost::lexical_cast<int32_t>(tokens[0]);
  A1_i = boost::lexical_cast<int32_t>(tokens[1]);
  Z2_i = boost::lexical_cast<int32_t>(tokens[2]);
  A2_i = boost::lexical_cast<int32_t>(tokens[3]);
  iform = boost::lexical_cast<int32_t>(tokens[4]);
  str = boost::lexical_cast<int32_t>(tokens[5]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  xmax = boost::lexical_cast<double>(tokens[0]);
  xmin = boost::lexical_cast<double>(tokens[1]);
  ymax = boost::lexical_cast<double>(tokens[2]);
  ymin = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  lam0 = boost::lexical_cast<double>(tokens[0]);
  w0x = boost::lexical_cast<double>(tokens[1]);
  w0y = boost::lexical_cast<double>(tokens[2]);
  chann_rad = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  a0 = boost::lexical_cast<double>(tokens[0]);
  lp_int = boost::lexical_cast<double>(tokens[1]);
  lp_pow = boost::lexical_cast<double>(tokens[2]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  targ_x1 = boost::lexical_cast<double>(tokens[0]);
  targ_x2 = boost::lexical_cast<double>(tokens[1]);
  n_over_nc = boost::lexical_cast<double>(tokens[2]);
  el_lp = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  np1 = boost::lexical_cast<double>(tokens[0]);
  lx1 = boost::lexical_cast<double>(tokens[1]);
  lx3 = boost::lexical_cast<double>(tokens[2]);
  np2 = boost::lexical_cast<double>(tokens[3]);
  lx5 = boost::lexical_cast<double>(tokens[4]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  ompe2 = boost::lexical_cast<double>(tokens[0]);
  nmacro = boost::lexical_cast<double>(tokens[1]);
  np_per_cell = boost::lexical_cast<double>(tokens[2]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  Nx = boost::lexical_cast<int32_t>(tokens[0]);
  Ny = boost::lexical_cast<int32_t>(tokens[1]);
  Nz = boost::lexical_cast<int32_t>(tokens[2]);
  n_cell = boost::lexical_cast<int32_t>(tokens[3]);
  Nsp = boost::lexical_cast<int32_t>(tokens[4]);
  Nsb = boost::lexical_cast<int32_t>(tokens[5]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  iter = boost::lexical_cast<int32_t>(tokens[0]);
  nst = boost::lexical_cast<int32_t>(tokens[1]);
  nvar = boost::lexical_cast<int32_t>(tokens[2]);
  npvar = boost::lexical_cast<int32_t>(tokens[3]);

  allocate_arrays();
  if (tipofile == 'd') decode_diag_v3(infile);
  if (tipofile == 's') decode_spec_v3(infile);
}

void Diag_data::read_header_v4(std::ifstream &infile) {
  std::string riga;
  std::vector<std::string> tokens;

  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  if (tokens[0] == "nfl") {
    multifile = true;
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    nfl = boost::lexical_cast<int32_t>(tokens[0]);
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    nstot = boost::lexical_cast<int32_t>(tokens[0]);
    riga.clear(), std::getline(infile, riga);
  }

  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  mod_id = boost::lexical_cast<int32_t>(tokens[0]);
  dmodel_id = boost::lexical_cast<int32_t>(tokens[1]);
  LP_ord = boost::lexical_cast<int32_t>(tokens[2]);
  der_ord = boost::lexical_cast<int32_t>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  Z1_i = boost::lexical_cast<int32_t>(tokens[0]);
  A1_i = boost::lexical_cast<int32_t>(tokens[1]);
  Z2_i = boost::lexical_cast<int32_t>(tokens[2]);
  A2_i = boost::lexical_cast<int32_t>(tokens[3]);
  iform = boost::lexical_cast<int32_t>(tokens[4]);
  str = boost::lexical_cast<int32_t>(tokens[5]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  xmax = boost::lexical_cast<double>(tokens[0]);
  xmin = boost::lexical_cast<double>(tokens[1]);
  ymax = boost::lexical_cast<double>(tokens[2]);
  ymin = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  lam0 = boost::lexical_cast<double>(tokens[0]);
  w0x = boost::lexical_cast<double>(tokens[1]);
  w0y = boost::lexical_cast<double>(tokens[2]);
  chann_rad = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  a0 = boost::lexical_cast<double>(tokens[0]);
  lp_int = boost::lexical_cast<double>(tokens[1]);
  lp_pow = boost::lexical_cast<double>(tokens[2]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  targ_x1 = boost::lexical_cast<double>(tokens[0]);
  targ_x2 = boost::lexical_cast<double>(tokens[1]);
  n_over_nc = boost::lexical_cast<double>(tokens[2]);
  el_lp = boost::lexical_cast<double>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  lx1 = boost::lexical_cast<double>(tokens[0]);
  lx2 = boost::lexical_cast<double>(tokens[1]);
  lx3 = boost::lexical_cast<double>(tokens[2]);
  lx4 = boost::lexical_cast<double>(tokens[3]);
  lx5 = boost::lexical_cast<double>(tokens[4]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  ompe2 = boost::lexical_cast<double>(tokens[0]);
  nmacro = boost::lexical_cast<double>(tokens[1]);
  np_per_cell = boost::lexical_cast<double>(tokens[2]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  Nx = boost::lexical_cast<int32_t>(tokens[0]);
  Ny = boost::lexical_cast<int32_t>(tokens[1]);
  Nz = boost::lexical_cast<int32_t>(tokens[2]);
  n_cell = boost::lexical_cast<int32_t>(tokens[3]);
  Nsp = boost::lexical_cast<int32_t>(tokens[4]);
  Nsb = boost::lexical_cast<int32_t>(tokens[5]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  iter = boost::lexical_cast<int32_t>(tokens[0]);
  nst = boost::lexical_cast<int32_t>(tokens[1]);
  nvar = boost::lexical_cast<int32_t>(tokens[2]);
  npvar = boost::lexical_cast<int32_t>(tokens[3]);

  allocate_arrays();
  if (tipofile == 'd') decode_diag_v4(infile);
  if (tipofile == 's') decode_spec_v4(infile);
}

void Diag_data::decode_diag_v1_v2(std::ifstream &infile) {
  std::string riga;

  std::getline(infile, riga);
  std::getline(infile, riga);


  for (int ik = 0; ik < nst; ik++)
  {
    infile >> timesteps[ik];
  }

  std::getline(infile, riga);
  std::getline(infile, riga);
  std::getline(infile, riga);



  for (int ik = 0; ik < nst; ik++)
  {
    infile >> Emean[0][ik] >> Emax[0][ik] >> Emean[1][ik] >> Emax[1][ik] >> Emean[2][ik] >> Emax[2][ik];
  }

  std::getline(infile, riga);
  std::getline(infile, riga);

  for (int ik = 0; ik < nst; ik++)
  {
    infile >> px[0][ik] >> py[0][ik] >> pz[0][ik] >> Jz[0][ik] >> charge_tot[ik];
  }

  std::getline(infile, riga);
  std::getline(infile, riga);
  std::getline(infile, riga);

  for (int ik = 0; ik < nst; ik++)
  {
    infile >> Ex2[ik] >> Ey2[ik] >> Ez2[ik] >> Ex_max[ik] >> Ey_max[ik] >> Ez_max[ik];
  }

  std::ostringstream nomefile_out;
  nomefile_out << std::string(nomefile) << ".txt";
  std::ofstream outfile;
  outfile.open(nomefile_out.str().c_str(), std::ofstream::out | std::ofstream::app);

  if (multifile)  outfile << "# " << nfl << "\t" << nstot << std::endl;

  outfile << "# " << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;

  outfile << "# " << Z_i << "\t" << A_i << "\t" << iform << "\t" << ibeam << std::endl;
  outfile << "# " << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
  outfile << "# " << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
  outfile << "# " << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
  outfile << "# " << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
  outfile << "# " << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << std::endl;
  outfile << "# " << ompe2 << "\t" << nmacro << "\t" << np_over_nmacro << std::endl;
  outfile << "# " << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
  outfile << "# " << iter << "\t" << nst << "\t" << sp_step << "\t" << nvar << "\t" << npvar << std::endl;

  for (int ik = 0; ik < nst; ik++)
  {
    if (btmax && (timesteps[ik] > tmax)) continue;
    outfile << timesteps[ik] << "\t" << Emean[0][ik] << "\t" << Emax[0][ik] << "\t" << Emean[1][ik] << "\t" << Emax[1][ik]
      << "\t" << Emean[2][ik] << "\t" << Emax[2][ik] << "\t" << px[0][ik] << "\t" << py[0][ik] << "\t" << pz[0][ik]
      << "\t" << Jz[0][ik] << "\t" << charge_tot[ik] << "\t" << Ex2[ik] << "\t" << Ey2[ik] << "\t" << Ez2[ik]
      << "\t" << Ex_max[ik] << "\t" << Ey_max[ik] << "\t" << Ez_max[ik] << std::endl;
  }

  outfile.close();
  deallocate_arrays();
  if (multifile && (++nfl_counter < nfl)) read_header_v1_v2(infile);
}

void Diag_data::decode_diag_v3(std::ifstream &infile) {
  std::string riga;
  std::vector<std::string> tokens;
  std::vector<std::string> tipo(Nsp);
  int ntimesteps = 0;

  riga.clear(), std::getline(infile, riga);
  riga.clear(), std::getline(infile, riga);

  while (ntimesteps < nst)
  {
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    for (size_t j = 0; j < tokens.size(); j++) timesteps[ntimesteps++] = boost::lexical_cast<double>(tokens[j]);
  }

  riga.clear(), std::getline(infile, riga);

  for (int i = 0; i < Nsp; i++)
  {
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    tipo[i] = tokens[0];

    riga.clear(), std::getline(infile, riga);
    for (int j = 0; j < nst; j++)
    {
      riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
      Etot[i][j] = boost::lexical_cast<double>(tokens[0]);
      Emax[i][j] = boost::lexical_cast<double>(tokens[1]);
      Jz[i][j] = boost::lexical_cast<double>(tokens[2]);
      px[i][j] = boost::lexical_cast<double>(tokens[3]);
      py[i][j] = boost::lexical_cast<double>(tokens[4]);
      pz[i][j] = boost::lexical_cast<double>(tokens[5]);
    }

    riga.clear(), std::getline(infile, riga);
    for (int j = 0; j < nst; j++)
    {
      riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
      sigma_px[i][j] = boost::lexical_cast<double>(tokens[0]);
      sigma_py[i][j] = boost::lexical_cast<double>(tokens[1]);
      sigma_pz[i][j] = boost::lexical_cast<double>(tokens[2]);
      mean_charge[i][j] = boost::lexical_cast<double>(tokens[3]);
      charge_per_cell[i][j] = boost::lexical_cast<double>(tokens[4]);
    }
  }

  riga.clear(), std::getline(infile, riga);
  riga.clear(), std::getline(infile, riga);

  if (Nz == 1) {
    for (int j = 0; j < nst; j++)
    {
      riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
      Ex2[j] = boost::lexical_cast<double>(tokens[0]);
      Ey2[j] = boost::lexical_cast<double>(tokens[1]);
      Ez2[j] = 0.0;
      Bx2[j] = 0.0;
      By2[j] = 0.0;
      Bz2[j] = boost::lexical_cast<double>(tokens[2]);
      Ex_max[j] = boost::lexical_cast<double>(tokens[3]);
      Ey_max[j] = boost::lexical_cast<double>(tokens[4]);
      Ez_max[j] = 0.0;
      Bx_max[j] = 0.0;
      By_max[j] = 0.0;
      Bz_max[j] = boost::lexical_cast<double>(tokens[5]);
    }
  }
  else {
    for (int j = 0; j < nst; j++)
    {
      riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
      Ex2[j] = boost::lexical_cast<double>(tokens[0]);
      Ey2[j] = boost::lexical_cast<double>(tokens[1]);
      Ez2[j] = boost::lexical_cast<double>(tokens[2]);
      Ex_max[j] = boost::lexical_cast<double>(tokens[3]);
      Ey_max[j] = boost::lexical_cast<double>(tokens[4]);
      Ez_max[j] = boost::lexical_cast<double>(tokens[5]);
    }
    riga.clear(), std::getline(infile, riga);
    for (int j = 0; j < nst; j++)
    {
      riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
      Bx2[j] = boost::lexical_cast<double>(tokens[0]);
      By2[j] = boost::lexical_cast<double>(tokens[1]);
      Bz2[j] = boost::lexical_cast<double>(tokens[2]);
      Bx_max[j] = boost::lexical_cast<double>(tokens[3]);
      By_max[j] = boost::lexical_cast<double>(tokens[4]);
      Bz_max[j] = boost::lexical_cast<double>(tokens[5]);
    }
  }

  riga.clear(), std::getline(infile, riga);
  riga.clear(), std::getline(infile, riga);

  for (int j = 0; j < nst; j++)
  {
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    Ex2_on_solid_target[j] = boost::lexical_cast<double>(tokens[0]);
    Ey2_on_solid_target[j] = boost::lexical_cast<double>(tokens[1]);
    if (Nz > 1) {
      Ez2_on_solid_target[j] = boost::lexical_cast<double>(tokens[2]);
      Bx2_on_solid_target[j] = boost::lexical_cast<double>(tokens[3]);
      By2_on_solid_target[j] = boost::lexical_cast<double>(tokens[4]);
      Bz2_on_solid_target[j] = boost::lexical_cast<double>(tokens[5]);
    }
    else {
      Ez2_on_solid_target[j] = 0.0;
      Bx2_on_solid_target[j] = 0.0;
      By2_on_solid_target[j] = 0.0;
      Bz2_on_solid_target[j] = boost::lexical_cast<double>(tokens[2]);
    }
  }

  std::ostringstream nomefile_out;
  nomefile_out << std::string(nomefile) << ".particles.txt";
  std::ofstream outfile;
  outfile.open(nomefile_out.str().c_str(), std::ofstream::out | std::ofstream::app);

  outfile << "# " << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;

  outfile << "# " << Z1_i << "\t" << A1_i << "\t" << Z2_i << "\t" << A2_i << "\t" << iform << "\t" << str << std::endl;
  outfile << "# " << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
  outfile << "# " << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
  outfile << "# " << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
  outfile << "# " << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
  outfile << "# " << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << std::endl;
  outfile << "# " << ompe2 << "\t" << nmacro << "\t" << np_per_cell << std::endl;
  outfile << "# " << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
  outfile << "# " << iter << "\t" << nst << "\t" << nvar << "\t" << npvar << std::endl;
  outfile << "# " << "timestep Efields2 Etot  Emax  Jz  px  py  pz  sigma_px  sigma_py  sigma_pz  mean_charge  charge_per_cell" << std::endl;
  outfile << "# " << "last line (except first two columns) is repeated horizontally for each species, so that each line is a full description of a timestep" << std::endl;

  for (int i = 0; i < nst; i++)
  {
    if (btmax && (timesteps[i] > tmax)) continue;
    outfile << timesteps[i] << "\t" << Ex2[i] + Ey2[i] + Ez2[i] + Bx2[i] + By2[i] + Bz2[i];
    for (int j = 0; j < Nsp; j++) outfile << "\t" << Etot[j][i] << "\t" << Emax[j][i] << "\t" << Jz[j][i] << "\t" << px[j][i] << "\t" << py[j][i] << "\t" << pz[j][i] << "\t" << sigma_px[j][i] << "\t" << sigma_py[j][i] << "\t" << sigma_pz[j][i] << "\t" << mean_charge[j][i] << "\t" << charge_per_cell[j][i];
    outfile << std::endl;
  }

  outfile.close();

  nomefile_out.str("\0");
  nomefile_out.seekp(0, std::ios::beg);
  nomefile_out << std::string(nomefile) << ".fields.txt";
  outfile.open(nomefile_out.str().c_str(), std::ofstream::out | std::ofstream::app);

  outfile << "# " << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;

  outfile << "# " << Z1_i << "\t" << A1_i << "\t" << Z2_i << "\t" << A2_i << "\t" << iform << "\t" << str << std::endl;
  outfile << "# " << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
  outfile << "# " << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
  outfile << "# " << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
  outfile << "# " << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
  outfile << "# " << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << std::endl;
  outfile << "# " << ompe2 << "\t" << nmacro << "\t" << np_per_cell << std::endl;
  outfile << "# " << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
  outfile << "# " << iter << "\t" << nst << "\t" << nvar << "\t" << npvar << std::endl;
  outfile << "# " << "timestep  Ex2  Ey2  Ez2  Ex_max  Ey_max  Ez_max  Bx2  By2  Bz2  Bx_max  By_max  Bz_max  Ex2_on_solid_target  Ey2_on_solid_target  Ez2_on_solid_target  Bx2_on_solid_target  By2_on_solid_target  Bz2_on_solid_target" << std::endl;

  for (int i = 0; i < nst; i++)
  {
    if (btmax && (timesteps[i] > tmax)) continue;
    outfile << timesteps[i] << "\t"
      << Ex2[i] << "\t" << Ey2[i] << "\t" << Ez2[i] << "\t" << Ex_max[i] << "\t" << Ey_max[i] << "\t" << Ez_max[i] << "\t" << Bx2[i] << "\t" << By2[i] << "\t" << Bz2[i] << "\t" << Bx_max[i] << "\t" << By_max[i] << "\t" << Bz_max[i] << "\t" << Ex2_on_solid_target[i] << "\t" << Ey2_on_solid_target[i] << "\t" << Ez2_on_solid_target[i] << "\t" << Bx2_on_solid_target[i] << "\t" << By2_on_solid_target[i] << "\t" << Bz2_on_solid_target[i]
      << std::endl;
  }

  outfile.close();
  deallocate_arrays();
  if (multifile && (++nfl_counter < nfl)) read_header_v3(infile);
}


void Diag_data::decode_diag_v4(std::ifstream &infile) {
  decode_diag_v3(infile);
}


void Diag_data::decode_spec_v1_v2(std::ifstream &infile) {
  int ntimestep = sp_step;
  std::string riga, tipo;
  double time, emax, energy, dE, *spectrum, *selected_spectrum;
  int nbin;
  std::getline(infile, riga);
  infile >> nbin;

  selected_spectrum = new double[nbin];
  spectrum = new double[nbin];


  std::getline(infile, riga);

  for (int i = 0; i < 3; i++) //la procedura va ripetuta 3 volte, per le tre specie: elettroni, protoni, ioni
  {
    infile >> tipo;
    std::getline(infile, riga);
    for (int j = 0; j < ntimestep; j++)
    {
      std::getline(infile, riga);
      infile >> time >> emax;
      energy = emax / nbin;
      dE = emax / nbin;
      std::getline(infile, riga);
      std::getline(infile, riga);
      for (int ik = 0; ik < nbin; ik++)
      {
        infile >> spectrum[ik];
      }
      std::getline(infile, riga);
      if (versione == 2) {
        std::getline(infile, riga);
        for (int ik = 0; ik < nbin; ik++)
        {
          infile >> selected_spectrum[ik];
        }
        std::getline(infile, riga);
      }

      std::ostringstream nomefile_out;
      nomefile_out.str("\0");
      nomefile_out.seekp(0, std::ios::beg);
      nomefile_out << std::string(nomefile) << "_" << tipo << "_" << time << ".txt";
      std::ofstream outfile;
      outfile.open(nomefile_out.str().c_str(), std::ofstream::out | std::ofstream::app);


      outfile << "# " << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;
      outfile << "# " << Z_i << "\t" << A_i << "\t" << iform << "\t" << ibeam << std::endl;
      outfile << "# " << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
      outfile << "# " << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
      outfile << "# " << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
      outfile << "# " << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
      outfile << "# " << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << std::endl;
      outfile << "# " << ompe2 << "\t" << nmacro << "\t" << np_over_nmacro << std::endl;
      outfile << "# " << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
      outfile << "# " << iter << "\t" << nst << "\t" << sp_step << "\t" << nvar << "\t" << npvar << std::endl;

      for (int ik = 0; ik < nbin; ik++) {
        if (versione == 2) {
          outfile << energy << "\t" << spectrum[ik] * np_over_nmacro << "\t" << selected_spectrum[ik] * np_over_nmacro << std::endl;
        }
        else
          outfile << energy << "\t" << spectrum[ik] * np_over_nmacro << std::endl;
        energy += dE;
      }

      outfile.close();
    }
  }
  deallocate_arrays();
}

void Diag_data::decode_spec_v3(std::ifstream &infile) {
  std::string riga;

  std::vector<std::string> tipo(Nsp);
  std::vector<std::string> tokens;

  double time, emax, energy, dE, *spectrum, *selected_spectrum;
  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  int32_t nbin = boost::lexical_cast<int32_t>(tokens[0]);
  int32_t contatore_bin;

  selected_spectrum = new double[nbin];
  spectrum = new double[nbin];

  for (int i = 0; i < Nsp; i++)
  {
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    tipo[i] = tokens[0];

    for (int j = 0; j < nst; j++)
    {
      riga.clear(), std::getline(infile, riga);
      riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
      time = boost::lexical_cast<double>(tokens[0]);
      emax = boost::lexical_cast<double>(tokens[1]);
      energy = emax / nbin;
      dE = emax / nbin;

      contatore_bin = 0;
      riga.clear(), std::getline(infile, riga);
      while (contatore_bin < nbin) {
        riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
        boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
        for (size_t j = 0; j < tokens.size(); j++) spectrum[contatore_bin++] = boost::lexical_cast<double>(tokens[j]);
      }

      contatore_bin = 0;
      riga.clear(), std::getline(infile, riga);
      while (contatore_bin < nbin) {
        riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
        boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
        for (size_t j = 0; j < tokens.size(); j++) selected_spectrum[contatore_bin++] = boost::lexical_cast<double>(tokens[j]);
      }

      std::ostringstream nomefile_out;
      nomefile_out.str("\0");
      nomefile_out.seekp(0, std::ios::beg);
      nomefile_out << std::string(nomefile) << "_" << tipo[i] << "_" << time << ".txt";
      std::ofstream outfile;
      outfile.open(nomefile_out.str().c_str(), std::ofstream::out | std::ofstream::app);


      outfile << "# " << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;
      outfile << "# " << Z1_i << "\t" << A1_i << "\t" << Z2_i << "\t" << A2_i << "\t" << iform << "\t" << str << std::endl;
      outfile << "# " << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
      outfile << "# " << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
      outfile << "# " << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
      outfile << "# " << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
      outfile << "# " << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << std::endl;
      outfile << "# " << ompe2 << "\t" << nmacro << "\t" << np_per_cell << std::endl;
      outfile << "# " << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
      outfile << "# " << iter << "\t" << nst << "\t" << nvar << "\t" << npvar << std::endl;

      for (int ik = 0; ik < nbin; ik++) {
        outfile << energy << "\t" << spectrum[ik] * np_per_cell << "\t" << selected_spectrum[ik] * np_per_cell << std::endl;
        energy += dE;
      }

      outfile.close();
    }
  }
  deallocate_arrays();
}

void Diag_data::decode_spec_v4(std::ifstream &infile) {
  std::string riga;

  std::vector<std::string> tipo(Nsp);
  std::vector<std::string> tokens;

  double time, emax, energy, dE, *spectrum, *selected_spectrum_1, *selected_spectrum_2;
  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  int32_t nbin = boost::lexical_cast<int32_t>(tokens[0]);
  int32_t contatore_bin;

  selected_spectrum_1 = new double[nbin];
  selected_spectrum_2 = new double[nbin];
  spectrum = new double[nbin];

  for (int i = 0; i < Nsp; i++)
  {
    riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
    boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
    tipo[i] = tokens[0];

    for (int j = 0; j < nst; j++)
    {
      riga.clear(), std::getline(infile, riga);
      riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
      time = boost::lexical_cast<double>(tokens[0]);
      emax = boost::lexical_cast<double>(tokens[1]);
      energy = emax / nbin;
      dE = emax / nbin;

      contatore_bin = 0;
      riga.clear(), std::getline(infile, riga);
      while (contatore_bin < nbin) {
        riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
        boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
        for (size_t j = 0; j < tokens.size(); j++) spectrum[contatore_bin++] = boost::lexical_cast<double>(tokens[j]);
      }

      contatore_bin = 0;
      riga.clear(), std::getline(infile, riga);
      while (contatore_bin < nbin) {
        riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
        boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
        for (size_t j = 0; j < tokens.size(); j++) selected_spectrum_1[contatore_bin++] = boost::lexical_cast<double>(tokens[j]);
      }

      contatore_bin = 0;
      riga.clear(), std::getline(infile, riga);
      while (contatore_bin < nbin) {
        riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
        boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
        for (size_t j = 0; j < tokens.size(); j++) selected_spectrum_2[contatore_bin++] = boost::lexical_cast<double>(tokens[j]);
      }

      std::ostringstream nomefile_out;
      nomefile_out.str("\0");
      nomefile_out.seekp(0, std::ios::beg);
      nomefile_out << std::string(nomefile) << "_" << tipo[i] << "_" << time << ".txt";
      std::ofstream outfile;
      outfile.open(nomefile_out.str().c_str(), std::ofstream::out | std::ofstream::app);


      outfile << "# " << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << std::endl;
      outfile << "# " << Z1_i << "\t" << A1_i << "\t" << Z2_i << "\t" << A2_i << "\t" << iform << "\t" << str << std::endl;
      outfile << "# " << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << std::endl;
      outfile << "# " << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << std::endl;
      outfile << "# " << a0 << "\t" << lp_int << "\t" << lp_pow << std::endl;
      outfile << "# " << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << std::endl;
      outfile << "# " << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << std::endl;
      outfile << "# " << ompe2 << "\t" << nmacro << "\t" << np_per_cell << std::endl;
      outfile << "# " << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << std::endl;
      outfile << "# " << iter << "\t" << nst << "\t" << nvar << "\t" << npvar << std::endl;

      for (int ik = 0; ik < nbin; ik++) {
        outfile << energy << "\t" << spectrum[ik] * np_per_cell << "\t" << selected_spectrum_1[ik] * np_per_cell << "\t" << selected_spectrum_2[ik] * np_per_cell << std::endl;
        energy += dE;
      }

      outfile.close();
    }
  }
  deallocate_arrays();
}

int main(int argc, const char* argv[]) {
  if (argc < 3) {
    std::cout << "Usage: ./leggi_diag filename v1/v2/v3/v4 [-tmax tt.ttt] " << std::endl;
    exit(1);
  }

  Diag_data oDiag_data;

  std::ifstream infile;

  infile.open(argv[1], std::ifstream::in);
  if (infile.fail()) {
    std::cout << "Unable to open" << argv[1] << std::endl;
    exit(2);
  }

  oDiag_data.versione = boost::lexical_cast<int32_t>(argv[2][1]);
  if ((argc == 5) && (std::string(argv[3]) == "-tmax")) {
    oDiag_data.btmax = true;
    oDiag_data.tmax = boost::lexical_cast<double>(argv[4]);
  }
  oDiag_data.tipofile = argv[1][0];
  oDiag_data.nomefile = argv[1];
  if (oDiag_data.tipofile != 's' && oDiag_data.tipofile != 'd')  std::cout << "File non riconosciuto" << std::endl;

  oDiag_data.start(infile);

  infile.close();

  return 0;
}

