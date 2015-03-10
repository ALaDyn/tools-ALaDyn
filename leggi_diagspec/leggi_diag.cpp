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

  int32_t mod_id, dmodel_id, LP_ord, der_ord;
  int32_t Z_i, A_i, Z1_i, A1_i, Z2_i, A2_i, iform, ibeam, str;
  double xmax, xmin, ymax, ymin;
  double lam0, w0x, w0y, chann_rad;
  double a0, lp_int, lp_pow;
  double targ_x1, targ_x2, n_over_nc, el_lp;
  double np1, lx1, lx3, np2, lx5;
  double ompe2, nmacro, np_over_nmacro, np_per_cell;
  int32_t Nx, Ny, Nz, n_cell, Nsp, Nsb;
  int32_t iter, nst, sp_step, nvar, npvar;

  float * timesteps;
  float * Emean_e;
  float * Emax_e;
  float * Emean_p;
  float * Emax_p;
  float * Emean_i;
  float * Emax_i;
  float * px;
  float * py;
  float * pz;
  float * pang;
  float * charge;
  float * Ex2;
  float * Ey2;
  float * Ez2;
  float * Ex_max;
  float * Ey_max;
  float * Ez_max;

  void allocate_arrays();

  void print_debug_v1_v2();
  void print_debug_v3();
  void read_header_v1_v2(std::ifstream &);
  void read_header_v3(std::ifstream &);
  void read_diag_v1_v2(std::ifstream &);
  void read_diag_v3(std::ifstream &);
  void read_spec_v1_v2(std::ifstream &);
  void read_spec_v3(std::ifstream &);
};

Diag_data::Diag_data() {
  timesteps = Emean_e = Emax_e = Emean_p = Emax_p = Emean_i = Emax_i = px = py = pz = pang = charge = Ex2 = Ey2 = Ez2 = Ex_max = Ey_max = Ez_max = NULL;
}


void Diag_data::allocate_arrays() {
  timesteps = new float[nst];
  Emean_e = new float[nst];
  Emax_e = new float[nst];
  Emean_p = new float[nst];
  Emax_p = new float[nst];
  Emean_i = new float[nst];
  Emax_i = new float[nst];
  px = new float[nst];
  py = new float[nst];
  pz = new float[nst];
  pang = new float[nst];
  charge = new float[nst];
  Ex2 = new float[nst];
  Ey2 = new float[nst];
  Ez2 = new float[nst];
  Ex_max = new float[nst];
  Ey_max = new float[nst];
  Ez_max = new float[nst];
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

void Diag_data::read_header_v1_v2(std::ifstream &infile) {
  std::string riga;
  std::vector<std::string> tokens;

  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  if (tokens[0][0] == 'n') {
    std::getline(infile, riga);
    std::getline(infile, riga);
    std::getline(infile, riga);
    std::getline(infile, riga);
  }

  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  mod_id = boost::lexical_cast<boost::int32_t>(tokens[0]);
  dmodel_id = boost::lexical_cast<boost::int32_t>(tokens[1]);
  LP_ord = boost::lexical_cast<boost::int32_t>(tokens[2]);
  der_ord = boost::lexical_cast<boost::int32_t>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  Z_i = boost::lexical_cast<boost::int32_t>(tokens[0]);
  A_i = boost::lexical_cast<boost::int32_t>(tokens[1]);
  iform = boost::lexical_cast<boost::int32_t>(tokens[2]);
  ibeam = boost::lexical_cast<boost::int32_t>(tokens[3]);

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
  Nx = boost::lexical_cast<boost::int32_t>(tokens[0]);
  Ny = boost::lexical_cast<boost::int32_t>(tokens[1]);
  Nz = boost::lexical_cast<boost::int32_t>(tokens[2]);
  n_cell = boost::lexical_cast<boost::int32_t>(tokens[3]);
  Nsp = boost::lexical_cast<boost::int32_t>(tokens[4]);
  Nsb = boost::lexical_cast<boost::int32_t>(tokens[5]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  iter = boost::lexical_cast<boost::int32_t>(tokens[0]);
  nst = boost::lexical_cast<boost::int32_t>(tokens[1]);
  sp_step = boost::lexical_cast<boost::int32_t>(tokens[2]);
  nvar = boost::lexical_cast<boost::int32_t>(tokens[3]);
  npvar = boost::lexical_cast<boost::int32_t>(tokens[4]);
}

void Diag_data::read_header_v3(std::ifstream &infile) {
  std::string riga;
  std::vector<std::string> tokens;

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  mod_id = boost::lexical_cast<boost::int32_t>(tokens[0]);
  dmodel_id = boost::lexical_cast<boost::int32_t>(tokens[1]);
  LP_ord = boost::lexical_cast<boost::int32_t>(tokens[2]);
  der_ord = boost::lexical_cast<boost::int32_t>(tokens[3]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  Z1_i = boost::lexical_cast<boost::int32_t>(tokens[0]);
  A1_i = boost::lexical_cast<boost::int32_t>(tokens[1]);
  Z2_i = boost::lexical_cast<boost::int32_t>(tokens[2]);
  A2_i = boost::lexical_cast<boost::int32_t>(tokens[3]);
  iform = boost::lexical_cast<boost::int32_t>(tokens[4]);
  str = boost::lexical_cast<boost::int32_t>(tokens[5]);

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
  Nx = boost::lexical_cast<boost::int32_t>(tokens[0]);
  Ny = boost::lexical_cast<boost::int32_t>(tokens[1]);
  Nz = boost::lexical_cast<boost::int32_t>(tokens[2]);
  n_cell = boost::lexical_cast<boost::int32_t>(tokens[3]);
  Nsp = boost::lexical_cast<boost::int32_t>(tokens[4]);
  Nsb = boost::lexical_cast<boost::int32_t>(tokens[5]);

  riga.clear(), std::getline(infile, riga);
  riga.clear(), tokens.clear(), std::getline(infile, riga), boost::algorithm::trim(riga);
  boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(": =\t"), boost::token_compress_on);
  iter = boost::lexical_cast<boost::int32_t>(tokens[0]);
  nst = boost::lexical_cast<boost::int32_t>(tokens[1]);
  nvar = boost::lexical_cast<boost::int32_t>(tokens[2]);
  npvar = boost::lexical_cast<boost::int32_t>(tokens[3]);
}

void Diag_data::read_diag_v1_v2(std::ifstream &infile) {
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
    infile >> Emean_e[ik] >> Emax_e[ik] >> Emean_p[ik] >> Emax_p[ik] >> Emean_i[ik] >> Emax_i[ik];
  }

  std::getline(infile, riga);
  std::getline(infile, riga);

  for (int ik = 0; ik < nst; ik++)
  {
    infile >> px[ik] >> py[ik] >> pz[ik] >> pang[ik] >> charge[ik];
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
  outfile.open(nomefile_out.str().c_str(), std::ifstream::out);

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
    outfile << timesteps[ik] << "\t" << Emean_e[ik] << "\t" << Emax_e[ik] << "\t" << Emean_p[ik] << "\t" << Emax_p[ik]
      << "\t" << Emean_i[ik] << "\t" << Emax_i[ik] << "\t" << px[ik] << "\t" << py[ik] << "\t" << pz[ik]
      << "\t" << pang[ik] << "\t" << charge[ik] << "\t" << Ex2[ik] << "\t" << Ey2[ik] << "\t" << Ez2[ik]
      << "\t" << Ex_max[ik] << "\t" << Ey_max[ik] << "\t" << Ez_max[ik] << std::endl;
  }

  outfile.close();

}

void Diag_data::read_diag_v3(std::ifstream &infile) {
  std::string riga;


  for (int ik = 0; ik < nst; ik++)
  {
    infile >> timesteps[ik];
  }

  for (int ik = 0; ik < nst; ik++)
  {
    infile >> Emean_e[ik] >> Emax_e[ik] >> Emean_p[ik] >> Emax_p[ik] >> Emean_i[ik] >> Emax_i[ik];
  }

  std::getline(infile, riga);
  std::getline(infile, riga);

  for (int ik = 0; ik < nst; ik++)
  {
    infile >> px[ik] >> py[ik] >> pz[ik] >> pang[ik] >> charge[ik];
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
  outfile.open(nomefile_out.str().c_str(), std::ifstream::out);

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
    outfile << timesteps[ik] << "\t" << Emean_e[ik] << "\t" << Emax_e[ik] << "\t" << Emean_p[ik] << "\t" << Emax_p[ik]
      << "\t" << Emean_i[ik] << "\t" << Emax_i[ik] << "\t" << px[ik] << "\t" << py[ik] << "\t" << pz[ik]
      << "\t" << pang[ik] << "\t" << charge[ik] << "\t" << Ex2[ik] << "\t" << Ey2[ik] << "\t" << Ez2[ik]
      << "\t" << Ex_max[ik] << "\t" << Ey_max[ik] << "\t" << Ez_max[ik] << std::endl;
  }

  outfile.close();

}

void Diag_data::read_spec_v1_v2(std::ifstream &infile) {
  int ntimestep = sp_step;
  std::string riga, tipo;
  float time, emax, energy, dE, *spectrum, *selected_spectrum;
  int nbin;
  std::getline(infile, riga);
  infile >> nbin;

  selected_spectrum = new float[nbin];
  spectrum = new float[nbin];


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
      outfile.open(nomefile_out.str().c_str(), std::ifstream::out);


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
          outfile << energy << "\t" << spectrum[ik] << "\t" << selected_spectrum[ik] << std::endl;
        }
        else
          outfile << energy << "\t" << spectrum[ik] << std::endl;
        energy += dE;
      }

      outfile.close();
    }
  }
}


int main(int argc, const char* argv[]) {
  if (argc < 3) {
    std::cout << "Usage: ./leggi_diag filename v1/v2/v3" << std::endl;
    exit(1);
  }

  Diag_data oDiag_data;

  std::ifstream infile;
  std::ostringstream nomefile_out;

  infile.open(argv[1], std::ifstream::in);
  if (infile.fail()) {
    std::cout << "Unable to open" << argv[1] << std::endl;
    exit(2);
  }

  oDiag_data.versione = boost::lexical_cast<boost::uint32_t>(argv[2][1]);
  oDiag_data.tipofile = argv[1][0];
  oDiag_data.nomefile = argv[1];
  if (oDiag_data.tipofile != 's' && oDiag_data.tipofile != 'd')  std::cout << "File non riconosciuto" << std::endl;

  if (oDiag_data.versione < 3) oDiag_data.read_header_v1_v2(infile);
  else oDiag_data.read_header_v3(infile);

  oDiag_data.allocate_arrays();

  if (oDiag_data.tipofile == 'd')
  {
    if (oDiag_data.versione < 3) oDiag_data.read_diag_v1_v2(infile);
    else oDiag_data.read_diag_v3(infile);
  }


  if (oDiag_data.tipofile == 's')
  {
    if (oDiag_data.versione < 3) oDiag_data.read_spec_v1_v2(infile);
    else oDiag_data.read_spec_v3(infile);
  }

  infile.close();

  return 0;
}


