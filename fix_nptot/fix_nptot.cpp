#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cstdlib>
#if defined (_MSC_VER)
#include<wchar.h>
#define fseeko _fseeki64
#define ftello _ftelli64
#endif

#if defined (__MINGW32__)
#define fseeko fseeko64
#define ftello ftello64
#endif

#define NUMERO_PARAMETRI_FILE_DAT	20

int main(int argc, char* argv[])
{
  int indice_multifile = 0;
  std::ostringstream nomefile_bin, nomefile_dat;
  nomefile_dat << std::string(argv[1]) << ".dat";


  std::ifstream datfile;
  std::vector<int> intpar(NUMERO_PARAMETRI_FILE_DAT,0);
  std::string riga_persa;
  int nptot_dat, ndv;
  unsigned long long int nptot_calculated = 0;

  datfile.open(nomefile_dat.str().c_str(), std::ios::in);
  if (datfile.fail())
  {
    std::cout << "Unable to find " << argv[1] << ".dat" << std::endl;
    exit(-1);
  }

  std::getline(datfile, riga_persa);	// per leggere la riga Integer parameters

  for (int i = 0; i < NUMERO_PARAMETRI_FILE_DAT; i++)
  {
    datfile >> intpar[i];
    if (datfile.fail())
    {
      datfile.clear();
      std::cout << "Unable to parse int_par #" << i + 1 << std::endl;
      datfile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    }
  }

  nptot_dat = intpar[16];
  ndv = intpar[17];

  unsigned long long int dim_file_in_bytes, num_of_floats_in_file, num_of_particles_in_file;
  FILE * binfile = NULL;

  while (true)
  {
    nomefile_bin.clear();
    nomefile_bin << std::string(argv[1]) << "_" << std::setfill('0') << std::setw(3) << indice_multifile << ".bin";
    binfile = fopen(nomefile_bin.str().c_str(), "rb");
    if (binfile == NULL)
    {
      std::cout << "Reading ended at file #" << indice_multifile << std::endl;
      break;
    }

    fseeko(binfile, 0, SEEK_END);
    dim_file_in_bytes = ftello(binfile);
    rewind(binfile);
    num_of_floats_in_file = (dim_file_in_bytes / sizeof(float));
    num_of_particles_in_file = (int)(num_of_floats_in_file / ndv);
    nptot_calculated += num_of_particles_in_file;
    indice_multifile++;
    fclose(binfile);
  }

  std::cout << "Read from dat this nptot value: " << nptot_dat << std::endl;
  std::cout << "Calculated from bin files this nptot value: " << nptot_calculated << std::endl;

  return 0;
}

