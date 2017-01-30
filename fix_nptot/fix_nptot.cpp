
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <climits>
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

#if defined (__CYGWIN__)
#define fseeko fseek
#define ftello ftell
#endif


#define NUMERO_PARAMETRI_FILE_DAT	20
#define ALADYN_VERSION 3


void fix_nptot_dat_file(char *, unsigned long long int, unsigned int);


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

  datfile.close();

  nptot_dat = intpar[16];
  ndv = intpar[17];

  unsigned long long int dim_file_in_bytes, num_of_floats_in_file, num_of_particles_in_file;
  FILE * binfile = NULL;

  while (true)
  {
    nomefile_bin.clear();
    nomefile_bin.str("");
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

  if (nptot_dat != nptot_calculated) fix_nptot_dat_file(argv[1], nptot_calculated, ALADYN_VERSION);

  return 0;
}



void fix_nptot_dat_file(char * nomefile, unsigned long long int nptot_corretto, unsigned int aladyn_version)
{
  std::ifstream datfile_in;
  std::ofstream datfile_out;
  std::ostringstream nomefile_in, nomefile_out, dati_out;

  nomefile_in << std::string(nomefile) << ".dat";
  nomefile_out << std::string(nomefile) << "_fix.dat";

  std::vector<unsigned int> intpar(NUMERO_PARAMETRI_FILE_DAT, 0);
  std::vector<std::string> datfile;
  std::string riga;

  datfile_in.open(nomefile_in.str().c_str(), std::ios::in);
  if (datfile_in.fail())
  {
    std::cout << "Unable to open " << nomefile << ".dat" << std::endl;
    exit(-1);
  }

  datfile_out.open(nomefile_out.str().c_str(), std::ios::out);
  if (datfile_out.fail())
  {
    std::cout << "Unable to open " << nomefile << "_fix.dat" << std::endl;
    exit(-1);
  }

  std::getline(datfile_in, riga);
  riga += '\n';
  datfile.push_back(riga);

  for (int i = 0; i < NUMERO_PARAMETRI_FILE_DAT; i++)
  {
    datfile_in >> intpar[i];
    if (datfile_in.fail())
    {
      datfile_in.clear();
      datfile_in.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    }
  }

  if (nptot_corretto < UINT_MAX) intpar[16] = (unsigned int)nptot_corretto;
  else std::cout << "nptot too large for a uint value!" << std::endl;
  if (intpar[18] != aladyn_version) intpar[18] = aladyn_version;

  for (int i = 0; i < NUMERO_PARAMETRI_FILE_DAT; i++)
  {
    dati_out << std::setw(14) << intpar[i];
    if (i>0 && !((i+1) % 4))
    {
      dati_out << std::endl;
      datfile.push_back(dati_out.str());
      dati_out.clear();
      dati_out.str("");
    }
  }
  if (dati_out.str().size())
  {
    dati_out << std::endl;
    datfile.push_back(dati_out.str());
  }

  std::getline(datfile_in, riga);	// per pulire i caratteri rimanenti sull'ultima riga degli interi, non salvata

  while (!datfile_in.eof())
  {
    std::getline(datfile_in, riga);
    riga += '\n';
    datfile.push_back(riga);
  }

//  for (auto i : datfile) datfile_out << i;
  for (unsigned int i = 0; i < datfile.size(); i++) datfile_out << datfile[i];


  datfile_in.close();
  datfile_out.close();
}


