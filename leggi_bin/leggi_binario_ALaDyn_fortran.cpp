
#include "leggi_binario_ALaDyn_fortran.h"
#include "leggi_griglia.h"
#include "leggi_particelle.h"


int main(const int argc, const char *argv[])
{
  Parameters params;
  bool testParameters = true;

  std::ostringstream bin_filename, dat_filename;
  std::string forget_this_line, endianness, columns;
  std::ifstream file_dat, file_bin;

  std::cout << std::endl << "ALaDyn output reader v" << MAJOR_RELEASE << "." << MINOR_RELEASE << "." << BUGFIX_RELEASE << std::endl;

  if (argc > 1)
  {
    if (argv[1] == "-help" || argv[1] == "/help" || argv[1] == "-h" || argv[1] == "/h")
    {
      params.man(argv[0]);
      exit(-1);
    }
  }

  if (argc == 1)
  {
    std::cout << "File basename (without extension): ";
    std::cin >> params.filebasename;
  }
  else params.filebasename = std::string(argv[1]);


  bin_filename << params.filebasename << ".bin";
  dat_filename << params.filebasename << ".dat";

  file_bin.open(bin_filename.str().c_str(), std::ios::binary | std::ios::in);
  if (file_bin.fail())
  {
    bin_filename.str("");
    bin_filename << params.filebasename << "_000.bin";
    file_bin.open(bin_filename.str().c_str(), std::ios::binary | std::ios::in);
    if (file_bin.fail())
    {
      std::cout << "Input file not found" << std::endl;
      exit(-2);
    }
    else params.multifile = true;
  }
  else params.multifile = false;
  if (params.multifile) std::cout << "Input files are " << params.filebasename << "_???.bin" << std::endl;
  else std::cout << "Input file is " << params.filebasename << ".bin" << std::endl;
  params.check_filename(params.filebasename.c_str());
  file_bin.close();
  params.check_forced_version(argc, argv);

  file_dat.open(dat_filename.str().c_str());
  if (file_dat.fail())
  {
    params.file_version = 1;
    std::cout << "Unable to find " << params.filebasename << ".dat, using routines for aladyn v1" << std::endl;

    if (params.we_dont_know_if_we_have_to_do_swap) params.ask_file_endianness();

    params.read_params_from_bin_file(params.filebasename.c_str());

    if (params.we_dont_know_if_sim_is_2d) params.ask_file_dims();
  }
  else
  {
    std::cout << "Found " << params.filebasename << ".dat, using new routines" << std::endl;
    params.read_params_from_dat_file(file_dat);
  }
  file_dat.close();

  params.check_swap();
  params.parse_command_line();
  testParameters = params.check_params();

  if (testParameters == false)
  {
    std::cout << "Incoherent parameters!" << std::endl;
    exit(-3);
  }

  if (params.grid_file)            read_grid_file(&params);
  else if (params.phasespace_file) read_phase_space_file(&params);
  else {
    std::cout << "Unknown file name, I don't know what to do with it! Rename it or change tool!" << std::endl;
    exit(-4);
  }

  return 0;
}


