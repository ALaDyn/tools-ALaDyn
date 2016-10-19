#include "binary_decoder.h"
#include "grid_file_decoder.h"
#include "phasespace_file_decoder.h"


int main(const int argc, const char *argv[])
{
  Parameters params;

  std::string bin_filename, dat_filename, json_filename;
  std::ifstream bin_file, dat_file, json_file;

  std::cout << "ALaDyn binary decoder v" << MAJOR_RELEASE << "." << MINOR_RELEASE << "." << BUGFIX_RELEASE << std::endl;

  if (argc > 1)
  {
    if (!(strcmp(argv[1], "-help") && strcmp(argv[1], "-h") && strcmp(argv[1], "/help") && strcmp(argv[1], "/h")))
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


  bin_filename = params.filebasename + ".bin";
  dat_filename = params.filebasename + ".dat";
  json_filename = params.filebasename + ".json";

  bin_file.open(bin_filename, std::ios::binary | std::ios::in);
  if (bin_file.fail())
  {
    bin_filename.clear();
    bin_filename = params.filebasename + "_000.bin";
    bin_file.open(bin_filename, std::ios::binary | std::ios::in);
    if (bin_file.fail())
    {
      std::cerr << "Input file not found" << std::endl;
      exit(-2);
    }
    else params.multifile = true;
  }
  else params.multifile = false;
  if (params.multifile) std::cout << "Input files are " << params.filebasename << "_???.bin" << std::endl;
  else std::cout << "Input file is " << params.filebasename << ".bin" << std::endl;
  params.check_filename(params.filebasename.c_str());
  bin_file.close();
  params.check_forced_version(argc, argv);

  dat_file.open(dat_filename);
  if (dat_file.fail())
  {
    std::cerr << "Unable to find " << params.filebasename << ".dat, using routines for aladyn v1" << std::endl;
    params.file_version = 1;

    if (params.we_dont_know_if_we_have_to_do_swap) {
      std::cerr << "Unable to determine if swap is required, please add command line flag" << std::endl;
      exit(-3);
    }

    params.read_params_from_bin_file(params.filebasename.c_str());

    if (params.we_dont_know_if_sim_is_2d) {
      std::cerr << "Unable to determine if sim is 2D or not, please add command line flag" << std::endl;
      exit(-4);
    }
  }
  else
  {
    std::cout << "Found " << params.filebasename << ".dat, using new routines" << std::endl;
    params.read_params_from_dat_file(dat_file);
  }
  dat_file.close();

  params.check_swap();
  params.parse_command_line();

  json_file.open(json_filename);
  if (json_file.fail()) {
    if (params.grid_file)            create_json_from_grid_file(&params);
    else if (params.phasespace_file) create_json_from_phasespace_file(&params);
  }
  else json_file.close();
  params.parse_json();

  if (params.grid_file)            read_grid_file(&params);
  else if (params.phasespace_file) read_phasespace_file(&params);
  else {
    std::cerr << "Unknown file name, I don't know what to do with it! Rename it or change tool!" << std::endl;
    exit(-5);
  }

  return 0;
}


