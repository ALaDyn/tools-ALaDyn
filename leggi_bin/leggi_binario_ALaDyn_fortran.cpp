
#include "leggi_binario_ALaDyn_fortran.h"



int main(const int argc, const char *argv[])
{
  Parametri parametri;
  bool testParametri = true;

  std::ostringstream nomefile_bin, nomefile_dat;
  std::string riga_persa, endianness, columns;
  std::ifstream file_dat, file_bin;

  std::cout << "ALaDyn output reader v" << MAJOR_RELEASE << "." << MINOR_RELEASE << "." << BUGFIX_RELEASE << std::endl;

  if (argc > 1) 
  {
    if (argv[1] == "-help" || argv[1] == "/help" || argv[1] == "-h" || argv[1] == "/h")
    {
      std::cout << "-interactive: ./" << argv[0] << std::endl;
      std::cout << "-batch:       ./" << argv[0] << " filebasename -arguments" << std::endl;

      std::cout << "----------Argument list------------------- " << std::endl;
      std::cout << "-params (write a .parameters file with params from .bin/.dat files)" << std::endl;
      std::cout << "-swap/-noswap (force endianess swap) -force_v2 -force_v3 (force new format)" << std::endl;
      std::cout << "-dump_vtk -dump_cutx #x -dump_cuty #y -dump_cutz #z  -dump_lineoutx -dump_gnuplot" << std::endl;
      std::cout << "-dump_vtk_nostretch (dumps in the vtk just the unstretched part of the grid)" << std::endl;
      std::cout << "(use -no_stretch_x if the grid is not stretched along x axis)" << std::endl;
      std::cout << "-dump_propaga -dump_csv -dump_clean -dump_xyzE -parameters -find_minxmax" << std::endl;
      std::cout << "-do_binning [REQUIRED TO ENABLE BINNING FOR PLOTTING]" << std::endl;
      std::cout << "-[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch]min/max #number" << std::endl;
      std::cout << "-plot_AB A,B={x,y,z,px,py,pz}" << std::endl;
      std::cout << "-plot_etheta -plot_ethetaT -plot_rfc -plot_espec -plot_thetaspec -plot_thetaTspec -plot_chspec" << std::endl;
      std::cout << "-nbin #num" << std::endl;
      std::cout << "-nbin[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch] #num" << std::endl;
      std::cout << "-dontask [TRY TO RUN IN NON-INTERACTIVE MODE]" << std::endl;
      std::cout << "Filters: \n +[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz,w,ch]min/max #num" << std::endl;
      std::cout << "----------Argument list------------------- " << std::endl;
      return -1;
    }
  }

  if (argc == 1)
  {
    std::cout << "File basename (without extension): ";
    std::cin >> parametri.filebasename;
  }
  else parametri.filebasename = std::string(argv[1]);

  nomefile_bin << parametri.filebasename << ".bin";
  nomefile_dat << parametri.filebasename << ".dat";

  /* Controllo file binario */
  file_bin.open(nomefile_bin.str().c_str(), std::ios::binary | std::ios::in);
  if (file_bin.fail())
  {
    nomefile_bin.str("");
    nomefile_bin << parametri.filebasename << "_000.bin";
    file_bin.open(nomefile_bin.str().c_str(), std::ios::binary | std::ios::in);
    if (file_bin.fail())
    {
      std::cout << "Input file non trovato" << std::endl;
      return -3;
    }
    else parametri.multifile = true;
  }
  else parametri.multifile = false;
  if (parametri.multifile) std::cout << "Input files are " << parametri.filebasename << "_???.bin" << std::endl;
  else std::cout << "Input file is " << parametri.filebasename << ".bin" << std::endl;
  parametri.check_filename(parametri.filebasename.c_str());
  file_bin.close();

  /* Controllo file ascii */
  file_dat.open(nomefile_dat.str().c_str());
  if (file_dat.fail())
  {
    parametri.aladyn_version = 1;
    std::cout << "Unable to find " << parametri.filebasename << ".dat, using routines for aladyn v1" << std::endl;

    if (parametri.p_b[SWAP]) parametri.chiedi_endian_file();

    parametri.leggi_parametri_da_file_bin(parametri.filebasename.c_str());

    if (parametri.file_spaziofasi && parametri.p_b[NCOLONNE]) parametri.chiedi_numero_colonne();
    if (parametri.file_griglia && parametri.p_b[NCOLONNE]) parametri.chiedi_2Do3D();
  }
  else
  {
    std::cout << "Found " << parametri.filebasename << ".dat, using new routines" << std::endl;
    parametri.leggi_file_dat(file_dat);
  }
  file_dat.close();


  if (parametri.endian_file == parametri.endian_machine)
  {
    parametri.p[SWAP] = 0;
    parametri.p_b[SWAP] = false;
  }
  else
  {
    parametri.p[SWAP] = 1;
    parametri.p_b[SWAP] = false;
  }

#ifdef ENABLE_DEBUG
  if (parametri.p[SWAP]) std::cout << "Swap is enabled" << std::endl;
  else std::cout << "Swap is disabled" << std::endl;
#endif

  parametri.parse_command_line(argc, argv);



#ifdef ENABLE_DEBUG
  for (int i = 0; i < NPARAMETRI; i++) std::cout << "p[" << i << "] = " << parametri.p[i] << std::endl;
#endif


  testParametri = parametri.check_parametri();

  if (testParametri == false)
  {
    std::cout << "Parametri non coerenti" << std::endl;
    return -2;
  }


  if (parametri.p[DO_BINNING])
  {
    parametri.organizza_minimi_massimi();
#ifdef ENABLE_DEBUG
    printf("Chiamata main parametri.organizza_minimi_massimi()\n");
    printf("Emin=%g    Emax=%g   dE=%g\n", parametri.minimi[8], parametri.massimi[8], parametri.dimmi_dim(8));
    fflush(stdout);
#endif
  }


#ifdef ENABLE_DEBUG
  printf("file_griglia? %i\n", parametri.file_griglia);
  printf("file_spaziofasi? %i\n", parametri.file_spaziofasi);

  printf("file_campi_Ex? %i\n", parametri.file_campi_Ex);
  printf("file_campi_Ey? %i\n", parametri.file_campi_Ey);
  printf("file_campi_Ez? %i\n", parametri.file_campi_Ez);
  printf("file_campi_Bx? %i\n", parametri.file_campi_Bx);
  printf("file_campi_By? %i\n", parametri.file_campi_By);
  printf("file_campi_Bz? %i\n", parametri.file_campi_Bz);

  printf("file_eden? %i\n", parametri.file_densita_elettroni);
  printf("file_pden? %i\n", parametri.file_densita_protoni);
  printf("file_hiden? %i\n", parametri.file_densita_HI);
  printf("file_liden? %i\n", parametri.file_densita_LI);
  printf("file_genionden? %i\n", parametri.file_densita_generic_ion);
  printf("file_driverden? %i\n", parametri.file_densita_driver);

  printf("file_Prpout? %i\n", parametri.file_particelle_P);
  printf("file_Elpout? %i\n", parametri.file_particelle_E);
  printf("file_Hipout? %i\n", parametri.file_particelle_HI);
  printf("file_Lipout? %i\n", parametri.file_particelle_LI);
  printf("file_genericIon? %i\n", parametri.file_particelle_generic_ion);

  printf("file_denen_el? %i\n", parametri.file_densita_energia_griglia_elettroni);
  printf("file_denen_pr? %i\n", parametri.file_densita_energia_griglia_protoni);
  printf("file_denen_Hi? %i\n", parametri.file_densita_energia_griglia_HI);
  printf("file_denen_Li? %i\n", parametri.file_densita_energia_griglia_LI);
  printf("file_denen_genericIon? %i\n", parametri.file_densita_energia_griglia_generic_ion);
  fflush(stdout);
#endif


  if (parametri.file_griglia)         leggi_campi(&parametri);
  else if (parametri.file_spaziofasi) leggi_particelle(&parametri);
  else printf("Unuseful file\n");
  fflush(stdout);

  return 0;
}


