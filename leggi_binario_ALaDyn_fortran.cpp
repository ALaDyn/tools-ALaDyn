
#include "leggi_binario_ALaDyn_fortran.h"



int main (const int argc, const char *argv[])
{

	Parametri parametri;
	bool testParametri = true;

	std::cout << "Binary file reader v" << MAJOR_RELEASE << "." << MINOR_RELEASE << "." << BUGFIX_RELEASE << std::endl;

	if (argc == 1)
	{
		std::cout << "-interactive: ./reader filebasename (so do not put file extension!)" << std::endl;
		std::cout << "-batch:       ./reader filebasename -arguments" << std::endl;

		std::cout <<"----------Argument list------------------- " << std::endl;
		std::cout << "-params (write a .parameters file with params from .bin/.dat files)" << std::endl;
		std::cout << "-swap/-noswap (force endianess swap) -force_new (force new format)" << std::endl;
		std::cout << "-dump_vtk -dump_cutx #x -dump_cuty #y -dump_cutz #z  -dump_lineoutx -dump_gnuplot" << std::endl;
		std::cout << "-dump_vtk_nostretch (dumps in the vtk just the unstretched part of the grid)" << std::endl;
		std::cout << "-do_binning [REQUIRED TO ENABLE BINNING FOR PLOTTING]" << std::endl;
		std::cout << "-[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz]min/max #number" << std::endl;
		std::cout << "-plot_AB A,B={x,y,z,px,py,pz}" << std::endl;
		std::cout << "-plot_etheta -plot_ethetaT -plot_rfc -plot_espec -plot_thetaspec -plot_thetaTspec" << std::endl;
		std::cout << "-nbin #num" << std::endl;
		std::cout << "-nbin[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz] #num" << std::endl;
		std::cout << "-dontask [TRIES TO RUN IN NON-INTERACTIVE MODE]" << std::endl;
		std::cout << "Filters: \n +[x,y,z,px,py,pz,theta,thetaT,gamma,E,ty,tz]min/max #num" << std::endl;
		std::cout <<"----------Argument list------------------- " << std::endl;
		return -1;
	}

	std::ostringstream nomefile_bin, nomefile_dat;
	nomefile_bin << std::string(argv[1]) << ".bin";
	nomefile_dat << std::string(argv[1]) << ".dat";
	std::string riga_persa, endianness, columns;
	std::ifstream file_dat, file_bin;

	/* Controllo file binario */
	file_bin.open(nomefile_bin.str().c_str(),std::ios::binary|std::ios::in);
	if ( file_bin.fail() )
	{
		nomefile_bin.str("");
		nomefile_bin << std::string(argv[1]) << "_000.bin";
		file_bin.open(nomefile_bin.str().c_str(),std::ios::binary|std::ios::in);
		if ( file_bin.fail() )
		{
			std::cout << "Input file non trovato" << std::endl;
			return -3;
		}
		else parametri.multifile = true;
	}
	else parametri.multifile = false;
	if(parametri.multifile) std::cout << "Input files are " << argv[1] << "_???.bin" << std::endl;
	else std::cout << "Input file is " << argv[1] << ".bin" << std::endl;
	parametri.check_filename(argv[1]);
	file_bin.close();

	/* Controllo file ascii */
	file_dat.open(nomefile_dat.str().c_str());
	if (file_dat.fail()) 
	{
		parametri.old_fortran_bin = true;
		std::cout << "Unable to find " << argv[1] << ".dat, using old routines" << std::endl;
		parametri.chiedi_endian_file();

		if (parametri.file_particelle_P || parametri.file_particelle_E || parametri.file_particelle_HI || parametri.file_particelle_LI) 

			parametri.chiedi_numero_colonne();



		if (parametri.file_campi_Ex || parametri.file_campi_Ey || parametri.file_campi_Ez 
			|| parametri.file_campi_Bx || parametri.file_campi_By || parametri.file_campi_Bz 
			|| parametri.file_densita_elettroni || parametri.file_densita_protoni 
			|| parametri.file_densita_HI || parametri.file_densita_LI
			|| parametri.file_densita_energia_griglia_elettroni || parametri.file_densita_energia_griglia_protoni 
			|| parametri.file_densita_energia_griglia_HI || parametri.file_densita_energia_griglia_LI)

			parametri.chiedi_2Do3D();
	}
	else
	{
		std::cout << "Found " << argv[1] << ".dat, using new routines" << std::endl;
		parametri.old_fortran_bin = false;
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

	parametri.parse_command_line(argc, argv);



#ifdef ENABLE_DEBUG
	for (int i = 0; i < NPARAMETRI; i++) std::cout << "p[" << i << "] = " << parametri.p[i] << std::endl;
#endif


	testParametri = parametri.check_parametri();

	if ( testParametri == false )
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
	printf("Tipo file: file_campi_Ex? %i\n", parametri.file_campi_Ex);
	printf("Tipo file: file_campi_Ey? %i\n", parametri.file_campi_Ey);
	printf("Tipo file: file_campi_Ez? %i\n", parametri.file_campi_Ez);
	printf("Tipo file: file_campi_Bx? %i\n", parametri.file_campi_Bx);
	printf("Tipo file: file_campi_By? %i\n", parametri.file_campi_By);
	printf("Tipo file: file_campi_Bz? %i\n", parametri.file_campi_Bz);
	printf("Tipo file: file_eden? %i\n", parametri.file_densita_elettroni);
	printf("Tipo file: file_pden? %i\n", parametri.file_densita_protoni);
	printf("Tipo file: file_hiden? %i\n", parametri.file_densita_HI);
	printf("Tipo file: file_liden? %i\n", parametri.file_densita_LI);
	printf("Tipo file: file_Prpout? %i\n", parametri.file_particelle_P);
	printf("Tipo file: file_Elpout? %i\n", parametri.file_particelle_E);
	printf("Tipo file: file_Hipout? %i\n", parametri.file_particelle_HI);
	printf("Tipo file: file_Lipout? %i\n", parametri.file_particelle_LI);
	fflush(stdout);
#endif


	if (parametri.file_campi_Ex || parametri.file_campi_Ey || parametri.file_campi_Ez 
		|| parametri.file_campi_Bx || parametri.file_campi_By || parametri.file_campi_Bz 
		|| parametri.file_densita_elettroni || parametri.file_densita_protoni 
		|| parametri.file_densita_HI || parametri.file_densita_LI
		|| parametri.file_densita_energia_griglia_elettroni || parametri.file_densita_energia_griglia_protoni 
		|| parametri.file_densita_energia_griglia_HI || parametri.file_densita_energia_griglia_LI)

		leggi_campi(argc, argv, &parametri);


	else if (parametri.file_particelle_P || parametri.file_particelle_E 
		|| parametri.file_particelle_HI || parametri.file_particelle_LI) 

		leggi_particelle(argc, argv, &parametri);

	return 0;
}


