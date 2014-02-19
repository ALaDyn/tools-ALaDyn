#define LUNGHEZZA_MAX_NOMEFILE 256
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

// #define DEBUGGING
#define NEW_SPEC


#include <iostream>
#include <cstdio>
#include <sstream>
#include <vector>
#include <fstream>

using namespace std;

int main()
{
	int tipofile = -1;		// 1 = diag, 2 = spec
	int mod_id, dmodel_id, LP_ord, der_ord;
	int Z_i, A_i, iform, ibeam;
	float xmax, xmin, ymax, ymin;
	float lam0, w0x, w0y, chann_rad;
	float a0, lp_int, lp_pow;
	float targ_x1, targ_x2, n_over_nc, el_lp;
	float np1, lx1, lx3, np2, lx5;
	float ompe2, nmacro, np_over_nmacro;
	int Nx, Ny, Nz, n_cell, Nsp, Nsb;
	int iter, nst, sp_step, nvar, npvar;

	string commenti[10];
	char * nomefile;
	nomefile = new char[LUNGHEZZA_MAX_NOMEFILE];
	ostringstream nomefile_in, nomefile_out;
	printf("dai file (no extension): ");
	scanf("%s",nomefile);
	if (nomefile[0] == 'd') tipofile = 1;
	else if (nomefile[0] == 's') tipofile = 2;
	else printf("File non riconosciuto\n");
	nomefile_in << string(nomefile) << ".dat";

	ifstream infile;
	infile.open(nomefile_in.str().c_str(),ifstream::in);

	infile >> commenti[0];
	if (commenti[0].c_str()[0] == 'n')
	{
		getline(infile,commenti[0]);
		getline(infile,commenti[0]);
		getline(infile,commenti[0]);
		getline(infile,commenti[0]);
		getline(infile,commenti[0]);
	}
	else getline(infile,commenti[0]);
	infile >> mod_id >> dmodel_id >> LP_ord >> der_ord;
	getline(infile,commenti[1]);
	getline(infile,commenti[1]);
	infile >> Z_i >> A_i >> iform >> ibeam;
	getline(infile,commenti[2]);
	getline(infile,commenti[2]);
	infile >> xmax >> xmin >> ymax >> ymin;
	getline(infile,commenti[3]);
	getline(infile,commenti[3]);
	infile >> lam0 >> w0x >> w0y >> chann_rad;
	getline(infile,commenti[4]);
	getline(infile,commenti[4]);
	infile >> a0 >> lp_int >> lp_pow;
	getline(infile,commenti[5]);
	getline(infile,commenti[5]);
	infile >> targ_x1 >> targ_x2 >> n_over_nc >> el_lp;
	getline(infile,commenti[6]);
	getline(infile,commenti[6]);
	infile >> np1 >> lx1 >> lx3 >> np2 >> lx5;
	getline(infile,commenti[7]);
	getline(infile,commenti[7]);
	infile >> ompe2 >> nmacro >> np_over_nmacro;
	getline(infile,commenti[8]);
	getline(infile,commenti[8]);
	infile >> Nx >> Ny >> Nz >> n_cell >> Nsp >> Nsb;
	getline(infile,commenti[9]);
	getline(infile,commenti[9]);
	infile >> iter >> nst >> sp_step >> nvar >> npvar;
	getline(infile,commenti[0]);

#ifdef DEBUGGING
	cout << commenti[0] << endl;
	cout << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << endl;
	cout << commenti[1] << endl;
	cout << Z_i << "\t" << A_i << "\t" << iform << "\t" << ibeam << endl;
	cout << commenti[2] << endl;
	cout << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << endl;
	cout << commenti[3] << endl;
	cout << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << endl;
	cout << commenti[4] << endl;
	cout << a0 << "\t" << lp_int << "\t" << lp_pow << endl;
	cout << commenti[5] << endl;
	cout << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << endl;
	cout << commenti[6] << endl;
	cout << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << endl;
	cout << commenti[7] << endl;
	cout << ompe2 << "\t" << nmacro << "\t" << np_over_nmacro << endl;
	cout << commenti[8] << endl;
	cout << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << endl;
	cout << commenti[9] << endl;
	cout << iter << "\t" << nst << "\t" << sp_step << "\t" << nvar << "\t" << npvar << endl;
#endif

	if (tipofile == 1)
	{
		getline(infile,commenti[0]);
		getline(infile,commenti[1]);

		float * timesteps = new float[nst];

		for (int ik=0; ik < nst; ik++)
		{
			infile >> timesteps[ik];
		}

		getline(infile,commenti[2]);
		getline(infile,commenti[3]);
		getline(infile,commenti[3]);

		float * Emean_e = new float[nst];
		float * Emax_e = new float[nst];
		float * Emean_p = new float[nst];
		float * Emax_p = new float[nst];
		float * Emean_i = new float[nst];
		float * Emax_i = new float[nst];
		float * px = new float[nst];
		float * py = new float[nst];
		float * pz = new float[nst];
		float * pang = new float[nst];
		float * charge = new float[nst];
		float * Ex2 = new float[nst];
		float * Ey2 = new float[nst];
		float * Ez2 = new float[nst];
		float * Ex_max = new float[nst];
		float * Ey_max = new float[nst];
		float * Ez_max = new float[nst];


		for (int ik=0; ik < nst; ik++)
		{
			infile >> Emean_e[ik] >> Emax_e[ik] >> Emean_p[ik] >> Emax_p[ik] >> Emean_i[ik] >> Emax_i[ik];
		}

		getline(infile,commenti[4]);
		getline(infile,commenti[4]);

		for (int ik=0; ik < nst; ik++)
		{
			infile >> px[ik] >> py[ik] >> pz[ik] >> pang[ik] >> charge[ik];
		}

		getline(infile,commenti[5]);
		getline(infile,commenti[5]);
		getline(infile,commenti[6]);

		for (int ik=0; ik < nst; ik++)
		{
			infile >> Ex2[ik] >> Ey2[ik] >> Ez2[ik] >> Ex_max[ik] >> Ey_max[ik] >> Ez_max[ik];
		}

		nomefile_out << string(nomefile) << ".txt";
		ofstream outfile;
		outfile.open(nomefile_out.str().c_str(),ifstream::out);

		outfile << "# " << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << endl;

		outfile << "# " << Z_i << "\t" << A_i << "\t" << iform << "\t" << ibeam << endl;
		outfile << "# " << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << endl;
		outfile << "# " << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << endl;
		outfile << "# " << a0 << "\t" << lp_int << "\t" << lp_pow << endl;
		outfile << "# " << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << endl;
		outfile << "# " << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << endl;
		outfile << "# " << ompe2 << "\t" << nmacro << "\t" << np_over_nmacro << endl;
		outfile << "# " << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << endl;
		outfile << "# " << iter << "\t" << nst << "\t" << sp_step << "\t" << nvar << "\t" << npvar << endl;

		for (int ik=0; ik < nst; ik++)
		{
			outfile << timesteps[ik] << "\t" << Emean_e[ik] << "\t" << Emax_e[ik] << "\t" << Emean_p[ik] << "\t" << Emax_p[ik] 
			<< "\t" << Emean_i[ik] << "\t" << Emax_i[ik] << "\t" << px[ik] << "\t" << py[ik] << "\t" << pz[ik] 
			<< "\t" << pang[ik] << "\t" << charge[ik] << "\t" << Ex2[ik] << "\t" << Ey2[ik] << "\t" << Ez2[ik] 
			<< "\t" << Ex_max[ik] << "\t" << Ey_max[ik] << "\t" << Ez_max[ik] << endl;
		}

		outfile.close();
	}



	if (tipofile == 2)
	{
		int ntimestep = sp_step;

		int nbin;
		float time, emax, energy, dE, *spectrum;
		string tipo;

		getline(infile,commenti[0]);
		infile >> nbin;

		spectrum = new float[nbin];
#ifdef NEW_SPEC
		float *selected_spectrum;
		selected_spectrum = new float[nbin];
#endif
		getline(infile,commenti[1]);

		for (int i = 0; i < 3; i++)	//la procedura va ripetuta 3 volte, per le tre specie: elettroni, protoni, ioni
		{
			infile >> tipo;
			getline(infile,commenti[1]);
			for (int j = 0; j < ntimestep; j++)
			{
				getline(infile,commenti[2]);
				infile >> time >> emax;
				energy = emax / nbin;
				dE = emax / nbin;
				getline(infile,commenti[3]);
				getline(infile,commenti[3]);
				for (int ik=0; ik < nbin; ik++)
				{
					infile >> spectrum[ik];
				}
				getline(infile,commenti[4]);
#ifdef NEW_SPEC
				getline(infile,commenti[4]);
				for (int ik=0; ik < nbin; ik++)
				{
					infile >> selected_spectrum[ik];
				}
				getline(infile,commenti[4]);
#endif

				nomefile_out.str("\0");
				nomefile_out.seekp(0, std::ios::beg);
				nomefile_out << string(nomefile) << "_" << tipo << "_" << time << ".txt";
				ofstream outfile;
				outfile.open(nomefile_out.str().c_str(),ifstream::out);


				outfile << "# " << mod_id << "\t" << dmodel_id << "\t" << LP_ord << "\t" << der_ord << endl;
				outfile << "# " << Z_i << "\t" << A_i << "\t" << iform << "\t" << ibeam << endl;
				outfile << "# " << xmax << "\t" << xmin << "\t" << ymax << "\t" << ymin << endl;
				outfile << "# " << lam0 << "\t" << w0x << "\t" << w0y << "\t" << chann_rad << endl;
				outfile << "# " << a0 << "\t" << lp_int << "\t" << lp_pow << endl;
				outfile << "# " << targ_x1 << "\t" << targ_x2 << "\t" << n_over_nc << "\t" << el_lp << endl;
				outfile << "# " << np1 << "\t" << lx1 << "\t" << lx3 << "\t" << np2 << "\t" << lx5 << endl;
				outfile << "# " << ompe2 << "\t" << nmacro << "\t" << np_over_nmacro << endl;
				outfile << "# " << Nx << "\t" << Ny << "\t" << Nz << "\t" << n_cell << "\t" << Nsp << "\t" << Nsb << endl;
				outfile << "# " << iter << "\t" << nst << "\t" << sp_step << "\t" << nvar << "\t" << npvar << endl;

				for (int ik=0; ik < nbin; ik++)
				{
#ifdef NEW_SPEC
					outfile << energy << "\t" << spectrum[ik] << "\t" << selected_spectrum[ik] << endl;
#else
					outfile << energy << "\t" << spectrum[ik] << endl;
#endif
					energy += dE;
				}

				outfile.close();
			}
		}
	}

	infile.close();

}
