#define LUNGHEZZA_MAX_NOMEFILE 256
#define BUFFER_SIZE 1024

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstring>

using namespace std;


int main (int argc, char *argv[])
{
	string commenti;
	string nomefile_in="part.dist.dat";
	string nomefile_out="part.dist.txt";

	bool fallita_lettura_file_in = true, fallita_apertura_file_out = true;
	istringstream input_buffer;
	char buffer[BUFFER_SIZE+1];
	ifstream input_file;
	ofstream output_file;
	input_file.open(nomefile_in.c_str(),ifstream::in);
	fallita_lettura_file_in = input_file.fail();
	output_file.open(nomefile_out.c_str(),ofstream::out);
	fallita_apertura_file_out = output_file.fail();


	if ((fallita_lettura_file_in || fallita_apertura_file_out)) 
	{
		cerr << "Unable to open input files" << endl;
		exit(253);
	}

	// lettura di tutto l'header fino alla stringa loc_ymin
	while(1)
	{
		input_file.getline(buffer, BUFFER_SIZE);
		if(strstr(buffer, "loc_ymin")) break;

		if(input_file.eof())
		{
			cout << "loc_ymin non trovato, file difettoso" << endl;
			exit(250);
		}
	}

	// salta tutti i punti di griglia alla ricerca di loc_ymax
	while(1)
	{
		input_file.getline(buffer, BUFFER_SIZE);
		if(strstr(buffer, "loc_ymax")) break;

		if(input_file.eof())
		{
			cout << "loc_ymax non trovato, file difettoso" << endl;
			exit(250);
		}
	}

	// salta tutti i punti di griglia alla ricerca di loc_ymax
	while(1)
	{
		input_file.getline(buffer, BUFFER_SIZE);
		if(strstr(buffer, "loc_zmin")) break;

		if(input_file.eof())
		{
			cout << "loc_zmin non trovato, file difettoso" << endl;
			exit(250);
		}
	}

	// salta tutti i punti di griglia alla ricerca di loc_zmax
	while(1)
	{
		input_file.getline(buffer, BUFFER_SIZE);
		if(strstr(buffer, "loc_zmax")) break;

		if(input_file.eof())
		{
			cout << "loc_zmax non trovato, file difettoso" << endl;
			exit(250);
		}
	}

	// salta tutti i punti di griglia alla ricerca di Node P-distribution
	while(1)
	{
		input_file.getline(buffer, BUFFER_SIZE);
		if(strstr(buffer, "Node P-distribution")) break;

		if(input_file.eof())
		{
			cout << "Node P-distribution non trovato, file difettoso" << endl;
			exit(250);
		}
	}

	int contaprocessori = 0;
	int contastep = 0;
	int contarighe_perstep = 0;
	int contarighe_primostep = 0;
	vector <int> vettore_numeri_particelle;
	int numero_particelle;


	//===========================
	//Lettura primo step
	//===========================
	// lettura prima riga ex: "at iter = 0 dcmp = 0"
	input_file.getline(buffer, BUFFER_SIZE);
	if(input_file.eof())
	{
		cout << "Impossibile trovare uno step" << endl;
		return -2;
	}
	// lettura seconda riga di ogni step ex: "npe_zloc = 0"
	input_file.getline(buffer, BUFFER_SIZE);

	while(!strstr(buffer, "at"))
	{
		input_file.getline(buffer, BUFFER_SIZE);
		if (!strstr(buffer, "at") && !strstr(buffer, "npe_zloc=")) 
		{
			contarighe_primostep++;
			stringstream ss(buffer);
			while (!ss.eof())
			{
				ss >> numero_particelle;
				contaprocessori++;
				vettore_numeri_particelle.push_back(numero_particelle);
			}
		}
	}

	//	cout << "debuginfo: contacpu = " << contaprocessori << ", size vettore = " << vettore_numeri_particelle.size() << endl;
	output_file << contastep << "\t";
	for (int i = 0; i < contaprocessori; i++) output_file << vettore_numeri_particelle[i] << "\t";
	output_file << endl;


	while(!input_file.eof())
	{
		contarighe_perstep = 0;
		contaprocessori = 0;
		vettore_numeri_particelle.clear();
		if (contastep>=1) 
		{
			input_file.getline(buffer, BUFFER_SIZE);
			if (input_file.eof()) break;
		}	
		input_file.getline(buffer, BUFFER_SIZE);

		while(contarighe_perstep < contarighe_primostep)
		{
			//			cout << "debug: " << contarighe_perstep << ", " << contarighe_primostep << endl;
			input_file.getline(buffer, BUFFER_SIZE);
			if (!strstr(buffer, "npe_zloc=")) 
			{
				contarighe_perstep++;
				stringstream ss(buffer);
				while (!ss.eof())
				{
					ss >> numero_particelle;
					contaprocessori++;
					//				cout << "cpu: " << contaprocessori << ", npart: " << numero_particelle << endl;
					vettore_numeri_particelle.push_back(numero_particelle);
				}
			}
		}

		contastep++;
		//		cout << "debuginfo: contacpu = " << contaprocessori << ", size vettore = " << vettore_numeri_particelle.size() << endl;
		output_file << contastep << "\t";
		for (int i = 0; i < contaprocessori; i++) output_file << vettore_numeri_particelle[i] << "\t";
		output_file << endl;
	}

	cout << "Finita lettura, #step: " << contastep << endl;

	input_file.close();
	output_file.close();
	return 0;

}

