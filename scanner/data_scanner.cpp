/******************************************************************************
Copyright 2014 Stefano Sinigardi
The program is distributed under the terms of the GNU General Public License
******************************************************************************/

/**************************************************************************

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

**************************************************************************/


#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 869)
#pragma warning(disable : 981)
#endif

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdint>
#include <string>
#include <cstring>
#include <iomanip>
#include <string>

#define RESERVE_SIZE_VECTOR 10000	
#define LINE_MAX_LENGTH     1024
#define EPSILON             1.0e-3




bool doubleEquality(double a, double b)
{
  return fabs(a - b) < EPSILON;
}


void seleziona(char* file_da_leggere, char* file_da_scrivere, int colonna_analizzata, double valore_riferimento)
{
  std::ifstream da_leggere, da_selezionare;
  std::ofstream da_scrivere;
  bool fallita_lettura_da_leggere = true, fallita_apertura_da_scrivere = true;

  da_leggere.open(file_da_leggere);
  fallita_lettura_da_leggere = da_leggere.fail();
  std::cout << "File in: " << std::string(file_da_leggere) << std::endl;

  if ((fallita_lettura_da_leggere))
  {
    std::cout << "Unable to open input file" << std::endl;
    return;
  }

  int contacolonne_da_leggere = 0;
  char str[LINE_MAX_LENGTH], backup_str[LINE_MAX_LENGTH];
  char * pch;
  da_leggere.getline(str, LINE_MAX_LENGTH);
  strcpy(backup_str, str);
  // non serve riportare lo stream all'inizio, anzi cosi' non dobbiamo piu' scartare la prima riga di commento!
  // da_leggere.clear();
  // da_leggere.seekg(0, std::ios::beg);

  pch = std::strtok(str, "\t");
  if (pch != NULL) contacolonne_da_leggere++;
  while (pch != NULL)
  {
    pch = std::strtok(NULL, "\t");
    if (pch != NULL) contacolonne_da_leggere++;
  }

  printf("%d columns were found in the input file\n", contacolonne_da_leggere);

  printf("Considering %d-th column as the one used to select data\n", colonna_analizzata);

  if (contacolonne_da_leggere <= colonna_analizzata)
  {
    printf("Not enough columns in the input file\n");
    return;
  }


  da_scrivere.open(file_da_scrivere);
  fallita_apertura_da_scrivere = da_scrivere.fail();
  std::cout << "File out: " << std::string(file_da_scrivere) << std::endl;

  if ((fallita_apertura_da_scrivere))
  {
    std::cout << "Unable to open output file" << std::endl;
    return;
  }


  //	size_t size_vector = RESERVE_SIZE_VECTOR;
  //	int multiplo = 1;

  std::vector<double> in_lettura;
  in_lettura.resize(contacolonne_da_leggere);
  std::vector< std::vector< double > > righe_salvate;

  //	righe_salvate.reserve(size_vector*multiplo);


  while (1)
  {
    for (std::vector<double>::size_type i = 0; i < in_lettura.size(); i++) da_leggere >> in_lettura[i];

    if (da_leggere.eof()) break;

    if (doubleEquality(in_lettura[colonna_analizzata - 1], valore_riferimento)) righe_salvate.push_back(in_lettura);
  }

  std::cout << righe_salvate.size() << " useful lines found" << std::endl;

  da_scrivere << backup_str << std::endl;

  for (std::vector<double>::size_type i = 0; i < righe_salvate.size(); i++)
  {
    for (std::vector<double>::size_type j = 0; j < in_lettura.size(); j++) da_scrivere << righe_salvate.at(i).at(j) << "\t";
    da_scrivere << std::endl;
  }


  da_leggere.close();
  da_scrivere.close();

}




int main(int argc, char *argv[])
{
  if (argc < 4)
  {
    std::cout << "Run as: ./a.out -in input_file -out output_file [-select column_number reference_value] [-average column_number]" << std::endl;
    return -254;
  }

  int file_da_leggere = -1, file_da_scrivere = -1;
  int colonna_da_selezionare = -1;
  //int colonna_da_mediare = -1;
  double valore_riferimento;
  bool do_select = false;
  //bool do_average = false;

  for (int i = 1; i < argc; i++)	// * We will iterate over argv[] to get the parameters stored inside.
  {								// * Note that we're starting on 1 because we don't need to know the path of the program, which is stored in argv[0]
    if (std::string(argv[i]) == "-in")
    {
      file_da_leggere = (i + 1);
      i++;		    			// so that we skip in the for cycle the parsing of the <da_leggere> file.
    }
    else if (std::string(argv[i]) == "-out")
    {
      file_da_scrivere = i + 1;
      i++;							// so that we skip in the for cycle the parsing of the <da_scrivere> file.
    }
    else if (std::string(argv[i]) == "-select")
    {
      do_select = true;
      colonna_da_selezionare = atoi(argv[i + 1]);
      valore_riferimento = atof(argv[i + 2]);
      i = i + 2;							// so that we skip in the for cycle the parsing of the selection properties.
    }
    //else if (std::string(argv[i]) == "-average")
    //{
    //  do_average = true;
    //  colonna_da_mediare = atoi(argv[i + 1]);
    //  i++;							// so that we skip in the for cycle the parsing of the average properties.
    //}
    else
    {
      std::cout << "Invalid argument: " << argv[i] << std::endl;
      return -243;
    }
  }

  //if (file_da_leggere < 1 || file_da_scrivere < 1 || (colonna_da_selezionare < 1 && do_select) || (colonna_da_mediare < 1 && do_average))
  if (file_da_leggere < 1 || file_da_scrivere < 1 || (colonna_da_selezionare < 1 && do_select))
  {
    printf("Something went wrong in the command line\n");
    return 200;
  }

  if (do_select) seleziona(argv[file_da_leggere], argv[file_da_scrivere], colonna_da_selezionare, valore_riferimento);
  //  if (do_average) media(argv[file_da_leggere], argv[file_da_scrivere], colonna_da_mediare); // bisogna in realta' passare anche altri dati, per ora non e' implementata perche' non interessante

  std::cout << "Done!" << std::endl;

  return 0;
}

