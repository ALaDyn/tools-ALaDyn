
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>

#if defined (USE_BOOST)
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#endif

int main(int argc, char *argv[])
{
  std::ifstream inputFile;
  inputFile.open(argv[1], std::fstream::in);
  bool stackFound = false;
  bool stackJustFound = false;
  bool frameJustFound = false;
  int contariga = 0;
#if !defined (USE_BOOST)
  char * pch;
  char * line_char;
  char * token;
#endif

  std::string riga;
  std::vector<std::string> tokens;
  std::vector<std::string> links;
  std::vector<std::string> octalLinks;

  if (inputFile.fail())
  {
    std::cerr << "Warning, " << argv[1] << " file not found!" << std::endl;
    return false;
  }
  else
  {
    while (!inputFile.eof())
    {
      contariga++;
      riga.clear();
      tokens.clear();
      std::getline(inputFile, riga);
#if defined (USE_BOOST)
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(" \t"), boost::token_compress_on);
#else
      token = new char[riga.size() * sizeof(char) + 1];
      line_char = new char[riga.size() * sizeof(char) + 1];

      //      strcpy(line_char, riga.c_str()/*, riga.size() * sizeof(char)*/);
      strncpy(line_char, riga.c_str(), riga.size() * sizeof(char) + 1);

      pch = strtok(line_char, " \t");
      while (pch != NULL)
      {
        sprintf(token, "%s", pch);
        tokens.push_back(token);
        pch = strtok(NULL, " \t");
      }
      delete[] line_char;
      line_char = NULL;
      delete[] token;
      token = NULL;
#endif
      if (!tokens.size()) continue;

      if (tokens[0] == "+++STACK")
      {
        stackJustFound = true;
        continue;
      }
      else if (tokens[0] == "---STACK") stackFound = false;

      if (stackJustFound)
      {
        frameJustFound = true;
        stackJustFound = false;
        continue;
      }

      if (frameJustFound)
      {
        frameJustFound = false;
        stackFound = true;
      }

      if (stackFound) links.push_back(tokens[1]);
    }

    std::cout << "Found " << links.size() << " links" << std::endl << std::flush;

    for (size_t i = 0; i < links.size(); i++)
    {
      std::string tempOctalLink = "0x";
      for (size_t j = 8; j < links[i].size(); j++) tempOctalLink += links[i][j];
      octalLinks.push_back(tempOctalLink);
    }

    if (octalLinks.size())
    {
      std::ofstream outputFile;
      std::string outputFileName = std::string(argv[1]) + ".txt";
      outputFile.open(outputFileName.c_str(), std::fstream::out);

#if defined (USE_CPP11)
      for (auto i : octalLinks) outputFile << i << std::endl;
#else
      for (size_t i = 0; i < octalLinks.size(); i++) outputFile << octalLinks[i] << std::endl;
#endif
      outputFile.close();
    }
    inputFile.close();
  }

  return 0;
}

