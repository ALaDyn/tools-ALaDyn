#define _SCL_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


int main(int argc, char *argv[])
{
  std::ifstream inputFile;
  inputFile.open(argv[1], std::ifstream::in);
  bool stackFound = false;
  bool stackJustFound = false;
  bool frameJustFound = false;
  int contariga = 0;

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
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(" \t"), boost::token_compress_on);

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

    std::cout << "Found " << links.size() << " links" << std::endl;

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
      outputFile.open(outputFileName, std::ofstream::out);

      for (auto i : octalLinks) outputFile << i << std::endl;
      outputFile.close();
    }
    inputFile.close();
  }

  return 0;
}