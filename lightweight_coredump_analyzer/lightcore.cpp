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
      riga.clear();
      tokens.clear();
      std::getline(inputFile, riga);
      boost::algorithm::split(tokens, riga, boost::algorithm::is_any_of(" \t"), boost::token_compress_on);
      boost::algorithm::to_lower(tokens[0]);

      if ((tokens[0] != "+++STACK" || stackFound) && tokens[0] != "---STACK") continue;
      else if (tokens[0] == "---STACK") stackFound = false;
      else stackFound = true;

      if (stackFound)
      {
        links.push_back(tokens[1]);
      }
    }

    for (size_t i = 0; i < links.size(); i++)
    {
      std::string tempOctalLink = "0x";
      for (size_t j = 0; j < links[i].size(); j++) tempOctalLink += links[i][j];
      octalLinks.push_back(tempOctalLink);
    }

    if (octalLinks.size())
    {
      std::ofstream outputFile;
      std::string outputFileName = std::string(argv[1]) + ".txt";
      outputFile.open(outputFileName, std::ofstream::out);

      for (size_t i = 0; i < octalLinks.size(); i++)
      {
        std::string tempOctalLink = "0x";
        for (size_t j = 0; j < links[i].size(); j++) tempOctalLink += links[i][j];
        octalLinks.push_back(tempOctalLink);
      }

      for (auto i : octalLinks) outputFile << i << std::endl;
      outputFile.close();
    }
    inputFile.close();
  }

  return 0;
}