#include "loadData.h"

std::vector<std::vector<double>> loadFromTxt(const std::string &fileName)
{
    std::ifstream inputFile(fileName);
    if (!inputFile.is_open()) // Check if file is open
    {
        std::cerr << "Can't open " << fileName << std::endl;
        return {{}}; // Return empty vector
    }
    std::vector<std::vector<double>> data2D;
    std::string line;
    while (std::getline(inputFile, line))
    {
        double value;
        std::istringstream iss(line);
        std::vector<double> row;
        while (iss >> value)
        {
            row.push_back(value);
        }
        data2D.push_back(row);
    };
    inputFile.close(); // Close the file
    return data2D;
}