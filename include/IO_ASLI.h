/* ==========================================================================
 *  This file is part of ASLI (A Simple Lattice Infiller)
 *  Copyright (C) KU Leuven, 2019-2022
 *
 *  ASLI is free software: you can redistribute it and/or modify it under the 
 *  terms of the GNU Affero General Public License as published by the Free 
 *  Software Foundation, either version 3 of the License, or (at your option) 
 *  any later version.
 *
 *  ASLI is distributed in the hope that it will be useful, but WITHOUT ANY 
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 *  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for 
 *  more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *  Please read the terms carefully and use this copy of ASLI only if you
 *  accept them.
 * ==========================================================================*/

#ifndef IO_ASLI_H
#define IO_ASLI_H

/* alglib headers */
#include "stdafx.h"
#include "alglibmisc.h"

/* Standard library headers */
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

/* Other headers */
#include "getline_multiOS.h"

struct unitCell {
	std::string type;
	std::vector<double> weights; //Parameter value
};

struct latticeTypeData {
	alglib::real_2d_array alglib_type_kdtree;	// Coordinates + IDX
	alglib::integer_1d_array IDX;
	std::vector<std::string> labels;
	std::vector<unitCell> unitCells;
};

struct latticeTypeTree {
	latticeTypeData alglib_data; // Unit cell typing data
	alglib::kdtree kdt;          // kd-tree
};

namespace ASLI_IO {
	// Declarations
	inline bool read_csv(std::string csvFile, char separator, 
	                     std::vector<std::vector<double>> &data);
	inline bool read_tap(std::string tapFile, char separator, 
	                     latticeTypeData &typeData);

	// Definitions
	bool read_csv(std::string csvFile, char separator, std::vector<std::vector<double>> &data) {
		/* Simple csv reader for numeric data
		* 
		*/
		
		// Get in file
		std::ifstream infile(csvFile.c_str());
		if(!infile) {
			std::cerr <<  "ERROR: Unable to open " << csvFile << std::endl;
			return(EXIT_FAILURE);
		} else {
			std::cout << "  Opening " << csvFile << std::endl;
		}

		// Read data line by line
		std::string currentLine;
		while(!getline_multiOS(infile, currentLine).eof()) {
			// Create a stringstream of the current line
			std::istringstream ss(currentLine);

			std::vector<double> record;
			// Extract each column
			for (size_t i = 0; i < 4; i++) {
				std::string currentCol;

				if (!getline(ss, currentCol, separator)) break;
				record.push_back(std::stod(currentCol));
			}
			data.push_back(record);
		}
		
		return EXIT_SUCCESS;
	}

	//
	bool read_tap(std::string tapFile, char separator, latticeTypeData &typeData) {
		/* Simple .tap file reader
		* 
		*/

		std::vector<std::string> types;
		std::vector<std::vector<std::string>> data;
		
		// Get in file
		std::ifstream infile(tapFile.c_str());
		if(!infile) {
			std::cerr <<  "ERROR: Unable to open " << tapFile << std::endl;
			return(EXIT_FAILURE);
		} else {
			std::cout << "  Opening " << tapFile << std::endl;
		}

		// Read data line by line
		std::string currentLine;
		while(!getline_multiOS(infile, currentLine).eof()) {
			// Create a stringstream of the current line
			std::istringstream ss(currentLine);

			std::vector<std::string> record;
			// Extract each column
			for (size_t i = 0; i < 4; i++) {
				std::string currentCol;

				if (!getline(ss, currentCol, separator)) break;
				record.push_back(currentCol);

				if (i == 3 && std::find(types.begin(), types.end(), currentCol) == types.end()) {// Keep track of lattice types
						types.push_back(currentCol);
				}
			}
			data.push_back(record);
		}

		// Store data in structure
		typeData.labels.resize(data.size());
		typeData.alglib_type_kdtree.setlength(data.size(), 4);
		
		typeData.IDX.setlength(data.size());
		//typeData.unitCells.resize(types.size(), unitCell());
		for (size_t i = 0; i < types.size(); i++) {
			typeData.unitCells.push_back(unitCell());
		}

		for (alglib::ae_int_t irow = 0; irow < data.size(); irow++) {
			//
			typeData.labels[irow] = data[irow][3];

			//
			typeData.alglib_type_kdtree[irow][0] = std::stod(data[irow][0]);
			typeData.alglib_type_kdtree[irow][1] = std::stod(data[irow][1]);
			typeData.alglib_type_kdtree[irow][2] = std::stod(data[irow][2]);
			typeData.alglib_type_kdtree[irow][3] = irow;

			typeData.IDX[irow] = irow;

			for (size_t i = 0; i < types.size(); i++) {
				//
				if (irow == 0) {
					typeData.unitCells[i].weights.resize(data.size());
				}
				
				// Lattice type
				typeData.unitCells[i].type = types[i];
				// Set weight to 1 at prescribed locations and to 0 elsewhere
				if (i == std::distance(types.begin(), std::find(types.begin(), types.end(), data[irow][3]))) {
					typeData.unitCells[i].weights[irow] = 1;
				}	else {
					typeData.unitCells[i].weights[irow] = 0;
				}
			}
		}	
		return EXIT_SUCCESS;
	}
}
#endif