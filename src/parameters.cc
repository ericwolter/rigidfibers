/*
 *  parameters.cc - provides services to parse the parameters and layout files
 *
 *  Copyright (C) 2014  Eric Wolter <eric.wolter@gmx.de>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
#include "parameters.h"

#include <sstream>
#include <limits>
#include <iomanip>

const Configuration Parameters::parseConfigurationFiles(const std::string layout_filename) {
    Configuration configuration;

    // we parse the layout file first so that we can insert the number of fibers
    // later directly into the parameters
    float4 *initial_positions = NULL;
    float4 *initial_orientations = NULL;
    int number_of_fibers;
    Parameters::parseInitialLayoutFile(layout_filename, &initial_positions, &initial_orientations, &number_of_fibers);

    configuration.initial_positions = initial_positions; 
    configuration.initial_orientations = initial_orientations;

    return configuration;
}

void Parameters::parseInitialLayoutFile(const std::string layout_filename, float4** initialPositions, float4** initialOrientations, int *number_of_fibers)
{
    std::cout << "Parsing layout file: " << layout_filename << std::endl;

    std::ifstream layout_file_stream;
    layout_file_stream.open(layout_filename.c_str());
    if (!layout_file_stream)
    {
        std::cerr << "Could not open layout file: " << layout_filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::getline(layout_file_stream, line);

    // @todo currently only old file format is used
    // // check for version in first line
    // if (line.find("#!version 2.0"))
    // {
    //     // this is a file format version 2.0
    // }
    // else
    // {

    // older file formats didn't specifc a version, so if we don't find it
    // in the first line we assume an old format
    // because the old format might already contain valid parameters in its
    // first line we backtrack to the beginning of the stream and let
    // the specialist file format take care of the file as a whole
    layout_file_stream.seekg(0, layout_file_stream.beg);
    Parameters::parseVersion1LayoutFile(layout_file_stream, initialPositions, initialOrientations, number_of_fibers);
}

void Parameters::parseVersion1LayoutFile(std::ifstream &layout_file_stream, float4 **initialPositions, float4 **initialOrientations, int *number_of_fibers)
{
    std::cout << "...detected file format version 1" << std::endl;

    // the first line contains the number of fibers to follow
    // after that for each fiber there are two lines, the first is the position
    // and the second the orientation

    std::string line;
    std::getline(layout_file_stream, line);

    std::istringstream parse_number_of_fibers(line);

    parse_number_of_fibers >> *number_of_fibers;

    *initialPositions = new float4[*number_of_fibers];
    *initialOrientations = new float4[*number_of_fibers];

    for (int fiber_index = 0; fiber_index < *number_of_fibers; ++fiber_index)
    {
        std::getline(layout_file_stream, line);
        float4 position;

        std::istringstream positionValues(line);
        positionValues >> position.x;
        positionValues >> position.y;
        positionValues >> position.z;
        position.w = 0;

        std::getline(layout_file_stream, line);
        float4 orientation;

        std::istringstream orientationValues(line);
        orientationValues >> orientation.x;
        orientationValues >> orientation.y;
        orientationValues >> orientation.z;
        orientation.w = 0;

        (*initialPositions)[fiber_index] = position;
        (*initialOrientations)[fiber_index] = orientation;
    }
}
