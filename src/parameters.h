#ifndef FIBERS_PARAMETERS_H_
#define FIBERS_PARAMETERS_H_
/*
 *  parameters.h - header for parameters.cc
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
#include "common.h"
#include <fstream>

// Simple container for the setup configuration which allows the simulation to
// set itself up correctly
typedef struct {
    float4 *initial_positions;
    float4 *initial_orientations;
} Configuration;

class Parameters
{
public:
    static const Configuration parseConfigurationFiles(const std::string layout_filename);
    static void parseInitialLayoutFile(const std::string layout_filename, float4** initialPositions, float4** initialOrientations, int *number_of_fibers);

private:
    static void parseVersion1LayoutFile(std::ifstream &layout_file_stream, float4** initialPositions, float4** initialOrientations, int *number_of_fibers);
};

#endif // FIBERS_PARAMETERS_H_
