/*
 *  resources.cc - handles loading of resources, i.e. kernel sources etc.
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
#include "resources.h"

#include <sstream>
#include <fstream>

// these includes are needed to determine the executable path on the various
// platforms
#if defined(__APPLE__)
    #include <mach-o/dyld.h>
    #include <libgen.h>
#elif defined(__WINDOWS__)
  // @todo
#else // Unix
	#include <linux/limits.h>
    #include <libgen.h>
    #include <unistd.h>
#endif

// see: http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
const std::string Resources::getExecutablePath()
{
    #if defined(__APPLE__)
        char path[PATH_MAX + 1];
        char absolute_path[PATH_MAX + 1];
        uint32_t size = sizeof(path) - 1;
        if (_NSGetExecutablePath(path, &size) == 0) {
            realpath(path, absolute_path);
        }
        
        return dirname(absolute_path);
    #elif defined(_WINDOWS)
        // @todo
    #else // Unix
		// see: https://www.securecoding.cert.org/confluence/display/seccode/POS30-C.+Use+the+readlink%28%29+function+properly
		char path[PATH_MAX + 1];
		char absolute_path[PATH_MAX + 1];
		ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
		if(len != 1) {
			path[len] = '\0';
			realpath(path, absolute_path);			
		}
		
		return dirname(absolute_path);
    #endif
}

const std::string Resources::getPathForKernel(const std::string kernel_filename)
{
    return Resources::getExecutablePath() + "/kernels/" + kernel_filename;
}

const std::string Resources::getKernelSource(const std::string kernel_filename)
{
    // Open file
    std::string kernel_fullpath = Resources::getPathForKernel(kernel_filename);
    std::ifstream ifs(kernel_fullpath.c_str());
    if ( !ifs.is_open() ) 
    {
        std::cerr << "Could not open: " << kernel_fullpath << std::endl;
        exit(EXIT_FAILURE);
    }

    // read content
    return std::string(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
}
