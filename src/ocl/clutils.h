#ifndef FIBERS_OCL_CLUTILS_H_
#define FIBERS_OCL_CLUTILS_H_
/*
 *  clutils.h - header for clutils.cc
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
#include <vector>
#include "../common.h"
#include "clplatform.h"
#include "cldevice.h"

class CLUtils
{
public:
    const static CLPlatform* selectPlatform();
    const static CLDevice* selectDevice(const CLPlatform *platform);
    static cl_context createContext(const CLPlatform *platform, const CLDevice *device);
};

#endif // FIBERS_OCL_CLUTILS_H_
