#ifndef FIBERS_OCL_CLDEVICE_H_
#define FIBERS_OCL_CLDEVICE_H_
/*
 *  cldevice.h - header for cldevice.cc
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

class CLDevice
{
public:
    CLDevice(cl_device_id id);
    ~CLDevice();

    cl_device_id id() const;
    const char* name() const;
private:
    DISALLOW_COPY_AND_ASSIGN(CLDevice);

    cl_device_id id_;

    const char* getInfo(cl_device_info param_name) const;
};

#endif // FIBERS_OCL_CLDEVICE_H_
