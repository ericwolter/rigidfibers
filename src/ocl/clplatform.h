#ifndef FIBERS_OCL_CLPLATFORM_H_
#define FIBERS_OCL_CLPLATFORM_H_
/*
 *  clplatform.h - header for clplatform.cc
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
#include "cldevice.h"

class CLPlatform
{
public:
    CLPlatform(cl_platform_id id);
    ~CLPlatform();

    const static std::vector<CLPlatform*> list();

    cl_platform_id id() const;
    const char* name() const;
    const char* vendor() const;

    const std::vector<CLDevice*> devices() const;

private:
    DISALLOW_COPY_AND_ASSIGN(CLPlatform);

    cl_platform_id id_;

    const char* getInfo(cl_platform_info param_name) const;
};

#endif // FIBERS_OCL_CLPLATFORM_H_
