/*
 *  cldevice.cc - wrapper for cl_device
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
#include "cldevice.h"

CLDevice::CLDevice(cl_device_id id) {
    id_ = id;
}

CLDevice::~CLDevice() {}

cl_device_id CLDevice::id() const {
    return id_;
}

const char* CLDevice::name() const {
    return getInfo(CL_DEVICE_NAME);
}

const char* CLDevice::getInfo(cl_device_info param_name) const {
    cl_int err;

    size_t size;
    err = clGetDeviceInfo(id_, param_name, 0, NULL, &size);
    clCheckError(err, "Could get size of device info");
    char *deviceInfo = new char[size];
    err = clGetDeviceInfo(id_, param_name, size, deviceInfo, NULL);
    clCheckError(err, "Could get device info");
    return deviceInfo;
}
