/*
 *  clplatform.cc - wrapper for cl_platform
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
#include "clplatform.h"

const std::vector<CLPlatform*> CLPlatform::list() {

    cl_int err;

    // each device vendor has its own platform so for instance on systems with
    // an Intel CPU and a Nvidia GPU this will be 2. On Mac's Apple is the sole
    // provide for OpenCL and thus this will be 1.
    cl_uint num_platforms;

    err = clGetPlatformIDs(0, NULL, &num_platforms);
    clCheckError(err, "Could not get number of platforms");

    // array to hold all OpenCL platforms available on the system
    cl_platform_id *platforms = new cl_platform_id[num_platforms];

    err = clGetPlatformIDs(num_platforms, platforms, NULL);
    clCheckError(err, "Could not get platforms");

    std::vector<CLPlatform*> result;
    for(cl_uint i = 0; i < num_platforms; i++) {
        CLPlatform *platform = new CLPlatform(platforms[i]);
        result.push_back(platform);
    }

    delete[] platforms;
    return result;
}

CLPlatform::CLPlatform(cl_platform_id id) {
    id_ = id;
}

CLPlatform::~CLPlatform() {}

cl_platform_id CLPlatform::id() const {
    return id_;
}

const char* CLPlatform::name() const {
    return getInfo(CL_PLATFORM_NAME);
}

const char* CLPlatform::vendor() const {
    return getInfo(CL_PLATFORM_VENDOR);
}

const char* CLPlatform::getInfo(cl_platform_info param_name) const {
    cl_int err;

    size_t size;
    err = clGetPlatformInfo(id_, param_name, 0, NULL, &size);
    clCheckError(err, "Could get size of platform info");
    char *platformInfo = new char[size];
    err = clGetPlatformInfo(id_, param_name, size, platformInfo, NULL);
    clCheckError(err, "Could get platform info");
    return platformInfo;
}

const std::vector<CLDevice*> CLPlatform::devices() const {
    cl_int err;

    cl_uint num_devices;
    err = clGetDeviceIDs(id_, CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);
    clCheckError(err, "Could not get number of devices");

    cl_device_id *devices = new cl_device_id[num_devices];

    err = clGetDeviceIDs(id_, CL_DEVICE_TYPE_ALL, num_devices, devices, NULL); 
    clCheckError(err, "Could not get devices");

    std::vector<CLDevice*> result;
    for(cl_uint i = 0; i < num_devices; i++) {
        CLDevice *device = new CLDevice(devices[i]);
        result.push_back(device);
    }

    delete[] devices;

    return result;    
}


