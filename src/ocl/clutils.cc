/*
 *  clutils.cc - various OpenCL utility functions
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
 #include "clutils.h"

const CLPlatform* CLUtils::selectPlatform() {

    std::vector<CLPlatform*> platforms = CLPlatform::list();

    if (platforms.size() <= 0) {
        std::cerr << "No platforms found" << std::endl;
        exit(EXIT_FAILURE);
    }

    cl_uint index;

    std::cout << "[?] What platform would you like to use?" << std::endl;
    std::vector<CLPlatform*>::const_iterator itPlatform;
    for (index = 0, itPlatform = platforms.begin(); itPlatform != platforms.end(); ++itPlatform, ++index) {
        CLPlatform *platform = *itPlatform;
        const char* platformName = platform->name();
        const char* platformVendor = platform->vendor();
        std::cout << index << ": " 
                  << platformName << " (" << platformVendor <<")" 
                  << std::endl;

        delete[] platformName;
        delete[] platformVendor;
    }

    cl_uint selectedPlatformIndex;

    std::cout << "> ";
    //std::cin >> selectedPlatformIndex;
    selectedPlatformIndex = 0;

    return platforms.at(selectedPlatformIndex);
}

const CLDevice* CLUtils::selectDevice(const CLPlatform *platform) {

    std::vector<CLDevice*> devices = platform->devices();

    if(devices.size() <= 0) {
        std::cerr << "No devices found" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "[?] Which device would you like to use?" << std::endl;
    std::vector<CLDevice*>::const_iterator itDevice;
    int index;
    
    for (index = 0, itDevice = devices.begin(); itDevice != devices.end(); ++itDevice, ++index) {
        CLDevice *device = *itDevice;
        const char* deviceName = device->name();
        std::cout << index << ": " << deviceName << std::endl;
        delete[] deviceName;
    }

    cl_uint selectedDeviceIndex;

    std::cout << "> ";
    //std::cin >> selectedDeviceIndex;
    selectedDeviceIndex = 0;

    return devices.at(selectedDeviceIndex);
}

cl_context CLUtils::createContext(const CLPlatform *platform, const CLDevice *device) {
    cl_context_properties properties[] =
    {
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform->id(),
        0
    };

    cl_int err;
    cl_device_id device_id = device->id();
    cl_context context = clCreateContext(properties, 1, &device_id, NULL, NULL, &err);
    clCheckError(err, "Could not create context");
    return context;
}
