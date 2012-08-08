
import os
import numpy
#import netCDF4
import sys
sys.path.append('/net2/nnz/opt/python/netCDF4-0.9.9/build/lib.linux-i686-2.4/')
import netCDF4

class ncFile(object):
    def __init__(self,path):
        self.path = path

        self.open()
        for dimname, dim in self.dims.iteritems():
            if dim.isunlimited():
                self.len = len(dim)

    def close(self):
        self.root.close()
        del self.root, self.dims, self.vars, self.gatts

    def delete(self):
        self.close()
        os.remove(self.path)

    def open(self):
        self.root  = netCDF4.Dataset(self.path,'r')
        self.dims  = self.root.dimensions
        self.vars  = self.root.variables
        self.gatts = self.root.ncattrs()
        
class newFile(object):
    def __init__(self,path,seed):
        self.path = path
        self.root = netCDF4.Dataset(self.path, 'w', format='NETCDF3_CLASSIC')
        self.vars = self.root.variables
        self.dims = self.root.dimensions

        for att in seed.root.ncattrs():
            self.root.setncattr(att,seed.root.getncattr(att))

    def close(self):
        self.root.close()
        del self.root, self.vars, self.dims


