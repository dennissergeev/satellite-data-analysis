# -*- coding: utf-8 -*-
"""
@author: D. E. Sergeev
"""
import datetime
import glob
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
from mpl_toolkits.basemap import Basemap
import numpy as np
import os
from osgeo import gdal
# My modules
import mypaths
import plot_params as pp 
from var_utils import str_input, int_input

pmc_id = int_input('PMC id')
outfname = str_input('Common output name','pmc_loc_time_ch4_20Mar-02Apr.txt')
file_name_sat = [] 
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130325*.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130321*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130322*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130323*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130324*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130325*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130326*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130327*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130328*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130329*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130330*ch4.r8.tif')
#file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130331*ch4.r8.tif')
file_name_sat += glob.glob(mypaths.satdir + os.sep + '20130402*ch4.r8.tif')
file_name_sat.sort()

lon_min = -20.
lon_max = 40.
lat_min = 65.
lat_max = 82.

class SatView(object):
    def __init__(self,file_names,out_file_name=None):
        self.file_names = file_names
        if out_file_name == None:
            self.out_file_name = 'pmc_loc_time_' + \
            os.path.basename(self.file_names[0]).split('.')[1] + '_' + \
            os.path.basename(self.file_names[0]).split('.')[0] + '-' + \
            os.path.basename(self.file_names[-1]).split('.')[0]+'.txt'
        else:
            self.out_file_name = out_file_name
        self.ind = 0
        self.track = []
        self.lonlist = []
        self.latlist = []

    def __call__(self,pmc_id,lon1=-180.,lon2=180.,lat1=-90.,lat2=90.,proj='mill',tick_incr=[1,1]):
        self.pmc_id = pmc_id
        self.fig, self.ax = plt.subplots()

        self.bm = Basemap(projection='lcc',\
                 llcrnrlon=-5,llcrnrlat=65,urcrnrlon=50.,urcrnrlat=80,\
                 lat_1=70.,lat_2=80.,lon_0=25.,\
                 resolution='l')
        #self.bm = Basemap(projection=proj,llcrnrlon=lon1,llcrnrlat=lat1, \
        #                  urcrnrlon=lon2,urcrnrlat=lat2, \
        #                  resolution='i')
        
        ticklon = np.array(tick_incr)[0]
        try:
            ticklat = np.array(tick_incr)[1]
        except IndexError:
            ticklat = ticklon

        #self.bm.drawmeridians(np.arange(round(lon1),lon2,ticklon),labels=[0,0,0,1],color='c')
        #self.bm.drawparallels(np.arange(round(lat1),lat2,ticklat),labels=[1,0,0,0],color='c') # Bug in drawparallels function
        self.bm.drawcoastlines(color='y')

        self.get_data() # get lons, lats, arr
        self.plot_data()
        self.fig.canvas.mpl_connect('button_press_event',self.on_click)
        self.cursor = Cursor(plt.gca(), useblit=True, color='red')

        pp.figview(False,maxfig=True)

    def on_click(self, event):
        if event.inaxes == self.ax:
        # Get nearest data
        #xpos = np.argmin(np.abs(event.xdata - self.x))
        #ypos = np.argmin(np.abs(event.ydata - self.y))

            if event.button == 1: # Left button
                # Convert (x, y) to (lon, lat) using Basemap
                ilon, ilat = self.bm(event.xdata, event.ydata,inverse=True)
                # Append coordinates to the lists
                self.lonlist.append(ilon)
                self.latlist.append(ilat)
                # Add a marker on the clicked spot
                self.markers.append(self.bm.plot(event.xdata,event.ydata,marker='o',mfc='r',mec='r',ms=10))
                plt.draw()
                # Write coordinates to text file
                self.write_xy()
        
            if event.button == 3: # Right button
                # Store data from the current frame
                self.track.append(dict(pmc_id=self.pmc_id,time=self.time_raw,lon=self.lonlist,lat=self.latlist))
                # Next frame
                self.ind += 1
                if self.ind < len(self.file_names):
                    self.pcolm.remove()
                    # Clear lon/lat lists
                    self.lonlist = []
                    self.latlist = []
                    # Clear the plot from the markers
                    if hasattr(self,'markers'):
                        for ms in self.markers:
                            self.ax.lines.remove(ms[0])
                    #    for ms in self.markers:
                    #        print ms
                    #        ms.lines = []
                    # Read data for the next frame and plot it
                    self.get_data() # get lons, lats, arr
                    self.plot_data()
                    plt.draw()
                else:
                    print 'Reached the end of the file list'
                    plt.close()

    def write_xy(self):
        outf = open(self.out_file_name,'a')
        outf.write(str(self.pmc_id)+\
        '\t'+self.time_raw+'\t'+'{0:.3f}\t{1:.3f}\n'.format(self.lonlist[-1],self.latlist[-1]))
        outf.close()

    def plot_data(self):
        xx, yy = self.bm(self.lons,self.lats)
        self.pcolm = self.bm.pcolormesh(xx,yy,self.arr,cmap=plt.get_cmap('gray'),alpha=1.)
        plt.title(self.time_string)
        self.markers = []

    def get_data(self):
        print('Reading: '+self.file_names[self.ind])
        self.get_time_str()
        gtif = gdal.Open(self.file_names[self.ind])
        trans = gtif.GetGeoTransform()
        self.arr = gtif.ReadAsArray()
        
        nx = gtif.RasterXSize*trans[1]
        ny = gtif.RasterYSize*trans[5]
        x = np.arange(0, gtif.RasterXSize*trans[1],trans[1])
        y = np.arange(0, gtif.RasterXSize*trans[5],trans[5])
        xx, yy = np.meshgrid(x, y)
        del gtif
        
        m = Basemap(projection='cass',lon_0=10.,lat_0=74, width=nx, height=ny)
        self.lons, self.lats = m(xx,yy,inverse=True) 

    def get_time_str(self):
        # Get file name from the path and parse it assuming '.' is a delimiter
        cur_file_name = os.path.basename(self.file_names[self.ind]).split('.')
        self.channel = cur_file_name[1]
        # Assume that the file name base consists of date and time
        cur_date_time = datetime.datetime.strptime(cur_file_name[0], '%Y%m%d%H%M')
        self.time_string = cur_date_time.strftime("%H:%M, %d %B %Y")
        self.time_raw = cur_file_name[0]

if __name__ == '__main__':
    SatView(file_name_sat,outfname)(pmc_id)
