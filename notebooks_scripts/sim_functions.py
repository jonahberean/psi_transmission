# -*- coding: utf-8 -*-
"""Functions for analysis of simulated data, relavant for the
 2017 PSI transmission experiment
"""

import shutil
# !pip install rootpy
from rootpy.io import root_open, DoesNotExist
from rootpy.plotting import Hist, Hist2D
from rootpy import testdata
from rootpy import asrootpy
import ROOT

# !pip install uncertainties
import uncertainties
from uncertainties import unumpy
from uncertainties import *

from rootpy.extern.six.moves import range
from rootpy.plotting import Hist, Hist2D, Hist3D, HistStack, Legend, Canvas
from rootpy.interactive import wait
import random

###############################################################################
###############################################################################

def hello_world():
    """test print function
    """    
    print('hello world!')
    
###############################################################################
###############################################################################
    
def sim_tau(config, run_type, diffuse_probability):
    """A function that fits the time constant of the detection peak's tail
    using ROOT. 
    
    Arguments:
        config {string} -- A string that determines the experimental 
            configuration of the data to be loaded. The options are:
                'NOMI' - nominal, guide-less
                'JPTI' - JP Ti guide with NiP
                'JPSU' - JP SUS guide with NiP
                'DISK' - SS Disk
                'GD01' - UGD01 guide
                'GD03' - UGD03 guide
                'EPSU' - EP SUS guide with NiP
        run_type {string} -- A string that determines what type of data will be 
            loaded. The options are:
                'shot' - direct shot measurements 
                's005' - 5 second storage
                's020' - 20 second storage
                's100' - 100 second storage
        diffuse_probability{int} -- in percentage, the diffuse reflection 
        probability value used in the simulation
    
    Returns:
        canvas -- the drawing canvas used for the PyROOT plot
    """
    # open the appropriate file
    filename = '../data_sim/' + config + '_' + run_type + \
                  '_dp' + str(diffuse_probability).zfill(2) + '.root'
    print(filename)
    f = root_open(filename)
    # print data information
    print('#####')
    print(config + ', ' + run_type + ', dp = {}'.format(diffuse_probability))
    print('#####\n')
    # get the neutronend tree
    end = f.Get('neutronend')

    # using rootpy to define the histogram
    h1 = Hist(4000, 0, 400)
    
    # a canvas for plotting
    canvas = Canvas(width=700, height=500)

    # fill a histogram with only the neutrons ending in the detector
    for evt in f.neutronend:
       if evt.solidend == 200: h1.Fill(evt.tend)
    
    # logarithmic y scale 
    canvas.SetLogy()

    # draw histogram to the canvas
    h1.Draw("hist")
    canvas.Draw()

    # fitting
    fit_start = 8 + 8.6 + int(run_type[1:4]) + 4
    fit_end   = fit_start + 20
    print('fitting range = [{}, {}]\n'.format(fit_start, fit_end))
    
#     TF1("f1","[0]*x*sin([1]*x)",-3,3);
    f1 = ROOT.TF1("m1","exp([0]-1/[1]*x)", fit_start, fit_end)
    f1.SetParameters(80,3)
    fit = h1.Fit(f1, 'SRL')
    redchi = fit.Chi2() / fit.Ndf()
#     redchi = 1
    
    # set the axis ranges for viewability
    h1.GetXaxis().SetRangeUser(0, fit_end + 20)
    
    # get the slope, or inverse tau, parameter
    slope = ufloat(f1.GetParameter(1), f1.GetParError(1))
    
    h1.GetXaxis().SetTitle('Time [s]')
    h1.GetYaxis().SetTitle('UCN Counts')
    h1.Draw()
#     canvas.Print('sim_tau.pdf')
    
    return canvas, slope, redchi