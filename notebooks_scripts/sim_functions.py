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

from rootpy.extern.six.moves import range
from rootpy.plotting import Hist, Hist2D, Hist3D, HistStack, Legend, Canvas
from rootpy.interactive import wait
import random

###############################################################################
###############################################################################

def hello_world():
    """test print function
    """    
    print('hello world')
    
###############################################################################
###############################################################################
    
def sim_tau():
    
    f = root_open('../data_sim/NOMI_005s_dp10.root')
    # end is a TNtupleD object, which is a type of TTree
    end = f.Get('neutronend')
    # from ROOT import gStyle
    # gStyle.SetPalette(51);
    # using rootpy to define the histogram
    h1 = Hist(4000, 0, 400)
    # TH1F h("h", "h", 64, -10, 10);
    # a = end.GetEntry(22)
    # for i in range(0, end.GetEntries()):
    #     entry = end.GetEnry(i)
    #     print(entry)
    #    auto xyzt = myntuple->GetArgs(); // Get a row
    #    if (xyzt[2] > 0) h.Fill(xyzt[0] * xyzt[1]);
    canvas = Canvas(width=700, height=500)
    # c = ROOT.TCanvas("myCanvasName","The Canvas Title",800,600)
    for evt in f.neutronend:
       if evt.solidend == 200: h1.Fill(evt.tend)

    # h1.GetYaxis().SetRangeUser(0, 50)
    h1.GetXaxis().SetRangeUser(0, 100)
    # h1.SetMarkerColorAlpha(kBlue)
    # h1.SetMarkerStyle(22)
    h1.Draw("hist")
    canvas.Draw()
    canvas.SetLogy()

    # fitting
    f1 = ROOT.TF1("m1","expo",25,80)
    f1.SetParameters(7,-0.2)
    fit = h1.Fit(f1, 'SR')
    h1.SetTitle('TR-norm: 5 second Pre-Storage, drp = 10')
    h1.GetXaxis().SetTitle('Time [s]')
    h1.GetYaxis().SetTitle('UCN Counts')
    h1.Draw()
    canvas.Print('TR_norm_5s_drp10.pdf')
    
    return 0