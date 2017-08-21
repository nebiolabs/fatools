#!/usr/bin/python
#

import os
import glob
import sys
from collections import namedtuple
import datetime

import matplotlib.pyplot as plt
import numpy as np

from get_peaks import MyPeak, get_peaks

from ebase import file_var_names

#from mpl_toolkits.mplot3d import Axes3D

def to_integer(dt_time):
    return 10000*dt_time.year + 100*dt_time.month + dt_time.day


def good_ladder(peak):
    return (peak.dye=='O' and peak.height>=100.)


def good_nonladder(peak):
    return (peak.dye=='B' and peak.height>=100.)

showplots = False

def main():

    plt.rcParams["patch.force_edgecolor"] = True
    plt.rcParams['figure.figsize'] = 12, 9

    os.chdir("../output/")
    
    dirnames = glob.glob("*")
    
    delta_hists_info = [ [ 'delta_s',      'Delta(size_s)',      8,   -4,    4],
                         [ 'delta_bp',     'Delta(size_bp)',     8,  -40,   40],
                         [ 'delta_h',      'Delta(height)',     20, -200, 1800],
                         [ 'delta_rel_h',  'Delta(rel height)', 10, -0.2,  0.8],
                         [ 'delta_area_s', 'Relative diff(area_s)', 20,-1, 1],
                         [ 'delta_area_bp','Relative diff(area_bp)', 20,-1, 1] ]
    
    single_alg_hists_info = [[ 's_fa',     'size_s FA',     60, 0,  6000 ],
                             [ 's_ps',     'size_s PS',     60, 0,  6000 ],
                             [ 'max_h_fa', 'max height FA', 50, 0, 50000 ],
                             [ 'max_h_ps', 'max height PS', 50, 0, 50000 ],
                             [ 'rel_h_fa', 'rel height FA', 10, 0,     1 ],
                             [ 'rel_h_ps', 'rel height PS', 10, 0,     1 ],
                             [ 's_fa_oor', 'size_s FA out-of-range',100,0,5000],
                             [ 'bp_fa_oor','size_bp FA out-of-range',100,-1000,5000],
                             [ 's_fa_matched',     'size_s FA (matched)',     60, 0,  6000 ],
                             [ 's_ps_matched',     'size_s PS (matched)',     60, 0,  6000 ],
                             [ 'bp_fa_matched',    'size_bp FA (matched)',    60, -100,  6000 ],
                             [ 'bp_ps_matched',    'size_bp PS (matched)',    60, -100,  6000 ],
                             [ 'area_s_fa_matched',  'area_s FA (matched)',     60, 0,  60000 ],
                             [ 'area_s_ps_matched',  'area_s PS (matched)',     60, 0,  60000 ],
                             [ 'area_bp_fa_matched', 'area_bp FA (matched)',    60, 0,  60000 ],
                             [ 'area_bp_ps_matched', 'area_bp PS (matched)',    60, 0,  60000 ],
                             [ 'h_fa_matched',       'height FA (matched)',    60, 0,  50000 ],
                             [ 'h_ps_matched',       'height PS (matched)',    60, 0,  50000 ],
                             [ 'rel_h_fa_matched', 'height rel. to highest peak - FA (matched)', 20, 0, 1 ],
                             [ 'rel_h_ps_matched', 'height rel. to highest peak - PS (matched)', 20, 0, 1 ] ]
    
    data = MyData()
    for h in delta_hists_info:
        data.add_array(h[0], h[1]+" (non-ladder)", h[2], h[3], h[4], key='B')
        data.add_array(h[0], h[1]+" (ladder)",     h[2], h[3], h[4], key='O')

    for h in single_alg_hists_info:
        data.add_array(h[0], h[1]+" (non-ladder)", h[2], h[3], h[4], key='B')
        data.add_array(h[0], h[1]+" (ladder)",     h[2], h[3], h[4], key='O')

    for i in range(9):
        data.add_array('h_peak'+str(i),  'height of peak '+str(i+1)+' (rfu)', 50, 0, 32000, key='O')
        data.add_array('fwhm_peak'+str(i),  'FWHM of peak '+str(i+1)+' (rfu)', 20, 0, 20, key='O')
        data.add_array('area_bp_peak'+str(i), 'area of peak '+str(i+1)+' (bp units)', 50, 0, 40000, key='O')
        data.add_array('area_bp_corr_peak'+str(i), 'area of peak '+str(i+1)+' corrected (bp units)', 50, 0, 40000, key='O')

        for j in range(i+1,9):
            data.add_array('area_bp_p'+str(j)+'_p'+str(i), 'ratio of areas of peak '+str(j+1)+' to peak '+str(i+1),
                           20, 0.5,2.5, key='O')
            #data.add_array('area_bp_corr_p'+str(j)+'_p'+str(i),
            #               'ratio of corrected areas of peak '+str(j+1)+' to peak '+str(i+1),
            #               20, 0.5,2.5, key='O')
            
    for idir in range(len(dirnames)):
        dir = dirnames[idir]

        print("dir: ", dir)
        
        os.chdir(dir)

        # get test info
        test_info = {}
        with open(dir+"_info.txt", 'r') as f:
            for line in f:
                line_info = line.split(':')
                test_info[line_info[0]] = line_info[1].strip()
        f.close()
        
        csv_files = glob.glob("*.csv")
        peaks_table_file_name = csv_files[0]

        # PeakScanner peaks
        peaks_peakscanner = get_peaks(peaks_table_file_name)
        if (len(peaks_peakscanner)==0):
            sys.stderr.write("\n  something wrong with PeakScanner output in %s\n\n" % dir)

        #FATools peaks
        peaks_fatools = get_peaks(dir + ".out")
        if (len(peaks_fatools)==0):
            sys.stderr.write("\n  something wrong with FATools output in %s\n\n" %dir)

        f_output = open("peak_performance_summary.out", "w")
        f_output.write('                             | size_s         | size_bp        | area                 | area_bp                | rfu         |\n')
        f_output.write('          FILE       DYE PEAK| PS   FA   DIFF | PS   FA   DIFF | PS     FA     %DIFF  | PS      FA    %DIFF    | PS     FA   |\n')
        for fsa_file in peaks_fatools.keys():

            try:
                peaks = peaks_peakscanner[fsa_file]
            except KeyError:
                sys.stderr.write("\nfsa_file mismatch between peakscanner and fatools: \n\t%s\n\n" % fsa_file)
                continue

            #print("\ndir: ", dir, ", fsa_file: ", fsa_file)
            peaks_ps_nl = [peak for peak in peaks_peakscanner[fsa_file] if good_nonladder(peak)]
            peaks_fa_nl = [peak for peak in peaks_fatools[fsa_file] if good_nonladder(peak)]
            peaks_ps_la = [peak for peak in peaks_peakscanner[fsa_file] if good_ladder(peak)]
            peaks_fa_la = [peak for peak in peaks_fatools[fsa_file] if good_ladder(peak)]

            # get peaks and write to file
            peak_info_nl = process_peaks(f_output, fsa_file, peaks_ps_nl, peaks_fa_nl)
            peak_info_la = process_peaks(f_output, fsa_file, peaks_ps_la, peaks_fa_la)
            
            # get data
            fill_data_arrays(peak_info_nl, peak_info_la, data, idir)
            fill_area_arrays(peak_info_nl, peak_info_la, data, idir)

        #print("output dir: ", os.getcwd())
        f_output.close()
        os.chdir('..')

    os.chdir('../scripts')
    
    histnames = [ x[0] for x in delta_hists_info ]

    data.plot_histograms(histnames, nonladder=True, dirs=[])
    data.plot_histograms(histnames, nonladder=False, dirs=[])

    data.plot_histograms(['h_fa_matched'], nonladder=True, dirs=[])
    data.plot_histograms(['h_fa_matched'], nonladder=False, dirs=[])

    data.plot_overlay_histograms( [ ('s_fa', 's_ps', 'peak size (s.t.u.)', 's_both.png'),
                                    ('rel_h_fa', 'rel_h_ps', 'relative height', 'rel_h_both.png') ] )
 
    data.plot_correlations( [ ('s_fa_matched', 's_ps_matched', 's_corr.png'),
                              ('bp_fa_matched', 'bp_ps_matched', 'bp_corr.png'),
                              ('area_s_fa_matched', 'area_s_ps_matched', 'area_s_corr.png'),
                              ('area_bp_fa_matched', 'area_bp_ps_matched', 'area_bp_corr.png'),
                              ('h_fa_matched', 'h_ps_matched', 'h_corr.png'),
                              ('rel_h_fa_matched', 'rel_h_ps_matched', 'rel_h_corr.png') ], ['L','L'], False)


    data.plot_histograms(['area_bp_peak0'], nonladder=False, dirs=[])
    data.plot_histograms(['area_bp_corr_peak0'], nonladder=False, dirs=[])
    #data.plot_histograms(['fwhm_peak0'], nonladder=False, dirs=[])
    
    for i in range(1,9):
        
        data.plot_histograms(['area_bp_peak'+str(i)], nonladder=False, dirs=[])
        data.plot_histograms(['area_bp_corr_peak'+str(i)], nonladder=False, dirs=[])

        #data.plot_histograms(['fwhm_peak'+str(i)], nonladder=False, dirs=[])
        
        data.plot_correlations([ ('h_peak'+str(i), 'h_peak0', 'peak'+str(i)+'_peak0.png')], ['L','L'], False)
        data.plot_correlations([ ('area_bp_peak'+str(i), 'area_bp_peak0','area_bp_peak'+str(i)+'_v_peak0.png')], ['L','L'], False)
        data.plot_correlations([ ('area_bp_corr_peak'+str(i), 'area_bp_corr_peak0','area_bp_corr_peak'+str(i)+'_v_peak0.png')], ['L','L'], False)
        #data.plot_correlations([ ('fwhm_peak'+str(i), 'fwhm_peak0', 'fwhm_peak'+str(i)+'_peak0.png')], ['L','L'], False)
        if i>1:
            data.plot_correlations([ ('area_bp_peak'+str(i), 'area_bp_peak'+str(i-1),'area_bp_peak'+str(i)+'_peak'+str(i-1)+'.png')], ['L','L'], False)
            data.plot_correlations([ ('area_bp_corr_peak'+str(i), 'area_bp_corr_peak'+str(i-1),'area_bp_corr_peak'+str(i)+'_peak'+str(i-1)+'.png')], ['L','L'], False)

    fig = plt.figure()
    full_ax = fig.add_subplot(111)
    for i in range(9):
        for j in range(i+1,9):
            ax = fig.add_subplot(9,9,j*9+i+1, sharex = full_ax, sharey = full_ax)
            data.plot_histograms_matrix(['area_bp_p'+str(j)+'_p'+str(i)], nonladder=False, dirs=[], axes=ax)
            ax.text(1.5,4,'p'+str(j+1)+'/p'+str(i+1))            
    full_ax.set_xlabel("Peak-to-peak ratio of areas")
    full_ax.set_ylabel("# entries")
    full_ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    full_ax.spines['top'].set_color('none')
    full_ax.spines['bottom'].set_color('none')
    full_ax.spines['left'].set_color('none')
    full_ax.spines['right'].set_color('none')
    if showplots: plt.show()
    plt.savefig('plots/area_ratio_matrix.png')
    plt.close()
    
    
def fill_data_arrays(peak_info_nl, peak_info_la, data, idir):

    max_ps = {}
    max_fa = {}
    max_ps['B'] = peak_info_nl.max_height_ps
    max_fa['B'] = peak_info_nl.max_height_fa
    max_ps['O'] = peak_info_la.max_height_ps
    max_fa['O'] = peak_info_la.max_height_fa

    for dye in ['B', 'O']:
        data.get_array('max_h_ps', dye).append(idir, max_ps[dye])
        data.get_array('max_h_fa', dye).append(idir, max_fa[dye])

    for pair in (peak_info_nl.matched_pairs + peak_info_la.matched_pairs):
        ps = pair.ps
        fa = pair.fa

        dye = ps.dye

        relheight_ps = float(ps.height)/max_ps[dye]
        relheight_fa = float(fa.height)/max_fa[dye]

        diff_size_s = fa.size_s - ps.size_s
        
        data.get_array('s_fa', dye).append(idir, fa.size_s)
        data.get_array('s_ps', dye).append(idir, ps.size_s)

        data.get_array('s_fa_matched', dye).append(idir, fa.size_s)
        data.get_array('s_ps_matched', dye).append(idir, ps.size_s)

        data.get_array('bp_fa_matched', dye).append(idir, fa.size_bp)
        data.get_array('bp_ps_matched', dye).append(idir, ps.size_bp)

        data.get_array('area_s_fa_matched', dye).append(idir, fa.area_s)
        data.get_array('area_s_ps_matched', dye).append(idir, ps.area_s)

        data.get_array('area_bp_fa_matched', dye).append(idir, fa.area_bp)
        data.get_array('area_bp_ps_matched', dye).append(idir, ps.area_bp)

        data.get_array('h_fa_matched', dye).append(idir, fa.height)
        data.get_array('h_ps_matched', dye).append(idir, ps.height)

        data.get_array('rel_h_fa_matched', dye).append(idir, relheight_fa)
        data.get_array('rel_h_ps_matched', dye).append(idir, relheight_ps)

        data.get_array('delta_s', dye).append(idir, fa.size_s - ps.size_s)
        data.get_array('delta_h', dye).append(idir, fa.height  - ps.height)
        data.get_array('rel_h_ps', dye).append(idir, relheight_ps)
        data.get_array('rel_h_fa', dye).append(idir, relheight_fa)
        data.get_array('delta_rel_h', dye).append(idir, relheight_fa - relheight_ps)

        if abs(fa.size_bp - ps.size_bp)<=40:
            data.get_array('delta_bp', dye).append(idir, int(fa.size_bp - ps.size_bp))
        else:
            data.get_array('bp_fa_oor', dye).append(idir, fa.size_bp)
            data.get_array('s_fa_oor', dye).append(idir, fa.size_s)

        data.get_array('delta_area_s', dye).append(idir, (fa.area_s - ps.area_s)/ps.area_s if ps.area_s!=0 else -999)
        data.get_array('delta_area_bp', dye).append(idir, (fa.area_bp - ps.area_bp)/ps.area_bp if ps.area_bp!=0 else -999)

    for peak in (peak_info_nl.extra_peaks_ps):
        data.get_array('s_ps', 'B').append(idir, peak.size_s)

    for peak in (peak_info_nl.extra_peaks_fa):
        data.get_array('s_fa', 'B').append(idir, peak.size_s)

    for peak in (peak_info_la.extra_peaks_ps):
        data.get_array('s_ps', 'O').append(idir, peak.size_s)

    for peak in (peak_info_la.extra_peaks_fa):
        data.get_array('s_fa', 'O').append(idir, peak.size_s)

def fill_area_arrays(peak_info_nl, peak_info_la, data, idir):

    if (len(peak_info_la.matched_pairs)==9 and len(peak_info_nl.matched_pairs)>0):

        heights_la, areas_la, areas_corr_la, fwhm_la = [], [], [], []
        for pair in peak_info_la.matched_pairs:
            heights_la.append(pair.fa.height)
            areas_la.append(pair.fa.area_bp)
            areas_corr_la.append(pair.fa.area_bp_corr)
            fwhm_la.append(pair.fa.fwhm)

        for i in range(9):
            data.get_array('h_peak'+str(i), 'O').append(idir, heights_la[i])
            data.get_array('area_bp_peak'+str(i), 'O').append(idir, areas_la[i])
            data.get_array('area_bp_corr_peak'+str(i), 'O').append(idir, areas_corr_la[i])
            data.get_array('fwhm_peak'+str(i), 'O').append(idir, fwhm_la[i])

            for j in range(i+1,9):
                if areas_la[j]>0:
                    data.get_array('area_bp_p'+str(j)+'_p'+str(i),'O').append(idir, areas_la[j]/areas_la[i])
                    #data.get_array('area_bp_corr_p'+str(j)+'_p'+str(i),'O').append(idir, areas_corr_la[j]/areas_corr_la[i])

def process_peaks(f, fsa_file, peaks_ps, peaks_fa):
    
    # match pairs of peaks
    MatchedPairs = namedtuple("MatchedPair", ["ps", "fa"])
    matched_pairs = []
    extra_peaks_ps = []
    for peak_ps in peaks_ps:

        best_size_diff = 999
        best_j = -1
        for j in range(len(peaks_fa)):
            peak_fa = peaks_fa[j]

            size_diff = abs(peak_ps.size_s - peak_fa.size_s)
            if size_diff < best_size_diff:
                best_j = j
                best_size_diff = size_diff

        if best_size_diff < 5:
            matched_pairs.append(MatchedPairs(ps=peak_ps, fa=peaks_fa[best_j]))
            del peaks_fa[best_j]
        else:
            extra_peaks_ps.append(peak_ps)

    npeaks = 0
    is_ladder = False
    for pair in matched_pairs:
        ps = pair.ps
        fa = pair.fa
        
        npeaks += 1
        is_ladder = ( ps.dye=='O' )

        f.write("%23s %s %2i | %4i  %4i  %2i | %4i %4i %4i | %6i %6i %6i | %6i  %6i  %6i | %5i %5i |\n" %
                (fsa_file if npeaks == 1 and not is_ladder else "", ps.dye, npeaks,
                 ps.size_s, fa.size_s, fa.size_s - ps.size_s,
                 ps.size_bp, fa.size_bp, fa.size_bp - ps.size_bp,
                 ps.area_s, fa.area_s, (fa.area_s - ps.area_s)/ps.area_s * 100 if ps.area_s!=0 else -999,
                 ps.area_bp, fa.area_bp, (fa.area_bp - ps.area_bp)/ps.area_bp * 100 if ps.area_bp!=0 else -999,
                 ps.height, fa.height))

    for ps in extra_peaks_ps:
        npeaks += 1
        is_ladder = ( ps.dye=='O' )
        f.write("%23s %s %2i | %4i   --   -- | %4i   --   -- | %6i     --     -- | %6i      --      -- | %5i    -- |\n" %
                (fsa_file if npeaks == 1 and not is_ladder else "", ps.dye, npeaks,
                 ps.size_s,
                 ps.size_bp,
                 ps.area_s,
                 ps.area_bp,
                 ps.height))

    for fa in peaks_fa:
        npeaks += 1
        is_ladder = ( fa.dye=='O' )
        f.write("%23s %s %2i |   --  %4i  -- |   -- %4i   -- |     -- %6i     -- |    --   %6i      -- |    -- %5i |\n" %
                (fsa_file if npeaks == 1 and not is_ladder else "", fa.dye, npeaks,
                 fa.size_s,
                 fa.size_bp,
                 fa.area_s,
                 fa.area_bp,
                 fa.height))
    
    if is_ladder:
        f.write("\n")

    #if not matched_pairs:
    #    print("problem with fsa_file: ", fsa_file)

    return PeakDiffInfo(fsa_file, matched_pairs, extra_peaks_ps, peaks_fa)


class MyData():

    def __init__(self):
        self.array = {}


    def add_array(self, name, title, nbins, low, high, key=None):

        if key=='B':
            title += " (nonladder channel)"
        elif key=='O':
            title += " (ladder channel)"
            
        arr = self.MyArray(nbins, low, high, title)
            
        if key:
            if name in self.array.keys():
                self.array[name][key] = arr
            else:
                self.array[name] = { key : arr }
        else:
            self.array[name] = self.MyArray(nbins, low, high, title)


    def append(self, name, key, value):
        self.array[name][key].append(value)

        
    def get_array(self, name, dye=None):
        if dye==None:
            return self.array[name]
        else:
            return self.get_array(name)[dye]

        
    def plot_histograms(self, histnames, nonladder=True, dirs=[]):

        dye = 'B' if nonladder else 'O'
        for histname in histnames:

            hist = self.get_array(histname, dye)
            hist.plot_hist(plt.figure().add_subplot(111), dirs=dirs) 

            plt.legend()
            
            pngfilename = histname + ".png"
            plt.savefig("plots/"+pngfilename)
            if showplots: plt.show()
            plt.close()

    def plot_histograms_matrix(self, histnames, nonladder=True, dirs=[], axes=None):

        dye = 'B' if nonladder else 'O'
        for histname in histnames:

            hist = self.get_array(histname, dye)
            if not axes:
                axes = plt.figure().add_subplot(111)

            hist.plot_hist(axes, dirs=dirs, labels=False) 

            
    def plot_overlay_histograms(self, names, nonladder=True):

        dye = 'B' if nonladder else 'O'
        for histname1, histname2, xtitle, pngfilename in names:

            hist1 = self.get_array(histname1, dye)
            hist2 = self.get_array(histname2, dye)
            
            ax = plt.figure().add_subplot(111)
            hist1.plot_hist(ax, color='r', histtype='step') 
            hist2.plot_hist(ax, color='b', histtype='step') 

            ax.set_xlabel(xtitle)
            ax.set_ylabel("# entries")
            
            plt.legend()
            plt.savefig("plots/"+pngfilename)
            if showplots: plt.show()
            plt.close()

            
    def plot_correlations(self, names, types=['L','L'], colors=False):
        
        dye1 = 'B' if types[0]=='NL' else 'O'
        dye2 = 'B' if types[1]=='NL' else 'O'
        dye1_title = ' (non-ladder channel)' if types[0]=='NL' else ' (ladder channel)'
        dye2_title = ' (non-ladder channel)' if types[1]=='NL' else ' (ladder channel)'
        
        for histname1, histname2, pngfilename in names:

            arr1 = self.get_array(histname1, dye1).array_
            arr2 = self.get_array(histname2, dye2).array_
            
            fullarr1 = [ val for sublist in arr1 for val in sublist ]
            fullarr2 = [ val for sublist in arr2 for val in sublist ]

            if colors:
                
                arr3 = self.get_array('date', dye).array_
                fullarr3 = [ val for sublist in arr3 for val in sublist ]
                max_arr3 = max(fullarr3)
                
            if not colors:
                ax = plt.figure().add_subplot(111)
                ax.scatter(fullarr2, fullarr1, s=1, label='data')
            else:
                ax = plt.figure().add_subplot(111)
                s = ax.scatter(fullarr2, fullarr1, c=fullarr3, s=0.25, label='data')
                cb = plt.colorbar(s)
                
            ax.set_title(pngfilename[:-4])
            ax.set_xlabel(self.get_array(histname2,dye2).get_title()+dye2_title)                      
            ax.set_ylabel(self.get_array(histname1,dye1).get_title()+dye1_title)

            xy_line = (max(min(fullarr1),min(fullarr2)), min(max(fullarr1),max(fullarr2)))
            
            ax.plot(xy_line, xy_line, 'r--', label='y=x', linewidth=0.5)

            # fit to line
            ax.plot(np.unique(fullarr2), np.poly1d(np.polyfit(fullarr2, fullarr1, 1))(np.unique(fullarr2)),
                    label='best fit', linewidth=0.5)
            
            if pngfilename=='bp_corr.png':
                x = [ 15, 20, 25, 35, 50, 62, 80, 110, 120 ]
                ax.plot(x,x, 'bx', label='ladder sizes', linewidth=.5, markersize=12)
            
            plt.legend()
            plt.savefig("plots/"+pngfilename)
            if showplots: plt.show()
            plt.close()

    class MyArray():

        def __init__(self, nbins, low, high, title):

            self.title = title

            step = (high-low)/nbins
            self.bins = np.arange(low, high+2*step, step)
            self.array_ = []
            self.narrays = 0

            
        def get_title(self):
            return self.title

        
        def append(self, idir, val):
            while idir >= self.narrays:
                self.array_.append([])
                self.narrays += 1
            self.array_[idir].append(val)

        
        def array(self, idir):
            return self.array_[idir]


        def plot_hist(self, axes, color=None, histtype='bar', dirs=[], labels=True):

            #self.fill_hist()

            width = self.bins[1] - self.bins[0]
            center = (self.bins[:-1] + self.bins[1:])/2

            if labels:
                axes.set_xlabel(self.title)
                axes.set_ylabel("log(# entries)")

            # array contains a list for each input directory
            # for now we combine into a single list and make a histogram
            # but we could in principle select files sharing experimental setups
            # and plot one histogram for each set of files
            if not dirs:
                fullarray = [ val for sublist in self.array_ for val in sublist]
            else:
                subarray = [ self.array(idir) for idir in dirs ]
                fullarray = [ val for sublist in subarray for val in sublist ]
                
            #axes.set_title(self.title)
            if labels:
                if len(fullarray)>1000:
                    axes.set_yscale('log')
                    axes.set_ylabel("log(# entries)")
                else:
                    axes.set_ylabel("# entries")
                    
                #if np.any(self.hist):
                #axes.bar(center, self.hist, align='center', width=width, fill=False,
                #         edgecolor='b', linewidth=.1)
                #axes.plot(center, self.hist, '_')

            if len(self.array_)>0:
                axes.hist(fullarray, self.bins, histtype=histtype, color=color,
                          label= self.title, align='left')

                
class PeakDiffInfo():

    def __init__(self, fsa_file, matched_pairs, extra_peaks_ps, extra_peaks_fa):
        
        self.fsa_file = fsa_file
        self.matched_pairs = matched_pairs
        self.extra_peaks_ps = extra_peaks_ps
        self.extra_peaks_fa = extra_peaks_fa

        # get max from matched pairs
        max_ps = max(self.matched_pairs, key=lambda x:x.ps.height).ps.height\
            if self.matched_pairs else -1
        
        max_fa = max(self.matched_pairs, key=lambda x:x.fa.height).fa.height\
            if self.matched_pairs else -1

        # now max from extra peaks
        max_extra_ps = max(self.extra_peaks_ps, key=lambda x:x.height).height\
            if self.extra_peaks_ps else -1

        max_extra_fa = max(self.extra_peaks_fa, key=lambda x:x.height).height\
            if self.extra_peaks_fa else -1

        self.max_height_ps = max(max_ps, max_extra_ps)
        self.max_height_fa = max(max_fa, max_extra_fa)


if __name__ == '__main__':

    main()
