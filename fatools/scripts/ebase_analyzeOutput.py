#!/usr/bin/python
#

import os
import glob
import sys
from collections import namedtuple
import datetime

import matplotlib.pyplot as plt
import numpy as np

from ebase import file_var_names

#from mpl_toolkits.mplot3d import Axes3D

def to_integer(dt_time):
    return 10000*dt_time.year + 100*dt_time.month + dt_time.day

class MyPeak():

    def __init__(self, dye, size_s, size_bp, area_s, area_bp, height):

        self.dye     = dye
        self.size_s  = size_s
        self.size_bp = size_bp
        self.area_s  = area_s
        self.area_bp = area_bp
        self.height  = height
        
    def __repr__(self):
        return "<P: dye %3s | size_s %3d | size_bp %4d | area_s %3d | area_bp %4d | height %5d >" % (
            self.dye, self.size_s, self.size_bp, self.area_s, self.area_bp, self.height )


def good_ladder(peak):
    return (peak.dye=='O' and peak.height>=100.)


def good_nonladder(peak):
    return (peak.dye=='B' and peak.height>=100.)


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
                             [ 'h_fa_matched',       'height FA (matched)',    60, 0,  200000 ],
                             [ 'h_ps_matched',       'height PS (matched)',    60, 0,  200000 ],
                             [ 'rel_h_fa_matched', 'height rel. to highest peak - FA (matched)', 20, 0, 1 ],
                             [ 'rel_h_ps_matched', 'height rel. to highest peak - PS (matched)', 20, 0, 1 ] ]
    
    data = MyData()
    for h in delta_hists_info:
        data.add_array(h[0], h[1]+" (non-ladder)", h[2], h[3], h[4], key='B')
        data.add_array(h[0], h[1]+" (ladder)",     h[2], h[3], h[4], key='O')

    for h in single_alg_hists_info:
        data.add_array(h[0], h[1]+" (non-ladder)", h[2], h[3], h[4], key='B')
        data.add_array(h[0], h[1]+" (ladder)",     h[2], h[3], h[4], key='O')

    for dir in dirnames:

        #print("dir: ", dir)
        os.chdir(dir)

        # get test info
        test_info = {}
        with open(dir+"_info.txt", 'r') as f:
            for line in f:
                line_info = line.split(':')
                test_info[line_info[0]] = line_info[1].strip()
        f.close()
        #print("test_info: ", test_info)
        print("excl: ", test_info['excluded'],", flag_th: ", test_info['flag_th'],
              ", op id: ", test_info['operator_id'], ", peak_th: ", test_info['peak_th'])
        
        if (test_info['operator_id']!='1182' or test_info['peak_th']!='150.0'):
            os.chdir("..")
            continue

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
            fill_data_arrays(peak_info_nl, peak_info_la, data)

        #print("output dir: ", os.getcwd())
        f_output.close()
        os.chdir('..')

    os.chdir('../scripts')
    
    histnames = [ x[0] for x in delta_hists_info ]

    data.plot_histograms(histnames, nonladder=True)
    data.plot_histograms(histnames, nonladder=False)
    
    data.plot_overlay_histograms( [ ('s_fa', 's_ps', 'peak size (s.t.u.)', 's_both.png'),
                                    ('rel_h_fa', 'rel_h_ps', 'relative height', 'rel_h_both.png') ] )

    data.plot_correlations( [ ('s_ps_matched', 's_fa_matched', 's_corr.png'),
                              ('bp_ps_matched', 'bp_fa_matched', 'bp_corr.png'),
                              ('area_s_ps_matched', 'area_s_fa_matched', 'area_s_corr.png'),
                              ('area_bp_ps_matched', 'area_bp_fa_matched', 'area_bp_corr.png'),
                              ('h_ps_matched', 'h_fa_matched', 'h_corr.png'),
                              ('rel_h_ps_matched', 'rel_h_fa_matched', 'rel_h_corr.png') ], False, False)

    
def fill_data_arrays(peak_info_nl, peak_info_la, data):

    max_ps = {}
    max_fa = {}
    max_ps['B'] = peak_info_nl.max_height_ps
    max_fa['B'] = peak_info_nl.max_height_fa
    max_ps['O'] = peak_info_la.max_height_ps
    max_fa['O'] = peak_info_la.max_height_fa

    for dye in ['B', 'O']:
        data.get_array('max_h_ps', dye).append(max_ps[dye])
        data.get_array('max_h_fa', dye).append(max_fa[dye])

    for pair in (peak_info_nl.matched_pairs + peak_info_la.matched_pairs):
        ps = pair.ps
        fa = pair.fa

        dye = ps.dye

        relheight_ps = float(ps.height)/max_ps[dye]
        relheight_fa = float(fa.height)/max_fa[dye]

        diff_size_s = fa.size_s - ps.size_s
        #if diff_size_s==-1 or diff_size_s==1:
        #    print("diff_size_s=", diff_size_s)

        data.get_array('s_fa', dye).append(fa.size_s)
        data.get_array('s_ps', dye).append(ps.size_s)

        data.get_array('s_fa_matched', dye).append(fa.size_s)
        data.get_array('s_ps_matched', dye).append(ps.size_s)

        data.get_array('bp_fa_matched', dye).append(fa.size_bp)
        data.get_array('bp_ps_matched', dye).append(ps.size_bp)

        data.get_array('area_s_fa_matched', dye).append(fa.area_s)
        data.get_array('area_s_ps_matched', dye).append(ps.area_s)

        data.get_array('area_bp_fa_matched', dye).append(fa.area_bp)
        data.get_array('area_bp_ps_matched', dye).append(ps.area_bp)

        data.get_array('h_fa_matched', dye).append(fa.height)
        data.get_array('h_ps_matched', dye).append(ps.height)

        data.get_array('rel_h_fa_matched', dye).append(relheight_fa)
        data.get_array('rel_h_ps_matched', dye).append(relheight_ps)

        data.get_array('delta_s', dye).append(fa.size_s - ps.size_s)
        data.get_array('delta_h', dye).append(fa.height  - ps.height)
        data.get_array('rel_h_ps', dye).append(relheight_ps)
        data.get_array('rel_h_fa', dye).append(relheight_fa)
        data.get_array('delta_rel_h', dye).append(relheight_fa - relheight_ps)

        if abs(fa.size_bp - ps.size_bp)<=40:
            data.get_array('delta_bp', dye).append(int(fa.size_bp - ps.size_bp))
        else:
            data.get_array('bp_fa_oor', dye).append(fa.size_bp)
            data.get_array('s_fa_oor', dye).append(fa.size_s)

        data.get_array('delta_area_s', dye).append((fa.area_s - ps.area_s)/ps.area_s if ps.area_s!=0 else -999)
        data.get_array('delta_area_bp', dye).append((fa.area_bp - ps.area_bp)/ps.area_bp if ps.area_bp!=0 else -999)

    for peak in (peak_info_nl.extra_peaks_ps):
        data.get_array('s_ps', 'B').append(peak.size_s)

    for peak in (peak_info_nl.extra_peaks_fa):
        data.get_array('s_fa', 'B').append(peak.size_s)

    for peak in (peak_info_la.extra_peaks_ps):
        data.get_array('s_ps', 'O').append(peak.size_s)

    for peak in (peak_info_la.extra_peaks_fa):
        data.get_array('s_fa', 'O').append(peak.size_s)

    
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


def get_peaks(filename):

    # get lines from output file about processed peaks
    f = open(filename)
    lines = f.readlines()
    f.close()

    all_peaks = {}
    peaks = []

    last_file=""

    # get first line and use this to see what type of output we have
    firstline = lines[0]
    names = firstline.split(',')

    ind = {}

    # first identify which item in names is Dye/Sample Peak. All items after this are shifted by one because Dye entries
    # have a comma
    try:
        ind['Dye/Sample Peak'] = names.index("Dye/Sample Peak")
    except KeyError:
        print("didn't find Dye/Sample Peak in csv file!")
        exit(5)

    # now find the rest
    selected_names = ['Sample File Name', 'Size', 'Height', 'Area in Point', 'Area in BP', 'Data Point']

    for name in selected_names:
        ind[name] = names.index(name)+1

    selected_lines = [ line for line in lines if '"O,' in line or '"B,' in line ]
    for line in selected_lines:

        items = line.split(',')

        # we assume the data lines will have one more field than the initial line... if this isn't the case,
        # stop and figure out what's going on!
        if len(items) != len(names)+1:
            print(len(items)," items in data rows, expect ",len(names)+1," items")
            exit(3)
            
        if items[ind['Sample File Name']] != last_file and last_file != "":
            all_peaks[last_file] = peaks
            peaks=[]
            last_file = items[ind['Sample File Name']]
        elif last_file=="":
            peakd=[]
            last_file = items[ind['Sample File Name']]
        
        size_bp  = float(items[ind['Size']]) if items[ind['Size']].strip()!="" else -1
        height   = int(items[ind['Height']])
        area_s   = int(items[ind['Area in Point']])
        area_bp  = float(items[ind['Area in BP']])
        size_s   = int(items[ind['Data Point']])
        dye      = items[ind['Dye/Sample Peak']][1]

        peak = MyPeak(dye, size_s, size_bp, area_s, area_bp, height)
        peaks.append(peak)

        
    if last_file != "" and len(peaks)>0:
        all_peaks[last_file] = peaks
    
    return all_peaks


class MyData():

    def __init__(self):
        self.array = {}

        
    def add_array(self, name, title, nbins, low, high, key=None):

        arr = self.MyArray(nbins, low, high, title)

        if key:
            if name in self.array.keys():
                self.array[name][key] = arr
            else:
                self.array[name] = { key : arr }
        else:
            self.array[name] = self.MyArray(nbins, low, high, title)


    def append(self, name, dye, value):
        self.array[name][dye].append(value)

        
    def get_array(self, name, dye=None):
        if dye==None:
            return self.array[name]
        else:
            return self.get_array(name)[dye]

        
    def plot_histograms(self, histnames, nonladder=True):

        dye = 'B' if nonladder else 'O'
        for histname in histnames:

            hist = self.get_array(histname, dye)
            hist.plot_hist(plt.figure().add_subplot(111)) 

            plt.legend()

            pngfilename = histname + ".png"
            plt.savefig(pngfilename)
            plt.show()

            
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
            plt.savefig(pngfilename)
            plt.show()

            
    def plot_correlations(self, names, nonladder=True, colors=False):

        dye = 'B' if nonladder else 'O'
        for histname1, histname2, pngfilename in names:

            arr1 = self.get_array(histname1, dye).array
            arr2 = self.get_array(histname2, dye).array

            if colors:
                arr3 = self.get_array('date', dye).array

                max_arr3 = max(arr3)
                
            if not colors:
                ax = plt.figure().add_subplot(111)
                ax.scatter(arr1, arr2, s=1)
            else:
                ax = plt.figure().add_subplot(111)
                s = ax.scatter(arr1, arr2, c=arr3, s=0.25)
                cb = plt.colorbar(s)
                
            ax.set_title(pngfilename[:-4])
            ax.set_xlabel(self.get_array(histname1,dye).get_title())                      
            ax.set_ylabel(self.get_array(histname2,dye).get_title())

            xy_line = (max(min(arr1),min(arr2)), min(max(arr1),max(arr2)))
            
            ax.plot(xy_line, xy_line, 'r-', label='y=x', linewidth=0.5)

            if pngfilename=='bp_corr.png':
                x = [ 15, 20, 25, 35, 50, 62, 80, 110, 120 ]
                ax.plot(x,x, 'bx', label='ladder sizes', linewidth=.5, markersize=12)
            
            #plt.legend()
            plt.savefig(pngfilename)
            plt.show()

    class MyArray():

        def __init__(self, nbins, low, high, title):

            self.title = title

            step = (high-low)/nbins
            self.bins = np.arange(low, high+2*step, step)
            self.array = []

            
        def get_title(self):
            return self.title

        
        def append(self, val):
            self.array.append(val)

            
        def array(self):
            return self.array

        
        def plot_hist(self, axes, color=None, histtype='bar'):

            #self.fill_hist()

            width = self.bins[1] - self.bins[0]
            center = (self.bins[:-1] + self.bins[1:])/2

            axes.set_xlabel(self.title)
            axes.set_ylabel("log(# entries)")
            
            #axes.set_title(self.title)
            if len(self.array)>1000:
                axes.set_yscale('log')
                axes.set_ylabel("log(# entries)")
            else:
                axes.set_ylabel("# entries")
                
                #if np.any(self.hist):
                #axes.bar(center, self.hist, align='center', width=width, fill=False,
                #         edgecolor='b', linewidth=.1)
                #axes.plot(center, self.hist, '_')
            if len(self.array)>0:
                axes.hist(self.array, self.bins, histtype=histtype, color=color,
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
