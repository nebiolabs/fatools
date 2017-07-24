#/usr/bin/python
#

from ebase import MyPeak
import os

def get_peaks(filename):

    # get lines from output file about processed peaks
    f = open(filename)
    lines = f.readlines()
    f.close()

    all_peaks = {}
    peaks = []
    
    last_file=""

    selected_lines = [ line for line in lines if '"G,' in line or '"B,' in line ]
    for line in selected_lines:

        items = line.split(',')
        if items[2] != last_file and last_file != "":
            all_peaks[last_file] = peaks
            peaks=[]
            last_file = items[2]
        elif last_file=="":
            peakd=[]
            last_file = items[2]
        
        size_bp  = float(items[3]) if items[3].strip()!="" else -1
        height   = int(items[4])
        area_s   = items[5]
        area_bp  = items[6]
        size_s   = int(items[7])
        dye      = items[0][1]

        peak = MyPeak(dye, size_s, size_bp, height)
        peaks.append(peak)

    if last_file != "" and len(peaks)>0:
        all_peaks[last_file] = peaks
    
    return all_peaks


fileroot = "trace-Jul-25-2013-Jul26-13-32-28_full"
peaks_table_file_name = "GL8-044.csv"

peaks_peakscanner = get_peaks(fileroot+"/GL8-044.csv")

peaks_fatools = get_peaks(fileroot+"/"+fileroot+".out")

def good_ladder(peak):
    return (peak.dye=='G' and peak.height>=100.)

def good_nonladder(peak):
    return (peak.dye=='B' and peak.height>=100.)

#for file in peaks_peakscanner.keys():
f = open("peak_performance_summary.out", "w")
f.write('                    \t     | \tsize_s\t   \t         | \tsize_bp\t   \t          |\trfu   \t      |\n')
f.write('          FILE      \tPEAK | \tPS    \tFA \t |DIFF|  | \tPS     \t FA\t |DIFF|   |\tPS    \tFA    |\n')
for file in peaks_fatools.keys():

    print("\nfile: ", file)
    peaks_ps = [peak for peak in peaks_peakscanner[file] if good_nonladder(peak) ]
    peaks_fa = [peak for peak in peaks_fatools[file] if good_nonladder(peak)]

    # match pairs of peaks
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
            matched_pairs.append((peak_ps, peaks_fa[best_j]))
            del peaks_fa[best_j]
        else:
            extra_peaks_ps.append(peak_ps)

    npeaks = 0
    for (ps,fa) in matched_pairs:
        npeaks += 1
        f.write("%20s\t%2i   |\t%4i\t%4i\t%2i   |\t%4i\t%4i\t%4i  |\t%5i\t%5i  |\n" %
                (file if npeaks==1 else "",npeaks,
                 ps.size_s,  fa.size_s,  abs(ps.size_s  - fa.size_s),
                 ps.size_bp, fa.size_bp, abs(ps.size_bp - fa.size_bp),
                 ps.height, fa.height))
        
    for ps in extra_peaks_ps:
        npeaks += 1
        f.write("%20s\t%2i   |\t%4i\t -- \t--   |\t%4i\t  -- \t --   |\t%5i\t  --   |\n" %
                (file if npeaks==1 else "",npeaks,ps.size_s,ps.size_bp,ps.height))

    for fa in peaks_fa:
        npeaks += 1
        f.write("%20s\t%2i   |\t -- \t%4i\t--   |\t  --\t%4i\t --   |\t  -- \t%5i  |\n" %
                (file if npeaks==1 else "",npeaks,fa.size_s,fa.size_bp,fa.height))
        
    f.write("\n")
    
f.close()
