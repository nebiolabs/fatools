#!/usr/bin/python
#


class MyPeak():

    def __init__(self, dye, size_s, size_bp, area_s, area_bp, height, area_bp_corr):

        self.dye     = dye
        self.size_s  = size_s
        self.size_bp = size_bp
        self.area_s  = area_s
        self.area_bp = area_bp
        self.height  = height
        self.area_bp_corr = area_bp_corr
        
    def __repr__(self):
        return "<P: dye %3s | size_s %3d | size_bp %4d | area_s %3d | area_bp %4d | height %5d >" % (
            self.dye, self.size_s, self.size_bp, self.area_s, self.area_bp, self.height )



def get_peaks(filename):

    print("in get_peaks, filename: ", filename)
    
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

    fatools_type = ("Type" in names)
    if fatools_type:
        ind['Type'] = names.index('Type')+1
        ind['Corrected Area in BP'] = names.index('Corrected Area in BP')+1
    else:
        ind['Corrected Area in BP'] = -1        
    # now find the rest
    selected_names = ['Sample File Name', 'Size', 'Height', 'Area in Point', 'Area in BP',
                      'Data Point']

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
        area_bp_corr  = float(items[ind['Corrected Area in BP']]) if fatools_type else -1
        size_s   = int(items[ind['Data Point']])
        dye      = items[ind['Dye/Sample Peak']][1]
        type     = items[ind['Type']].strip() if fatools_type else 'broad'
        
        if type != 'broad': continue
        #if fatools_type and items[ind['Type']].strip() != 'broad':
        #    continue

        peak = MyPeak(dye, size_s, size_bp, area_s, area_bp, height, area_bp_corr)
        peaks.append(peak)
        
    if last_file != "" and len(peaks)>0:
        all_peaks[last_file] = peaks

    return all_peaks
