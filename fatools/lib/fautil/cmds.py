# provide commands for Fragment Analysis (FA)

from fatools.lib import params
from fatools.lib.utils import cerr, cout, cverr, cexit, tokenize, detect_buffer, set_verbosity, is_verbosity
from fatools.lib.fautil.algo2 import LadderMismatchException

import sys, argparse, csv, os
from io import StringIO


def init_argparser(parser=None):

    p = parser if parser else argparse.ArgumentParser('facmd')

    p.add_argument('--sqldb', default=False,
            help = 'SQLite3 database filename')

    p.add_argument('--file', default=None,
            help = "Comma-separated FSA filenames (optional)")

    p.add_argument('--indir', default=False,
            help = 'input directory (eg. containing FSA files)')

    p.add_argument('--fsdb', default=None,
            help = 'Filesystem-based database')

    p.add_argument('--infile', default=None,
            help = 'Tab-delimited or CSV manifest file')

    p.add_argument('--outfile', default='-',
            help = 'output filename')

    # command in ascending order

    p.add_argument('--clear', default=False, action='store_true',
            help = 'clear (or remove) all peaks from FSA file')

    p.add_argument('--align', default=False, action='store_true',
            help = 'scan ladder channel, preannotate and align with size standards')

    p.add_argument('--call', default=False, action='store_true',
            help = 'scan non-ladder channels, preannotate peaks and determine their sizes')

    p.add_argument('--bin', default=False, action='store_true',
            help = 'bin non-ladder peaks')

    p.add_argument('--annotate', default=False, action='store_true',
            help = 'annotate non-ladder peaks')

    p.add_argument('--split-plot', action='store_true',
            help='plot dye separately')

    p.add_argument('--dendogram', default=False, action='store_true',
            help = 'plot dendograms of ladders and alleles')

    p.add_argument('--ladder-plot', action='store_true',
            help='report and plot ladder alignment for assessment purposes')

    p.add_argument('--normalize', default=False, action='store_true',
            help = 'calculate normalized areas for all peaks')

    p.add_argument('--listpeaks', default=False, action='store_true',
            help = 'list all peaks')

    p.add_argument('--listrawdata', default=False, action='store_true',
            help = 'list all data (base pair units and relative frequency units)')

    p.add_argument('--peaks_format', default="standard",
                   help = "format for peaks output file (standard, peakscanner)")

    p.add_argument('--plot', default="", type=str,
            help = 'plot trace (delimited list input)')

    p.add_argument('--range', default="", type=str,
            help = 'range for plotting trace (delimited list input)')

    # semi-mandatory

    p.add_argument('--panel', default="",
            help = 'comma-separated panel code(s) (prebuilt options: GS600LIZ, GS500LIZ, GS120LIZ')

    p.add_argument('--marker', default="",
            help = 'comma-separated marker code(s)')

    p.add_argument('--panelfile', default="",
            help = 'YAML panel file')

    p.add_argument('--markerfile', default="",
            help = "YAML marker file")

    # options

    p.add_argument('--score', default=1.0, type=float,
            help = 'minimum alignment score threshold to be plotted')

    p.add_argument('--rss', default=-1, type=float,
            help = 'maximum rss threshold to be plotted')

    p.add_argument('--cluster', default=0, type=int,
            help = 'number of cluster for hierarchical clustering alignment')

    p.add_argument('--merge', default=False, action='store_true',
            help = 'merge smeared peaks into single peaks and write to output file')

    p.add_argument('--verbose', default=0, type=int,
            help = 'show verbosity')

    p.add_argument('--cache-path',
            help='store cache in other location (defaults to home)')

    p.add_argument('--no-cache', default=False, action='store_true',
            help = 'do not use caches')

    p.add_argument('--plot-file',
            help='save --plot or --split-plot result into a file')

    p.add_argument('--commit', default=False, action='store_true',
            help = 'commit to database')

    p.add_argument('--plot_merged_peaks', default=False, action='store_true',
                   help = 'make plots of merged peaks and save to output directory')

    ## Override params

    p.add_argument('--ladder_rfu_threshold', default=-1, type=float,
                   help='ladder rfu threshold')

    p.add_argument('--ladder_rfu_ratio_threshold', default=-1, type=float,
                   help='ladder rfu ratio threshold')

    p.add_argument('--nonladder_rfu_threshold', default=-1, type=float,
                   help='nonladder rfu threshold')

    p.add_argument('--nonladder_peak_window', default=-1, type=float,
                   help='size of window (in scan time units) used to calculate first derivatives')

    p.add_argument('--nonladder_smoothing_window', default=-1, type=int,
                   help='size of window (in scan time units) used to smooth non-ladder data (no smoothing if -1)')

    p.add_argument('--nonladder_smoothing_order', default=-1, type=int,
                   help='order of polynomial used to smooth non-ladder data (no smoothing if -1)')

    p.add_argument('--allelemethod', default='', type=str,
                   help='allele method (leastsquare, cubicspline, localsouthern)')

    p.add_argument('--baselinemethod', default='median', type=str,
                   help='baseline method (none, median, minimum)')

    p.add_argument('--baselinewindow', default=399, type=int,
                   help='size of running window for baseline determination (default 399)')

    return p


def main(args):

    if args.verbose != 0:
        set_verbosity(args.verbose)

    dbh = None

    # set parameter for baseline correction and allelemethod
    from fatools.lib.const import allelemethod, baselinemethod
    _params = params.Params()

    _params.baselinewindow = args.baselinewindow

    if args.baselinemethod !="":
        if args.baselinemethod=='none':
            _params.baselinemethod = baselinemethod.none
        elif args.baselinemethod=='median':
            _params.baselinemethod = baselinemethod.median
        elif args.baselinemethod=='minimum':
            _params.baselinemethod = baselinemethod.minimum
        else:
            raise NotImplementedError()

    if args.allelemethod !="":
        if args.allelemethod=='leastsquare':
            _params.allelemethod = allelemethod.leastsquare
        elif args.allelemethod=='cubicspline':
            _params.allelemethod = allelemethod.cubicspline
        elif args.allelemethod=='localsouthern':
            _params.allelemethod = allelemethod.localsouthern
        else:
            raise NotImplementedError()

    if args.nonladder_smoothing_window > 0:
        _params.nonladder.smoothing_window = args.nonladder_smoothing_window
        _params.nonladder.smoothing_order = args.nonladder_smoothing_order

    cerr('I: Aligning size standards...')
    if args.file or args.infile or args.indir:
        cverr(4, 'D: opening FSA file(s)')
        fsa_list = open_fsa(args, _params)
    elif dbh is None:
        cverr(4, 'D: connecting to database')
        dbh = get_dbhandler(args)
        fsa_list = get_fsa_list(args, dbh)

    cerr('I: obtained %d FSA' % len(fsa_list))

    if args.commit:
        with transaction.manager:
            do_facmd(args, fsa_list, dbh)
            cerr('** COMMIT to database **')
    elif dbh:
        cerr('WARNING ** running without database COMMIT! All changes will be discarded!')
        if not ( args.test or args.y ):
            keys = input('Do you want to continue [y/n]? ')
            if not keys.lower().strip().startswith('y'):
                sys.exit(1)
        do_facmds(args, fsa_list, _params, dbh)
    else:
        do_facmds(args, fsa_list, _params)


def do_facmds(args, fsa_list, _params, dbh=None):

    if args.ladder_rfu_threshold >= 0:
        _params.ladder.min_rfu = args.ladder_rfu_threshold

    if args.ladder_rfu_ratio_threshold >= 0:
        _params.ladder.min_rfu_ratio = args.ladder_rfu_ratio_threshold

    if args.nonladder_rfu_threshold >= 0:
        _params.nonladder.min_rfu = args.nonladder_rfu_threshold

    if args.nonladder_peak_window >0 :
        _params.nonladder.peakwindow = args.nonladder_peak_window


    bad_files_filename = args.indir + "/badfiles.out"

    f_bad_files = open(bad_files_filename,'w')

    args.plotrange = []

    aligned_fsa_list = fsa_list

    executed = 0
    if args.clear:
        do_clear( args, fsa_list, dbh )
        executed += 1
    if args.align:
        aligned_fsa_list = do_align( args, fsa_list, _params, f_bad_files, dbh )
        executed += 1
    if args.call:
        do_call( args, aligned_fsa_list, _params, dbh )
        executed += 1
    if args.plot or args.split_plot or args.ladder_plot:
        args.plotlist = [item for item in args.plot.split(',')]
        if args.range is not "":
            args.plotrange = [int(item) for item in args.range.split(',')]
        do_plot( args, fsa_list, dbh )
        executed += 1
    if args.normalize is not False:
        do_normalize( args, aligned_fsa_list, _params)
        executed += 1
    if args.dendogram:
        do_dendogram( args, fsa_list, dbh)
        executed += 1
    if args.merge:
        do_merge( args, aligned_fsa_list, _params )
        executed += 1
    if args.listpeaks:
        do_listpeaks( args, aligned_fsa_list, dbh )
        executed += 1
    if args.listrawdata:
        do_listrawdata( args, fsa_list, dbh )
        executed += 1
    if executed == 0:
        cerr('W: please provide a relevant command')
    else:
        cerr('I: executed %d command(s)' % executed)

    f_bad_files.close()


def do_clear( args, fsa_list, dbh ):
    pass


def do_align( args, fsa_list, _params, f_bad_files, dbh ):
    """
    This takes an input list of FSA instances , calls FSA.align for
    FSA in the list, and returns a list of good FSAs.
    """

    cerr('I: Aligning size standards...')

    good_fsa = []
    for (fsa, fsa_index) in fsa_list:
        cverr(3, 'D: aligning FSA %s' % fsa.filename)
        try:
            fsa.align(_params)
            good_fsa.append( (fsa, fsa_index) )
        except LadderMismatchException:
            f_bad_files.write(("LadderMismatch: %s\n") % fsa.filename)

    return good_fsa


def do_call( args, fsa_list, params, dbh ):

    cerr('I: Calling non-ladder peaks...')

    for (fsa, fsa_index) in fsa_list:
        cverr(3, 'D: calling FSA %s' % fsa.filename)
        fsa.call(params)


def do_merge( args, fsa_list, params ):

    cerr('I: merging smeared peaks...')

    for (fsa, fsa_index) in fsa_list:
        print("fsa_index: ", fsa_index)
        cverr(3, 'D: calling merge for FSA %s' % fsa.filename)
        fsa.merge(params, args.plot_merged_peaks)


def do_normalize( args, fsa_list, params ):

    cerr('I: Normalizing all peaks...')

    # use panel method to set scale factors for all FSA
    from fatools.lib.fileio.models import Panel
    panel = Panel.get_panel(args.panel)

    ladder_means = panel.get_ladder_area_means(fsa_list)

    # normalize areas for each FSA
    for (fsa, fsa_index) in fsa_list:
        cverr(3, 'D: calling normalize for %s' % fsa.filename)
        fsa.normalize(params, ladder_means)


def do_plot(args, fsa_list, dbh):

    cerr('I: Creating plot...')

    from fatools.lib.fautil import plot

    plot.plot(args, fsa_list, dbh)


def do_dendogram( args, fsa_list, dbh ):

    from fatools.lib.fautil import hcalign
    from matplotlib import pyplot as plt

    for (fsa, fsa_index) in fsa_list:

        c = fsa.get_ladder_channel()
        c.scan(params.Params()) # scan first if necessary

        ladder = fsa.panel.get_ladder()
        peaks = c.get_alleles(False)

        #initial_pair, P, L = hclustalign.hclust_align(peaks, ladder)
        P = hcalign.generate_tree( [ (n.rtime, 0) for n in peaks ] )
        L = hcalign.generate_tree( [ (e, 0) for e in ladder['sizes'] ] )

        clusters = hcalign.fcluster(L.z, args.cluster or ladder['k'], criterion="maxclust")
        print(clusters)

        clusters = hcalign.fcluster(P.z, args.cluster or ladder['k'], criterion="maxclust")
        print(clusters)

        plt.figure()
        plt.subplot(121)
        hcalign.dendrogram(L.z, leaf_rotation=90, leaf_font_size=8,
                labels = [ x[0] for x in L.p ])
        plt.subplot(122)
        hcalign.dendrogram(P.z, leaf_rotation=90, leaf_font_size=8,
                labels = [ x[0] for x in P.p ])
        plt.show()

def do_listpeaks( args, fsa_list, dbh ):

    if args.outfile != '-':
        out_stream = open(args.outfile, 'w')
    else:
        out_stream = sys.stdout

    if args.peaks_format=='standard':
        out_stream.write('SAMPLE\tFILENAME   \tDYE\tRTIME\tSIZE\tHEIGHT\tAREA\tSCORE\n')
    elif args.peaks_format == 'peakscanner':
        out_stream.write("Dye/Sample Peak,Sample File Name,Type,Size,Height,Area in Point,Area in BP,Corrected Area in BP,Data Point,Begin Point,")
        if args.merge:
            out_stream.write("Begin BP,End Point,End BP,Width in Point,Width in BP,Score,Peak Group,User Comments,User Edit\n")
        else:
            out_stream.write("Begin BP,End Point,End BP,Width in Point,Width in BP,Score,User Comments,User Edit\n")

    else:
        raise RuntimeError("Unknown value for args.peaks_format")
    out_stream.close()

    for (fsa, fsa_index) in fsa_list:
        cverr(3, 'D: calling FSA %s' % fsa.filename)

        markers = fsa.panel.data['markers']

        if args.outfile != '-':
            out_stream = open(args.outfile, 'a')
        else:
            out_stream = sys.stdout

        for channel in fsa.channels:
            if channel.is_ladder():
                color = markers['x/ladder']['filter']
            else:
                color = markers['x/'+channel.dye]['filter']

            alleles = channel.get_alleles(broad_peaks_only=False)

            if is_verbosity(4):
                cout('Marker => %s | %s [%d]' % (channel.marker.code, channel.dye,
                                                 len(alleles)))
                cout("channel has alleles :",len(alleles))

            i=1

            smeared_alleles = channel.smeared_alleles
            if (not args.merge) or channel.is_ladder():
                for p in alleles:
                    if args.peaks_format=='standard':
                        out_stream.write('%6s\t%10s\t%3s\t%d\t%d\t%5i\t%3.2f\t%3.2f\n' %
                                         (fsa_index, fsa.filename[:-4], color, p.rtime, p.size, p.height, p.area, p.qscore))
                    else:
                        out_stream.write('"%s, %i",%s, %s, %f, %i, %i, %i, %i, %i, %i, %f, %i, %f, %i, %f, %f,,\n' %
                                         (color, i, fsa.filename, p.type, p.size, p.height, p.area, p.area_bp, p.area_bp_corr, p.rtime, p.brtime,p.begin_bp,p.ertime,p.end_bp,p.wrtime,p.width_bp,p.qscore))
                    i = i+1


            else:
                if is_verbosity(4):
                    cout('Marker => %s | %s [%d]' % (channel.marker.code, channel.dye,
                                                     len(smeared_alleles)))
                    cout("channel has smeared alleles :",len(smeared_alleles))
                i=1
                for p in smeared_alleles:
                    out_stream.write('"%s, %i", %s, %s, %f, %i, %i, %i, %i, %i, %i, %f, %i, %f, %i, %f, %f, %i,,\n' %
                                     (color, i, fsa.filename, p.type, p.size, p.height, p.area, p.area_bp, p.area_bp_corr, p.rtime, p.brtime,p.begin_bp,p.ertime,p.end_bp,p.wrtime,p.width_bp,p.qscore, p.group))
                    i = i+1


        out_stream.close()

def do_listrawdata( args, fsa_list, dbh ):

    outfile = '-'
    if args.outfile != '-':
        print("outfile: ", args.outfile)
        outfile = args.outfile.rsplit('.',1)[0]
        outfile += "_rawdata."
        outfile += args.outfile.rsplit('.',1)[1]

        out_stream = open(outfile, 'w')
    else:
        out_stream = sys.stdout

    out_stream.write('SAMPLE NAME,WELL ID,TRACE DYE,RAW DATA\n')
    out_stream.close()

    for (fsa, fsa_index) in fsa_list:
        cverr(3, 'D: calling FSA %s' % fsa.filename)

        if outfile != '-':
            out_stream = open(outfile, 'a')
        else:
            out_stream = sys.stdout

        # sample name
        sample_name = fsa.filename.rsplit('.',1)[0]

        # iterate through channels
        markers = fsa.panel.data['markers']
        trace = fsa.get_trace()

        # get well ID
        well_id = fsa.get_well_id()

        for channel in fsa.channels:

            # get trace dye
            if channel.is_ladder():
                trace_dye = markers['x/ladder']['filter']
            else:
                trace_dye = markers['x/'+channel.dye]['filter']

            # get raw data
            data = channel.data
            datastring = "["

            data = channel.data
            basepairs = channel.get_basepairs()
            #channel.set_basepairs(fsa.allele_fit_func)
            for i in range(len(data)):
                rfu = data[i]
                bp  = basepairs[i] if basepairs else -999
                if bp>-999:
                    datastring+="[%i,%.2f,%i]," % (i,bp,rfu)
                else:
                    datastring+="[%i,null,%i]," % (i,rfu)
            datastring = datastring[:-1] + "]"

            out_stream.write("\"%10s\",\"%s\",\"%s\",\"%s\"\n" % (sample_name, well_id, trace_dye, datastring))

        out_stream.close()

def open_fsa( args, _params ):
    """ open FSA file(s) and prepare fsa instances
        requires: args.file, args.panel, args.panelfile
    """

    from fatools.lib.fileio.models import Marker, Panel, FSA

    if not args.panel:
        cexit('ERR: using FSA file(s) requires --panel argument!')

    if not args.panelfile:
        cerr('WARN: using default built-in panels')
        Panel.upload(params.default_panels)
    else:
        with open(args.panelfile) as f:
            # open a YAML file that describe panel sets
            import yaml
            Panel.upload(yaml.load(f))

    if not args.markerfile:
        Marker.upload(params.default_markers)
    else:
        raise NotImplementedError()

    panel = Panel.get_panel(args.panel)
    fsa_list = []
    index = 1

    # prepare caching
    cache_path = None
    if not args.no_cache:
        cache_path = os.path.join(os.path.expanduser('~'), '.fatools_caches', 'channels')
        if args.cache_path is not None:
            cache_path = os.path.join(args.cache_path, '.fatools_caches', 'channels')
        if not os.path.exists(cache_path):
            os.makedirs(cache_path)

    if args.file:
        for fsa_filename in args.file.split(','):
            fsa_filename = fsa_filename.strip()

            if args.indir != "":
                filename = args.indir + "/" + fsa_filename
            else:
                filename = fsa_filename

            fsa = FSA.from_file(filename, panel, _params,
                                cache=not args.no_cache,
                                cache_path=cache_path)
            # yield (fsa, str(i))
            fsa_list.append( (fsa, str(index)) )
            index += 1

    elif args.infile:

        with open(args.infile) as f:
            buf, delim = detect_buffer( f.read() )
        inrows = csv.DictReader( StringIO(buf), delimiter=delim )
        line = 1
        index = 1

        for r in inrows:

            line += 1

            fsa_filename = r['FILENAME'].strip()
            if fsa_filename.startswith('#'):
                continue

            if r.get('OPTIONS', None):
                options = tokenize( r['OPTIONS'] )
            else:
                options = None

            panel_code = r.get('PANEL', None) or args.panel
            panel = Panel.get_panel(panel_code)

            fsa = FSA.from_file(fsa_filename, panel, options,
                                cache=not args.no_cache,
                                cache_path=cache_path)
            if 'SAMPLE' in inrows.fieldnames:

                # yield (fsa, r['SAMPLE'])
                fsa_list.append( (fsa, r['SAMPLE']) )
            else:

                # yield (fsa, str(index))
                fsa_list.append( (fsa, str(index)) )
                index += 1

    elif args.indir:
        import glob
        for fsa_filename in sorted(glob.glob(args.indir+"/*.fsa")):

            fsa_filename = fsa_filename.strip()
            fsa = FSA.from_file(fsa_filename, panel, _params, cache = not args.no_cache)
            # yield (fsa, str(i))
            fsa_list.append( (fsa, str(index)) )
            index += 1

    return fsa_list


def get_fsa_list( args, dbh ):
    """
    get fsa instance from database based on parameters in args
    """

    if not args.batch:
        cexit('ERR: using database requires --batch argument!', 1)

    batch = dbh.get_batch( args.batch )
    if not batch:
        cexit('ERR: batch %s not found!' % args.batch, 1)

    samples = []
    if args.sample:
        samples = args.sample.split(',')

    fsas = []
    if args.fsa:
        fsas = args.assay.split(',')

    panels = []
    if args.panel:
        panels = args.panel.split(',')

    markers = []
    if args.marker:
        markers = dbh.get_markers(args.panel.split(','))

    fsa_list = []
    for sample in batch.samples:
        if samples and sample.code not in samples: continue
        for assay in sample.assays:
            if assays and assay.filename not in assays: continue
            if panels and assay.panel.code not in panels: continue
            fsa_list.append( (assay, sample.code) )

    cerr('I: number of assays to be processed: %d' % len(assay_list))
    return fsa_list
