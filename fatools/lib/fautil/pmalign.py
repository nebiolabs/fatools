"""
pair minimization algorithm
"""


import numpy as np
import itertools
from scipy.optimize import minimize

from fatools.lib.utils import cout, cerr, cverr, is_verbosity
from fatools.lib.fautil.alignutils import (estimate_z, pair_f, align_dp,
        pair_sized_peaks, DPResult, AlignResult, plot)
from fatools.lib.fautil.gmalign import ZFunc, align_gm
from fatools.lib import const



def align_pm(peaks, ladder, anchor_pairs=None):

    if not anchor_pairs:
        anchor_peaks = [ p for p in peaks if 1500 < p.rtime < 5000 ]

        # this finds the pair of peaks that best match to the 2nd and next-to-last ladder steps, and
        # does a linear fit to the rest of peaks to find the peaks matched to ladder steps
        anchor_pairs, initial_z = estimate_pm( anchor_peaks, ladder['signature'] )

    else:
        rtimes, bpsizes = zip( *anchor_pairs )
        initial_z = estimate_z(rtimes, bpsizes, 1)

    anchor_pairs.sort()

    # if the number of anchor pairs equals the number of ladder steps, no need to do pair matching
    if len(anchor_pairs)==len(ladder['sizes']):

        f = ZFunc(peaks, ladder['sizes'], anchor_pairs, estimate=True)

        anchor_rtimes, anchor_bpsizes = zip( *anchor_pairs )
        zres = estimate_z(anchor_rtimes, anchor_bpsizes, 2)
        score, z = minimize_score(f, zres.z, 2)        
        pairs, rss = f.get_pairs(z)

    else:
        pairs, z, rss, f = align_upper_pm(peaks, ladder, anchor_pairs, initial_z)
        pairs, z, rss, f = align_lower_pm(peaks, ladder, pairs, z)

    #print(rss)
    #plot(f.rtimes, f.sizes, z, pairs)
    # last dp
    dp_result = align_dp(f.rtimes, f.sizes, f.similarity, z, rss)

    if is_verbosity(4):
        import pprint; pprint.pprint(dp_result.sized_peaks)
        plot(f.rtimes, f.sizes, dp_result.z, [(x[1], x[0]) for x in dp_result.sized_peaks])

    dp_result.sized_peaks = f.get_sized_peaks(dp_result.sized_peaks)

    score, msg = ladder['qcfunc'](dp_result, method='strict')
    if score > 0.9:
        return AlignResult(score, msg, dp_result, const.alignmethod.pm_strict)

    score, msg = ladder['qcfunc'](dp_result, method='relax')
    return AlignResult(score, msg, dp_result, const.alignmethod.pm_relax)


    f = ZFunc(peaks, ladder['sizes'], anchor_pairs)

    z = initial_z
    score = last_score = 0
    last_z = None

    for order in [1, 2, 3]:

        last_rss = -1
        rss = 0

        niter = 0
        while abs(rss - last_rss) > 1e-3:

            niter += 1
            print('Iter: %d' % niter)

            print(z)
            score = f(z)
            if last_score and last_score < score:
                # score does not converge; just exit
                print('does not converge!')
                break

            pairs, cur_rss = f.get_pairs(z)
            rtimes, bpsizes = zip( *pairs )
            zres = estimate_z(rtimes, bpsizes, order)

            last_z = z
            z = zres.z
            last_rss = rss
            rss = zres.rss
            print(rss)

    dp_result = align_dp(f.rtimes, f.sizes, last_z, last_rss)

    return align_gm2(peaks, ladder, anchor_pairs, dp_result.z)



    new_anchor_pairs = []
    zf = np.poly1d(dp_result.z)
    for p in dp_result.sized_peaks:
        if (p[0] - zf(p[1]))**2 < 2:
            new_anchor_pairs.append( (p[1], p[0]) )
    if is_verbosity(4):
        #import pprint; pprint.pprint(dp_result.sized_peaks)
        plot(f.rtimes, f.sizes, dp_result.z, [(x[1], x[0]) for x in dp_result.sized_peaks])

    return align_gm(peaks, ladder, anchor_pairs, dp_result.z)


def align_lower_pm(peaks, ladder, anchor_pairs, anchor_z):

    # anchor pairs must be in asceding order


    last_rtime = anchor_pairs[-1][0]
    last_size = anchor_pairs[-1][1]
    anchor_rtimes, anchor_bpsizes = zip( *anchor_pairs )
    anchor_rtimes = list(anchor_rtimes)
    anchor_bpsizes = list(anchor_bpsizes)

    lower_peaks= [ p for p in peaks if p.rtime <= last_rtime ]
    lower_sizes = [ s for s in ladder['sizes'] if s <= last_size]

    # we try to pair-minimize lower_peaks and lower_sizes

    f = ZFunc(peaks, ladder['sizes'], anchor_pairs, estimate=True)

    scores = []

    # check the first
    est_first_bpsize = np.poly1d(anchor_z)(lower_peaks[0].rtime)
    remaining_sizes = [ s for s in lower_sizes if s >= est_first_bpsize ]

    if remaining_sizes:
        first_bpsize = remaining_sizes[0]

        for first_peak in lower_peaks[:-2]:
            if first_peak.rtime >= anchor_pairs[0][0]:
                break

            for first_bpsize in ladder['sizes'][:2]:
                zres = estimate_z( [ first_peak.rtime ] + anchor_rtimes, [ first_bpsize ] + anchor_bpsizes, 3 )
                #print('rss:', zres.rss)
                #plot(f.rtimes, f.sizes, zres.z, [ (first_peak.rtime, first_bpsize), ] )
                score, z = minimize_score(f, zres.z, 3)
                
                scores.append( (score, z) )
                #plot(f.rtimes, f.sizes, z, [ (first_peak.rtime, first_bpsize), ] )
                
        scores.sort( key = lambda x: x[0] )
        #import pprint; pprint.pprint( scores[:10] )

    if scores:
        z = scores[0][1]        
        pairs, rss = f.get_pairs(z)
    else:
        z = anchor_z
        pairs = anchor_pairs
        rss = None
        
    #plot(f.rtimes, f.sizes, z, pairs )

    return pairs, z, rss, f


    raise RuntimeError



def align_upper_pm(peaks, ladder, anchor_pairs, anchor_z):

    # anchor pairs must be in ascending order
    anchor_rtimes, anchor_bpsizes = zip( *anchor_pairs )
    anchor_rtimes = list(anchor_rtimes)
    anchor_bpsizes = list(anchor_bpsizes)

    # we try to pair-minimize higher peaks and sizes
    first_rtime = anchor_rtimes[0]
    first_bpsize = anchor_bpsizes[0]
    peaks = [ p for p in peaks if p.rtime >= first_rtime]
    sizes = [ s for s in ladder['sizes'] if s >= first_bpsize]
    remaining_sizes = [ s for s in ladder['sizes'] if s > anchor_bpsizes[-1] ]

    scores = []

    f = ZFunc(peaks, sizes, anchor_pairs, estimate=True)

    if remaining_sizes:
        
        #sizes = ladder['sizes']

        # check the first
        est_last_bpsize = np.poly1d(anchor_z)(peaks[-1].rtime)

        last_bpsize = max( remaining_sizes[1] if remaining_sizes else 0, [ s for s in sizes if s < est_last_bpsize ][-3] )
        
        for last_peak in reversed(peaks[-14:]):
            if last_peak.rtime <= anchor_pairs[-1][0]:
                break

            zres = estimate_z(anchor_rtimes + [last_peak.rtime], anchor_bpsizes + [last_bpsize], 2)
            #plot(f.rtimes, f.sizes, zres.z, [ (last_peak.rtime, last_bpsize)] )
            score, z = minimize_score(f, zres.z, 2)
            #print(score)
            #plot(f.rtimes, f.sizes, z, [] )
            
            scores.append( (score, z) )

        scores.sort( key = lambda x: x[0] )
        #import pprint; pprint.pprint( scores[:10] )

    if scores:
        z = scores[0][1]
        pairs, rss = f.get_pairs(z)        
    else:
        z = anchor_z
        pairs = anchor_pairs
        rss = None

    #print(rss)
    #plot(f.rtimes, f.sizes, z, pairs )

    return pairs, z, rss, f


def minimize_score( f, z, order ):

    last_score = score = 0

    niter = 1
    while niter  < 50:

        score = f(z)
        #print(score)

        if last_score and abs(last_score - score) < 1e-6:
            break

        pairs, rss = f.get_pairs(z)
        rtimes, bpsizes = zip( *pairs )
        zres = estimate_z(rtimes, bpsizes, order)

        z = zres.z
        last_score = score
        niter += 1

    return last_score, z



def estimate_pm(peaks, bpsizes):
    """
    returns sorted list of bp sizes matched to peaks and fits
    """
    
    rtimes = [ p.rtime for p in peaks ]

    rtime_points = prepare_rtimes( rtimes )
    bpsize_pair = [ bpsizes[1], bpsizes[-2]]

    f = ZFunc(peaks, bpsizes, [], estimate = True)


    # find linear fit for pair of peaks best matched to 2nd and next-to-last ladder step
    scores = []
    for rtime_pair in rtime_points:
        
        if rtime_pair[0] >= rtime_pair[1]:
            continue

        # y = ax + b
        # y1 = ax1 + b
        # y2 = ax2 + b
        # ------------ -
        # y1 - y2 = a(x1 - x2)
        # a = (y1 - y2)/(x1 - x2)
        # b = y1 - ax1

        #slope = (bpsize_pair[1] - bpsize_pair[0]) / (rtime_pair[1] - rtime_pair[0])
        #intercept = bpsize_pair[0] - slope * rtime_pair[0]
        #z = [ slope intercept ]

        # get the linear fit to this pair of peaks
        zres = estimate_z(rtime_pair, bpsize_pair, 1)

        # see how well all ladder peaks fit to this linear fit
        score = f(zres.z)

        scores.append( (score, zres) )
        #plot(f.rtimes, f.sizes, zres.z, [] )

    scores.sort( key = lambda x: x[0] )
    
    #import pprint; pprint.pprint(scores[:5])
    zresult = scores[0][1]

    # check this linear fit
    dp_result = align_dp(f.rtimes, f.sizes, f.similarity, zresult.z, zresult.rss)
    
    #import pprint; pprint.pprint(dp_result.sized_peaks)
    #plot(f.rtimes, f.sizes, dp_result.z, [(x[1], x[0]) for x in dp_result.sized_peaks])

    return ( [(x[1], x[0]) for x in dp_result.sized_peaks], dp_result.z )


def prepare_rtimes(rtimes):
    # prepare combination of begin and end rtimes

    mid_size = round(len(rtimes)/2)
    return list( itertools.product(rtimes[:mid_size], rtimes[mid_size-2:]) )

