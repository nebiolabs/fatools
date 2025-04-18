import numpy as np
import math

from fatools.lib.utils import cverr, is_verbosity
from fatools.lib import const
from fatools.lib.fautil.hcalign import align_hc
from fatools.lib.fautil.gmalign import align_gm, align_sh, align_de
from fatools.lib.fautil.pmalign import align_pm

from sortedcontainers import SortedListWithKey

from scipy import signal, ndimage
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy import asarray as ar,exp
import scipy.integrate as integrate
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline

import pandas as pd
import matplotlib.pyplot as plt

from peakutils import indexes

import attr

import logging, os
LOGLEVEL = os.environ.get("LOGLEVEL", "WARNING").upper()
logging.basicConfig(level=LOGLEVEL)

class LadderMismatchException(Exception):
    """Raised when number of peaks in ladder channel not equal to number of ladder steps."""
    pass

@attr.s(repr=False)
class Peak(object):
    rtime = attr.ib(default=-1)
    rfu = attr.ib(default=-1)
    area = attr.ib(default=-1)
    brtime = attr.ib(default=-1)
    ertime = attr.ib(default=-1)
    srtime = attr.ib(default=-1)
    beta = attr.ib(default=-1)
    theta = attr.ib(default=-1)
    omega = attr.ib(default=-1)

    size = attr.ib(default=-1)
    bin = attr.ib(default=-1)

    def __repr__(self):
        return "<P: %4d | %4d | %5d | %2d | %+3.2f | b %4.1f | t %4.2f | o %3d>" % (
            self.rtime, self.rfu, self.area, self.ertime - self.brtime, self.srtime,
            self.beta, self.theta, self.omega)

def scan_peaks(channel, parameters, offset=0):
    """
    """

    params = parameters.ladder if channel.is_ladder() else parameters.nonladder

    # check if channel is ladder channel, and adjust expected_peak_number accordingly
    expected_peak_number = params.expected_peak_number
    if channel.is_ladder():
        expected_peak_number = len(channel.fsa.panel.get_ladder()['sizes'])

    initial_peaks = find_peaks(channel, params, offset, expected_peak_number)

    # For LIZ120 ladder only! This will need to be modified for other dyes.
    if len(initial_peaks) > 9:
        # Find the largest gap between 2 peaks
        max_gap = 0
        extra_peak_idx = 0

        # Compare consecutive peaks
        for i in range(1, len(initial_peaks)):
            delta = initial_peaks[i].rtime - initial_peaks[i - 1].rtime

            if delta > max_gap:
                max_gap = delta
                extra_peak_idx = i  # Potential extra peak index

        # If gap was between peak indices 0 and 1 it implies it was first peak
        if extra_peak_idx == 1:
            extra_peak_idx = 0

        # Find peaks with min/max RFU
        min_rfu = float('inf')
        max_rfu = float('-inf')
        min_idx = max_idx = 0

        for idx, peak in enumerate(initial_peaks):
            if peak.rfu > max_rfu:
                max_rfu = peak.rfu
                max_idx = idx
            if peak.rfu < min_rfu:
                min_rfu = peak.rfu
                min_idx = idx

        # Remove extra peak if it's both min/max RFU and largest gap
        if extra_peak_idx in (max_idx, min_idx):
            del initial_peaks[extra_peak_idx]

    # create alleles based on these peaks
    alleles = []
    for p in initial_peaks:
        allele = channel.Allele(
            rtime = p.rtime,
            rfu = p.rfu,
            rfu_uncorr = p.rfu_uncorr,
            area = p.area,
            brtime = p.brtime,
            ertime = p.ertime,
            wrtime = p.wrtime,
            srtime = p.srtime,
            beta = p.beta,
            theta = p.theta,
            omega = p.omega,
        )
        allele.type = const.peaktype.scanned
        allele.method = const.binningmethod.notavailable
        allele.marker = channel.marker
        channel.add_allele( allele )
        alleles.append( allele )

    channel.status = const.channelstatus.scanned
    return alleles

def preannotate_peaks(channel, params):

    """
    pre-annotate peaks as peak-scanned / peak-broad / peak-stutter/ peak-overlap
    based on criteria defined in params
    """

    # peak_1 is overlap of peak_2 if
    #   brtime2 < ertime1 and ertime1 < ertime2
    #   and height1 at rtime1 is a_fraction of height2 at rtime1 and
    #   height1 at rtime2 is a_fraction of height2 at rtime2.

    # peak is broad if beta > beta_broad_threshold

    channel_peaks = [ (list(channel.get_alleles(broad_peaks_only = False)),
                       np.median(channel.data)) ]

    # reset all peak type, score the peaks and set the peak type to peak-noise,
    # peak-broad

    # collect beta * theta first, and used beta * theta as descriptor for noise
    # also if height at either brtime or ertime is higher than 50% at rtime, it is
    # likely a noise

    for (peaks, med_baseline) in channel_peaks:

        if len(peaks) == 0: continue

        beta_theta = sorted([ p.beta * p.theta for p in peaks ])
        sampled_beta_theta = beta_theta[2:len(beta_theta)-2]
        if len(sampled_beta_theta) == 0: sampled_beta_theta = beta_theta
        avg_beta_theta = sum(sampled_beta_theta) / len(sampled_beta_theta)

        for p in peaks: # these are actually Alleles, not Peaks

            p.type = const.peaktype.scanned
            p.size = -1
            p.bin = -1

            peak_beta_theta = p.beta * p.theta
            score = 1.0

            # extreme noise

            if p.height < 2 * med_baseline:
                p.qscore = 0.25
                p.type = const.peaktype.noise
                continue

            if p.wrtime < 6 or (p.wrtime < 10 and peak_beta_theta < 0.275 * avg_beta_theta):
                p.qscore = 0.25
                p.type = const.peaktype.noise
                continue

            # moderately noise

            if peak_beta_theta < 0.33 * avg_beta_theta:
                if (    channel.data[p.brtime] > 0.5 * p.height or
                        channel.data[p.ertime] > 0.5 * p.height ):
                    p.qscore = 0.25
                    p.type = const.peaktype.noise
                    continue
                score -= 0.15

            score = 1.0
            if p.beta > params.max_beta:
                p.type = const.peaktype.broad
                score -= 0.20
            elif p.beta < 5:
                score -= 0.20

            # check theta
            if p.theta < 4:
                # perhaps an artifact
                score -= 0.20

            # penalty by height
            if p.height < 75:
                # decrease the score
                score -= 0.1
            if p.height < 50:
                # decrease even further
                score -= 0.1

            # penalty by symmetrics
            if not ( -1.32 < p.srtime < 1.32 ):
                score -= 0.1

            p.qscore = score
            if p.qscore < 0.5 and p.type == const.peaktype.scanned:
                p.type = const.peaktype.noise

            if p.qscore < 0:
                p.qscore = 0.0  # reset to zero

    # checking for stutter peaks based on minimum rtime & rfu

    for (peaks, med_baseline) in channel_peaks:
        alleles = sorted( [ p for p in peaks ],
                        key = lambda x: x.rtime )

        for idx in range( len(alleles) ):
            allele = alleles[idx]            
            if idx > 0:
                allele_0 = alleles[idx-1]
                if allele.rtime - allele_0.rtime < params.stutter_rtime_threshold:
                    if allele_0.height * params.stutter_height_threshold > allele.height:
                        allele.type = const.peaktype.stutter
                        allele.qscore -= 0.2
            if idx < len(alleles) - 1:
                allele_1 = alleles[idx+1]
                if allele_1.rtime - allele.rtime < params.stutter_rtime_threshold:
                    if allele_1.height * params.stutter_height_threshold > allele.height:
                        allele.type = const.peaktype.stutter
                        allele.qscore -= 0.2


def call_peaks( channel, params, func, min_rtime, max_rtime ):
    """
    call (determine size) each of peaks with type peak-scanned, and annotate as either
    peak-called or peak-unassigned
    """

    
    for allele in channel.get_alleles(broad_peaks_only=False):

        if not min_rtime < allele.rtime < max_rtime:
            if allele.type == const.peaktype.scanned:
                allele.type = const.peaktype.unassigned
            if is_verbosity(2):
                print("allele at ", allele.rtime," not called... outside range [", min_rtime, ", ", max_rtime, "] (height=", allele.height,") !")
            continue
        ladder = channel.is_ladder()
        size, deviation, qcall, method = func(allele.rtime, ladder)
        allele.size = size
        allele.bin = round(size)
        allele.deviation = deviation
        allele.qcall = qcall
        if allele.type == const.peaktype.scanned:
            allele.type = const.peaktype.called
        allele.method = const.binningmethod.notavailable

        # set values in base pair units
        allele.begin_bp = func(allele.brtime, ladder)[0]
        allele.end_bp = func(allele.ertime, ladder)[0]
        allele.width_bp = allele.end_bp - allele.begin_bp
        size_delta = func(allele.rtime+1, ladder)[0] - size
        allele.area_bp = float(allele.area) * size_delta

        if is_verbosity(4):
            print(allele)

def align_peaks(channel, params, ladder, anchor_pairs=None, saturated_peak_rtimes=[]):
    """
    returns (score, rss, dp, aligned_peak_number)
    """

    alleles = channel.get_alleles()

    if (len(alleles) != len(ladder['sizes'])):
        # attempt to filter out cross talk peaks from saturated signal in other channels
        logging.debug(f"observed rtimes:{[a.rtime for a in alleles]}, saturated_rtimes:{saturated_peak_rtimes}")
        alleles = [ a for a in alleles if not any(np.isclose(a.rtime,s, atol= 2) for s in saturated_peak_rtimes)] 
        if (len(alleles) != len(ladder['sizes'])):  
            raise LadderMismatchException( ("alleles not same length as ladder for file: %s!") % channel.fsa.filename)

    if anchor_pairs:
        return align_pm( alleles, ladder, anchor_pairs)

    if len(alleles) <= len(ladder['sizes']) + 5:
        result = align_hc( alleles, ladder )

        if result.score > 0.9:
            return result

    return align_pm( alleles, ladder )

    if result.initial_pairs:
        result = align_gm( alleles, ladder, result.initial_pairs )
        if result.score > 0.75:
            return result

    result = align_sh( alleles, ladder )
    if result.score > 0.75:
            return result

    # perform differential evolution
    return align_de( alleles, ladder )

    raise RuntimeError


    result = hclust_align( alleles, ladder )

    # add relevant info to peaks
    aligned_peaks = result[2][3]
    f = np.poly1d( result[2][2] )
    for (size, p) in aligned_peaks:
        p.dev = abs( f(p.rtime) - size)
        p.size = size
    z, rss = estimate_z( [p[1].rtime for p in aligned_peaks],
                          [p[0] for p in aligned_peaks], 3)
    print('>>> RSS:', rss)
    #import pprint; pprint.pprint( aligned_peaks )
    return result

def normalize_peaks(channel, params):

    #scale_factor = channel.fsa.area_scale_factor
    sf_poly1d = np.poly1d(channel.fsa.area_scale_factor_params)
    
    for allele in channel.get_alleles(False):
        #allele.area_bp_corr = allele.area_bp * scale_factor
        allele.area_bp_corr = allele.area_bp * sf_poly1d(allele.size)

def merge_peaks( channel, params, func, plot=False ):
    """
    Group peaks and fit to a series of Gaussians on top of a polynomial
    """

    smeared_peaks = []
    def gaus(x, *p0):
        ngauss = int(len(p0)/3+.01)
        val = 0.
        for i in range(ngauss):
            a  = p0[3*i]
            x0 = p0[3*i+1]
            sigma = p0[3*i+2]
            val += a*exp(-(x-x0)**2/(2*sigma**2))
        return val

    def poly(x, xbar, *p0):
        val = p0[0]+(x-xbar)*p0[1]+(x-xbar)*(x-xbar)*p0[2]
        return val
    
    def gauspoly(x, xbar, *p0):
        return gaus(x,*p0[0:-3]) + poly(x,xbar,*p0[-3:])

    broad_alleles = channel.get_alleles(True)[:]
    all_alleles = channel.get_alleles()[:]
    basepairs = channel.get_basepairs()

    nmergedpeaks = 0
    
    while len(broad_alleles)>0:

        # find tallest peak
        peak = max(broad_alleles, key=lambda p: p.rfu)
        
        if peak.rfu<params.min_rfu or peak.size<-999:
            break

        # rtime is the location of the current peak
        rtime = peak.rtime

        # array containing peak centers and heights (starting with tallest peak)
        peak_centers = [basepairs[peak.rtime]]
        peak_heights = [channel.data[peak.rtime]]

        # get left and right bounds
        left  = find_smeared_peak_bound(peak, channel, peak_centers, peak_heights,
                                        move_to_right=False)
        right = find_smeared_peak_bound(peak, channel, peak_centers, peak_heights,
                                        move_to_right=True)

        # remove peaks from broad_alleles between left and right
        broad_alleles = [ allele for allele in broad_alleles
                          if (allele.rtime<left or allele.rtime>right) ]

        #
        # do an initial fit to make sure we have a real peak with a sensible width
        #
        height = channel.data[peak.rtime]
        mean = basepairs[peak.rtime]
        p0 = [ height, mean, 0.25]
        fitbounds = [[0, mean-1., 0.], [height+2.*np.sqrt(height), mean+1., .6]]

        first_peak_right = peak.rtime
        while basepairs[first_peak_right+1]>-999 and \
              abs(basepairs[first_peak_right]-basepairs[peak.rtime])<0.5:
            first_peak_right += 1
        first_peak_left = peak.rtime
        while basepairs[first_peak_left-1]>-999 and \
              abs(basepairs[first_peak_left]-basepairs[peak.rtime])<0.5:
            first_peak_left -= 1

        x = ar(basepairs[first_peak_left:first_peak_right+1])
        y = ar(channel.data[first_peak_left:first_peak_right+1])

        try:
            popt, pcov = curve_fit(gaus, x, y, p0=p0, bounds=fitbounds)
            perr = np.sqrt(np.diag(pcov))

        except RuntimeError:
            #print("initial fit failed, will skip this peak")
            perr = [ 100. ]

        # use uncertainty in width of 1st peak to see how good fit is
        if np.sqrt(pcov[2,2])>0.05:
            continue

        #
        # We have a good peak, so fit to series of gaussians
        #
        p0=[]
        fitbounds=[[],[]]
        for peak_center, peak_height in zip(peak_centers, peak_heights):

            height = peak_height
            mean = peak_center
            sigma = 0.25
            p0.extend([height, mean, sigma])
            fitbounds[0].extend([0,mean-.5,.2])
            fitbounds[1].extend([height+2.*np.sqrt(height), mean+.5, .6])

        # x and y arrays for defining background and total functions 
        x = ar(basepairs[left:right+1])
        npoints = int(200 * (basepairs[right]-basepairs[left])/5)
        xfine = np.linspace(basepairs[left],basepairs[right],npoints)
        y = ar(channel.data[left:right+1])
        
        # background and total functions
        xbar = basepairs[peak.rtime]
        f_poly = lambda x, *p: poly(x,xbar,*p[-3:])
        f = lambda x, *p: gaus(x,*p[0:-3]) + f_poly(x,*p[-3:])

        # background parameters
        p0.extend([0,0,0])
        fitbounds[0].extend([0, -np.inf, -np.inf])
        fitbounds[1].extend([np.inf, np.inf, np.inf])
        
        try:
            popt, pcov = curve_fit(f, x, y, p0=p0, bounds=fitbounds)
            p0 = popt.tolist()
            
        except RuntimeError:
            #print("runtime error")
            perr = [ 100. ]

        # remove any degenerate or too narrow peaks (always keep first peak)
        new_p0 = p0[0:3]
        new_bounds = [fitbounds[0][0:3],fitbounds[1][0:3]]
        ngaus = int(int(len(p0)-3)/3)

        for i in range(ngaus-1, 0, -1): # for each gauss in p0, add to new_p0 if good

            # see if gauss is too narrow
            good_gauss = p0[3*i+2]>.1
            if not good_gauss:
                continue
            
            # see if any other peaks are too close
            peakrtime_i = p0[3*i+1]
            for j in range(i-1, -1, -1):
                peakrtime_j = p0[3*j+1]
                if np.abs(peakrtime_j - peakrtime_i)<.55:
                    good_gauss = False
                    break
            if not good_gauss:
                continue

            new_p0.extend(p0[3*i:3*i+3])
            new_bounds[0].extend(fitbounds[0][3*i:3*i+3])
            new_bounds[1].extend(fitbounds[1][3*i:3*i+3])

        new_p0.extend(p0[-3:])
        new_bounds[0].extend(fitbounds[0][-3:])
        new_bounds[1].extend(fitbounds[1][-3:])

        p0 = new_p0
        fitbounds = new_bounds
        
        # fit again if number parameters changed
        ngaus_new = int(int(len(p0)-3)/3)
            
        if  ngaus_new<ngaus:                    
            try:
                popt, pcov = curve_fit(f, x, y, p0=p0, bounds=fitbounds)
                p0 = popt.tolist()
                
            except RuntimeError:
                #print("runtime error")
                perr = [ 100. ]

        if plot:
            fig = plt.figure()
            plt.plot(x,y,'b+:',label='data')
            plt.plot(xfine,gauspoly(xfine,xbar,*p0),'ro:',label='fit', markersize=2)
            plt.plot(xfine, poly(xfine,xbar,*p0[-3:]),'r--', label='bkg', markersize=2)
            for i in range(ngaus_new):
                plt.axvline(p0[3*i+1])
            plt.xlabel("peak size (base pair units]")
            plt.ylabel("intensity (relative freq. units)")

        # append individual peaks based on gaussians
        means = [ p0[3*i+1] for i in range(ngaus_new) ]
        sorted_means = sorted(means)
        
        for i in range(ngaus_new):

            height = p0[3*i]
            mean = p0[3*i+1]
            sigma = p0[3*i+2]

            # get rtime and peak boundaries in scan time units
            sorted_i = sorted_means.index(mean)
            rtime = left
            while basepairs[rtime]<mean:
                rtime+=1
            rtime-=1

            if sorted_i==0:
                brtime = left
            else:
                left_bp = (mean + sorted_means[sorted_i-1])/2
                brtime = rtime
                while basepairs[brtime]>left_bp:
                    brtime-=1

            if sorted_i==ngaus_new-1:
                ertime=right
            else:
                right_bp = (mean + sorted_means[sorted_i+1])/2
                ertime = rtime
                while basepairs[ertime]<right_bp:
                    ertime+=1
                ertime-=1

            begin_bp = basepairs[brtime]
            end_bp = basepairs[ertime]
            area = integrate.quad(gaus, begin_bp, end_bp, args=tuple(p0))[0]

            l_area = integrate.quad(gaus, begin_bp, mean, args=tuple(p0))[0]
            r_area = integrate.quad(gaus, mean, end_bp, args=tuple(p0))[0]
            
            srtime = 0
            if r_area>0 and l_area>0:
                srtime = math.log2(r_area / l_area)

            # use Allele to store info for the smeared allele
            allele = channel.Allele(
                rtime = rtime,
                rfu = height,
                rfu_uncorr = height,
                area = area, #np.sqrt(2.*np.pi) * height * sigma,
                brtime = brtime,
                ertime = ertime,
                wrtime = ertime - brtime,
                srtime = srtime,
                beta = 0, theta = 0, omega = 0
            )
            allele.size = mean
            allele.type = const.peaktype.smeared
            allele.method = const.binningmethod.notavailable
            allele.marker = channel.marker
            allele.group = nmergedpeaks

            # set values in base pair units
            allele.begin_bp = begin_bp
            allele.end_bp = end_bp
            allele.width_bp = end_bp - begin_bp

            delta = (basepairs[right]-basepairs[left])/float(right-left)            
            allele.area_bp = float(allele.area) * delta
            
            smeared_peaks.append(allele)

            if plot:
                # add dotted lines for peak boundaries to plot
                plt.axvline(basepairs[brtime], linestyle='--')
                plt.axvline(basepairs[ertime], linestyle='--')

        if plot:
            #plt.show()
            figname = "fit_"+channel.fsa.filename[:-4]+"_"+str(left)+"_"+str(right)+".png"
            fig.savefig(figname)
    
        nmergedpeaks += 1

    return smeared_peaks

def find_smeared_peak_bound( peak, channel, peak_centers, peak_heights, move_to_right ):

    basepairs = channel.basepairs
    rtime = peak.rtime

    if move_to_right:
        step = 1
    else:
        step = -1

    total_width = 0.5

    # go half a basepair unit to get to end of main peak
    while basepairs[rtime+step]>-999 and abs(basepairs[rtime] - basepairs[peak.rtime]) < total_width:
        rtime += step # move one scan-time unit

    while basepairs[rtime+step]>-999:
        max_val = 0.
        min_val = 1.e5
        rtime_max_val = 0                   
        last_rtime = rtime
        while abs(basepairs[rtime] - basepairs[peak.rtime]) < (1.+total_width) and basepairs[rtime+step]>-999:

            # get max over the range
            val = channel.data[rtime]
            if val > max_val:
                max_val = val
                rtime_max_val = rtime
            if val < min_val:
                min_val = val
                
            # move one scan-time unit
            rtime += step

        if last_rtime == rtime:
            break
        
        if max_val>0.:
            peak_centers.append(basepairs[rtime_max_val])
            peak_heights.append(max_val)

        # check max over range... if below threshold then stop
        if (max_val-min_val)<0.05*channel.data[peak.rtime]:
            break
        
        total_width += 1.

    return rtime

    
# helper functions


def find_raw_peaks(channel, params, offset, expected_peak_number=0):
    
    """
    params.min_dist
    params.norm_thres
    params.min_rfu
    params.min_rfu_ratio
    params.max_peak_number
    """

    data = channel.data
    
    #print("expected:", expected_peak_number)
    # cut and pad data to overcome peaks at the end of array
    obs_data = np.append(data[offset:], [0,0,0])
    if False: #expected_peak_number:
        min_dist = params.min_dist
        indices = []
        norm_threshold = params.norm_thres
        expected_peak_number = expected_peak_number * 1.8
        while len(indices) <= expected_peak_number and norm_threshold > 1e-7:
            indices = indexes( obs_data, norm_threshold, min_dist)
            norm_threshold *= 0.5
    elif False:
        indices = indexes( obs_data, params.norm_thres, params.min_dist)

    if params.peakwindow==0 or params.peakdegree<1:
        indices = indexes( obs_data, 1e-7, params.min_dist)
        channel.firstderiv = np.hstack([np.diff(obs_data),0.]).tolist()

    else:

        # do a fit in a sliding window of size given by peakwindow
        # to a polynomial of degree given by peakdegree
        firstderiv_ = []
        winsize = int(params.peakwindow/2)*2+1 # make sure window is odd
        degree = params.peakdegree

        data_np = np.asarray(data)


        firstderiv = signal.savgol_filter(data_np, winsize, degree, deriv=1)
        channel.firstderiv = firstderiv.tolist()

        thres = 1e-7*(np.max(data)-np.min(data))+np.min(data)

        indices = []
        last = False
        for i in range(len(firstderiv)):
            current = (firstderiv[i]>0)
            if last==True and current==False:
                if data[i]>last_data:
                    indices.append(i)
                else:
                    indices.append(i-1)
            last_data = data[i]
            last = current
        indices = np.array(indices)

        """import matplotlib.pyplot as plt
        plt.plot(data, label="data")
        plt.plot(channel.firstderiv, '-',label="first deriv")
        plt.plot((0.,6000.),(0.,0.),'--')
        plt.legend()
        plt.show()"""

    cverr(5, '## indices: %s' % str(indices))
    cverr(3, '## raw indices: %d' % len(indices))

    if len(indices) == 0:
        return []

    # normalize indices
    if offset > 0:
        indices = indices + offset

    # filter peaks by minimum rfu and minimum rfu ratio

    def pass_threshold(h, params, maxheight):
        return ((h >= params.min_rfu) and
                (h >= params.min_rfu_ratio * maxheight))

    max_rfu = max(data)
    peaks = [ Peak( int(i), int(data[i]) ) for i in indices
              if ( pass_threshold(data[i], params, max_rfu) and
                   params.min_rtime < i < params.max_rtime) ]

    if not peaks: return peaks
    
    # call measure_peaks here to adjust the rfu by the baseline
    measure_peaks(peaks, channel, params.baseline_correct, offset)

    # check the thresholds again
    max_rfu = max([peak.rfu for peak in peaks])
    peaks = [ peak for peak in peaks if pass_threshold(peak.rfu,params,max_rfu) ]
    
    # filter peaks by maximum peak number after sorted by rfu
    #peaks = sorted( peaks, key = lambda x: x.rfu )[:params.max_peak_number * 2]
    #import pprint; pprint.pprint(peaks)
    #print('======')

    if expected_peak_number:
        peaks.sort( key = lambda x: x.rfu, reverse = True )
        peaks = peaks[: round(expected_peak_number * 2)]
        peaks.sort( key = lambda x: x.rtime )

    cverr(3, '## peak above min rfu: %d' % len(peaks))

    return peaks


def find_peaks(channel, params, offset=0, expected_peak_number=0):

    data = channel.data
    
    peaks = find_raw_peaks(channel, params, offset, expected_peak_number)

    # check for any peaks
    if not peaks:
        return peaks

    # measure peaks parameters
    #measure_peaks(peaks, channel, offset)

    #import pprint; pprint.pprint(peaks)

    # filter artefact peaks
    non_artifact_peaks = peaks
    if not params.keep_artifacts:
        non_artifact_peaks = filter_for_artifact(peaks, params, expected_peak_number)

    # for ladder, special filtering is applied
    if params.expected_peak_number:
        peaks = filter_for_ladder(non_artifact_peaks, params)
    else:
        peaks = non_artifact_peaks

    return peaks


def measure_peaks(peaks, channel, baseline_correct = True, offset=0):

    data = channel.data
    firstderiv = channel.firstderiv
    
    (q50, q70) = np.percentile( data[offset:], [50, 75] )
    for p in peaks:
        if False:
            p.area, p.brtime, p.ertime, p.srtime, ls, rs = \
                calculate_area( data, p.rtime, 5e-2, q50 )
        else:
            p.area, p.brtime, p.ertime, p.srtime, ls, rs = \
                calculate_area_firstderiv( data, firstderiv, p.rtime )

        # correct for baseline
        baseline = min(data[p.brtime],data[p.ertime])
        p.wrtime = p.ertime - p.brtime

        p.rfu_uncorr = p.rfu
        if baseline_correct:            
            p.area -= baseline * p.wrtime
            p.rfu -= baseline
            p.rfu = max(p.rfu, 0.)

        if p.rfu==0.: continue
        
        p.beta = p.area / p.rfu
        if p.wrtime == 0:
            p.theta = 0
            p.omega = 0
        else:
            p.theta = p.rfu / p.wrtime
            p.omega = p.area / p.wrtime


def calculate_area(y, t, threshold, baseline):
    """ return (area, brtime, ertime, srtime)
        area: area
        brtime: begin rtime
        ertime: end rtime
    """

    # right area
    data = y[t:]
    r_area, ertime, r_shared = half_area(data, threshold, baseline)

    # left area
    data = y[:t+1][::-1]
    l_area, brtime, l_shared = half_area(data, threshold, baseline)


    return ( l_area + r_area - y[t], t - brtime, ertime + t, math.log2(r_area / l_area),
                l_shared, r_shared )


def half_area(y, threshold, baseline):
    """ return (area, ertime, shared_status)
    """

    winsize = 3
    threshold = threshold/2
    shared = False
    area = y[0]
    edge = float(np.sum(y[0:winsize]))/winsize
    old_edge = 2 * edge

    index = 1
    limit = len(y)

    while ( edge > area * threshold and edge < old_edge and
            index < limit and y[index] >= baseline ):
        old_edge = edge
        area += y[index]
        edge = float(np.sum(y[index:index+winsize]))/winsize
        index += 1
    if edge >= old_edge:
        shared = True
    index -= 1

    return area, index, shared

def calculate_area_firstderiv(y, dy, t):
    """ return (area, brtime, ertime, srtime)
        area: area
        brtime: begin rtime
        ertime: end rtime
    """
    
    # right area
    data = y[t:]
    firstderiv = dy[t:]
    r_area, ertime, r_shared = half_area_firstderiv(data, firstderiv, True)

    # left area
    data = y[:t+1][::-1]
    firstderiv = dy[:t+1][::-1]
    
    l_area, brtime, l_shared = half_area_firstderiv(data, firstderiv, False)

    log2areas = 0
    if r_area>0 and l_area>0:
        log2areas = math.log2(r_area / l_area)
        
    return ( l_area + r_area - y[t], t - brtime, ertime + t, log2areas,
             l_shared, r_shared )


def half_area_firstderiv(y, firstderiv, decreasing):
    """ return (area, ertime, shared_status)
    """

    shared = False
    area = y[0]
    dy = firstderiv[0]

    limit = len(firstderiv)

    index = 1

    if index<limit:
        pos = ( firstderiv[index] < 0. ) if decreasing else ( firstderiv[index] > 0. ) 
    else:
        pos = False
        
    while ( pos and index < limit-1 ):
        area += y[index]
        index += 1
        pos = ( firstderiv[index] < 0. ) if decreasing else ( firstderiv[index] > 0. )
    index -= 1

    if area<0:
        area=0
        
    return area, index, shared


def math_func(x, a, b):
    #return a*np.exp(x*b)
    return a*x + b

def quadratic_math_func(x, a, b, c):
    return a*x**2 + b*x + c


def filter_for_artifact(peaks, params, expected_peak_number = 0):
    """
    params.max_peak_number
    params.artifact_ratio
    params.artifact_dist ~ 5
    """

    # the following code in this function performs the necessary acrobatic act
    # to select the most likely peaks that can be considered as true signals,
    # which is especially necessary for ladder - size assignment

    if len(peaks) == expected_peak_number:
        return peaks

    # we need to adapt to the noise level of current channel
    if expected_peak_number > 0:
        epn = expected_peak_number
        theta_peaks = sorted(peaks, key = lambda x: x.theta, reverse=True)[round(epn/2)+3:epn-1]
        #theta_peaks = theta_peaks[2:4] + theta_peaks[round(epn/2):epn-1]
        omega_peaks = sorted(peaks, key = lambda x: x.omega, reverse=True)
        omega_peaks = omega_peaks[2:4] + omega_peaks[round(epn/2):epn-1]
        rfu_peaks = sorted(peaks, key = lambda x: x.rfu, reverse=True)[:epn-1]

        if theta_peaks[-1].theta < 8:
            theta_peaks.sort()
            thetas = np.array([ p.theta for p in theta_peaks ])
            rtimes = [ p.rtime for p in theta_peaks ]

            #plt.scatter(rtimes, thetas)
            #plt.show()
            popt, pcov = curve_fit( math_func, rtimes, 0.5 * thetas, p0 = [ -1, 1 ])

            if is_verbosity(4):
                import matplotlib.pyplot as plt
                xx = np.linspace( rtimes[0], rtimes[-1]+2000, 100 )
                yy = math_func(xx, *popt)
                plt.plot(xx, yy)
                plt.scatter( [p.rtime for p in peaks], [p.theta for p in peaks])
                plt.show()

            q_theta = lambda x: x.theta >= math_func(x.rtime, *popt) or x.theta > 100

        else:
            q_theta = lambda x: x.theta >= min(theta_peaks[-1].theta, params.min_theta)


        if omega_peaks[-1].omega < 200:
            omega_peaks.sort()
            omegas = np.array([ p.omega for p in omega_peaks ])
            rtimes = np.array([ p.rtime for p in omega_peaks ])

            # generate a quadratic threshold for omega

            # generate a quadratic ratio series first
            popt, pcov = curve_fit( quadratic_math_func,
                    [rtimes[0], (rtimes[0] + rtimes[-1])/2, rtimes[-1]],
                    [0.1, 0.3, 0.1])
            ratios = quadratic_math_func(rtimes, *popt)
            if is_verbosity(4):
                plt.plot(rtimes, ratios)
                plt.show()

            # use the ratios to enforce quadratic threshold
            popt, pcov = curve_fit( quadratic_math_func, rtimes, ratios * omegas,
                                        p0 = [ -1, 1, 0 ])
            if popt[0] > 0:
                # enforce small flat ratio
                popt, pcov = curve_fit( math_func, rtimes, 0.25 * omegas, p0 = [ 1, 0 ])
                popt = np.insert(popt, 0, 0.0)  # convert to 3 params
            if is_verbosity(4):
                plt.scatter(rtimes, omegas)
                xx = np.linspace( rtimes[0], rtimes[-1]+2000, 100 )
                yy = quadratic_math_func(xx, *popt)
                plt.plot(xx, yy)
                plt.scatter( [p.rtime for p in peaks], [p.omega for p in peaks])
                plt.show()

            q_omega = lambda x: (   x.omega >= 100 or
                                    x.omega >= quadratic_math_func(x.rtime, *popt) )

        else:

            q_omega = lambda x: x.omega >= min(omega_peaks[-1].omega, 125)


        min_rfu = rfu_peaks[-1].rfu * 0.125

    else:
        min_theta = 0
        min_omega = params.min_omega
        min_theta_omega = 0
        min_rfu = 2
        q_omega = lambda x: x.omega >= min_omega


    # filter for too sharp/thin peaks
    filtered_peaks = []
    for p in peaks:
        #filtered_peaks.append(p); continue\

        if len(filtered_peaks) < 2 and p.area > 50:
            # first two real peaks might be a bit lower
            filtered_peaks.append(p)
            continue

        if not q_omega(p):
            print('! q_omega')
            continue
        #if not q_theta(p):
        #    print('! q_theta')
        #    continue

        #if min_theta and min_omega and p.omega < min_omega and p.theta < min_theta:
        #    print('! omega & theta')
        #    continue
        #if min_theta_omega and p.theta * p.omega < min_theta_omega:
        #    print('! theta_omega')
        #    continue
        if p.theta < 1.0 and p.area < 25 and p.omega < 5:
            print('! extreme theta & area & omega')
            continue
        if p.rfu < min_rfu:
            print('! extreme min_rfu')
            continue
        if p.beta > 25 and p.theta < 0.5:
            print('! extreme beta')
            continue
        if p.wrtime < 3:
            continue
        if p.rfu >= 25 and p.beta * p.theta < 6:
            continue
        if p.rfu < 25 and p.beta * p.theta < 3:
            continue
        #if p.omega < 50:
        #    continue
        #if p.omega < 100 and p.theta < 5:
        #    continue
        #if ( params.max_beta and min_theta and
        #        (p.beta > params.max_beta and p.theta < min_theta) ):
        #    print('! max_beta')
        #    continue
        filtered_peaks.append(p)

    #import pprint; pprint.pprint(filtered_peaks)

    # filter for distance between peaks and their rfu ratio
    peaks = sorted(filtered_peaks, key = lambda x: x.rtime)
    non_artifact_peaks = []
    for idx in range(len(peaks)):
        p = peaks[idx]

        if idx > 0:
            prev_p = peaks[idx-1]
            if ( p.brtime - prev_p.ertime < params.artifact_dist
                    and p.rfu < params.artifact_ratio * prev_p.rfu ):
                # we are artifact, just skip
                print("artifact: ",p)
                continue

        if idx < len(peaks)-1:
            next_p = peaks[idx+1]
            if ( next_p.brtime - p.ertime < params.artifact_dist
                    and p.rfu < params.artifact_ratio * next_p.rfu ):
                # we are artifact, just skip
                print("artifact: ",p)
                continue

        non_artifact_peaks.append( p )

    #import pprint; pprint.pprint(non_artifact_peaks)
    #print(len(non_artifact_peaks))

    peaks = non_artifact_peaks

    cverr(3, '## non artifact peaks: %d' % len(peaks))

    return peaks


def filter_for_ladder(peaks, params):
    """
    we need to obtaine enough peaks for ladder alignment purpose, but not too much to avoid
    excessive alignment process and potentially incorrect alignment

    peaks must in rtime ascending order
    """

    epn = params.expected_peak_number   # this is the number of ladder peaks

    #
    return peaks


def baseline_als(y, lam, p, niter=10):
    pass


@attr.s
class NormalizedTrace(object):
    signal = attr.ib()
    baseline = attr.ib()

    def get_qc(self):
        """ return tuple of qcfunc
        """
        return tuple()


def normalize_baseline( raw, params, is_ladder):
    """
    params.medwin_size
    params.savgol_order
    params.savgol_size
    """

    medwinsize = params.baselinewindow

    if is_ladder:
        savgol_size = params.ladder.smoothing_window
        savgol_order = params.ladder.smoothing_order
    else:
        savgol_size = params.nonladder.smoothing_window
        savgol_order = params.nonladder.smoothing_order
        
    if params.baselinemethod == const.baselinemethod.median:
        baseline_raw = signal.medfilt(raw, [medwinsize])

    elif params.baselinemethod == const.baselinemethod.minimum:
        df = pd.Series(raw)
        baseline_df = df.rolling(medwinsize,center=True).min()
        baseline_raw = baseline_df.tolist()

        # correct for NaNs in beginning and end of list
        halfwin = (int)(medwinsize/2)
        baseline_raw[:halfwin] = [baseline_raw[halfwin+1]]*halfwin
        baseline_raw[-halfwin:] = [baseline_raw[-halfwin-1]]*halfwin

    elif params.baselinemethod == const.baselinemethod.none:
        baseline_raw = raw

    else:
        raise RuntimeError("invalid option for baseline method")

    baseline = signal.savgol_filter( baseline_raw, medwinsize, 5)
    corrected_baseline = raw - baseline
    np.maximum(corrected_baseline, 0, out=corrected_baseline)

    #print("savgol_size: ", savgol_size,", savgol_order: ", savgol_order)

    """
    plt.figure()
    plt.plot(corrected_baseline, color='red')

    plt.plot(signal.savgol_filter(corrected_baseline, 3, 1), color='blue')
    plt.plot(signal.savgol_filter(corrected_baseline, 5, 1), color='orange')
    plt.plot(signal.savgol_filter(corrected_baseline, 7, 1), color='violet')
    plt.show()
    """
    
    if savgol_size>-1 and savgol_order>-1:
        corrected_baseline = signal.savgol_filter(corrected_baseline, savgol_size, savgol_order)

    #smooth = ndimage.white_tophat(savgol, None,
    #                np.repeat([1], int(round(raw.size * tophat_factor))))    

    return NormalizedTrace( signal=corrected_baseline, baseline = baseline )


@attr.s
class TraceChannel(object):
    dye_name = attr.ib()
    dye_wavelength = attr.ib()
    raw_channel = attr.ib()
    smooth_channel = attr.ib()


def ladder_area_means(ladders, fsa_list):

    # get mean size for each ladder step
    areas = [ [] for i in range(len(ladders)) ]
    for (fsa, fsa_index) in fsa_list:

        ladder = fsa.get_ladder_channel()

        for i in range(len(ladders)):

            alleles = ladder.get_alleles()
            min_index = min(range(len(alleles)),
                            key = lambda j: abs(ladders[i] - alleles[j].size))
            areas[i].append(alleles[min_index].area_bp)

            """
            # iterate through alleles and find the closest one            
            best_area = -1
            best_val = 99999
            for allele in ladder.get_alleles():
                val = abs(allele.size - ladders[i])
                if val < best_val:
                    best_val = val
                    best_area = allele.area_bp
            areas[i].append(best_area)
            """
            
    means = []
    for arr in areas:
        mu, sigma = norm.fit(arr)
        means.append(mu)

        #axes = plt.figure().add_subplot(111)
        #axes.hist(arr)
        #plt.show()

        #print("means: ", means)

    return means

def set_scale_factor(ladder, ladder_area_means):

    ladder_sizes = ladder.fsa.panel.get_ladder()['sizes']
    alleles = ladder.get_alleles()

    # calculate scale factor for each fsa
    sum_area_mean = 0.
    sum_area2 = 0.

    size_bp = []
    scale_factors = []

    for allele in alleles:
        # get ladder step closest to this allele
        min_index = min(range(len(ladder_sizes)),
                        key = lambda i: abs(ladder_sizes[i] - allele.size))
        
        sum_area_mean  += allele.area_bp * ladder_area_means[min_index] 
        sum_area2 += allele.area_bp * allele.area_bp

        size_bp.append(allele.size)
        scale_factors.append(ladder_area_means[min_index]/allele.area_bp)

    ladder.fsa.area_scale_factor = sum_area_mean / sum_area2 # this is a single scale factor for all ladders
    ladder.fsa.area_scale_factor_params = np.polyfit(size_bp, scale_factors, 1) # parameterization of s.f.

    #plt.figure()
    #plt.plot(size_bp, scale_factors)
    #plt.plot(size_bp, np.poly1d(ladder.fsa.area_scale_factor_params)(size_bp))
    #plt.axis([10,125,.625,2.55])
    #plt.xlabel('ladder size (base pair units)')
    #plt.ylabel('scale factor required to scale to mean ladder area')
    #plt.savefig('areasf_v_sizebp_'+ladder.fsa.filename[:-4]+'.png')
    #plt.show()
    #plt.close()


def b(txt):
    """ return a binary string aka bytes """
    return txt.encode('UTF-8')

def get_well_id( trace ):

    well_id = '00'

    try:
        well_id = trace.get_data(b('TUBE1')).decode('ASCII')

    except KeyError:
        pass

    return well_id

def separate_channels( trace, params):
    # return a list of [ 'dye name', dye_wavelength, numpy_array, numpy_smooth_baseline ]

    results = []
    for (idx, data_idx) in [ (1,1), (2,2), (3,3), (4,4), (5,105) ]:
        try:
            dye_name = trace.get_data(b('DyeN%d' % idx)).decode('UTF-8')
            dye_wavelength = trace.get_data(b('DyeW%d' % idx))
            raw_channel = np.array( trace.get_data(b('DATA%d' % data_idx)) )
            is_ladder = (idx==5)
            nt = normalize_baseline( raw_channel, params, is_ladder )

            results.append(
                TraceChannel(dye_name, dye_wavelength, raw_channel, nt.signal)
            )
        except KeyError:
            pass

    return results


def generate_scoring_function( strict_params, relax_params ):

    def _scoring_func( dp_result, method ):
        # alignment_result is (dp_score, dp_rss, dp_z, dp_peaks)
        dp_score = dp_result.dpscore
        dp_rss = dp_result.rss
        dp_peaks = dp_result.sized_peaks

        if method == 'strict':
            if ( dp_score >= strict_params['min_dpscore'] and
                    dp_rss <= strict_params['max_rss'] and
                    len(dp_peaks) >= strict_params['min_sizes'] ):
                return (1, None)
            return (0, None)
        elif method == 'relax':
            msg = []
            # scoring based on parts of results

            # score based on DP score compared to minimum DP score
            delta_score = relax_params['min_dpscore'] - dp_score
            if delta_score <= 0:
                dp_score_part = 1
            else:
                dp_score_part = 1e-2 ** (1e-2 * delta_score)

            # score based on RSS compared to the maximum allowed RSS
            delta_rss = dp_rss - relax_params['max_rss']
            if delta_rss <= 0:
                dp_rss_part = 1
            else:
                dp_rss_part = 1e-2 ** ( 1e-3 * delta_rss )
                msg.append( 'RSS > %d' % ( relax_params['max_rss'] ) )

            # score based on how many peaks we might miss compared to minimum number of peaks
            delta_peaks = relax_params['min_sizes'] - len(dp_peaks)
            if delta_peaks <= 0:
                dp_peaks_part = 1
            else:
                dp_peaks_part = max( 0, - delta_peaks / 0.5 * relax_params['min_sizes'] - 1)
                msg.append( 'Missing peaks = %d' % delta_peaks )

            # total overall score
            score = 0.3 * dp_score_part + 0.5 * dp_rss_part + 0.2 * dp_peaks_part
            return (score, msg)

        raise RuntimeError("Shouldn't be here!")

    return _scoring_func

def least_square( ladder_alleles, z ):

    """ 3rd order polynomial resolver
    """

    ladder_allele_sorted = SortedListWithKey( ladder_alleles, key = lambda k: k.rtime )
    x = [ p.rtime for p in ladder_allele_sorted ]
    y = [ p.size for p in ladder_allele_sorted ]

    for p in ladder_alleles:
        if p.qscore>=99:
            print("missing qscore for allele: ", p,"!")
            exit(6)
            
    f = np.poly1d(z)

    def _f( rtime ):
        size = f(rtime)
        # get the left-closest and right-closest ladder

        #left_idx = ladder_allele_sorted.bisect_key_left( rtime )
        right_idx = ladder_allele_sorted.bisect_key_right( rtime )
        left_idx = right_idx - 1

        # check if in range
        if left_idx<0:
            # do linear extrapolation based on 1st 3 ladder points
            z1 = np.polyfit(x[0:3], y[0:3], 1)
            left_ladder = np.poly1d(z1)(rtime)
            left_deviation = max((z.size - f(z.rtime)) for z in ladder_allele_sorted[0:3])**2
            left_qscore = min(z.qscore for z in ladder_allele_sorted[0:3])
        else:
            left_ladder = ladder_allele_sorted[left_idx]
            left_deviation = (left_ladder.size - f(left_ladder.rtime))**2
            left_qscore = left_ladder.qscore
            
        if right_idx >= len(ladder_alleles):
            # do linear extrapolation based on last 3 points
            z2 = np.polyfit(x[-3:], y[-3:], 1)
            left_ladder = np.poly1d(z2)(rtime)
            right_deviation = max((z.size - f(z.rtime)) for z in ladder_allele_sorted[-3:])**2
            right_qscore = min(z.qscore for z in ladder_allele_sorted[-3:])
        else:
            right_ladder = ladder_allele_sorted[right_idx]
            right_deviation = (right_ladder.size - f(right_ladder.rtime))**2
            right_qscore = right_ladder.qscore
            
        #cerr(' ==> rtime: %d/%4.2f  [ %d/%4.2f | %d/%4.2f ]' % ( rtime, size,
        #            left_ladder.rtime, left_ladder.size,
        #            right_ladder.rtime, right_ladder.size))

        return (size, (left_deviation + right_deviation) / 2,
                        min( left_qscore, right_qscore ),
                        const.allelemethod.leastsquare)

    return _f


def cubic_spline( ladder_alleles ):
    """ cubic spline interpolation
        x is peaks, y is standard size
    """

    ladder_allele_sorted = SortedListWithKey( ladder_alleles, key = lambda k: k.rtime )

    ladder_peaks = []
    ladder_sizes = []
    for ladder_allele in ladder_allele_sorted:
        ladder_peaks.append( ladder_allele.rtime )
        ladder_sizes.append( ladder_allele.size )
    f = UnivariateSpline(ladder_peaks, ladder_sizes, k=3, s=0)

    def _f( rtime ):
        size = f(rtime)

        right_idx = ladder_allele_sorted.bisect_key_right( rtime )
        left_idx = right_idx - 1

        if left_idx<0:
            # do linear extrapolation based on 1st 3 ladder points
            z1 = np.polyfit(ladder_peaks[0:3], ladder_sizes[0:3], 1)
            left_ladder = np.poly1d(z1)(rtime)
            left_deviation = max((z.size - f(z.rtime)) for z in ladder_allele_sorted[0:3])**2
            left_qscore = min(z.qscore for z in ladder_allele_sorted[0:3])
        else:
            left_ladder = ladder_allele_sorted[left_idx]
            left_deviation = (left_ladder.size - f(left_ladder.rtime))**2
            left_qscore = left_ladder.qscore
            
        if right_idx >= len(ladder_alleles):
            # do linear extrapolation based on last 3 points
            z2 = np.polyfit(ladder_peaks[-3:], ladder_sizes[-3:], 1)
            left_ladder = np.poly1d(z2)(rtime)
            right_deviation = max((z.size - f(z.rtime)) for z in ladder_allele_sorted[-3:])**2
            right_qscore = min(z.qscore for z in ladder_allele_sorted[-3:])
        else:
            right_ladder = ladder_allele_sorted[right_idx]
            right_deviation = (right_ladder.size - f(right_ladder.rtime))**2
            right_qscore = right_ladder.qscore

        return (size, (left_deviation + right_deviation) / 2,
                        min( left_qscore, right_qscore),
                        const.allelemethod.cubicspline)

    return _f


def local_southern( ladder_alleles ):
    """ southern local interpolation """

    ladder_allele_sorted = SortedListWithKey( ladder_alleles, key = lambda k: k.rtime )
    x = [ p.rtime for p in ladder_allele_sorted ]
    y = [ p.size for p in ladder_allele_sorted ]

    def southern3(xp, yp, rtime):

        x1, x2, x3 = xp[0], xp[1], xp[2]
        y1, y2, y3 = yp[0], yp[1], yp[2]

        L_denom = x1*(y2-y3) - x2*(y1-y3) + x3*(y1-y2)
        if abs(L_denom)<1.e-7:
            # return linear fit instead of local southern
            m = (y1-y2)/(x1-x2)
            b = y1 - m*x1
            size = m*rtime + b

        else:
            L_num = x1*y1*(y2-y3) - x2*y2*(y1-y3) + x3*y3*(y1-y2)
            L0 = L_num/L_denom
            M0 = (x1*(y1-L0) - x3*(y3-L0))/(y1-y3)
            c  = (y1 - L0) * (x1 - M0)
            if (rtime==M0):
                print("There's a problem... M0=",M0,", rtime=",rtime)
                return -9999
            size = c/(rtime - M0) + L0

        return size
    
    def _f( rtime, is_ladder ):
        """ return (size, deviation)
            deviation is calculated as delta square between curve1 and curve2
        """

        idx = ladder_allele_sorted.bisect_key_right( rtime )

        # special case for points to left or right of ladder steps
        if (idx==0):
            if is_ladder and (x[0] - rtime) > 100:
                return (-9999, 0, 0, const.allelemethod.localsouthern)
            else:
                z = np.polyfit( x[0:4], y[0:4], 2)
                min_score = .5 * min( z.qscore for z in ladder_allele_sorted[0:3] )
                return ( np.poly1d(z)(rtime), 0, min_score, const.allelemethod.localsouthern)
        if (idx==len(x)):
            if is_ladder and (rtime - x[-1]) > 100:
                return (-9999, 0, 0, const.allelemethod.localsouthern)
            else:
                z = np.polyfit( x[-4:], y[-4:], 2)
                min_score = .5 * min( z.qscore for z in ladder_allele_sorted[-3:] )
                return ( np.poly1d(z)(rtime), 0, min_score, const.allelemethod.localsouthern)
            
        # left curve
        if (idx>1 and idx<len(x)):
            size1 = southern3(x[idx-2:idx+1], y[idx-2:idx+1], rtime)
            min_score1 = min( z.qscore for z in ladder_allele_sorted[idx-2:idx+1] )
        else: # these should be only idx==1, so we do left curve from 1st 3 points
            size1 = southern3(x[0:3], y[0:3], rtime)
            min_score1 = .5 * min( z.qscore for z in ladder_allele_sorted[0:3] )

        # right curve
        if (idx<len(x)-2 and idx>0):
            size2 = southern3(x[idx:idx+3], y[idx:idx+3], rtime)
            min_score2 = min( z.qscore for z in ladder_allele_sorted[idx-1:idx+2] )
        else: # these should be only idx==len(x)-1, so we do right curve from last 3 points
            size2 = southern3(x[-3:], y[-3:], rtime)
            min_score2 = .5 * min( z.qscore for z in ladder_allele_sorted[-3:] )

        return ( (size1 + size2)/2, (size1 - size2) ** 2, (min_score1 + min_score2)/2,
                const.allelemethod.localsouthern)

    return _f


def mark_overlap_peaks(channels, params):

    # checking overlaps against channel !
    for channel in channels:
        for channel_r in channels:
            if channel == channel_r or channel.is_ladder():
                continue
            
            for p in channel.get_alleles(broad_peaks_only=False):

                if p.type == const.peaktype.noise:
                    continue
                
                if p.ertime - p.brtime < 3:
                    brtime = p.brtime
                    ertime = p.ertime
                elif p.ertime - p.brtime < 6:
                    brtime = p.brtime + 1
                    ertime = p.ertime - 1
                else:
                    brtime = p.brtime + 3
                    ertime = p.ertime - 3

                if brtime > p.rtime: brtime = p.rtime
                if ertime < p.rtime: ertime = p.rtime

                brtime = max(0, brtime)
                ertime = min(len(channel.data), len(channel_r.data), ertime)

                #cerr('checking %d | %s with channel %s' % (p.rtime, channel.dye,
                #            channel_r.dye))

                if ( channel.data[brtime] < channel_r.data[brtime] and
                     channel.data[ertime] < channel_r.data[ertime] and
                     p.height < channel_r.data[p.rtime] ):

                    # check how much is the relative height of this particular peak
                    rel_height = p.height / channel_r.data[p.rtime]
                    if rel_height > 1.0:
                        continue

                    (o_state, o_ratio, o_sym) = calc_overlap_ratio( channel.data,
                                                                    channel_r.data, p.rtime,
                                                                    brtime, ertime )

                    # if not really overlap, just continue reiteration
                    if not o_state:
                        continue

                    if is_verbosity(4):
                        print('peak: %d | %s | %s <> %f | %f | %f' % (p.rtime, channel.dye, p.type, rel_height, o_ratio, o_sym))
                    if rel_height < 0.15:
                        if p.type != const.peaktype.noise:
                            p.type = const.peaktype.overlap
                            if is_verbosity(4):
                                print('peak: %d | %s -> overlap' % (p.rtime, channel.dye))
                        p.qscore -= 0.10
                        continue

                    if ((rel_height < params.overlap_height_threshold and -0.5 < o_sym < 0.5) or
                        (o_ratio < 0.25 and -1.5 < o_sym < 1.5 ) or
                        (o_ratio < 0.75 and -0.5 < o_sym < 0.5 )):
                        if p.type != const.peaktype.noise:
                            p.type = const.peaktype.overlap
                            if is_verbosity(4):
                                print('peak: %d | %s -> overlap' % (p.rtime, channel.dye))
                        p.qscore -= 0.10
                        continue

def calc_overlap_ratio(data, data_r, rtime, brtime, ertime):
    """ calculate the difference ratio of overlapping area to the left and right area
        of rtime
        return ( boolean, total_ratio, log(right_ratio/left_ratio) )
    """

    lr = rr = 0.0
    lc = rc = 0
    for x in range(brtime, rtime+1):
        if data[x] > data_r[x]:
            return (False, 0, 0)
        lr += data[x]/data_r[x]
        lc += 1
    for x in range(rtime, ertime+1):
        if data[x] > data_r[x]:
            return (False, 0, 0)
        rr += data[x]/data_r[x]
        rc += 1
    lrc = lr / lc
    rrc = rr / rc

    log2vals = 0
    if lrc/rrc > 0:
        log2vals = math.log2(rrc/lrc)
        
    return (True, (lrc + rrc)/2, log2vals)

## this is a new algorithm and steps to perform peak analysis
##
## fsa = import_fsa()
## ladder_channel = fsa.ladder_channel()
## alleles = scan_peaks(ladder_channel, params)
## alleles = preannotate_peaks(ladder_channel, params)
## result = align_ladder(ladder_channel, params, size_standards)
##
## for channel in fsa.non_ladder_channel():
##     scan_peaks(channel, params)
##     preannotate_peaks(channel, params)
##     call_peaks(channel, params)
##     bin_peaks(channel, params)
##     postannotate_peaks(channel, params)

## the high level methods
##
##  fsa = import_fsa()
##  fsa.align_ladder(params.ladder)
##  fsa.scan_peaks(params.nonladder, marker=None)
##  fsa.preannotate_peaks(params.nonladder, marker=None)
##  fsa.call_peaks(params.nonladder, marker=None)
##  fsa.bin_peaks(params.nonladder, marker=None)

