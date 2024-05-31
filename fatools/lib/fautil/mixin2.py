
import itertools
from fatools.lib.fautil import algo2 as algo
from fatools.lib.utils import cout, cerr, cexit, is_verbosity
from fatools.lib import const

from functools import lru_cache

import attr
import time

# FA

class AlleleMixIn(object):

    def __repr__(self):
        return "<A: %7s | %3d | %4d | %5d | %2d | %+3.1f | %4.1f | %5.1f | %6d | %4.2f >" % (
                    self.type, self.size, self.rtime, self.rfu, self.wrtime,
                    self.srtime, self.beta, self.theta, self.omega, self.dev
        )

    @property
    def height(self):
        return self.rfu

    @property
    def height_uncorr(self):
        return self.rfu_uncorr
        
    def __lt__(self, other):
        return self.rtime < other.rtime


class ChannelMixIn(object):
    """
    attrs: alleles, fsa, status
    """

    __slots__ = [   'data', 'dye', 'wavelen', 'alleles', 'fsa', 'status', 'marker',
                    'mma', 'mmb', 'p80', 'basepairs', 'smeared_alleles'
                ]

    
    def __init__(self):
        self.basepairs = []
        self.smeared_alleles = []
        
    def add_allele(self, allele):
        """ add this allele to channel """
        raise NotImplementedError()


    def is_ladder(self):
        return self.marker.code == 'ladder'


    def assign(self):
        marker = self.fsa.panel.get_marker(self.dye)
        if marker.label.lower() in self.fsa.excluded_markers:
            self.marker = self.fsa.panel.get_marker('undefined')
        else:
            self.marker = marker


    def get_alleles(self, broad_peaks_only=True):
        if self.status == const.channelstatus.reseted:
            # create alleles first
            raise RuntimeError('E: channel needs to be scanned first')
        if broad_peaks_only:
            return [ allele for allele in self.alleles if allele.type=='broad' ]
        else:
            return self.alleles

    # ChannelMixIn scan method
    def scan(self, parameters):

        if self.status != const.channelstatus.reseted:
            return
        
        alleles = algo.scan_peaks(self, parameters)


    def preannotate(self, parameters):

        params = parameters.ladder if self.is_ladder() else parameters.nonladder

        # determine qscore and mark peak type for alleles
        algo.preannotate_peaks(self, params)

    # ChannelMixIn align method
    def align(self, parameters, ladder=None, anchor_pairs=None, saturated_peak_rtimes=[]):
        # sanity checks
        if self.marker.code != 'ladder':
            raise RuntimeError('E: align() must be performed on ladder channel!')

        ladder = self.fsa.panel.get_ladder()

        # prepare ladder qcfunc
        if 'qcfunc' not in ladder:
            ladder['qcfunc'] =  algo.generate_scoring_function(
                                            ladder['strict'], ladder['relax'] )

        start_time = time.process_time()
        result = algo.align_peaks(self, parameters, ladder, anchor_pairs, saturated_peak_rtimes)
        dpresult = result.dpresult
        fsa = self.fsa
        fsa.z = dpresult.z
        fsa.rss = dpresult.rss
        fsa.nladder = len(dpresult.sized_peaks)
        fsa.score = result.score
        fsa.duration = time.process_time() - start_time

        # set allele sizes from ladder steps
        alleles = self.get_alleles()
        alleles.sort(key = lambda x: x.rtime)

        ladder_sizes = ladder['sizes']
        ladder_sizes.sort()

        for allele, ladder_size in zip(alleles, ladder_sizes):
            allele.size = ladder_size

        # check the allele method
        method = parameters.allelemethod

        if method == const.allelemethod.leastsquare:
            fsa.allele_fit_func = algo.least_square( alleles, self.fsa.z)
        elif method == const.allelemethod.cubicspline:
            fsa.allele_fit_func = algo.cubic_spline( alleles )
        elif method == const.allelemethod.localsouthern:
            fsa.allele_fit_func = algo.local_southern( alleles )
        else:
            raise RuntimeError
        
        #min_rtime = ladders[1].rtime
        #max_rtime = ladders[-2].rtime
        fsa.min_rtime = parameters.ladder.min_rtime
        fsa.max_rtime = parameters.ladder.max_rtime
    
        #import pprint; pprint.pprint(dpresult.sized_peaks)
        #print(fsa.z)
        if is_verbosity(4):
            cout('O: Score %3.2f | %5.2f | %d/%d | %s | %5.1f | %s' %
                 (fsa.score, fsa.rss, fsa.nladder, len(ladder['sizes']),
                  result.method, fsa.duration, fsa.filename) )


    # ChannelMixIn call method
    def call(self, parameters, ladder):

        params = parameters.ladder if self.is_ladder() else parameters.nonladder

        # skip if ladder
        if self.marker.code == 'ladder':
            pass
        
        algo.call_peaks(self, params, self.fsa.allele_fit_func,
                        self.fsa.min_rtime, self.fsa.max_rtime)
            
        #algo.bin_peaks(self, params, self.marker)
        #algo.postannotate_peaks(self, params)
        
        #import pprint; pprint.pprint(alleles)


    # ChannelMixIn merge method
    def merge(self, parameters, ladder, plot=False):

        if self.is_ladder():
            return

        params = parameters.nonladder

        #print("calling algo.merge_peaks for ",self.dye)
        
        self.smeared_alleles = algo.merge_peaks(self, params, self.fsa.allele_fit_func, plot)


    # ChannelMixIn normalize method
    def normalize(self, parameters):

        params = parameters.ladder if self.is_ladder() else parameters.nonladder

        algo.normalize_peaks(self, params)


    def get_basepairs(self):

        if not self.fsa.allele_fit_func:
            return []
        
        if not self.basepairs:
            #min=0
            #max=0
            for st in range(len(self.data)):
                basepair = self.fsa.allele_fit_func(st)[0]
                # get max/min for plotting
                #if basepair>-999 and min==0:
                #    min=st
                #if min>0 and max==0 and basepair<-999:
                #    max=st-1
                self.basepairs.append(basepair)
        """
        import matplotlib.pyplot as plt
        import numpy as np    
        plt.plot(np.arange(0,len(self.data))[min:max+1], self.basepairs[min:max+1])
        plt.xlabel('scan times')
        plt.ylabel('base pairs')
        plt.show()
        """
        return self.basepairs
    
class FSAMixIn(object):
    """
    attrs: channels
    """

    __slots__ = [   'panel', 'channels', 'excluded_markers', 'filename',
                    'date', 'rss', 'z', 'score', 'nladder', 'duration', 'status',
                    'allele_fit_func', 'area_scale_factor_params', 'area_scale_factor',
                    'min_rtime', 'max_rtime', 'scan_done', 'align_done', 'call_done'
                ]

    def __init__(self):
        self.area_scale_factor_params = [] # array of parameters for polynomial fit of scale factors vs. ladder sie
        self.area_scale_factor = -1 # single scale factor for all ladders
        self.scan_done = False
        self.call_done = False
        self.allele_fit_func = None
        
    def get_data_stream(self):
        """ return stream of data """
        raise NotImplementedError()


    def add_channel(self, trace_channel):
        """ add this channel to fsa """
        raise NotImplementedError()

    def get_well_id(self):
        return algo.get_well_id(self.get_trace())

    def get_trace(self):
        if not hasattr(self, '_trace'):
            from fatools.lib.fautil import traceio
            self._trace = traceio.read_abif_stream( self.get_data_stream() )
            self.close_file()
            
        return self._trace


    def set_panel(self, panel, options=None):
        """ set panel and add excluded markers """
        self.panel = panel

        if not options:
            return

        excluded_markers = [ x.strip() for x in options['exclude'].split(',') ]

        for marker_code in excluded_markers:
            if panel.has_marker(marker_code):
                self.excluded_markers.append(marker_code.lower())


    def create_channels(self, params):
        if is_verbosity(4): cerr('I: Generating channels for %s' % self.filename)
        trace = self.get_trace()
        trace_channels = algo.separate_channels(trace, params)
        for tc in trace_channels:
            channel = self.Channel(data=tc.smooth_channel, dye=tc.dye_name,
                        wavelen=tc.dye_wavelength,
                        status = const.channelstatus.reseted,
                        fsa=self)
            self.add_channel(channel)
            #channel.status = const.channelstatus.reseted

    # FSAMixIn scan method
    def scan(self, parameters):

        if self.scan_done: return

        for c in self.channels:
            c.scan(parameters)
            c.preannotate(parameters)

        algo.mark_overlap_peaks(self.channels, parameters.nonladder)

        self.scan_done = True


    # FSAMixIn align method
    def align(self, parameters=None):

        self.scan(parameters)

        ladder = self.get_ladder_channel()
        ladder.align(parameters, None, None, self.get_saturated_rtimes())

    # FSAMixIn call method
    def call(self, parameters):

        if self.call_done: return
        
        ladder = self.get_ladder_channel()

        self.scan(parameters)
        
        for c in self.channels:
            c.call(parameters, ladder)

        self.call_done = True

    # FSAMixIn merge method
    def merge(self, parameters, plot=False):

        ladder = self.get_ladder_channel()
        
        for c in self.channels:
            #print("calling merge for channel: ", c.dye)
            c.merge(parameters, ladder, plot)

    # FSAMixIn call method
    def normalize(self, parameters, ladder_means):

        self.call(parameters)

        # get scale factor for this FSA using ladder channel
        algo.set_scale_factor(self.get_ladder_channel(), ladder_means)
        
        for c in self.channels:
            c.normalize(parameters)

    def get_ladder_channel(self):

        for c in self.channels:
            if c.marker.code == 'ladder':
                return c
        raise RuntimeError('E: ladder channel not found')
    
    def get_saturated_rtimes(self):
        saturation_threshold = 29000
        #makes a flat list of alleles from all channels
        alleles = itertools.chain(*[c.get_alleles() for c in self.channels])
        saturated_rtimes = [a.rtime for a in alleles if a.rfu >= saturation_threshold]
        return saturated_rtimes



class MarkerMixIn(object):
    """
    attrs: id, code, species, min_size, max_size, repeats, z_params
    """

    def update(self, obj):

        if isinstance(obj, dict):
            if hasattr(self, 'id') and self.id is not None:
                raise RuntimeError('ERR: can only update from dictionary on new instance')
            if 'code' in obj:
                self.code = obj['code']
            if 'species' in obj:
                self.species = obj['species']
            if 'min_size' in obj:
                self.min_size = obj['min_size']
            if 'max_size' in obj:
                self.max_size = obj['max_size']
            if 'repeats' in obj:
                self.repeats = obj['repeats']
            if 'z_params' in obj:
                self.z_params = obj['z_params']

            # field 'related_to' must be handled by implementation object

        else:
            if self.code != obj.code:
                raise RuntimeError('ERR: attempting to update Marker with different code')
            if obj.min_size is not None:
                self.min_size = obj.min_size
            if obj.max_size is not None:
                self.max_size = obj.max_size
            if obj.repeats is not None:
                self.repeats = obj.repeats
            if obj.related_to is not None:
                self.related_to = obj.related_to
            if obj.z_params is not None:
                self.z_params = obj.z_params

    @property
    def label(self):
        return '%s/%s' % (self.species, self.code)

    @classmethod
    def from_dict(cls, d):
        marker = cls()
        marker.update(d)
        return marker

    @lru_cache(maxsize=32)
    def get_sortedbins(self, batch):
        # get Bin from this batch
        bin = self.get_bin(batch)
        return bin.sortedbins


class PanelMixIn(object):
    """
    attrs: id, code, data, dyes, Marker
    """

    def __init__(self):
        self.ladder_area_means = []


    def set_ladder_dye(self, ladder):
        raise NotImplementedError()


    def update(self, obj):

        if isinstance(obj, dict):
            self.code = obj['code']
            self.data = obj['data']

        else:
            raise NotImplementedError()


    def get_markers(self):
        """ return all instance of markers """
        return [ self.Marker.get_marker(code) for code in self.data['markers'] ]


    def get_marker(self, dye):
        if not hasattr(self, '_dyes') or not self._dyes:
            self._dyes = {}
            for m in self.get_markers():
                self._dyes[ self.data['markers'][m.label]['dye'] ] = m
            ladder = self.get_ladder()
            self._dyes[ ladder['dye'] ] = self.Marker.get_marker('ladder')
        return self._dyes[dye]


    def get_ladder(self):
        return const.ladders[ self.data['ladder']]


    def get_ladder_area_means(self, fsa_list):

        ladders = self.get_ladder()['sizes']

        if not self.ladder_area_means:
            self.ladder_area_means = algo.ladder_area_means(ladders, fsa_list)

        return self.ladder_area_means

    @classmethod
    def from_dict(cls, d):
        panel = cls()
        panel.update(d)
        return panel


# Sample


class SampleMixIn(object):
    pass


class BatchMixIn(object):
    pass
