
from fatools.lib.fautil import algo2 as algo
from fatools.lib.utils import cout, cerr, cexit
from fatools.lib import const

from functools import lru_cache

import attr
import time

# FA


class AlleleMixIn(object):

    def __repr__(self):
        return "<A: %3d | %4d | %5d | %2d | %+3.1f | %4.1f | %5.1f | %6d | %4.2f >" % (
                    self.size, self.rtime, self.rfu, self.wrtime,
                    self.srtime, self.beta, self.theta, self.omega, self.dev
        )

    @property
    def height(self):
        return self.rfu

    def __lt__(self, other):
        return self.rtime < other.rtime


class ChannelMixIn(object):
    """
    attrs: alleles, fsa, status
    """

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


    def get_alleles(self):       
        if self.status == const.channelstatus.reseted:
            # create alleles first
            raise RuntimeError('E: channel needs to be scanned first')
        return self.alleles

    # ChannelMixIn scan method
    def scan(self, parameters, ladder=None):
        
        if self.status != const.channelstatus.reseted:
            return
        
        params = parameters.ladder if self.is_ladder() else parameters.nonladder
        
        alleles = algo.scan_peaks(self, params)

        if len(alleles)==0: return
        
        # determine qscore and mark peak type for alleles
        algo.preannotate_peaks(self, params)
        
        if ladder==None:
            return
        
        # get function for call_peaks        
        ladders = ladder.alleles

        # check the allele method
        method = parameters.allelemethod
        
        if method == const.allelemethod.leastsquare:
            func = algo.least_square( ladders, self.fsa.z)
        elif method == const.allelemethod.cubicspline:
            func = algo.cubic_spline( ladders )
        elif method == const.allelemethod.localsouthern:
            func = algo.local_southern( ladders )
        else:
            raise RuntimeError
        
        #min_rtime = ladders[1].rtime
        #max_rtime = ladders[-2].rtime
        
        min_rtime = params.min_rtime
        max_rtime = ladders[-1].rtime
        print("min/max rtime: ",min_rtime,", ", max_rtime)
        
        algo.call_peaks(self, params, func, min_rtime, max_rtime)
        #algo.bin_peaks(self, params, self.marker)
        #algo.postannotate_peaks(self, params)
        
        #import pprint; pprint.pprint(alleles)


    def preannotate(self, parameters):
        pass

    # ChannelMixIn align method
    def align(self, parameters, anchor_pairs=None):

        # sanity checks
        if self.marker.code != 'ladder':
            raise RuntimeError('E: align() must be performed on ladder channel!')

        self.scan( parameters )         # in case this channel hasn't been scanned

        ladder = self.fsa.panel.get_ladder()

        # prepare ladder qcfunc
        if 'qcfunc' not in ladder:
            ladder['qcfunc'] =  algo.generate_scoring_function(
                                            ladder['strict'], ladder['relax'] )

        start_time = time.process_time()
        result = algo.align_peaks(self, parameters, ladder, anchor_pairs)
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

        if (len(alleles) != len(ladder_sizes)):
            cerr("alleles not same length as ladder!")
            import sys
            sys.exit() # exit for now because this code can't handle the wrong number of peaks yet

        for allele, ladder_size in zip(alleles, ladder_sizes):
            allele.size = ladder_size
            
        #import pprint; pprint.pprint(dpresult.sized_peaks)
        #print(fsa.z)
        cout('O: Score %3.2f | %5.2f | %d/%d | %s | %5.1f | %s' %
            (fsa.score, fsa.rss, fsa.nladder, len(ladder['sizes']),
             result.method, fsa.duration, fsa.filename) )


    # ChannelMixIn call method
    def call(self, ladder, parameters=None):

        # skip if ladder
        if self.marker.code == 'ladder':
            pass
        
        if parameters:
            self.scan( parameters, ladder)  # in case this channel hasn't been scanned



class FSAMixIn(object):
    """
    attrs: channels
    """

    __slots__ = [   'panel', 'channels', 'excluded_markers', 'filename',
                    'rss', 'z', 'score', 'nladder', 'duration',
                ]

    def get_data_stream(self):
        """ return stream of data """
        raise NotImplementedError()


    def add_channel(self, trace_channel):
        """ add this channel to fsa """
        raise NotImplementedError()


    def get_trace(self):
        if not hasattr(self, '_trace'):
            from fatools.lib.fautil import traceio
            self._trace = traceio.read_abif_stream( self.get_data_stream() )
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


    def create_channels(self):
        cerr('I: Generating channels for %s' % self.filename)
        trace = self.get_trace()
        trace_channels = algo.separate_channels(trace)
        for tc in trace_channels:
            channel = self.Channel(data=tc.smooth_channel, dye=tc.dye_name,
                        wavelen=tc.dye_wavelength,
                        status = const.channelstatus.reseted,
                        fsa=self)
            self.add_channel(channel)
            #channel.status = const.channelstatus.reseted


    def align(self, parameters=None):

        c = self.get_ladder_channel()
        c.align( parameters )

    # FSAMixIn call method
    def call(self, parameters):

        ladder = self.get_ladder_channel()
        
        for c in self.channels:
            if c.marker.code != 'ladder':
                c.scan( parameters, ladder)

    def get_ladder_channel(self):

        for c in self.channels:
            if c.marker.code == 'ladder':
                return c
        raise RuntimeError('E: ladder channel not found')



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