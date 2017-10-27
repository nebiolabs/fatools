#!/usr/bin/python
#

import argparse
from fatools.lib import params

def main():

    p = argparse.ArgumentParser('test_fatools')
    p.add_argument('--type', default='', help = "type of test")

    args = p.parse_args()

    print("args: ", args)
    
    trace_dir = "116"
    file_list = "05-M13ii-polD-5min.fsa"
    
    files = {} # dictionary containing directories and names of files within directories
    files['2'] = trace_dir

    # get FSA
    from fatools.lib.fileio.models import Marker, Panel, FSA

    Panel.upload(params.default_panels)
    Marker.upload(params.default_markers)
    
    panel = Panel.get_panel("GS120LIZ")
    fsa_list = []
    index = 1

    # set parameters for baseline correction
    from fatools.lib.const import allelemethod, baselinemethod
    
    _params = params.Params()
    _params.baselinewindow = 51
    _params.baselinemethod = baselinemethod.minimum
    _params.ladder.min_rfu = 500
    _params.ladder.min_rfu_ratio = 0.2
    
    for fsa_filename in file_list.split(','):

        fsa_filename = fsa_filename.strip()
        filename = trace_dir + "/" + fsa_filename

        fsa = FSA.from_file(filename, panel, _params, cache = False)
        fsa_list.append( (fsa, str(index)) )
        index += 1

    if args.type == 'allelemethods':
        
        import matplotlib.pyplot as plt
        import numpy as  np

        fig = plt.figure()
        
        (fsa, fsa_index) = fsa_list[0]

        print('D: aligning FSA %s' % fsa.filename)
        try:
            fsa.align(_params)
        except LadderMismatchException:
            print(("LadderMismatch: %s\n") % fsa.filename)
            
        c = fsa.get_ladder_channel()
        
        # get ladder and times for peaks fit to ladder
        ladder_sizes = fsa.panel.get_ladder()['sizes']
        alleles = c.get_alleles()
        allele_sizes = [allele.rtime for allele in alleles]
        
        plt.plot(allele_sizes, ladder_sizes, 'p',
                 label='peaks matched to ladder steps')

        for method in [ allelemethod.leastsquare, allelemethod.cubicspline, allelemethod.localsouthern ]:
        #for method in [ allelemethod.localsouthern ]:
            
            print("\nmethod: ", method)
            
            _params.allelemethod = method

            # call align again just to set the allelemethod
            print('D: aligning FSA %s' % fsa.filename)
            try:
                fsa.align(_params)
            except LadderMismatchException:
                print(("LadderMismatch: %s\n") % fsa.filename)

            func = fsa.allele_fit_func 

            # plot fit of ladder scan times to base pairs
            fit = np.poly1d(c.fsa.z)
            #x = np.arange(allele_sizes[0] - 150, allele_sizes[-1] + 100)  # len(c.data))
            x = np.arange(800, allele_sizes[-1] + 100)  # len(c.data))
            vecfunc = np.vectorize(func)
            print("vecfunc([1,2,3])=", vecfunc([1,2,3]))
            y_all  = vecfunc(x)
            plt.plot(x, vecfunc(x)[0], label=method)
            
        plt.legend()
        plt.xlabel("peak scan times")
        plt.ylabel("# base pairs")
        
        plt.show()


if __name__ == '__main__':

    main()
