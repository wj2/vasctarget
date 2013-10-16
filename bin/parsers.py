import argparse

def get_damage_quant_parser():
    parser = argparse.ArgumentParser(description='provide quantification of '
                                     'damage between pre and post datasets')
    parser.add_argument('pre', help='path to your 2photon *pre* data')
    parser.add_argument('post', help='path to your 2photon *post* data')
    parser.add_argument('probesize', help='length (first) and thickness '
                        '(second) of probe in microns', type=int, nargs=2)
    parser.add_argument('-i', '--interval', help='number of microns between '
                        'line profiles (default 50)', type=int, default=50)
    parser.add_argument('-n', '--n-profiles', help='number of profiles '
                        '(default 3)', type=int, default=3)
    parser.add_argument('-t', '--z-thickness', help='designate thickness (in '
                        'microns) of z-collapse in the sub-z stack (default 50)',
                        type=int, default=50)
    parser.add_argument('-d', '--debug', help='print debug pdfs of each mask on '
                        'rotated layers (default off)', action='store_true')
    parser.add_argument('-s', '--smooth', help='turn on extra smoothing for '
                        'data that has bleeds, may miss some smaller vessels '
                        'but will hopefully reduce false positives '
                        '(default off)', action='store_true')
    parser.add_argument('-c', '--channel', help='indicate a channel different '
                        'than Ch1 to use for analysis', nargs=1, default='Ch1')
    parser.add_argument('-r', '--room', help='inverse proportion of room for '
                        'cross correlation', default=8, type=int)

    return parser

def get_vasc_target_parser():
    parser = argparse.ArgumentParser(description='provide targeting for probe '
                                     'insertion')
    parser.add_argument('path', help='path to your 2photon data')
    parser.add_argument('probesize', help='length (first) and thickness '
                        '(second) of probe in microns', type=int, nargs=2)
    parser.add_argument('-g', '--gaussian', help='smooth line profiles with a '
                        'gaussian filter (default on)', action='store_false', 
                        default=True)
    parser.add_argument('-t', '--z-thickness', help='designate thickness (in '
                        'microns) of z-collapse in the sub-z stack (default 50)', 
                        type=int, default=50)
    parser.add_argument('-a', '--angles', help='number of angles to use '
                        '(default 6, max 180)', type=int, default=6)
    parser.add_argument('-d', '--downsample', help='number by which to '
                        'downsample line profile analysis (default 3)', 
                        type=int, default=3)
    parser.add_argument('-s', '--regionsize', help='size of region for best '
                        'options marks (default 50)', type=int, default=50)
    parser.add_argument('-b', '--best-in-region', help='toggle best-in-region '
                        'results rather than all results (default on)', 
                        action='store_false', default=True)
    parser.add_argument('-c', '--channel', help='indicate a channel different '
                        'than Ch1 to use for analysis', nargs=1, default='Ch1')

    return parser

def get_vasc_stats_parser():
    parser = argparse.ArgumentParser(description='provide targeting for probe '
                                     'insertion')
    parser.add_argument('paths', help='paths to your 2photon data', nargs='+')
    parser.add_argument('-g', '--gaussian', help='smooth line profiles with a '
                        'gaussian filter (default on)', action='store_false', 
                        default=True)
    parser.add_argument('-l', '--lengths', help='list lengths to use for stats',
                        nargs='+', type=int)
    parser.add_argument('-w', '--widths', help='list widths to use for stats',
                        nargs='+', type=int)
    parser.add_argument('-u', '--luminance', help='use luminance measure '
                        '(default off)', action='store_true')
    parser.add_argument('-L', '--luminance-only', help='use ONLY luminance '
                        '(default off)', action='store_true')
    parser.add_argument('-t', '--z-thickness', help='designate thickness (in '
                        'microns) of z-collapse in the sub-z stack (default 50)', 
                        type=int, default=50)
    parser.add_argument('-a', '--angles', help='number of angles to use '
                        '(default 6, max 180)', type=int, default=6)
    parser.add_argument('-d', '--downsample', help='number by which to '
                        'downsample line profile analysis (default 3)', type=int, 
                        default=3)
    parser.add_argument('-c', '--channel', help='indicate a channel different '
                        'than Ch1 to use for analysis', nargs=1, default='Ch1')


    return parser

    
