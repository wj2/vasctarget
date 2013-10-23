
import sys
import proc, analysis

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec
import matplotlib.colors as color

from matplotlib.widgets import Slider
from scipy.misc import imrotate
from random import choice

def make_redscale_cm():
    reddict = {'red' : [(0.0, 0.0, 0.0),
                        (1.0, 1.0, 1.0)],
               'green' : [(0.0, 0.0, 0.0),
                          (1.0, 0.0, 0.0)],
               'blue' : [(0.0, 0.0, 0.0),
                         (1.0, 0.0, 0.0)]}
    redscale = color.LinearSegmentedColormap('redscale', reddict)
    return redscale    

def make_red_alpha_scale_cm():
    redalphadict = {'red' : [(0.0, 0.0, 0.0),
                        (1.0, 1.0, 1.0)],
               'green' : [(0.0, 0.0, 0.0),
                          (1.0, 0.0, 0.0)],
               'blue' : [(0.0, 0.0, 0.0),
                         (1.0, 0.0, 0.0)],
               'alpha' : [(0.0, 0.0, 0.0),
                          (0.5, 0.0, 0.0),
                          (1.0, 1.0, 1.0)]}
    redalphascale = color.LinearSegmentedColormap('redalphascale', redalphadict)
    return redalphascale    

def make_greenscale_cm():
    reddict = {'red' : [(0.0, 0.0, 0.0),
                        (1.0, 0.0, 0.0)],
               'green' : [(0.0, 0.0, 0.0),
                          (1.0, 1.0, 1.0)],
               'blue' : [(0.0, 0.0, 0.0),
                         (1.0, 0.0, 0.0)]}
    greenscale = color.LinearSegmentedColormap('greenscale', reddict)
    return greenscale

def make_disp_probes(probesize, thetas):
    probe = proc.make_probe(probesize)
    return [imrotate(probe, t) for t in thetas]

def get_all_locs(lines, down):
    locs = {}
    for i in xrange(int(lines.min()), int(lines.max() + 1)):
        temp = np.where(lines == i)
        locs[i] = (temp[0], temp[1] * down, temp[2] * down)

    return locs

def get_best_locs(lines, down, regionsize, ori):
    regionsize = regionsize / down
    yregions = int(np.ceil(lines.shape[1] / regionsize))
    xregions = int(np.ceil(lines.shape[2] / regionsize))
    best = {}
    aos, axs, ays, mins = [], [], [], []
    for y in [y*regionsize for y in xrange(yregions)]:
        for x in [x*regionsize for x in xrange(xregions)]:
            # print y, y+regionsize, x, x+regionsize
            region = lines[:, y:y+regionsize, x:x+regionsize]
            rmin = region.min()
            regionmin = np.where(region == rmin)
            # get absolute index 
            ys = (regionmin[1] + y) * down; xs = (regionmin[2] + x) * down
            i = choice(xrange(ys.size))
            if rmin in best.keys():
                # print rmin, ys[i], xs[i]
                pr = best[rmin]
                pr[0].append(regionmin[0][i])
                pr[1].append(ys[i]); pr[2].append(xs[i])

            else:
                best[rmin] = ([regionmin[0][i]],
                               [ys[i]], [xs[i]])
    return best

def next_loc(curr, locs, mod):
    curr += mod
    while curr not in locs.keys():
        curr += mod
    return curr

def make_targeting_gui(stack, thetas, probesize, views, args):
    # find data to use
    if 'linegauss' in views.keys():
        dat = views['linegauss']
    else:
        dat = views['line']

    down = args.downsample
    
    # interaction functs!
    def slider_moved(newLoc):
        currloc[0] = int(np.floor(newLoc))
        index = (newLoc - np.floor(newLoc)) * len(locs[currloc[0]])
        currindex = [int(np.floor(index)) ]
        update_display(True)

    def key_modify(mod):
        currsize = len(locs[currloc[0]][0])
        locmax = max(locs)
        locmin = min(locs)
        print currloc, currindex, currsize, mod
        if currindex[0] + mod < 0:
            if currloc[0] + mod < locmin:
                currloc[0] = locmax
            else:
                currloc[0] = next_loc(currloc[0], locs, mod)
            currindex[0] = len(locs[currloc[0]][0]) - 1
        elif currindex[0] + mod > currsize - 1:
            if currloc[0] + mod > locmax:
                currloc[0] = locmin
            else:
                currloc[0] = next_loc(currloc[0], locs, mod)
            currindex[0] = 0
        else:
            currindex[0] += mod
        print currloc, currindex, len(locs[currloc[0]])
        update_display()

    def print_current():
        printfig = plt.figure()
        ax = printfig.add_subplot(111)
        axim = ax.imshow(stack[0], 'gray', interpolation='none')
        
        info, probe = get_curr_info()
        probeim = ax.imshow(probe, extent=[0,0,0,0])
        
        probeim.set_extent([info[2], info[2] + probe.shape[1],
                            info[1], info[1] + probe.shape[0]])
        ax.set_ylim(stack[0].shape[0], 0)
        ax.set_xlim(0, stack[0].shape[1])

        name = 'dam'+str(currloc[0])+'at'+str(info[1])+','+str(info[2])+'.pdf'
        printfig.savefig(name, format='pdf')

        plt.figure(1)

    def get_curr_info():
        choice = locs[currloc[0]]
        info = (choice[0][currindex[0]], choice[1][currindex[0]], 
                choice[2][currindex[0]])
        probe = dprobes[info[0]]

        return info, probe

    def update_display(fromSliderMove=False):
        info, p = get_curr_info()
        pim.set_data(p)
        pim.set_extent([info[2], info[2] + p.shape[1],
                        info[1], info[1] + p.shape[0]])
        meanvasMap.set_ylim(meanvas.shape[0], 0)
        meanvasMap.set_xlim(0, meanvas.shape[1])

        dotim.set_extent([(info[2] / down) - (dot.shape[1] / 2), 
                          (info[2] / down) + (dot.shape[1] / 2),
                          (info[1] / down) - (dot.shape[0] / 2), 
                          (info[1] / down) + (dot.shape[0] / 2)])

        lpdamMap.set_ylim(dat.shape[1], 0)
        lpdamMap.set_xlim(0, dat.shape[2])

        if not fromSliderMove:
            # sliderItem.set_val(currloc[0] + currindex[0] 
            #                    / len(locs[currloc[0]]))
            pass
        meanvasMap.set_title('damage: '+str(currloc[0]))

        plt.draw()

    def key_press(event):
        print event.key
        if event.key == 'left':
            key_modify(-1)
        elif event.key == 'right':
            key_modify(1)
        elif event.key == 'enter':
            print_current()
        elif event.key == 'escape':
            plt.close(fig)

    meanvas = stack.mean(axis=0)
    dprobes = make_disp_probes(probesize, thetas)
    
    # get location data
    if args.best_in_region:
        locs = get_best_locs(dat, down, args.regionsize, thetas)
    else:
        locs = get_all_locs(dat, down)

    currindex = [0]
    currloc = [int(min(locs.keys()))]

    # set up figure
    fig = plt.figure(1)
    fig.canvas.mpl_connect('key_press_event', key_press)

    # create grid to hold all items
    gs = gspec.GridSpec(6, 6)

    # line profile damage
    lpmean = round(dat.mean(), 2)
    lpstd = round(dat.std(), 2)
    lpmin = dat.min()
    lpmax = dat.max()
    lpdamMap = fig.add_subplot(gs[:3, :3])
    lpdamMap.set_title('mean: '+str(lpmean)+', std: '+str(lpstd)+
                       ', range: ['+str(lpmin)+','+str(lpmax)+']')
    lpdamMap.imshow(dat.min(axis=0), interpolation='none')
    
    # create dot
    dot = np.ones((5,5))
    dotim = lpdamMap.imshow(dot, 'spring',  extent=[0,0,0,0])

    # mean vasculature
    meanvasMap = fig.add_subplot(gs[:3, 3:])
    meanvasMap.imshow(meanvas, 'gray', interpolation='none')
    # create probe image
    pim = meanvasMap.imshow(dprobes[0], extent=[0,0,0,0])
    
    # histogram
    histoPlot = fig.add_subplot(gs[3:5, :])
    histoDat = plt.hist(dat.flatten(), 
                        bins=xrange(int(lpmin), 
                                    int(lpmax + 2)))
    
    # sliderbar
    sliderBar = fig.add_subplot(gs[5, :])
    sliderItem = Slider(sliderBar, '', 
                        lpmin, lpmax + 1,
                        valinit=currloc[0], closedmax=False)
    sliderItem.on_changed(slider_moved)

    # best location (maybe)
    update_display()

    # make it known!
    plt.draw()
    plt.show()
    plt.ioff()

    
def make_damage_gui(postack, psize):

    def probe_abs_move(loc):
        if locked[0]:
            locs[:] = np.array(loc)
        else:
            locs[i[0]] = np.array(loc)
        
        update_display()

    def probe_rel_move(trans):
        new = locs[i[0]] + trans
        if (new[0] >= 0 and new[0] < postack[i[0]].shape[0]
            and new[1] >= 0 and new[1] < postack[i[0]].shape[1]):
            if locked[0]:
                locs[:] += trans
            else:
                locs[i[0]] += trans

        update_display()

    def probe_rel_rotate(direc):
        if locked[0]:
            rots[:] += direc
        else:
            rots[i[0]] += direc
        
        update_display()

    def subz_rel_change(direc):
        if i[0] + direc < postack.shape[0] and i[0] + direc >= 0:
            i[0] += direc

        update_display()

    def coord_to_extent(yx):
        return [yx[1] - (p.shape[1] / 2),
                yx[1] + (p.shape[1] / 2),
                yx[0] - (p.shape[0] / 2),
                yx[0] + (p.shape[0] / 2)]

    def update_display():

        vasIm.set_data(postack[i[0]])
        pim.set_data(proc.clean_rotate(p, rots[i[0]]))
        pim.set_extent(coord_to_extent(locs[i[0]]))
        y = vasMap.set_ylim(postack[i[0]].shape[0], 0)
        x = vasMap.set_xlim(0, postack[i[0]].shape[1])
        vasMap.set_title('locked: '+str(locked[0])+' (toggle with r)')

        plt.draw()

    def mouse_click(event):
        probe_abs_move((event.ydata, event.xdata))

    def key_press(event):
        if event.key == 'j' or event.key == 'J' or event.key == 'left':
            probe_rel_move(np.array([0, -1]))
        elif event.key == 'l' or event.key == 'L' or event.key == 'right':
            probe_rel_move(np.array([0, 1]))
        elif event.key == 'i' or event.key == 'I' or event.key == 'up':
            probe_rel_move(np.array([-1, 0]))
        elif event.key == 'k' or event.key == 'K' or event.key == 'down':
            probe_rel_move(np.array([1, 0]))
        elif event.key == 'q' or event.key == 'Q':
            probe_rel_rotate(-1)
        elif event.key == 'e' or event.key == 'E':
            probe_rel_rotate(1)
        elif (event.key == 't' or event.key == 'ctrl+i' 
              or event.key == 'ctrl+up'):
            subz_rel_change(-1)
        elif (event.key == 'g' or event.key == 'ctrl+j' 
              or event.key == 'ctrl+down'):
            subz_rel_change(1)
        elif event.key == 'r' or event.key == 'R':
            locked[0] = not locked[0]
            update_display()
        elif event.key == 'enter':
            plt.close(fig)

        
    # if locked, then all layers move together
    locked = [True]
    # display top level initially
    i = [0]

    # initial probe location and rotation on each level
    locs = np.ones((postack.shape[0], 2))
    rots = np.zeros((postack.shape[0]))

    # create probe
    p = proc.make_probe(psize)

    # set up figure
    fig = plt.figure(1)
    fig.canvas.mpl_connect('key_press_event', key_press)
    fig.canvas.mpl_connect('button_press_event', mouse_click)

    # main display
    vasMap = fig.add_subplot(111)
    vasIm = vasMap.imshow(postack[i[0]], 'gray', interpolation='none')
    
    # probe display
    pim = vasMap.imshow(p, extent=[0,0,0,0], alpha=0.5)

    # make it so!
    update_display()
    plt.show()
    plt.ioff()
    plt.close()

    return locs, rots


def make_locate_gui(pre, post, offset):

    def coord_to_extent(p, y, x):
        return [x, x + p.shape[1],
                y + p.shape[0], y]

    def top_rel_move(adjust):
        offset[0] = offset[0] + adjust
        update_display()

    def on_press(event):
        press[0] = True 
        clicked[0] = np.array([event.ydata, event.xdata])
        orig_offset[0] = np.array([offset[0][0], offset[0][1]])

        update_display()

    def on_move(event):
        if press[0]:
            offset[0] = orig_offset[0] - (clicked[0] - np.array([event.ydata, 
                                                                 event.xdata]))

            update_display()

    def on_release(event):
        press[0] = False

    def update_display():

        base_im.set_data(post_mean)
        # top_im.set_data()
        top_im.set_extent(coord_to_extent(pre_mean, *offset[0]))

        y = base_map.set_ylim(0, post_mean.shape[0])
        x = base_map.set_xlim(0, post_mean.shape[1])

        plt.draw()

    def key_press(event):
        if event.key == 'j' or event.key == 'J' or event.key == 'left':
            top_rel_move(np.array([0, -1]))
        elif event.key == 'l' or event.key == 'L' or event.key == 'right':
            top_rel_move(np.array([0, 1]))
        elif event.key == 'i' or event.key == 'I' or event.key == 'up':
            top_rel_move(np.array([-1, 0]))
        elif event.key == 'k' or event.key == 'K' or event.key == 'down':
            top_rel_move(np.array([1, 0]))
        elif event.key == 'enter':
            plt.close(fig)


    press = [False] 
    offset = [offset]
    orig_offset = [np.array([offset[0][0], offset[0][1]])]
    clicked = [np.array([0,0])]

    fig = plt.figure(1)
    fig.canvas.mpl_connect('key_press_event', key_press)
    fig.canvas.mpl_connect('button_press_event', on_press)
    fig.canvas.mpl_connect('button_release_event', on_release)
    fig.canvas.mpl_connect('motion_notify_event', on_move)
    
    pre_mean = np.mean(pre, axis=0)
    post_mean = np.mean(post, axis=0)

    base_map = fig.add_subplot(111)
    redsc = make_redscale_cm()
    greensc = make_greenscale_cm()
    base_im = base_map.imshow(post_mean, redsc, interpolation='none', 
                              alpha=1.0)
    top_im = base_map.imshow(pre_mean, greensc, interpolation='none', alpha=0.5,
                             extent=coord_to_extent(pre_mean, *offset[0]))

    update_display()
    plt.show()
    plt.ioff()
    plt.close()

    return offset[0]

def make_thresh_setter(prof, low, up, mask, section, buff, smooth):

    sect_mean = section.mean()
    sect_std = section.std()

    renormed = [False]
    recalc = [True]
    on = [True]

    def coord_to_extent(p, y, x):
        return [x, x + p.shape[1],
                y + p.shape[0], y]

    def renorm_data():
        if renormed[0]:
            prof[0] = ((prof[0] * sect_std) / sect_mean)
            renormed[0] = False
        else:
            prof[0] = ((prof[0] * sect_mean) / sect_std)
            renormed[0] = True

        update_display(force_recalc=True)

    def recalc_prof():
        print 'recalc'
        dam_re, line_re, smooth_re = analysis.profile(prof[0], thresh=low[0], 
                                                      upper=up[0], smooth=smooth)
        dam[0] = dam_re
        smoothed[0] = smooth_re
        line[0] = line_re
        last_low[0] = low[0]; last_up[0] = up[0]

    def adjust_threshs(low_delt, up_delt):
        print 'adjust'
        print low_delt, up_delt
        low[0] += low_delt
        up[0] += up_delt
        
        update_display()

    def update_display(force_recalc=False): 
        print 'update'
        print last_low, low, last_up, up
        diff = (last_low[0] != low[0] or last_up[0] != up[0])

        diffd = '*' if diff else ''

        if (recalc[0] and diffd) or force_recalc:
            recalc_prof()
        
        print last_low, low, last_up, up
        print 'finish update'
        profile_graph.clear()
        profile_graph.plot(smoothed[0])
        profile_graph.axhline(low[0], color='r')
        profile_graph.axhline(up[0], color='g')
        
        profile_graph.set_title('vessels: '+str(dam[0]))
        section_map.set_title('lower thresh: '+str(low[0])+' | '
                              'upper thresh: '+str(up[0])+' | '
                              'renormed: '+str(renormed[0])+' | '
                              'constant recalc: '+str(recalc[0])+diffd)
        print 'reset titles'
        plt.draw()
                              
    def toggle_recalc():
        recalc[0] = not recalc[0]
        update_display()    

    def key_press(event):
        print 'key press'
        print event.key
        if event.key == 'j' or event.key == 'J' or event.key == 'left':
            adjust_threshs(0, -.1)
        elif event.key == 'l' or event.key == 'L' or event.key == 'right':
            adjust_threshs(0, .1)
        elif event.key == 'i' or event.key == 'I' or event.key == 'up':
            adjust_threshs(.1, 0)
        elif event.key == 'k' or event.key == 'K' or event.key == 'down':
            adjust_threshs(-.1, 0)
        elif event.key == 'r' or event.key == 'R':
            renorm_data()
        elif event.key in ['c', 'C']:
            toggle_recalc()
        elif event.key in ['space']:
            update_display(force_recalc=True)
        elif event.key == 'enter':
            plt.close(fig)
        elif event.key in ['escape']:
            on[0] = False
            plt.close(fig)       


    prof = [prof]

    dam, line, smoothed = analysis.profile(prof[0], thresh=low, upper=up, 
                                           smooth=smooth)

    low = [low]; up = [up]; smoothed = [smoothed]; dam = [dam]; line = [line]
    last_low = [low[0]]; last_up = [up[0]]

    fig = plt.figure(1)
    fig.canvas.mpl_connect('key_press_event', key_press)
    
    gs = gspec.GridSpec(6, 6)

    section_map = fig.add_subplot(gs[:4, :])
    section_map.imshow(section, 'gray', interpolation='none')

    r_alpha_cm = make_red_alpha_scale_cm()
    if min(mask.shape) > 1:
        section_map.imshow(mask, r_alpha_cm, interpolation='none', 
                           extent=coord_to_extent(mask, buff, buff),
                           alpha=0.5)
    else:
        section_map.imshow(mask,  interpolation='none', 
                           extent=coord_to_extent(mask, buff, buff),
                           alpha=0.5)
        
    section_map.set_xlim(0, section.shape[1])
    section_map.set_ylim(0, section.shape[0])
    
    profile_graph = fig.add_subplot(gs[5:, :])
    profile_graph.plot(smoothed[0])
    profile_graph.axhline(low[0], color='r')
    profile_graph.axhline(up[0], color='g')

    update_display()
    plt.show()
    plt.ioff()
    plt.close()
    
    return dam[0], line[0], on[0]

    
