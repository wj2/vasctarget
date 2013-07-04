
import sys
import proc

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec

from matplotlib.widgets import Slider
from scipy.misc import imrotate
from random import choice

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
            print y, y+regionsize, x, x+regionsize
            region = lines[:, y:y+regionsize, x:x+regionsize]
            rmin = region.min()
            regionmin = np.where(region == rmin)
            # get absolute index 
            ys = (regionmin[1] + y) * down; xs = (regionmin[2] + x) * down
            i = choice(xrange(ys.size))
            if rmin in best.keys():
                print rmin, ys[i], xs[i]
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
