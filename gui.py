
import sys
import proc

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspec

from matplotlib.widgets import Slider
from scipy.misc import imrotate

def make_disp_probes(probesize, thetas):
    probe = proc.make_probe(probesize)
    return [imrotate(probe, t) for t in thetas]

def get_all_locs(lines, down):
    locs = {}
    print lines.min(), lines.max()
    sys.stdout.write('constructing range...\n')
    for i in xrange(int(lines.min()), int(lines.max() + 1)):
        temp = np.where(lines == i)
        locs[i] = (temp[0], temp[1] * down, temp[2] * down)

    return locs


def make_gui(stack, thetas, probesize, views, args):
    currindex = [0]
    currloc = [int(views['line'].min())]
    down = args.downsample
    
    # interaction functs, using closure!
    def slider_moved(newLoc):
        print 'moved'
        currindex = [0]
        currloc[0] = int(round(newLoc, 0))
        update_display()

    def key_modify(mod):
        currsize = len(allLocs[currloc[0]][0])
        print currloc, currindex, currsize, mod
        if currindex[0] + mod < 0 or currindex[0] + mod > currsize:
            currloc[0] += mod
            currindex[0] = 0
            if len(allLocs[currloc[0]][0]) < 0:
                key_modify(mod)
        else:
            currindex[0] += mod
        print currloc, currindex, len(allLocs[currloc[0]])
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
        choice = allLocs[currloc[0]]
        info = (choice[0][currindex[0]], choice[1][currindex[0]], 
                choice[2][currindex[0]])
        probe = dprobes[info[0]]

        return info, probe

    def update_display():
        info, p = get_curr_info()
        pim.set_data(p)
        pim.set_extent([info[2], info[2] + p.shape[1],
                        info[1], info[1] + p.shape[0]])
        meanvasMap.set_ylim(meanvas.shape[0], 0)
        meanvasMap.set_xlim(0, meanvas.shape[1], 0)

        sys.stdout.write('damage: '+str(currloc[0])+'\n')

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
            quit_targetting()

    meanvas = stack.mean(axis=0)
    dprobes = make_disp_probes(probesize, thetas)
    
    # get histogram-rel data
    allLocs = get_all_locs(views['line'], down)

    # set up figure
    fig = plt.figure(1)
    fig.canvas.mpl_connect('key_press_event', key_press)

    # create grid to hold all items
    gs = gspec.GridSpec(6, 6)
    
    # line profile damage
    lpdamMap = fig.add_subplot(gs[:3, :3])
    lpdamMap.imshow(views['line'].min(axis=0), interpolation='none')

    # mean vasculature
    meanvasMap = fig.add_subplot(gs[:3, 3:])
    meanvasMap.imshow(meanvas, 'gray', interpolation='none')
    # create probe image
    pim = meanvasMap.imshow(dprobes[0], extent=[0,0,0,0])
    
    # histogram
    histoPlot = fig.add_subplot(gs[3:5, :])
    histoDat = plt.hist(views['line'].flatten(), bins=len(allLocs))
    
    # sliderbar
    sliderBar = fig.add_subplot(gs[5, :])
    sliderItem = Slider(sliderBar, '', 
                        views['line'].min(), views['line'].max(),
                        valinit=currloc[0])
    sliderItem.on_changed(slider_moved)

    # best location (maybe)
    update_display()

    # make it known!
    plt.draw()
    plt.show()
    plt.ioff()
    
        
