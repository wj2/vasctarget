
def mouse_move(event):

    
    # plot best orientation at a given location
    if event.inaxes == ax1 or event.inaxes == ax2:
        i = np.argmin(damage[:,event.ydata,event.xdata])
        probe1.set_data(probe[i])
        extent = [event.xdata - kxy/2, 
                  event.xdata + kxy/2, 
                  event.ydata - kxy/2, 
                  event.ydata + kxy/2]
        probe1.set_extent(extent)
        ax2.set_ylim(512,0)
        ax2.set_xlim(0,512)
        print('damage=%1d'%np.min(damage[:,event.ydata,event.xdata]))
        
        ax3.clear()
        ax3.plot(lineprof[0, i, event.ydata, event.xdata])

        plt.draw()
    elif event.inaxes == ax4:
        i = np.argmin(damage[:,event.ydata,event.xdata])
        probe1.set_data(probe[i])
        extent = [event.xdata - kxy/2, 
                  event.xdata + kxy/2, 
                  event.ydata - kxy/2, 
                  event.ydata + kxy/2]
        probe1.set_extent(extent)
        ax2.set_ylim(512,0)
        ax2.set_xlim(0,512)
        print('damage=%1d'%np.min(lumdam[:,event.ydata,event.xdata]))

        plt.draw()
