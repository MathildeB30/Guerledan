import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

def print_(message, type_="info"):
    print(message)
    if type_ == "error":
        raise ValueError

def plot_map(map, longitude_bounds, latitude_bounds,
             title, nro_fig=None):
    if nro_fig is None:
        ff, ax = plt.subplots(clear=True)
    else:
        ff, ax = plt.subplots(num=nro_fig, clear=True)

    ax.imshow(map, extent=[longitude_bounds[0], longitude_bounds[1],
                           latitude_bounds[0], latitude_bounds[1]],
              aspect="equal")

    ax.set_xlabel("longitude (°)")
    ax.set_ylabel("latitude (°)")
    ax.set_title(title)
    ax.grid(True)
    return ff, ax

def plot_xy(data_x, data_y, x_label, y_label, title, symbol=None, label=None,
            nro_fig=None):
    if nro_fig is None:
        ff, ax = plt.subplots(clear=True)
    else:
        ff, ax = plt.subplots(num=nro_fig, clear=True)

    if symbol is None:
        ax.plot(data_x, data_y, label=label)
    else:
        ax.plot(data_x, data_y, symbol, label=label)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True)
    return ff, ax

def plot_xy_add(ax, data_x, data_y, symbol=None, label=None):
    if symbol is None:
        ax.plot(data_x, data_y, label=label)
    else:
        ax.plot(data_x, data_y, symbol, label=label)


def plot_logxlogy(data_x, data_y, x_label, y_label, title, symbol=None, nro_fig=None):
    if nro_fig is None:
        ff, ax = plt.subplots(clear=True)
    else:
        ff, ax = plt.subplots(num=nro_fig, clear=True)

    if symbol is None:
        ax.loglog(data_x, data_y)
    else:
        ax.loglog(data_x, data_y, symbol)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True)
    #plt.pause(0.01)
    return ff, ax

def plot_logxlogy_add(ax, data_x, data_y, symbol=None):
    if symbol is None:
        ax.loglog(data_x, data_y)
    else:
        ax.loglog(data_x, data_y, symbol)
    #plt.pause(0.01)

def plot_logxy(data_x, data_y, x_label, y_label, title, symbol=None, nro_fig=None):
    if nro_fig is None:
        ff, ax = plt.subplots(clear=True)
    else:
        ff, ax = plt.subplots(num=nro_fig, clear=True)

    ax = ff.add_subplot(111)
    if symbol is None:
        ax.semilogx(data_x, data_y)
    else:
        ax.semilogx(data_x, data_y, symbol)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True)
    #plt.pause(0.01)
    return ff, ax

def plot_logxy_add(ax, data_x, data_y, symbol=None):
    if symbol is None:
        ax.semilogx(data_x, data_y)
    else:
        ax.semilogx(data_x, data_y, symbol)
    #plt.pause(0.01)

    
def plot_image(img, x_label, y_label, title, nro_fig=None,
               xlim=None, ylim=None, zlim=None, q_origin="lower",
               colormap=None, q_grid=True, q_interpolation="linear",
               aspect="auto", q_colorbar=False):
    if nro_fig is None:
        ff, ax = plt.subplots(clear=True)
    else:
        ff, ax = plt.subplots(num=nro_fig, clear=True)

    if zlim is None:
        v_min = np.min(img[np.isfinite(img)])
        v_max = np.max(img[np.isfinite(img)])
    else:
        v_min = zlim[0]
        v_max = zlim[1]
    
    if xlim is None:
        x0 = -0.5
        x1 = img.shape[1] - 0.5
    else:
        x0 = xlim[0]
        x1 = xlim[1]

    if ylim is None:
        y0 = -0.5
        y1 = img.shape[0] - 0.5
    else:
        y0 = ylim[0]
        y1 = ylim[1]

    if colormap is None:
        cmap = None
    else:
        cmap = colormap
        
    y = ax.imshow(img, vmin=v_min, vmax=v_max, extent=[x0, x1, y0, y1],
                  aspect=aspect, cmap=colormap, origin=q_origin)
    
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    if q_grid == True:
        ax.grid(True)
    #plt.pause(0.01)
    if q_colorbar == True:
        plt.colorbar(y)
    return ff, ax

def plot_mesh2d(m, title, line_symbol=None, pnt_symbol=None, nro_fig=None,
                aspect="auto"):
    """ Tracé d'un maillage 2D """
    if nro_fig is None:
        ff, ax = plt.subplots(clear=True)
    else:
        ff, ax = plt.subplots(num=nro_fig, clear=True)

    ax = ff.add_subplot(111)
    ax.triplot(m.pnt[:,0], m.pnt[:,1], m.sP, marker=pnt_symbol,
               linestyle=line_symbol)
    ax.set_aspect(aspect)
    if title is not None:
        ax.set_title(title)
    return ff, ax

def get_Simrad_colorbar(q_reverse=False):
                
    color = np.array((
        (( 1,1,1), ( 1,1,1)),
        ((0.179831936955452,0.586554646492004,0.811764717102051),
         (0.179831936955452,0.586554646492004,0.811764717102051)),
        ((0.0431372560560703,0.517647087574005,0.780392169952393),
         (0.0431372560560703,0.517647087574005,0.780392169952393)),
        ((0.0215686280280352,0.507843136787415,0.390196084976196),
         (0.0215686280280352,0.507843136787415,0.390196084976196)),
        ((0,0.498039215803146,0), (0,0.498039215803146,0)),
        ((0.500000000000000,0.749019622802734,0),
         (0.500000000000000,0.749019622802734,0)),
        ((1,1,0), (1,1,0)),
        ((1,0.500000000000000,0), (1,0.500000000000000,0)),
        ((1,0,0), (1,0,0)),
        ((0.500000000000000,0,0), (0.500000000000000,0,0))), np.float32)

    n = color.shape[0]
    r = np.arange(n) / float(n-1)

    if q_reverse == True:
        red_data = color[::-1,:,0]
        green_data = color[::-1,:,1]
        blue_data = color[::-1,:,2]
    else:
        red_data = color[:,:,0]
        green_data = color[:,:,1]
        blue_data = color[:,:,2]

    color_dict = {}
    color_dict['red'] = np.hstack((r[:,None], red_data))
    color_dict['green'] = np.hstack((r[:,None], green_data))
    color_dict['blue'] = np.hstack((r[:,None], blue_data))
    
    return matplotlib.colors.LinearSegmentedColormap\
        ("simrad", color_dict, 1024)
