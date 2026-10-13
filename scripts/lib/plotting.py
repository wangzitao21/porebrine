from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .fields import COLUMNS, MINERALS, frames, read_frame


def plot_result(run_directory):
    run_directory = Path(run_directory)
    stamps = frames(run_directory)
    first, last = (read_frame(run_directory, stamp) for stamp in (stamps[0], stamps[-1]))
    nx, ny = len(np.unique(first[:, 0])), len(np.unique(first[:, 1]))
    if nx*ny != len(first) or not np.array_equal(first[:, :2], last[:, :2]):
        raise ValueError('Plotting requires one complete, unchanged rectangular grid')
    x, y = np.unique(first[:, 0]), np.unique(first[:, 1])
    dx, dy = np.median(np.diff(x)), np.median(np.diff(y))
    extent = np.array([x[0]-dx/2, x[-1]+dx/2, y[0]-dy/2, y[-1]+dy/2])*1e6
    changes = np.array([last[:, COLUMNS.index(name)]-first[:, COLUMNS.index(name)] for name in MINERALS])
    active = np.argmax(np.max(np.abs(changes), axis=1))
    mineral = MINERALS[active] if np.any(changes) else 'Total mineral'
    change = changes[active] if np.any(changes) else changes.sum(axis=0)
    speed = np.hypot(last[:, COLUMNS.index('Ux_m_s')], last[:, COLUMNS.index('Uy_m_s')])
    images = [first[:, 2], last[:, 2], speed, change]
    titles = ['Initial porosity', 'Final porosity', 'Final speed (m/s)', f'{mineral}: final minus initial fraction']
    fig, axes = plt.subplots(2, 2, figsize=(9, 5.6), constrained_layout=True)
    for i, (ax, values, title) in enumerate(zip(axes.flat, images, titles)):
        kw = {'vmin': 0., 'vmax': 1.} if i < 2 else {}
        if i == 3:
            bound = max(np.max(np.abs(change)), 1e-12)
            kw = {'vmin': -bound, 'vmax': bound}
        image = ax.imshow(values.reshape(ny, nx), origin='lower', extent=extent,
                          interpolation='nearest', cmap='RdBu_r' if i == 3 else 'viridis', **kw)
        ax.set(title=title, xlabel='x (µm)', ylabel='y (µm)')
        fig.colorbar(image, ax=ax, shrink=.8)
    fig.suptitle(f'{float(stamps[0]):g}–{float(stamps[-1]):g} s · stored cell averages')
    target = run_directory/'figures/overview.png'
    target.parent.mkdir(exist_ok=True)
    fig.savefig(target, dpi=140)
    plt.close(fig)
    return target
