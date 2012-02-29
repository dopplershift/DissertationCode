from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from .basic import Plot
from disser.datatypes import TypePlotInfo


class PPIPlot(Plot):

    def __init__(self, data, var, x='x', y='y', ax=None, rings=None,
        labels=None):
        if ax is None:
            self._ax = plt.gca()
        else:
            self._ax = ax

        if labels is None:
            labels = ['S', 'W']

        if rings is None:
            import quantities as pq
            rings = range(10, 70, 10) * pq.kilometer

        # Rely on TypePlotInfo being always returning something
        typeInfo = TypePlotInfo[var.type]
        cmap = typeInfo['cmap']
        norm = typeInfo['norm']

        self._mesh = ax.pcolormesh(data[x], data[y],
            data[var].magnitude, cmap=cmap, norm=norm)

        self._ring_patches = []
        for rng in rings:
            c = Circle(xy=(0, 0), radius=rng, fill=False)
            ax.add_patch(c)
            self._ring_patches.append(c)

        ax.set_xlabel('X (%s)' % data[x].dimensionality.latex)
        ax.set_ylabel('Y (%s)' % data[y].dimensionality.latex)
