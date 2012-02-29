import numpy as np
from functools import partial
import quantities as pq
pq.markup.format_units_latex = partial(pq.markup.format_units_latex,
    font='mathsf')

rings = np.arange(10, 60, 10) * pq.km

from .ctables import get_cmap
from disser import datatypes
import matplotlib.pyplot as plt
datatypes.TypePlotInfo.set_defaults(cmap=get_cmap('Carbone42'))
datatypes.TypePlotInfo[datatypes.Reflectivity].update(
    norm=plt.Normalize(-20, 70))
datatypes.TypePlotInfo[datatypes.DopplerVelocity].update(
    norm=plt.Normalize(-30, 30))
datatypes.TypePlotInfo[datatypes.SpectrumWidth].update(
    norm=plt.Normalize(0, 15))
