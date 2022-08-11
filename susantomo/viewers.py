# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from matplotlib.figure import Figure

from pwem.emlib.image import ImageHandler
from tomo.viewers.viewers_data import CtfEstimationTomoViewer

from .protocols import ProtSusanEstimateCtf


def getPlotSubtitle(ctf):
    """ Create plot subtitle using CTF values. """
    ang = u"\u212B"
    deg = u"\u00b0"
    def1, def2, angle = ctf.getDefocus()
    phSh = ctf.getPhaseShift()
    score = ctf.getFitQuality()
    res = ctf.getResolution()

    title = "Def1: %d %s | Def2: %d %s | Angle: %0.1f%s | " % (
        def1, ang, def2, ang, angle, deg)

    if phSh is not None:
        title += "Phase shift: %0.2f %s | " % (phSh, deg)

    title += "Fit: %0.1f %s | Score: %0.3f" % (res, ang, score)

    return title


class CtfEstimationTomoViewerSusan(CtfEstimationTomoViewer):
    _targets = [ProtSusanEstimateCtf]

    def plot2D(self, ctfSet, ctfId):
        ctfModel = ctfSet[ctfId]
        psdFn = ctfModel.getPsdFile()
        img = ImageHandler().read(psdFn)
        fig = Figure(figsize=(7, 7), dpi=100)
        psdPlot = fig.add_subplot(111)
        psdPlot.get_xaxis().set_visible(False)
        psdPlot.get_yaxis().set_visible(False)
        psdPlot.set_title('%s # %d\n' % (ctfSet.getTsId(), ctfId) + getPlotSubtitle(ctfModel))
        psdPlot.imshow(img.getData(), cmap='gray')

        return fig
