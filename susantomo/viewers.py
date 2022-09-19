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
import pickle
import os
from matplotlib.figure import Figure

import pyworkflow.protocol.params as params
from pyworkflow.viewer import DESKTOP_TKINTER
from pwem.emlib.image import ImageHandler
from pwem.objects import FSC
from pwem.viewers import EmPlotter
from pwem.viewers import ObjectView, ChimeraView, EmProtocolViewer, FscViewer
from tomo.viewers.viewers_data import CtfEstimationTomoViewer

from .constants import *
from .protocols import ProtSusanEstimateCtf, ProtSusanMRA


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
        index, psdFn = ctfModel.getPsdFile().split("@")
        if not os.path.exists(psdFn):
            return None
        img = ImageHandler().read((int(index), psdFn))
        fig = Figure(figsize=(7, 7), dpi=100)
        psdPlot = fig.add_subplot(111)
        psdPlot.get_xaxis().set_visible(False)
        psdPlot.get_yaxis().set_visible(False)
        psdPlot.set_title('%s # %d\n' % (ctfSet.getTsId(), ctfId) + getPlotSubtitle(ctfModel))
        psdPlot.imshow(img.getData(), cmap='gray')

        return fig


class MRAViewer(EmProtocolViewer):
    """ Visualization of SUSAN MRA results. """

    _environments = [DESKTOP_TKINTER]
    _targets = [ProtSusanMRA]
    _label = 'analyze results'

    def __init__(self, **kwargs):
        EmProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('viewIter', params.EnumParam,
                      choices=['last', 'selection'], default=ITER_LAST,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize",
                      help="""
        *last*: only the last iteration will be visualized.
        *selection*: you may specify a range of iterations.
        Examples:
        "1,5-8,10" -> [1,5,6,7,8,10]
        "2,6,9-11" -> [2,6,9,10,11]
        "2 5, 6-8" -> [2,5,6,7,8]""")
        form.addParam('iterSelection', params.NumericRangeParam,
                      condition='viewIter==%d' % ITER_SELECTION,
                      label="Iterations list",
                      help="Write the iteration list to visualize.")

        group = form.addGroup('Volumes')
        group.addParam('showClasses3D', params.EnumParam,
                       default=CLS_ALL,
                       choices=['all', 'selection'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='3D Class to visualize')
        group.addParam('class3DSelection', params.NumericRangeParam,
                       default='1',
                       condition=f'showClasses3D=={CLS_SELECTION}',
                       label='Classes list')
        group.addParam('showHalves', params.EnumParam, default=3,
                       choices=['half1', 'half2', 'both', 'final'],
                       label='Volume to visualize',
                       help='Select which half do you want to visualize.')
        group.addParam('displayVol', params.EnumParam,
                       choices=['slices', 'chimera'],
                       default=VOLUME_SLICES,
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')

        group = form.addGroup('Resolution')
        group.addParam('figure', params.EnumParam, default=0,
                       choices=['new', 'active'],
                       label='Figure',
                       display=params.EnumParam.DISPLAY_HLIST)
        group.addParam('resolutionPlotsFSC', params.LabelParam,
                       default=True,
                       label='Display resolution plots (FSC)')
        group.addParam('resolutionThresholdFSC', params.FloatParam,
                       default=0.143,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Threshold in resolution plots')

        form.addParam('showCC', params.LabelParam,
                      label="Show CC histogram for the last iteration")
        form.addParam('nBins', params.IntParam,
                      default=100,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Number of bins for histogram plot')

    def _getVisualizeDict(self):
        self._load()
        visDict = {
            'displayVol': self._showVolumes,
            'resolutionPlotsFSC': self._showFSC,
            'showCC': self._showCC,
            }

        # If the is some error during the load,
        # just show that instead of any viewer
        if self._errors:
            for k in visDict.keys():
                visDict[k] = self._showErrors

        return visDict

    def _showErrors(self, param=None):
        views = []
        self.errorList(self._errors, views)
        return views

    # --------- Show volumes --------------------------------------------------
    def _showVolumes(self, paramName=None):
        try:
            if self.displayVol == VOLUME_CHIMERA:
                return self._showVolumesChimera()
            elif self.displayVol == VOLUME_SLICES:
                return self._createVolumesSqlite()
        except Exception as e:
            self.showError(str(e))

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        prot = self.protocol
        volumes = self._getVolumeNames()
        extra = prot._getExtraPath()
        cmdFile = prot._getExtraPath('chimera_volumes.cxc')

        with open(cmdFile, 'w+') as f:
            for vol in volumes:
                localVol = os.path.relpath(vol, extra)
                if os.path.exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')

        view = ChimeraView(cmdFile)

        return [view]

    def _createVolumesSqlite(self):
        """ Write a sqlite with all volumes selected for visualization. """
        prot = self.protocol  # shortcut
        path = prot._getExtraPath('viewer_volumes.sqlite')
        samplingRate = prot._getInputTs().getSamplingRate()
        files = self._getVolumeNames()
        self.createVolumesSqlite(files, path, samplingRate)

        return [ObjectView(self._project, prot.strId(), path)]

    # -------------------------- show FSC -------------------------------------
    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()
        fscViewer = FscViewer(project=self.protocol.getProject(),
                              threshold=threshold,
                              protocol=self.protocol,
                              figure=self._getFigure(),
                              addButton=True)
        fscSet = self.protocol._createSetOfFSCs()
        pixSize = self.protocol._getInputTs().getSamplingRate()
        with open(self.protocol._getFileName("info"), "rb") as fn:
            frc = pickle.load(fn)[0]
            for it in self._iterations:
                fsc = self._plotFSC(frc, it, pix=pixSize, ref3d=1)
                fscSet.append(fsc)
        fscViewer.visualize(fscSet)
        return [fscViewer]

    def _plotFSC(self, frc, iter, pix=1.0, ref3d=1):
        frc = frc[ref3d][iter-1].tolist()  # iters are 0-indexed
        nqfp = len(frc)  # half-box size in px
        resolution_inv = [fpix/(2*pix*(nqfp-1)) for fpix in range(1, nqfp+1)]
        fsc = FSC(objLabel=f"iter {iter}")
        fsc.setData(resolution_inv, frc)

        return fsc

    # ----------------------------- Show CC -----------------------------------
    def _showCC(self, param=None):
        fn = self.protocol._getFileName('info')
        result = []
        with open(fn, "rb") as fn:
            cc = pickle.load(fn)[1]
            for ref3d in self._refsList:
                lastIter = len(cc[ref3d])
                result.append(cc[ref3d][lastIter-1].tolist())

        numberOfBins = self.nBins.get()
        plotter = EmPlotter()
        plotter.createSubPlot(f"Cross-correlation (iter {lastIter})",
                              "CC", "Number of particles")
        for r in result:
            plotter.plotHist(r, nbins=numberOfBins)
        plotter.show()

        return [plotter]

    # ------------------------------ Utils funcs ------------------------------
    def _getFigure(self):
        return None if self.figure == 0 else 'active'

    def _getVolumeNames(self):
        vols = []
        halves = self.showHalves.get()
        if halves in [0, 1]:
            keysFn = ["outvol_" + self.getEnumText("showHalves")]
        elif halves == 2:
            keysFn = ["outvol_half1", "outvol_half2"]
        else:
            keysFn = ["outavg"]

        for it in self._iterations:
            for ref3d in self._refsList:
                for key in keysFn:
                    volFn = self.protocol._getFileName(key,
                                                       iter=it, ref3d=ref3d)
                    if os.path.exists(volFn):
                        vols.append(volFn)
                    else:
                        raise Exception(f"Volume {volFn} does not exists.\n"
                                        "Please select a valid iteration "
                                        "number.")
        return vols

    def _load(self):
        """ Load selected iterations and classes for visualization. """
        self._refsList = [1]
        self._errors = []
        prot = self.protocol

        if self.showClasses3D == CLS_ALL:
            self._refsList = range(1, prot.getNumRefs() + 1)
        else:
            self._refsList = self._getRange(self.class3DSelection, 'classes 3d')

        prot._initialize()  # Load filename templates

        if self.viewIter.get() == ITER_LAST or self.showHalves == 3:
            self._iterations = [prot._lastIter()]
        else:
            self._iterations = self._getRange(self.iterSelection, 'iterations')

    def _getRange(self, var, label):
        """ Check if the range is not empty.
        :param var: The variable to retrieve the value
        :param label: the label used for the message string
        :return: the list with the range of values, empty
        """
        value = var.get()
        if value is None or not value.strip():
            self._errors.append(f'Provide {label} selection.')
            result = []
        else:
            result = self._getListFromRangeString(value)

        return result
