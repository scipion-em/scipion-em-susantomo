# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es) [2]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# * [2] BCU, Centro Nacional de Biotecnologia, CSIC
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

import os

import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA
import pyworkflow.utils as pwutils
from pyworkflow.plugin import Domain
from pwem import emlib
from pwem.protocols import EMProtocol

from tomo.objects import SetOfCTFTomoSeries

from ..convert import writeDynTable


class ProtSusanBase(EMProtocol):
    """ Base protocol for SUSAN. """
    _label = None
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        (self.stacks, self.tilts, self.ids,
         self.refs, self.masks) = [], [], [], [], []

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0', label="Choose GPU IDs")

        form.addSection(label='Input')

        self._defineContinueParams(form)

        form.addParam('inputSetOfCoords3D', params.PointerParam,
                      pointerClass="SetOfCoordinates3D",
                      condition="not doContinue",
                      important=True,
                      label='Input 3D coordinates',
                      help="Set of 3D coordinates defining "
                           "subtomograms positions.")
        form.addParam('inputTiltSeries',
                      params.PointerParam,
                      condition="not doContinue",
                      pointerClass='SetOfCTFTomoSeries, SetOfTiltSeries',
                      important=True,
                      label='CTF tomo series or tilt-series (aligned)',
                      help='Set of tilt-series that correspond to the '
                           'input above. The matching is done using tsId.')
        form.addParam('tomoSize', params.IntParam,
                      default=800, label='Tomogram thickness (px)',
                      help='Z height of a tomogram volume in '
                           'pixels. Required for tilt series stack.')
        form.addParam('boxSize', params.IntParam,
                      default=32, label='Output box size (voxels)',
                      help='Size of the reconstructed average volume in '
                           'voxels. Pixel size will be the same as input '
                           'tilt series.')
        form.addParam('sym', params.StringParam,
                      default='c1',
                      label='Symmetry group',
                      help="Specify the particle symmetry.")

        self._defineProcessParams(form)

        form.addSection(label='Averaging')
        form.addParam('ctfCorrAvg', params.EnumParam,
                      choices=['none', 'phase_flip', 'wiener', 'wiener_ssnr'],
                      default=2,
                      label="CTF correction method (averager)")
        form.addParam('norm', params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['none', 'zero_mean', 'zero_mean_one_std',
                               'zero_mean_proj_weight'],
                      default=2,
                      label="Normalization type")
        form.addParam('padding', params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['zero', 'noise'], default=1,
                      label="Padding type")
        form.addParam('doHalfSets', params.BooleanParam, default=False,
                      label="Reconstruct half-sets?")

        form.addParallelSection(threads=1, mpi=1)

    def _defineContinueParams(self, form):
        """ Can be re-defined in subclasses. """
        form.addParam('doContinue', params.BooleanParam, default=False,
                      label="Continue previous run?")
        form.addParam('previousRun', params.PointerParam,
                      pointerClass=self.getClassName(),
                      important=True,
                      allowsNull=True,
                      condition="doContinue",
                      label="Select previous run")

    def _defineProcessParams(self, form):
        """ Should be implemented in subclasses. """
        pass

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        if self.isContinue():
            self._insertFunctionStep(self.continueStep)
        else:
            self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.runSusanStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Prepare input files. """
        inputCoords = self.inputSetOfCoords3D.get()
        fnTable = self._getTmpPath("input_particles.tbl")
        factor = self.getScaleFactor()
        if abs(factor - 1.0 > 0.00001):
            self.info(f"Scaling coordinates by a factor of {factor:0.2f}")

        tsSet = self._getInputTs()
        tsIds_from_ts = set(item.getTsId() for item in tsSet)
        tsIds_from_coords = set(inputCoords.getPrecedentsInvolved().keys())
        if not tsIds_from_ts.issubset(tsIds_from_coords):
            self.warning("Found coords with tsId that did not match "
                         "provided tilt-series: "
                         f"{set.difference(tsIds_from_coords, tsIds_from_ts)}")

        angleMax = tsSet.getAcquisition().getAngleMax() or 0
        angleMin = tsSet.getAcquisition().getAngleMin() or 0

        with open(fnTable, 'w') as fn:
            writeDynTable(fn, inputCoords, angleMin, angleMax, scaleFactor=factor)

        if self.hasCtf():
            # generate defocus files
            setOfCtfTomoSeries = self.inputTiltSeries.get()
            imodUtils = Domain.importFromPlugin('imod.utils')
            for ctf in setOfCtfTomoSeries:
                tsId = ctf.getTsId()
                defocusFilePath = self._getTmpPath(tsId + ".defocus")
                imodUtils.generateDefocusIMODFileFromObject(ctf, defocusFilePath)

        for ts in tsSet:
            tsId = ts.getTsId()
            tsInputFn = ts.getFirstItem().getFileName()
            tsFn = self._getTmpPath(tsId + ".mrc")
            tiltFn = self._getTmpPath(tsId + ".tlt")

            # has to be float32
            ih = emlib.image.ImageHandler()
            if ih.getDataType(tsInputFn) != emlib.DT_FLOAT:
                ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)
            elif pwutils.getExt(tsInputFn) in ['.mrc', '.st', '.mrcs', '.ali']:
                pwutils.createAbsLink(os.path.abspath(tsInputFn), tsFn)
            else:
                ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)

            ts.generateTltFile(tiltFn)
            self.stacks.append(os.path.abspath(tsFn))
            self.tilts.append(os.path.abspath(tiltFn))
            self.ids.append(ts.getObjId())

        pwutils.makePath(self._getExtraPath("input"))
        self.convertInputRefs()

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        """ Should be re-defined in subclasses. """
        return []

    def _summary(self):
        """ Should be re-defined in subclasses. """
        return []

    # --------------------------- UTILS functions -----------------------------
    def convertInputRefs(self):
        """ Should be defined in subclasses. """
        pass

    def _getInputTs(self, pointer=False):
        if isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries):
            return self.inputTiltSeries.get().getSetOfTiltSeries(pointer=pointer)
        return self.inputTiltSeries.get() if not pointer else self.inputTiltSeries

    def getScaleFactor(self):
        samplingRateCoords = self.inputSetOfCoords3D.get().getSamplingRate()
        samplingRateTS = self._getInputTs().getSamplingRate()
        return samplingRateCoords / samplingRateTS

    def hasCtf(self):
        """ Should be re-defined in subclasses. """
        return isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries)

    def isContinue(self):
        """ Should be re-defined in subclasses. """
        return False