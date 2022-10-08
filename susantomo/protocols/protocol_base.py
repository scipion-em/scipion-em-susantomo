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
import re
from glob import glob

import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA
import pyworkflow.utils as pwutils
from pyworkflow.plugin import Domain
from pwem import emlib
from pwem.protocols import EMProtocol

from tomo.objects import SetOfCTFTomoSeries, SetOfTiltSeries

from ..convert import writeDynTable


class ProtSusanBase(EMProtocol):
    """ Base protocol for SUSAN. """
    _label = None
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _initialize(self):
        """ This function is meant to be called after the
        working dir for the protocol have been set
        (maybe after recovery from mapper)
        """
        (self.stacks, self.tilts, self.ids,
         self.refs, self.masks) = [], [], [], [], []
        self._createFilenameTemplates()
        self._createIterTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        iterDir = self._getExtraPath("mra/ite_%(iter)04d")
        inputDir = self._getExtraPath("input")

        self._updateFilenamesDict({
            'input_ptcls': inputDir + '/input_particles.ptclsraw',
            'input_refs': inputDir + '/input_refs.refstxt',
            'input_tomos': inputDir + '/input_tomos.tomostxt',
            'ptcls': iterDir + '/particles.ptclsraw',
            'refs': iterDir + '/reference.refstxt',
            'outvol': iterDir + '/map_class%(ref3d)03d.mrc',
            'outvol_half1': iterDir + '/map_class%(ref3d)03d_half1.mrc',
            'outvol_half2': iterDir + '/map_class%(ref3d)03d_half2.mrc',
            'outavg': self._getExtraPath('average_class%(ref3d)03d.mrc'),
            'outavg_half1': self._getExtraPath('average_class%(ref3d)03d_half1.mrc'),
            'outavg_half2': self._getExtraPath('average_class%(ref3d)03d_half2.mrc'),
            'info': self._getExtraPath('mra/info.pkl')
        })

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('ptcls', iter=0).replace('0000', '????')
        self._iterRegex = re.compile('ite_(\d{4})')

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0', label="Choose GPU IDs")

        form.addSection(label='Input')

        self._defineContinueParams(form)

        form.addParam('inputSetOfSubTomograms', params.PointerParam,
                      pointerClass="SetOfSubTomograms",
                      condition="not doContinue",
                      important=True,
                      label='Input subtomograms',
                      help="Set of subtomograms that will be used to fetch "
                           "coordinates information. *Input alignment is "
                           "ignored at the moment.*")
        form.addParam('inputTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfCTFTomoSeries, SetOfTiltSeries',
                      important=True,
                      label='CTF tomo series or tilt-series (aligned)',
                      help='Set of tilt-series that correspond to the '
                           'input above. The matching is done using tsId.')
        form.addParam('tomoSize', params.IntParam,
                      default=110, label='Tomogram thickness (px)',
                      help='Z height of a tomogram volume in '
                           'pixels. Required for tilt series stack. '
                           'Pixel size will be the same as input '
                           'tilt series.')
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
        form.addParam('doHalfSets', params.BooleanParam, default=True,
                      label="Reconstruct half-sets?")

        if self.getClassName() == "ProtSusanMRA":
            form.addParam('threshold', params.FloatParam,
                          label='CC threshold', default=0.9,
                          help='Threshold value for CC (cross-correlation). '
                               'Determines the fraction of particles to be used '
                               'for reconstruction.')

        form.addParallelSection(threads=1, mpi=1)

    def _defineContinueParams(self, form):
        """ Can be re-defined in subclasses. """
        form.addParam('doContinue', params.BooleanParam, default=False,
                      label="Continue previous run?",
                      help="Particles metadata from the previous "
                           "protocol run will be used. You still have "
                           "to provide tilt-series, references and masks. "
                           "They can have a different binning compared to "
                           "the previous run.")
        form.addParam('inputSubstacks', params.PointerParam,
                      pointerClass="TomoSubStacks",
                      important=True,
                      allowsNull=True,
                      condition="doContinue",
                      label="Substacks from the previous run")

    def _defineProcessParams(self, form):
        """ Should be implemented in subclasses. """
        pass

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        pwutils.makePath(self._getExtraPath("input"))
        if self.isContinue():
            self._insertFunctionStep(self.continueStep)
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.runSusanStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def continueStep(self):
        """ Copy the ptclsraw from the previous run. """
        prevPtcls = self.inputSubstacks.get().getFileName()
        self.info(f"Copying particles from the previous run: {prevPtcls}")
        pwutils.copyFile(prevPtcls, self._getExtraPath("input/input_particles.ptclsraw"))

        if getattr(self, "reuseRefs", False):
            # copy refstxt file
            prevRefs = prevPtcls.replace("particles.ptclsraw", "reference.refstxt")
            self.info(f"Copying references from the previous run: {prevRefs}")
            pwutils.copyFile(prevRefs, self._getExtraPath("input/input_refs.refstxt"))

    def convertInputStep(self):
        """ Prepare input files. """
        tsSet = self._getInputTs()

        if not self.isContinue():
            inputSubTomos = self.inputSetOfSubTomograms.get()
            scaleCoords = self.getScaleCoords()
            if abs(scaleCoords - 1.0) > 0.00001:
                self.info(f"Scaling coordinates by a factor of {scaleCoords:0.2f}")

            tsSet = self._getInputTs()
            tsIds_from_ts = set(item.getTsId() for item in tsSet)
            tsIds_from_subtomos = set(inputSubTomos.getTomograms().keys())
            if not tsIds_from_ts.issubset(tsIds_from_subtomos):
                self.warning("Found subtomos with tsId that did not match "
                             "provided tilt-series: "
                             f"{set.difference(tsIds_from_subtomos, tsIds_from_ts)}")

            angleMax = tsSet.getAcquisition().getAngleMax() or 0
            angleMin = tsSet.getAcquisition().getAngleMin() or 0

            fnTable = self._getTmpPath("input_particles.tbl")
            with open(fnTable, 'w') as fn:
                writeDynTable(fn, inputSubTomos, angleMin, angleMax,
                              scaleCoords=scaleCoords,
                              scaleShifts=self.getScaleShifts())

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

        self.convertInputRefs()

    def runSusanStep(self):
        """ Should be implemented in sub-classes. """
        raise NotImplementedError

    def createOutputStep(self):
        """ Should be implemented in sub-classes. """
        raise NotImplementedError

    # --------------------------- INFO functions ------------------------------
    def _validateBase(self):
        """ Should be re-defined in subclasses. """
        errors = []
        if self.doContinue and not self.inputSubstacks.hasValue():
            errors.append("Please input the tomo substacks")

        if self.doCtf() and isinstance(self.inputTiltSeries.get(),
                                       SetOfTiltSeries):
            errors.append("CTF correction requires that you provide "
                          "CTFTomoSeries as input")

        return errors

    def _summary(self):
        summary = []

        output = self.getOutputName(self.getNumRefs())
        if hasattr(self, output):
            msg = "Computed average subtomogram(s) using the tilt-series "
            msg += "substacks" if self.doContinue else "stacks"
            summary.append(msg)
        else:
            summary.append("Output is not ready")

        return summary

    # --------------------------- UTILS functions -----------------------------
    def convertInputRefs(self):
        """ Should be defined in subclasses. """
        pass

    def _getInputTs(self, pointer=False):
        if isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries):
            return self.inputTiltSeries.get().getSetOfTiltSeries(pointer=pointer)
        return self.inputTiltSeries.get() if not pointer else self.inputTiltSeries

    def getScaleCoords(self):
        samplingRateCoords = self.inputSetOfSubTomograms.get().getCoordinates3D().getSamplingRate()
        samplingRateTS = self._getInputTs().getSamplingRate()
        return samplingRateCoords / samplingRateTS

    def getScaleShifts(self):
        samplingRateSubtomos = self.inputSetOfSubTomograms.get().getSamplingRate()
        samplingRateTS = self._getInputTs().getSamplingRate()
        return samplingRateSubtomos / samplingRateTS

    def hasCtf(self):
        return isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries)

    def doCtf(self):
        """ Should be re-defined in subclasses. """
        return False

    def isContinue(self):
        """ Should be re-defined in subclasses. """
        return self.doContinue

    def getNumRefs(self):
        return int(self.inputSubstacks.get().getNumRefs()) if self.doContinue else 1

    def getOutputName(self, nrefs):
        return f"outputAverage{'s' if nrefs > 1 else ''}"

    def _getIterNumber(self, index):
        """ Return the last iteration number. """
        result = None
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 4 digit iteration number
        return result

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0)
