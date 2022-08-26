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

import os
import json
import math
from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import NEW
import pyworkflow.utils as pwutils

from tomo.objects import SetOfSubTomograms
from tomo.utils import getNonInterpolatedTsFromRelations
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging

from .. import Plugin
from ..convert import writeDynTable


class outputs(Enum):
    outputSubtomograms = SetOfSubTomograms


class ProtSusanAverage(ProtTomoSubtomogramAveraging):
    """ Average and reconstruct a 3D volume (subtomogram). """
    _label = 'average and reconstruct 3D'
    _devStatus = NEW
    _possibleOutputs = outputs

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0', label="Choose GPU IDs")

        form.addSection(label='Input')
        form.addParam('inputSetOfSubTomograms', params.PointerParam,
                      pointerClass="SetOfSubTomograms",
                      important=True,
                      label='Input set of Subtomograms',
                      help="Set of subtomograms that will be used to fetch "
                           "alignment and coordinates information.")
        form.addParam('inputSetOfTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Input set of Tilt-Series',
                      help='Set of tilt-series that correspond to subtomograms. '
                           'The matching is checked using tsId.')

        form.addParam('sym', params.StringParam,
                      default='c1',
                      label='Symmetry group',
                      help="Specify the particle symmetry.")
        form.addParam('ctfCorr', params.EnumParam,
                      choices=['none', 'phase_flip', 'wiener', 'wiener_ssnr'],
                      default=2,
                      label="CTF correction method")
        form.addParam('norm', params.EnumParam,
                      choices=['none', 'zero_mean', 'zero_mean_one_std',
                               'zero_mean_proj_weight'],
                      default=2,
                      label="Normalization type")
        form.addParam('padding', params.EnumParam,
                      choices=['zero', 'noise'], default=1,
                      label="Padding type")
        form.addParam('doHalfSets', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Reconstruct half-sets?")

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('averageStep')
        self._insertFunctionStep('createOutputStep')
        self._insertFunctionStep('closeSetsStep')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """Run susan_reconstruct program"""

        # create ptcls tbl
        inputSubTomos = self.inputSetOfSubTomograms.get()
        fnTable = self._getExtraPath("input_particles.tbl")
        with open(fnTable, 'w') as fn:
            writeDynTable(fn, inputSubTomos)

        # create tomos file
        tsSet = self.inputSetOfTiltSeries.get()
        tsIds_from_ts = set(item.getTsId() for item in tsSet)
        tsIds_from_subtomos = set(inputSubTomos.getTomograms().keys())
        if tsIds_from_subtomos != tsIds_from_ts:
            self.warning("Found tsId that did not match between "
                         "provided tilt-series and subtomograms!")

        for ts in tsSet:
            tsId = ts.getTsId()
            ts.applyTransform(self._getTmpPath(tsId + ".mrc"))
            ts.generateTltFile(self._getTmpPath(tsId + ".tlt"))




        tsId = ts.getTsId()
        ts.generateTltFile(self.getFilePath(tsObjId, tmpPrefix, ".tlt"))

        paramDict = self.getCtfParamsDict()
        tomo_size = [ts.getDim()[0], ts.getDim()[1], self.tomoSize.get()]
        ts_num = ts.getObjId()

        params = {
            'ts_num': ts_num,
            'inputStack': self.getFilePath(tsObjId, tmpPrefix, ".mrc"),
            'inputAngles': self.getFilePath(tsObjId, tmpPrefix, ".tlt"),
            'output_dir': extraPrefix,
            'num_tilts': ts.getDim()[-1],
            'pix_size': ts.getSamplingRate(),
            'tomo_size': tomo_size,
            'sampling': self.gridSampling.get(),
            'binning': int(math.log2(self.ctfDownFactor.get())),
            'gpus': self.getGpuList(),
            'min_res': paramDict['lowRes'],
            'max_res': paramDict['highRes'],
            'def_min': paramDict['minDefocus'],
            'def_max': paramDict['maxDefocus'],
            'patch_size': paramDict['windowSize'],
            'voltage': paramDict['voltage'],
            'sph_aber': paramDict['sphericalAberration'],
            'amp_cont': paramDict['ampContrast'],
            'thr_per_gpu': self.numberOfThreads.get()
        }

        jsonFn = self.getFilePath(tsObjId, tmpPrefix, ".json")
        with open(jsonFn, "w") as fn:
            json.dump(params, fn, indent=4)

        self.runJob(Plugin.getProgram("estimate_ctf.py"), jsonFn,
                    env=Plugin.getEnviron())

        ctfs, psds = self.getCtf(tsObjId, ts_num)
        for i, ti in enumerate(self._tsDict.getTiList(tsId)):
            ti.setCTF(self.getCtfTi(ctfs, i, psds))

        self._tsDict.setFinished(tsId)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    # --------------------------- UTILS functions -----------------------------
    def getFilePath(self, tsObjId, prefix, ext=None):
        ts = self._getTiltSeries(tsObjId)
        return os.path.join(prefix,
                            ts.getFirstItem().parseFileName(extension=ext))

    def getOutputPath(self, tsObjId, tsNum):
        return os.path.join(self._getExtraPath(tsObjId),
                            f'ctf_grid/Tomo{tsNum:03d}')
