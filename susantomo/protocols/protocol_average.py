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
from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA
import pyworkflow.utils as pwutils
from pyworkflow.plugin import Domain
from pwem import emlib

from tomo.objects import SetOfCTFTomoSeries, SetOfTiltSeries, AverageSubTomogram
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging

from .. import Plugin
from ..convert import writeDynTable


class outputs(Enum):
    outputAverage = AverageSubTomogram


class ProtSusanAverage(ProtTomoSubtomogramAveraging):
    """ Average and reconstruct a 3D volume (subtomogram average).
        This protocol is used to check subtomograms positions import.
    """
    _label = 'average and reconstruct'
    _devStatus = BETA
    _possibleOutputs = outputs

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)
        self.stacks, self.tilts, self.ids = [], [], []

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0', label="Choose GPU IDs")

        form.addSection(label='Input')
        form.addParam('inputSetOfCoords3D', params.PointerParam,
                      pointerClass="SetOfCoordinates3D",
                      important=True,
                      label='Input 3D coordinates',
                      help="Set of 3D coordinates defining subtomograms positions.")
        form.addParam('inputTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfCTFTomoSeries, SetOfTiltSeries',
                      important=True,
                      label='CTF tomo series or tilt-series (aligned)',
                      help='Set of tilt-series that correspond to subtomograms. '
                           'The matching is done using tsId.')
        form.addParam('tomoSize', params.IntParam,
                      default=800, label='Tomogram thickness (px)',
                      help='Z height of a tomogram volume in '
                           'pixels. Required for tilt series stack.')
        form.addParam('boxSize', params.IntParam,
                      default=32, label='Output box size (voxels)',
                      help='Size of the reconstructed average volume in '
                           'voxels. Pixel size will be the same as input '
                           'tilt series.')

        form.addSection(label='Averaging')
        form.addParam('sym', params.StringParam,
                      default='c1',
                      label='Symmetry group',
                      help="Specify the particle symmetry.")
        form.addParam('ctfCorr', params.EnumParam,
                      choices=['none', 'phase_flip', 'wiener', 'wiener_ssnr'],
                      default=2,
                      label="CTF correction method")
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

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.averageStep)
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

        if self.ctfCorr.get():
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
            self.stacks.append(tsFn)
            self.tilts.append(tiltFn)
            self.ids.append(ts.getObjId())

    def averageStep(self):
        """ Run susan_reconstruct program. """
        tsSet = self._getInputTs()
        tomo_size = [tsSet.getDim()[0], tsSet.getDim()[1], self.tomoSize.get()]

        params = {
            'ts_nums': self.ids,
            'inputStacks': self.stacks,
            'inputAngles': self.tilts,
            'output_dir': self._getExtraPath(),
            'num_tilts': max([ts.getSize() for ts in tsSet.iterItems()]),
            'pix_size': tsSet.getSamplingRate(),
            'tomo_size': tomo_size,
            'box_size': self.boxSize.get(),
            'gpus': self.getGpuList(),
            'voltage': tsSet.getAcquisition().getVoltage(),
            'sph_aber': tsSet.getAcquisition().getSphericalAberration(),
            'amp_cont': tsSet.getAcquisition().getAmplitudeContrast(),
            'thr_per_gpu': self.numberOfThreads.get(),
            'ctf_corr': self.getEnumText('ctfCorr'),
            'do_halfsets': self.doHalfSets.get(),
            'symmetry': self.sym.get(),
            'padding': self.getEnumText('padding')
        }

        jsonFn = self._getTmpPath("params.json")
        with open(jsonFn, "w") as fn:
            json.dump(params, fn, indent=4)

        self.runJob(Plugin.getProgram("average.py"), jsonFn,
                    env=Plugin.getEnviron())

    def createOutputStep(self):
        imgSet = self._getInputTs()
        volume = AverageSubTomogram()
        volumeFile = self._getExtraPath("average_class001.mrc")
        volume.setFileName(volumeFile)
        volume.setSamplingRate(imgSet.getSamplingRate())
        if self.doHalfSets:
            volume.setHalfMaps([self._getExtraPath("average_class001_half1.mrc"),
                                self._getExtraPath("average_class001_half2.mrc")])

        self._defineOutputs(**{outputs.outputAverage.name: volume})
        self._defineSourceRelation(self._getInputTs(pointer=True), volume)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.ctfCorr.get() and isinstance(self.inputTiltSeries.get(),
                                             SetOfTiltSeries):
            errors.append("CTF correction requires that you provide "
                          "CTFTomoSeries as input")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getInputTs(self, pointer=False):
        if isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries):
            return self.inputTiltSeries.get().getSetOfTiltSeries(pointer=pointer)
        return self.inputTiltSeries.get() if not pointer else self.inputTiltSeries

    def getScaleFactor(self):
        samplingRateCoords = self.inputSetOfCoords3D.get().getSamplingRate()
        samplingRateTS = self._getInputTs().getSamplingRate()
        return samplingRateCoords / samplingRateTS
