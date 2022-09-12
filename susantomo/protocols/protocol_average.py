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

import json
import os.path
from enum import Enum

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params

from tomo.objects import SetOfTiltSeries, AverageSubTomogram
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging

from .. import Plugin
from .protocol_base import ProtSusanBase


class outputs(Enum):
    outputAverage = AverageSubTomogram


class ProtSusanAverage(ProtSusanBase, ProtTomoSubtomogramAveraging):
    """ Average and reconstruct a 3D volume (subtomogram average).
        This protocol is used to check subtomograms positions import.
    """
    _label = 'average and reconstruct'
    _devStatus = BETA
    _possibleOutputs = outputs

    def __init__(self, **args):
        ProtSusanBase.__init__(self, **args)

    # --------------------------- STEPS functions -----------------------------
    def runSusanStep(self):
        """ Run susan_reconstruct program. """
        tsSet = self._getInputTs()
        tomo_size = [tsSet.getDim()[0], tsSet.getDim()[1], self.tomoSize.get()]

        params = {
            'continue': bool(self.doContinue),
            'ts_nums': self.ids,
            'inputStacks': self.stacks,
            'inputAngles': self.tilts,
            'num_tilts': max([ts.getSize() for ts in tsSet.iterItems()]),
            'pix_size': tsSet.getSamplingRate(),
            'tomo_size': tomo_size,
            'box_size': self.boxSize.get(),
            'gpus': self.getGpuList(),
            'voltage': tsSet.getAcquisition().getVoltage(),
            'sph_aber': tsSet.getAcquisition().getSphericalAberration(),
            'amp_cont': tsSet.getAcquisition().getAmplitudeContrast(),
            'thr_per_gpu': self.numberOfThreads.get(),
            'has_ctf': self.hasCtf(),
            'ctf_corr_avg': self.getEnumText('ctfCorrAvg'),
            'do_halfsets': self.doHalfSets.get(),
            'symmetry': self.sym.get(),
            'padding': self.getEnumText('padding')
        }

        jsonFn = self._getTmpPath("params.json")
        with open(jsonFn, "w") as fn:
            json.dump(params, fn, indent=4)

        self.runJob(Plugin.getProgram("average.py"),
                    os.path.abspath(jsonFn),
                    env=Plugin.getEnviron(),
                    cwd=self._getExtraPath())

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
        errors = self._validateBase()

        return errors

    def _summary(self):
        summary = []

        if hasattr(self, outputs.outputAverage.name):
            summary.append("Computed a 3D average using input subtomograms.")
        else:
            summary.append("Output is not ready")

    # --------------------------- UTILS functions -----------------------------
    def doCtf(self):
        return self.ctfCorrAvg.get()
