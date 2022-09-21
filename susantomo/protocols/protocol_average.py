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

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params

from tomo.objects import AverageSubTomogram, SetOfAverageSubTomograms
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging

from .. import Plugin
from .protocol_base import ProtSusanBase


class ProtSusanAverage(ProtSusanBase, ProtTomoSubtomogramAveraging):
    """ Average and reconstruct a 3D volume (subtomogram average). """
    _label = 'average and reconstruct'
    _devStatus = BETA
    _possibleOutputs = {'outputAverage': AverageSubTomogram,
                        'outputAverages': SetOfAverageSubTomograms}

    def __init__(self, **args):
        ProtSusanBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineContinueParams(self, form):
        form.addParam('doContinue', params.BooleanParam, default=False,
                      label="Use MRA output?",
                      help="Particles metadata from the previous MRA "
                           "protocol run will be used. You still have "
                           "to provide tilt-series. "
                           "They can have a different binning compared to "
                           "the previous run.")
        form.addParam('inputSubstacks', params.PointerParam,
                      pointerClass="TomoSubStacks",
                      important=True,
                      allowsNull=True,
                      condition="doContinue",
                      label="Substacks from the previous run")

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
            'do_halfsets': bool(self.doHalfSets.get()),
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
        pixSize = self._getInputTs().getSamplingRate()
        nrefs = self.getNumRefs()
        volume = AverageSubTomogram()

        def _createVolume(ind):
            volumeFile = self._getFileName("outavg", ref3d=ind)
            volume.setObjId(None)
            volume.setFileName(volumeFile)
            volume.setSamplingRate(pixSize)
            if self.doHalfSets:
                volume.setHalfMaps([self._getFileName("outavg_half1", ref3d=ind),
                                    self._getFileName("outavg_half2", ref3d=ind)])
            return volume

        if nrefs > 1:
            volumes = self._createSetOfAverageSubTomograms()
            volumes.setSamplingRate(pixSize)
            for i in range(nrefs):
                volume = _createVolume(i + 1)
                volume.setClassId(i + 1)
                volumes.append(volume)
        else:
            volumes = _createVolume(1)

        self._defineOutputs(**{self.getOutputName(nrefs): volumes})
        self._defineSourceRelation(self._getInputTs(pointer=True), volumes)

        if self.doContinue:
            self._defineSourceRelation(self.inputSubstacks, volumes)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = self._validateBase()

        return errors

    # --------------------------- UTILS functions -----------------------------
    def doCtf(self):
        return self.ctfCorrAvg.get()
