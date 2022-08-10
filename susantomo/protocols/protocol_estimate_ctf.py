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
from pyworkflow.constants import NEW
import pyworkflow.utils as pwutils

from tomo.objects import CTFTomo, SetOfCTFTomoSeries
from tomo.protocols.protocol_ts_estimate_ctf import ProtTsEstimateCTF

from .. import Plugin
from ..convert import parseCtf


class outputs(Enum):
    outputCTFs = SetOfCTFTomoSeries


class ProtSusanEstimateCtf(ProtTsEstimateCTF):
    """ CTF estimation on a set of tilt series using SUSAN. """
    _label = 'ctf estimation'
    _devStatus = NEW
    _possibleOutputs = outputs

    # -------------------------- DEFINE param functions -----------------------
    def _defineProcessParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0', label="Choose GPU IDs")

        form.addParam('tomoSize', params.IntParam,
                      default=800, label='Tomogram thickness (voxels)',
                      help='Z height of the reconstructed volume in '
                           '*unbinned* voxels.')

        form.addParam('gridSampling', params.IntParam, default=180,
                      label="CTF grid spacing (px)",
                      help="Spacing between grid points, in pixels.")

        group = form.addGroup('Search limits')
        line = group.addLine('Resolution (A)',
                             help='The CTF model will be fit to regions '
                                  'of the amplitude spectrum corresponding '
                                  'to this range of resolution.')
        line.addParam('lowRes', params.FloatParam, default=30., label='Min')
        line.addParam('highRes', params.FloatParam, default=5., label='Max')

        line = group.addLine('Defocus search range (A)',
                             help='Select _minimum_ and _maximum_ values for '
                                  'defocus search range (in A). Underfocus'
                                  ' is represented by a positive number.')
        line.addParam('minDefocus', params.FloatParam, default=5000.,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=50000.,
                      label='Max')

        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      help='')

        form.addParam('windowSize', params.IntParam, default=512,
                      label='FFT box size (px)',
                      help='')

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- STEPS functions -----------------------------
    def processTiltSeriesStep(self, tsObjId):
        """Run susan_estimate_ctf program"""
        ts = self._getTiltSeries(tsObjId)
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        pwutils.makePath(tmpPrefix)
        pwutils.makePath(extraPrefix)

        ts.applyTransform(self.getFilePath(tsObjId, tmpPrefix, ".mrc"))
        ts.generateTltFile(self.getFilePath(tsObjId, tmpPrefix, ".tlt"))

        paramDict = self.getCtfParamsDict()
        tomo_size = [ts.getDim()[0], ts.getDim()[1], self.tomoSize.get()]

        params = {
            'ts_num': ts.getObjId(),
            'inputStack': self.getFilePath(tsObjId, tmpPrefix, ".mrc"),
            'inputAngles': self.getFilePath(tsObjId, tmpPrefix, ".tlt"),
            'output_dir': extraPrefix,
            'num_tilts': ts.getDim()[-1],
            'pix_size': ts.getSamplingRate(),
            'tomo_size': tomo_size,
            'sampling': self.gridSampling.get(),
            'binning': self.ctfDownFactor.get(),
            'gpus': self.getGpuList(),
            'min_res': paramDict['lowRes'],
            'max_res': paramDict['highRes'],
            'def_min': paramDict['minDefocus'],
            'def_max': paramDict['maxDefocus'],
            'patch_size': paramDict['windowSize'],
            'voltage': paramDict['voltage'],
            'sph_aber': paramDict['sphericalAberration'],
            'amp_cont': paramDict['ampContrast']
        }

        jsonFn = self.getFilePath(tsObjId, tmpPrefix, ".json")
        with open(jsonFn, "w") as fn:
            json.dump(params, fn, indent=4)

        self.runJob(Plugin.getProgram("estimate_ctf.py"), jsonFn,
                    env=Plugin.getEnviron())

        for i, ti in enumerate(self._tsDict.getTiList(tsId)):
            ti.setCTF(self.getCtf(tsObjId, ti))

        self._tsDict.setFinished(tsId)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _doInsertTiltImageSteps(self):
        """ Default True, but if return False, the steps for each
        TiltImage will not be inserted. """
        return False

    def getFilePath(self, tsObjId, prefix, ext=None):
        ts = self._getTiltSeries(tsObjId)
        return os.path.join(prefix,
                            ts.getFirstItem().parseFileName(extension=ext))

    def getPsdName(self, tsObjId, ti):
        return os.path.join(self._getExtraPath(tsObjId),
                            f'ctf_grid/Tomo{tsObjId:03g}/???')

    def getCtfName(self, tsObjId):
        return os.path.join(self._getExtraPath(tsObjId),
                            f'ctf_grid/Tomo{tsObjId:03g}/defocus.txt')

    def getCtf(self, tsObjId, ti):
        """ Parse the CTF object estimated for this Tilt-Image. """
        ctfModel = parseCtf(self.getCtfName(tsObjId), ti,
                            psdFile=self.getPsdName(tsObjId, ti))
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctfModel)

        return ctfTomo
