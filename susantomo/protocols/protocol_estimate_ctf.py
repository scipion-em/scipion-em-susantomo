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
from pyworkflow.constants import BETA, SCIPION_DEBUG_NOCLEAN
import pyworkflow.utils as pwutils
from pwem.objects import CTFModel
from pwem import emlib

from tomo.objects import CTFTomo, SetOfCTFTomoSeries
from tomo.protocols.protocol_ts_estimate_ctf import ProtTsEstimateCTF

from .. import Plugin
from ..convert import parseImodCtf, readCtfModelStack


class outputs(Enum):
    outputCTFs = SetOfCTFTomoSeries


class ProtSusanEstimateCtf(ProtTsEstimateCTF):
    """ CTF estimation on a set of tilt series using SUSAN. """
    _label = 'ctf estimation'
    _devStatus = BETA
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
                      help='Set to 1 for no downsampling. '
                           'This downsampling is only used '
                           'for estimating the CTF and it does not affect any '
                           'further calculation. Ideally the estimation of the '
                           'CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) '
                           'and not occupying the whole power spectrum (since '
                           'this downsampling might entail aliasing).')

        form.addParam('windowSize', params.IntParam, default=512,
                      label='FFT box size (px)',
                      help='Boxsize in pixels to be used for FFT')

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- STEPS functions -----------------------------
    def processTiltSeriesStep(self, tsId):
        """Run susan_estimate_ctf program"""
        ts = self._getTiltSeries(tsId)
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)
        pwutils.makePath(tmpPrefix)
        pwutils.makePath(extraPrefix)
        tsInputFn = ts.getFirstItem().getFileName()
        tsFn = self.getFilePath(tsId, tmpPrefix, ".mrc")
        tiltFn = self.getFilePath(tsId, tmpPrefix, ".tlt")

        # has to be float32
        ih = emlib.image.ImageHandler()
        if ih.getDataType(tsInputFn) != emlib.DT_FLOAT:
            ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)
        elif pwutils.getExt(tsInputFn) in ['.mrc', '.st', '.mrcs']:
            pwutils.createAbsLink(os.path.abspath(tsInputFn), tsFn)
        else:
            ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)

        ts.generateTltFile(tiltFn)

        paramDict = self.getCtfParamsDict()
        tomo_size = [ts.getDim()[0], ts.getDim()[1], self.tomoSize.get()]
        ts_num = ts.getObjId()

        params = {
            'ts_nums': [ts_num],
            'inputStacks': [os.path.abspath(tsFn)],
            'inputAngles': [os.path.abspath(tiltFn)],
            'num_tilts': ts.getSize(),
            'pix_size': ts.getSamplingRate(),
            'tomo_size': tomo_size,
            'sampling': self.gridSampling.get(),
            'binning': int(math.log2(self.ctfDownFactor.get())),
            'has_ctf': self.hasCtf(),
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

        jsonFn = self.getFilePath(tsId, tmpPrefix, ".json")
        with open(jsonFn, "w") as fn:
            json.dump(params, fn, indent=4)

        try:
            self.runJob(Plugin.getProgram("estimate_ctf.py"),
                        os.path.abspath(jsonFn),
                        env=Plugin.getEnviron(),
                        cwd=extraPrefix,
                        numberOfThreads=1)

            outputLog = self.getOutputPath(tsId, ts_num) + "/defocus.txt"
            ctfResult = parseImodCtf(outputLog)
            psds = self.getOutputPath(tsId, ts_num) + "/ctf_fitting_result.mrc"
            ctf = CTFModel()

            for i, ti in enumerate(self._tsDict.getTiList(tsId)):
                ti.setCTF(self.getCtfTi(ctf, ctfResult, i, psds))

            if not pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN):
                pwutils.cleanPath(tsFn)

            self._tsDict.setFinished(tsId)
        except Exception as e:
            self.error(f"ERROR: SUSAN ctf estimation has failed for {tsFn}: {e}")

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.numberOfMpi.get() > 1:
            errors.append("MPI not supported, use threads (per GPU) instead.")

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

    def getOutputPath(self, tsObjId, tsNum):
        return os.path.join(self._getExtraPath(tsObjId),
                            f'ctf_grid/Tomo{tsNum:03d}')

    def getCtfTi(self, ctf, ctfArray, tiIndex, psdStack):
        """ Parse the CTF object estimated for this Tilt-Image. """
        readCtfModelStack(ctf, ctfArray, item=tiIndex)
        ctf.setPsdFile(f"{tiIndex + 1}@" + psdStack)
        ctfTomo = CTFTomo.ctfModelToCtfTomo(ctf)

        return ctfTomo

    def hasCtf(self):
        return isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries)
