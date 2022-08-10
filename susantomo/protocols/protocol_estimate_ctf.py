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
from tomo.objects import CTFTomoSeries, SetOfCTFTomoSeries

from imod.protocols.protocol_base import ProtImodBase
from .. import Plugin


class outputs(Enum):
    outputCTFs = SetOfCTFTomoSeries


class ProtSusanEstimateCtf(ProtImodBase):
    """ CTF estimation on a set of tilt series using SUSAN. """
    _label = 'ctf estimation'
    _devStatus = NEW
    _possibleOutputs = outputs

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputSet',
                      params.PointerParam,
                      pointerClass='SetOfTiltSeries, SetOfCTFTomoSeries',
                      label='Input set of tilt-series',
                      help='')

        form.addParam('sampling', params.IntParam, default=180,
                      label="CTF grid sampling (px)",
                      help="spacing between particles, in pixels")

        group = form.addGroup('Search limits')
        line = group.addLine('Resolution (A)', condition='not recalculate',
                             help='The CTF model will be fit to regions '
                                  'of the amplitude spectrum corresponding '
                                  'to this range of resolution.')
        line.addParam('lowRes', params.FloatParam, default=30., label='Min')
        line.addParam('highRes', params.FloatParam, default=5., label='Max')

        line = group.addLine('Defocus search range (A)',
                             condition='not recalculate',
                             help='Select _minimum_ and _maximum_ values for '
                                  'defocus search range (in A). Underfocus'
                                  ' is represented by a positive number.')
        line.addParam('minDefocus', params.FloatParam, default=5000.,
                      label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=50000.,
                      label='Max')

        form.addParam('ctfDownFactor', params.FloatParam, default=1.,
                      label='CTF Downsampling factor',
                      help='Set to 1 for no downsampling. Non-integer downsample '
                           'factors are possible. This downsampling is only used '
                           'for estimating the CTF and it does not affect any '
                           'further calculation. Ideally the estimation of the '
                           'CTF is optimal when the Thon rings are not too '
                           'concentrated at the origin (too small to be seen) '
                           'and not occupying the whole power spectrum (since '
                           'this downsampling might entail aliasing).')

        form.addParam('windowSize', params.IntParam, default=512,
                      label='FFT box size (px)', condition='not recalculate',
                      help='The dimensions (in pixels) of the amplitude '
                           'spectrum CTFfind will compute. Smaller box '
                           'sizes make the fitting process significantly '
                           'faster, but sometimes at the expense of '
                           'fitting accuracy. If you see warnings '
                           'regarding CTF aliasing, consider '
                           'increasing this parameter.')

        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0', label="Choose GPU IDs")

        form.addParallelSection(threads=1, mpi=0)

    # -------------------------- INSERT steps functions ---------------------
    def _insertAllSteps(self):
        # This assignment is needed to use methods from base class
        self.inputSetOfTiltSeries = self._getSetOfTiltSeries()

        for item in self.inputSet.get():
            self._insertFunctionStep(self.ctfEstimationStep, item.getObjId(),
                                     imodInterpolation=False)
            self._insertFunctionStep(self.createOutputStep, item.getObjId())
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def ctfEstimationStep(self, tsObjId):
        """Run susan_estimate_ctf program"""
        ts = self._getTiltSeries(tsObjId)
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)
        tmpPrefix = self._getTmpPath(tsId)

        params = {
            'ts_num': 1,
            'inputStack': self.getFilePath(tsObjId, tmpPrefix),
            'inputAngles': self.getFilePath(tsObjId, tmpPrefix, ".tlt"),
            'output_dir': extraPrefix,
            'num_tilts': ts.getAnglesCount(),
            'pix_size': ts.getSamplingRate(),
            'tomo_size': ts.getDim(),
            'sampling': self.sampling.get(),
            'binning': self.ctfDownFactor.get(),
            'gpus': '%(GPU)s'.split(","),
            'min_res': self.lowRes.get(),
            'max_res': self.highRes.get(),
            'def_min': self.minDefocus.get(),
            'def_max': self.maxDefocus.get(),
            'patch_size': self.windowSize.get()
        }

        jsonFn = self.getFilePath(tsObjId, extraPrefix, ".json")
        with open(jsonFn, "w") as fn:
            json.dump(params, fn)

        self.runJob(Plugin.getProgram("protocol_estimate_ctf.py", jsonFn),
                    env=Plugin.getEnviron())

    def createOutputStep(self, tsObjId):
        ts = self._getTiltSeries(tsObjId)
        tsId = ts.getTsId()
        extraPrefix = self._getExtraPath(tsId)

        self.outputSetName = outputs.outputCTFs.name

        # f'ctf_grid/Tomo{ts_id:03g}/defocus.txt'
        defocusFn = self.getFilePath(tsObjId, extraPrefix, ".txt")
        if os.path.exists(defocusFn):
            self.getOutputSetOfCTFTomoSeries(self.outputSetName)
            self.addCTFTomoSeriesToSetFromDefocusFile(ts, defocusFn)

    def closeOutputSetsStep(self):
        output = getattr(self, self.outputSetName, None)
        if output is not None:
            output.setStreamState(output.STREAM_CLOSED)
            output.write()
        self._store()

    # --------------------------- INFO functions ------------------------------
    
    def _summary(self):
        summary = []

        return summary
    
    def _validate(self):
        errors = []

        return errors
    
    # --------------------------- UTILS functions -----------------------------
    def _getSetOfTiltSeries(self, pointer=False):
        if isinstance(self.inputSet.get(), SetOfCTFTomoSeries):
            return self.inputSet.get().getSetOfTiltSeries(pointer=pointer)

        return self.inputSet.get() if not pointer else self.inputSet

    def _getTiltSeries(self, itemId):
        obj = None
        inputSetOfTiltseries = self._getSetOfTiltSeries()
        for item in inputSetOfTiltseries.iterItems(iterate=False):
            if item.getTsId() == itemId:
                obj = item
                if isinstance(obj, CTFTomoSeries):
                    obj = item.getTiltSeries()
                break

        if obj is None:
            raise ("Could not find tilt-series with tsId = %s" % itemId)

        return obj

    def getFilePath(self, tsObjId, prefix, ext=None):
        ts = self._getSetOfTiltSeries()[tsObjId]
        return os.path.join(prefix,
                            ts.getFirstItem().parseFileName(extension=ext))
