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
import json

import pyworkflow.protocol.params as params
from pyworkflow.utils import getListFromValues
from pyworkflow.constants import BETA

from .. import Plugin
from ..objects import TomoSubStacks
from .protocol_base import ProtSusanBase


class ProtSusanSubset(ProtSusanBase):
    """ Create a subset of tomo substacks. """
    _label = 'multi-reference alignment'
    _devStatus = BETA
    _possibleOutputs = {'outputSubstacks': TomoSubStacks}

    def __init__(self, **args):
        ProtSusanBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addHidden('doContinue', params.BooleanParam, default=False)
        form.addParam('inputSubstacks', params.PointerParam,
                      pointerClass="TomoSubStacks",
                      important=True,
                      label="Substacks from the previous MRA run")

        form.addParam('selectRefs', params.BooleanParam, default=False,
                      label="Filter by MRA references")
        form.addParam('refsList', params.NumericRangeParam,
                      condition="selectRef",
                      label="References numbers list",
                      help="Provide a list of references (separated by space) "
                           "you would like to keep. "
                           "The first reference number is 1.")

        form.addParam('doThreshold', params.BooleanParam,
                      default=False)
        line = form.addLine("CC range (0-1)", condition="doThreshold")
        line.addParam('lowCC', params.IntParam,
                      condition="doThreshold",
                      label='Low limit', default=0)
        line.addParam('highCC', params.IntParam,
                      condition="doThreshold",
                      label='High limit', default=1)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        self.params = {
            'input_parts': os.path.abspath(self.inputSubstacks.get().getFileName()),
            'cc_min': self.lowCC.get(),
            'cc_max': self.highCC.get(),
            'select_refs': bool(self.selectRefs),
            'do_thr_cc': bool(self.doThreshold),
            'refs_list': getListFromValues(self.refsList.get())
        }

        jsonFn = self._getTmpPath("params.json")
        with open(jsonFn, "w") as fn:
            json.dump(self.params, fn, indent=4)

        self.runJob(Plugin.getProgram("subsets.py"),
                    os.path.abspath(jsonFn),
                    env=Plugin.getEnviron(), cwd=self._getExtraPath())

        nRefs = self.getNumRefs()
        substacks = TomoSubStacks(filename=self._getExtraPath("particles.ptclsraw"),
                                  n_ptcl=self.getNumParts(),
                                  n_refs=nRefs)
        self._defineOutputs(**{"outputSubstacks": substacks})

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if not self.selectRefs and not self.doThreshold:
            errors.append("You did not choose any operation!")

        if self.selectRefs:
            origRefs = range(1, self.inputSubstacks.get().getNumRefs()+1)
            keepRefs = getListFromValues(self.refsList.get())

        return errors

    # --------------------------- UTILS functions -----------------------------
    def getNumRefs(self):
        if not self.selectRefs:
            return int(self.inputSubstacks.get().getNumRefs())
        else:
            return len(getListFromValues(self.refsList.get()))
