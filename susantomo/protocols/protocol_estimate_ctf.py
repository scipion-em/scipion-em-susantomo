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
import sys

import pyworkflow.protocol.params as params
from pyworkflow.constants import NEW
from pwem.protocols import ProtAnalysis3D
from pwem.objects import Volume
from pwem.emlib.image import ImageHandler

from susantomo import Plugin


class ProtSusanEstimateCtf(ProtAnalysis3D):
    """
    Protocol for ...
    """
    _label = 'ctf estimation'
    _devStatus = NEW
    _possibleOutputs = {

    }

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
                  }

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    
    def _insertAllSteps(self):
        #self._createFilenameTemplates()
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------

    def createOutputStep(self):
        sys.path.insert(0, Plugin.getHome())
        import susan as SUSAN

        os.chdir("/lmb/home/gsharov/cephfs/data/susan_tutorial")
        N = 4  # number of tilt series (TS) to process together
        K = 59 # maximal number of projections in processed TS
        tomos = SUSAN.data.Tomograms(None, N, K)
        apix = 2.62            # angstroms per pixel
        tsiz = [3710, 3710, 1400] # pixels (should be the same for all processed TS)

        for i in range(N):
            tomo_base = './mixedCTEM_tomo%d%s'
            tomos.tomo_id[i] = i+1
            tomos.set_stack(i, tomo_base % (i+1, '.b1.ali.mrc'))
            tomos.set_angles(i, tomo_base % (i+1, '.tlt'))
            tomos.pix_size[i] = apix
            tomos.tomo_size[i] = tsiz

        tomos.save('tomos_raw_b1.tomostxt')

        raise Exception("lala")
        test = ""
        self.runJob(test, "",
                    cwd=self._getExtraPath(),
                    env=Plugin.getEnviron())

        #outputs = {'outputVolume1': vol,
        #           'outputVolume2': vol2}
        #self._defineOutputs(**outputs)
        #self._defineSourceRelation(inputVol, vol)
        #self._defineSourceRelation(inputVol, vol2)

    # --------------------------- INFO functions ------------------------------
    
    def _summary(self):
        summary = []

        return summary
    
    def _validate(self):
        errors = []

        return errors
    
    # --------------------------- UTILS functions -----------------------------
