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

from pyworkflow.utils import magentaStr
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from tomo.protocols import ProtImportTs

from ..protocols import ProtSusanEstimateCtf


class TestSusanBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('tomo-em')
        cls.inputSoTS = cls.inputDataSet.getFile('ts1')

    @classmethod
    def _runImportTiltSeries(cls):
        cls.protImportTS = cls.newProtocol(ProtImportTs,
                                           filesPath=os.path.split(cls.inputSoTS)[0],
                                           filesPattern="BB{TS}.st",
                                           anglesFrom=0,
                                           voltage=300,
                                           magnification=105000,
                                           sphericalAberration=2.7,
                                           amplitudeContrast=0.1,
                                           samplingRate=20.2,
                                           doseInitial=0,
                                           dosePerFrame=0.3,
                                           minAngle=-55,
                                           maxAngle=65.0,
                                           stepAngle=2.0,
                                           tiltAxisAngle=-12.5)
        cls.launchProtocol(cls.protImportTS)
        return cls.protImportTS


class TestSusan(TestSusanBase):
    def test_estimateCtf(self):
        print(magentaStr("\n==> Importing data - TiltSeries:"))
        protImport = self._runImportTiltSeries()

        print(magentaStr("\n==> Testing SUSAN: ctf estimation"))
        prot = self.newProtocol(ProtSusanEstimateCtf,
                                inputSetOfTiltSeries=protImport.outputTiltSeries,
                                )
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.outputCTFs,
                             "SetOfCTFTomoSeries has not been produced.")
