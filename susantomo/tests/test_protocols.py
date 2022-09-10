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
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr

from tomo.protocols import ProtImportTs, ProtImportCoordinates3D
from imod.protocols import ProtImodImportSetOfCtfTomoSeries, ProtImodTomoReconstruction
from ..protocols import ProtSusanEstimateCtf, ProtSusanMRA, ProtSusanAverage


class TestBase(BaseTest):
    @classmethod
    def runImportTiltSeries(cls, **kwargs):
        cls.protImportTS = cls.newProtocol(ProtImportTs, **kwargs)
        cls.launchProtocol(cls.protImportTS)
        cls.assertIsNotNone(cls.protImportTS.outputTiltSeries, "SetOfTiltSeries has not been produced.")
        return cls.protImportTS

    @classmethod
    def runImportCtf(cls, **kwargs):
        cls.protImportCtf = cls.newProtocol(ProtImodImportSetOfCtfTomoSeries, **kwargs)
        cls.launchProtocol(cls.protImportCtf)
        cls.assertIsNotNone(cls.protImportCtf.outputSetOfCTFTomoSeries, "SetOfCTFTomoSeries has not been produced.")
        return cls.protImportCtf

    @classmethod
    def runImportCoords(cls, **kwargs):
        cls.protImportCoords = cls.newProtocol(ProtImportCoordinates3D, **kwargs)
        cls.launchProtocol(cls.protImportCoords)
        cls.assertIsNotNone(cls.protImportCoords.outputCoordinates, "SetOfCoordinates3D has not been produced.")
        return cls.protImportCoords


class TestSusanCtf(TestBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('tutorialDataImodCTF')
        cls.inputSoTS = cls.inputDataSet.getFile('tsCtf1')

        print(magentaStr("\n==> Importing data - tilt series:"))
        cls.protImportTS = cls.runImportTiltSeries(filesPath=os.path.dirname(cls.inputSoTS),
                                                   filesPattern="WTI042413_1series4.mdoc",
                                                   voltage=300,
                                                   sphericalAberration=2.7,
                                                   amplitudeContrast=0.07,
                                                   anglesFrom=2)

    def testCtfEstimation(self):
        print(magentaStr("\n==> Testing susan - ctf estimation:"))
        protCTF = ProtSusanEstimateCtf()
        protCTF.inputTiltSeries.set(self.protImportTS.outputTiltSeries)
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputSetOfCTFTomoSeries, "SetOfCTFTomoSeries has not been produced.")


class TestSusanMRA(TestBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.inputDataSet = DataSet.getDataSet('empiar10064')
        cls.path = cls.inputDataSet.getPath()

        print(magentaStr("\n==> Importing data - tilt series:"))
        cls.protImportTS = cls.runImportTiltSeries(filesPath=cls.path,
                                                   filesPattern="mixedCTEM_{TS}.mrcs",
                                                   anglesFrom=2,  # tlt file
                                                   voltage=300,
                                                   sphericalAberration=2.7,
                                                   amplitudeContrast=0.07,
                                                   magnification=50000,
                                                   samplingRate=20.96,
                                                   tiltAxisAngle=85.0,
                                                   dosePerFrame=1.67)

        print(magentaStr("\n==> Importing data - tomo CTFs:"))
        cls.protImportCtf = cls.runImportCtf(filesPath=cls.path,
                                             filesPattern="mixedCTEM_tomo*.defocus",
                                             inputSetOfTiltSeries=cls.protImportTS.outputTiltSeries)

        print(magentaStr("\n==> Reconstructing tomograms:"))
        protRecon = ProtImodTomoReconstruction(tomoThickness=110.0)
        protRecon.inputSetOfTiltSeries.set(cls.protImportTS.outputTiltSeries)
        cls.launchProtocol(protRecon)
        cls.assertIsNotNone(protRecon.outputSetOfTomograms, "SetOfTomograms has not been produced.")

        print(magentaStr("\n==> Importing data - coordinates 3D:"))
        cls.protImportCoords = cls.runImportCoords(filesPath=cls.path,
                                                   filesPattern="mixed*.tbl",
                                                   importFrom=3,  # dynamo
                                                   samplingRate=20.96,
                                                   boxSize=32,
                                                   importTomograms=protRecon.outputSetOfTomograms)

        print(magentaStr("\n==> Testing susan - average and reconstruct:"))
        cls.protAvg = ProtSusanAverage(tomoSize=110, boxSize=32)
        cls.protAvg.inputSetOfCoords3D.set(cls.protImportCoords.outputCoordinates)
        cls.protAvg.inputTiltSeries.set(cls.protImportCtf.outputSetOfCTFTomoSeries)
        cls.launchProtocol(cls.protAvg)
        cls.assertIsNotNone(cls.protAvg.outputAverage, "AverageSubTomogram has not been produced.")

    def testMRA(self):
        print(magentaStr("\n==> Testing susan - multi-reference alignment:"))
        protMRA = ProtSusanMRA(tomoSize=110, boxSize=32, numberOfIters=2,
                               coneRange=0, coneSampling=1, inplaneRange=0, inplaneSampling=1,
                               refine=2, refineFactor=1, allowDrift=True)
        protMRA.inputSetOfCoords3D.set(self.protImportCoords.outputCoordinates)
        protMRA.inputTiltSeries.set(self.protImportCtf.outputSetOfCTFTomoSeries)
        protMRA.inputRefs.set(self.protAvg.outputAverage)
        protMRA.inputMasks.set(self.protAvg.outputAverage)
        self.launchProtocol(protMRA)
        self.assertIsNotNone(protMRA.outputAverage, "AverageSubtomogram has not been produced.")
