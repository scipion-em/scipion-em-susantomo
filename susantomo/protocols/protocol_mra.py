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
from pyworkflow.constants import BETA
import pyworkflow.utils as pwutils
from pwem.objects import Volume, VolumeMask

from tomo.objects import (SetOfTiltSeries, AverageSubTomogram,
                          SetOfAverageSubTomograms)
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging

from .. import Plugin
from .protocol_base import ProtSusanBase


class ProtSusanMRA(ProtSusanBase, ProtTomoSubtomogramAveraging):
    """ Align subtomograms using multiple references. """
    _label = 'multi-reference alignment'
    _devStatus = BETA
    _possibleOutputs = {'outputAverage': AverageSubTomogram,
                        'outputAverages': SetOfAverageSubTomograms}

    def __init__(self, **args):
        ProtSusanBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineProcessParams(self, form):
        form.addParam('numberOfIters', params.IntParam,
                      label='Iterations', default=3,
                      help="Number of iterations to run.")
        form.addParam('ctfCorrAln', params.EnumParam,
                      choices=['none', 'on_reference', 'on_substack',
                               'wiener_ssnr', 'wiener_white', 'cfsc'],
                      default=1,
                      label="CTF correction method (aligner)")

        form.addSection(label='References')
        form.addParam('contMsg', params.LabelParam,
                      condition="doContinue",
                      label="In continue mode these options are not available.")
        form.addParam('generateRefs', params.BooleanParam,
                      condition="not doContinue",
                      default=False,
                      label='Generate reference(s) and mask(s)',
                      help="Generate a sphere template(s) based "
                           "on radius provided.")
        form.addParam('nref', params.IntParam,
                      label='Number of references for MRA',
                      default=2, condition="generateRefs and not doContinue")
        form.addParam('refRadius', params.StringParam, default='',
                      condition="generateRefs and not doContinue",
                      label="Radii for generated template(s)",
                      help="Input values separated by a space.")
        form.addParam('maskRadius', params.StringParam, default='',
                      condition="generateRefs and not doContinue",
                      label="Radii for generated mask(s)",
                      help="Input values separated by a space. Values bigger "
                           "than reference radii are expected.")

        form.addParam('inputRefs', params.PointerParam,
                      pointerClass='AverageSubTomogram, '
                                   'SetOfAverageSubTomograms,'
                                   'Volume, SetOfVolumes',
                      allowsNull=True,
                      label="Input reference(s)",
                      condition="not generateRefs and not doContinue")
        form.addParam('inputMasks', params.PointerParam,
                      condition="not generateRefs and not doContinue",
                      label="Alignment mask(s)",
                      pointerClass='VolumeMask, SetOfVolumes',
                      allowsNull=True,
                      help='Need to have the same dimensions as the template. '
                           'Masks order should match the references.')

        form.addSection(label='Angular scanning')
        form.addParam('randomizeAngles', params.BooleanParam,
                      label='Rotate subtomograms randomly at the start?',
                      condition="not doContinue",
                      default=False)
        form.addParam('coneRange', params.IntParam,
                      label='Cone range', default=360,
                      help="The first two Euler angles are used to define "
                           "the orientation of the vertical axis of the "
                           "protein. First Euler angle (tdrot) rotates the "
                           "template around its z axis. Second Euler "
                           "angle (tilt) rotates the template around its "
                           "x axis. Susan scans for this axis inside a "
                           "cone: The 'cone range' parameter defines the "
                           "angular aperture of this cone. 360 degrees is "
                           "thus the value for a global scan. To skip the "
                           "part of the angular search that looks for "
                           "orientations, you have to set 'cone range' "
                           "to 0, and 'cone sampling' to 1.")
        form.addParam('coneSampling', params.IntParam,
                      label='Cone sampling', default=60,
                      help="This parameter expresses the discretization "
                           "inside the cone. The sampling is given in "
                           "degrees, and corresponds to a representative "
                           "angular distance between two neighboring "
                           "orientations inside the cone.")
        form.addParam('inplaneRange', params.IntParam,
                      label='In-plane rotation range', default=360,
                      help="The third Euler angle (narot) defines rotations "
                           "about the new axis of the reference. 360 "
                           "degrees is the value for a global scan. To skip "
                           "the part of the angular search that looks "
                           "for azimuthal rotations, you have to set "
                           "'in-plane range' to 0 and "
                           "'in-plane sampling' to 1")
        form.addParam('inplaneSampling', params.IntParam,
                      label='In-plane rotation sampling', default=45,
                      help="Defines the angular interval (in degrees) for "
                           "in-plane rotations")
        form.addParam('refine', params.IntParam,
                      label='Refine', default=5,
                      help="When Dynamo finds the triplet [tdrot,tilt,narot] that "
                           "maximizes the similarity between rotated template and "
                           "particle, a new grid is generated around this angles. "
                           "The 'refine' parameter defines how many times this "
                           "process will be repeated for each particle.")
        form.addParam('refineFactor', params.IntParam,
                      label='Refine factor', default=2,
                      help="This defines the range of the new set of angles. "
                           "A value of 2, for instance, means that the range "
                           "of the new grid should twice as large as the step "
                           "of the old grid.")

        form.addSection(label='Shifts & thresholds')
        form.addParam('allowDrift', params.BooleanParam, default=True,
                      label="Allow drift?",
                      help="Take into account shifts from previous iteration. "
                           "Equivalent to DYNAMO's area search modus 2.")
        form.addParam('offsetRange', params.IntParam, default=20,
                      label="Offset range (px)")
        form.addParam('offsetStep', params.IntParam, default=2,
                      label="Offset step (px)")
        form.addParam('threshold', params.FloatParam,
                      label='CC threshold', default=0.8,
                      help='Threshold value for CC (cross-correlation). '
                           'Determines the fraction of particles to be used '
                           'for reconstruction.')
        line = form.addLine("Band-pass filter (freq.)")
        line.addParam('low', params.IntParam, label='Low frequency', default=32)
        line.addParam('high', params.IntParam, label='High frequency', default=2)
        form.addParam('incLowpass', params.BooleanParam, default=False,
                      label="Increase the lowpass filter on each iteration?")

    # --------------------------- STEPS functions -----------------------------
    def continueStep(self):
        """ Copy project from the previous run. """
        mraDir = self._getExtraPath("mra")
        pwutils.makePath(mraDir)
        prevRun = self.previousRun.get()
        prevProj = prevRun._getExtraPath("mra")

        pwutils.copyTree(prevProj, mraDir)
        pwutils.copyPattern(prevRun._getExtraPath("input_tomos.tomostxt"),
                            self._getExtraPath())

    def convertInputRefs(self):
        if self.generateRefs:
            self.refs = self.refRadius.get().split()
            self.masks = self.maskRadius.get().split()
        else:
            if self.isSingleObject(self.inputRefs.get()):
                self.refs = [os.path.abspath(self.inputRefs.get().getFileName())]
                self.masks = [os.path.abspath(self.inputMasks.get().getFileName())]
            else:
                self.refs = [os.path.abspath(i.getFileName()) for i in self.inputRefs.get()]
                self.masks = [os.path.abspath(i.getFileName()) for i in self.inputMasks.get()]

    def runSusanStep(self):
        """ Run susan_aligner and susan_reconstruct programs. """
        tsSet = self._getInputTs()
        tomo_size = [tsSet.getDim()[0], tsSet.getDim()[1], self.tomoSize.get()]

        self.params = {
            'continue': bool(self.doContinue),
            'ts_nums': self.ids,
            'refs_nums': self.nref.get() if self.generateRefs else len(self.refs),
            'generate_refs': bool(self.generateRefs),
            'inputStacks': self.stacks,
            'inputAngles': self.tilts,
            'inputRefs': self.refs,
            'inputMasks': self.masks,
            'num_tilts': max([ts.getSize() for ts in tsSet.iterItems()]),
            'pix_size': tsSet.getSamplingRate(),
            'tomo_size': tomo_size,
            'box_size': self.boxSize.get(),
            'gpus': self.getGpuList(),
            'voltage': tsSet.getAcquisition().getVoltage(),
            'sph_aber': tsSet.getAcquisition().getSphericalAberration(),
            'amp_cont': tsSet.getAcquisition().getAmplitudeContrast(),
            'thr_per_gpu': self.numberOfThreads.get(),
            'ctf_corr_avg': self.getEnumText('ctfCorrAvg'),
            'ctf_corr_aln': self.getEnumText('ctfCorrAln'),
            'do_halfsets': self.doHalfSets.get(),
            'symmetry': self.sym.get(),
            'padding': self.getEnumText('padding'),
            'iter': self.numberOfIters.get(),
            'allow_drift': bool(self.allowDrift),
            'cc': self.threshold.get(),
            'low': self.low.get(),
            'high': self.high.get(),
            'refine': self.refine.get(),
            'refine_factor': self.refineFactor.get(),
            'inc_lowpass': bool(self.incLowpass),
            'angles': [self.coneRange.get(), self.coneSampling.get(),
                       self.inplaneRange.get(), self.inplaneSampling.get()],
            'offsets': [self.offsetRange.get(), self.offsetStep.get()],
            'randomize': bool(self.randomizeAngles)
        }

        jsonFn = self._getTmpPath("params.json")
        with open(jsonFn, "w") as fn:
            json.dump(self.params, fn, indent=4)

        self.runJob(Plugin.getProgram("mra.py"),
                    os.path.abspath(jsonFn),
                    env=Plugin.getEnviron(), cwd=self._getExtraPath())

    def createOutputStep(self):
        imgSet = self._getInputTs()
        pixSize = imgSet.getSamplingRate()
        nRefs = self.params['refs_nums']
        volume = AverageSubTomogram()

        def _createVolume(ind):
            volumeFile = self._getExtraPath(f"average_class{ind:03d}.mrc")
            volume.setObjId(None)
            volume.setFileName(volumeFile)
            volume.setSamplingRate(pixSize)
            if self.doHalfSets:
                volume.setHalfMaps([self._getExtraPath(f"average_class{ind:03d}_half1.mrc"),
                                    self._getExtraPath(f"average_class{ind:03d}_half2.mrc")])
            return volume

        if nRefs > 1:
            volumes = self._createSetOfAverageSubTomograms()
            volumes.setSamplingRate(pixSize)
            for i in range(nRefs):
                volume = _createVolume(i+1)
                volumes.append(volume)
        else:
            volumes = _createVolume(1)

        self._defineOutputs(**{f"outputAverage{'s' if nRefs > 1 else ''}": volumes})
        self._defineSourceRelation(self._getInputTs(pointer=True), volumes)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.doCtf() and isinstance(self.inputTiltSeries.get(),
                                       SetOfTiltSeries):
            errors.append("CTF correction requires that you provide "
                          "CTFTomoSeries as input")

        if self.numberOfIters.get() == 1 and self.incLowpass:
            errors.append("You cannot increase lowpass when doing only 1 iteration.")

        if (self.isSingleObject(self.inputRefs.get()) and
                not self.isSingleObject(self.inputMasks.get())):
            errors.append("Number of references and masks must be the same.")

        if self.doContinue and not self.previousRun.hasValue():
            errors.append("Please input the previous protocol run.")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def doCtf(self):
        return self.ctfCorrAvg.get() or self.ctfCorrAln.get()

    def isSingleObject(self, obj):
        return (isinstance(obj, Volume) or isinstance(obj, AverageSubTomogram) or
                isinstance(obj, VolumeMask))

    def isContinue(self):
        return self.doContinue
