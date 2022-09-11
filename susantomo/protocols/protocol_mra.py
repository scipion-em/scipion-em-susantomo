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
import sys
import fileinput
import json

import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA
import pyworkflow.utils as pwutils
from pwem import emlib
from pwem.objects import EMSet

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
        form.addParam('inputRefs', params.PointerParam,
                      important=True,
                      allowsNull=True,
                      pointerClass='AverageSubTomogram, '
                                   'SetOfAverageSubTomograms,'
                                   'Volume, SetOfVolumes',
                      label="Input reference(s)",
                      condition="not doContinue")
        form.addParam('inputMasks', params.PointerParam,
                      important=True,
                      allowsNull=True,
                      condition="not doContinue",
                      label="Alignment mask(s)",
                      pointerClass='Volume, VolumeMask, SetOfVolumes',
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
        """ Copy the project and inputs from the previous run. """
        prevRun = self.previousRun.get()
        for p in ["mra", "input"]:
            prevProj = prevRun._getExtraPath(p)
            pwutils.copyTree(prevProj, self._getExtraPath(p))

        # convert previous TS stacks again
        tsSet = prevRun._getInputTs()
        for ts in tsSet:
            tsId = ts.getTsId()
            tsInputFn = ts.getFirstItem().getFileName()
            tsFn = self._getTmpPath(tsId + ".mrc")

            # has to be float32
            ih = emlib.image.ImageHandler()
            if ih.getDataType(tsInputFn) != emlib.DT_FLOAT:
                ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)
            elif pwutils.getExt(tsInputFn) in ['.mrc', '.st', '.mrcs', '.ali']:
                pwutils.createAbsLink(os.path.abspath(tsInputFn), tsFn)
            else:
                ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)

        # update stack paths in tomostxt
        tomosFn = self._getExtraPath("input/input_tomos.tomostxt")
        for line in fileinput.input(tomosFn, inplace=True):
            if "stack_file:" in line:
                tomo = os.path.basename(line.split(":")[-1])
                line = f"stack_file:{os.path.abspath(self._getTmpPath(tomo))}"
            sys.stdout.write(line)

        # count previous refs
        refs = prevRun.inputRefs.get()
        self.refs_nums = refs.getSize() if self.isSet(refs) else 1

    def convertInputRefs(self):
        if not self.isSet(self.inputRefs.get()):
            self.refs.append(os.path.abspath(self.inputRefs.get().getFileName()))
            self.masks.append(os.path.abspath(self.inputMasks.get().getFileName()))
        else:
            self.refs.extend([os.path.abspath(i.getFileName()) for i in self.inputRefs.get()])
            self.masks.extend([os.path.abspath(i.getFileName()) for i in self.inputMasks.get()])

    def runSusanStep(self):
        """ Run susan_aligner and susan_reconstruct programs. """
        self.params = {
            'continue': bool(self.doContinue),
            'refs_nums': self.refs_nums if self.doContinue else len(self.refs),
            'box_size': self.boxSize.get(),
            'gpus': self.getGpuList(),
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

        if not self.doContinue:
            tsSet = self._getInputTs()
            tomo_size = [tsSet.getDim()[0], tsSet.getDim()[1], self.tomoSize.get()]
            self.params.update({
                'ts_nums': self.ids,
                'num_tilts': max([ts.getSize() for ts in tsSet.iterItems()]),
                'pix_size': tsSet.getSamplingRate(),
                'tomo_size': tomo_size,
                'voltage': tsSet.getAcquisition().getVoltage(),
                'sph_aber': tsSet.getAcquisition().getSphericalAberration(),
                'amp_cont': tsSet.getAcquisition().getAmplitudeContrast(),
                'inputStacks': self.stacks,
                'inputAngles': self.tilts,
                'inputRefs': self.refs,
                'inputMasks': self.masks,
                'has_ctf': self.hasCtf()
            })

        jsonFn = self._getTmpPath("params.json")
        with open(jsonFn, "w") as fn:
            json.dump(self.params, fn, indent=4)

        self.runJob(Plugin.getProgram("mra.py"),
                    os.path.abspath(jsonFn),
                    env=Plugin.getEnviron(), cwd=self._getExtraPath())

    def createOutputStep(self):
        if self.doContinue:
            imgSet = self.previousRun.get()._getInputTs()
            imgSetPointer = self.previousRun.get()._getInputTs(pointer=True)
        else:
            imgSet = self._getInputTs()
            imgSetPointer = self._getInputTs(pointer=True)

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
        self._defineSourceRelation(imgSetPointer, volumes)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        refs = self.inputRefs
        masks = self.inputMasks

        if not self.doContinue:
            if self.doCtf() and isinstance(self.inputTiltSeries.get(),
                                           SetOfTiltSeries):
                errors.append("CTF correction requires that you provide "
                              "CTFTomoSeries as input")
            if any([not refs.hasValue(), not masks.hasValue()]):
                errors.append("Input references and masks are required!")

        if self.numberOfIters.get() == 1 and self.incLowpass:
            errors.append("You cannot increase lowpass when doing only 1 iteration.")

        msg = "Number of references and masks must be the same."
        if self.isSet(refs.get()) and self.isSet(masks.get()):
            if refs.get().getSize() != masks.get().getSize():
                errors.append(msg)
        elif self.isSet(refs.get()) and not self.isSet(masks.get()):
            errors.append(msg)
        elif not self.isSet(refs.get()) and self.isSet(masks.get()):
            errors.append(msg)

        if self.doContinue and not self.previousRun.hasValue():
            errors.append("Please input the previous protocol run.")

        return errors

    def _summary(self):
        summary = []

        return summary

    # --------------------------- UTILS functions -----------------------------
    def doCtf(self):
        return self.ctfCorrAvg.get() or self.ctfCorrAln.get()

    def isContinue(self):
        return self.doContinue

    def isSet(self, obj):
        return isinstance(obj, EMSet)
