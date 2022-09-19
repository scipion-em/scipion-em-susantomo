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
import fileinput
import sys

import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA

from tomo.objects import AverageSubTomogram, SetOfAverageSubTomograms
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging

from .. import Plugin
from ..objects import TomoSubStacks
from .protocol_base import ProtSusanBase


class ProtSusanMRA(ProtSusanBase, ProtTomoSubtomogramAveraging):
    """ Align subtomograms using multiple references. """
    _label = 'multi-reference alignment'
    _devStatus = BETA
    _possibleOutputs = {'outputAverage': AverageSubTomogram,
                        'outputAverages': SetOfAverageSubTomograms,
                        'outputSubstacks': TomoSubStacks}

    def __init__(self, **args):
        ProtSusanBase.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineProcessParams(self, form):
        form.addParam('numberOfIters', params.IntParam,
                      label='Iterations', default=3,
                      help="Number of iterations to run in a given round.")
        form.addParam('ctfCorrAln', params.EnumParam,
                      choices=['none', 'on_reference', 'on_substack',
                               'wiener_ssnr', 'wiener_white', 'cfsc'],
                      default=1,
                      label="CTF correction method (aligner)")

        form.addSection(label='References')
        form.addParam('reuseRefs', params.BooleanParam, default=False,
                      condition='doContinue',
                      label="Re-use references and masks from the previous run?")
        form.addParam('inputRefs', params.MultiPointerParam,
                      condition="not reuseRefs",
                      important=True,
                      allowsNull=True,
                      pointerClass='AverageSubTomogram, Volume',
                      label="Input references")
        form.addParam('inputMasks', params.MultiPointerParam,
                      condition="not reuseRefs",
                      important=True,
                      allowsNull=True,
                      label="Alignment masks",
                      pointerClass='Volume, VolumeMask',
                      help='Need to have the same dimensions as the template. '
                           'Masks order should match the references.')

        form.addSection(label='Angular scanning')
        form.addParam('randomizeAngles', params.BooleanParam,
                      label='Randomize input orientations?',
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
                      label='In-plane rotation sampling', default=60,
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
                      help="This controls the size of the angular neighbourhood "
                           "during the local refinement of the angular grid. "
                           "In particular, the scanning range at the refinement "
                           "level [i+1] will be determined by the scanning step "
                           "at level [i].\n\n Example:\n"
                           "cone_range[i+1] = refine_factor * cone_sampling[i];\n"
                           "inplane_range[i+1] = refine_factor * inplane_sampling[i];")

        form.addSection(label='Shifts & thresholds')
        form.addParam('allowDrift', params.BooleanParam, default=True,
                      label="Allow drift?",
                      help="Take into account shifts from previous iteration. "
                           "Equivalent to DYNAMO's area search modus 2.")
        form.addParam('offsetRange', params.IntParam, default=12,
                      label="Offset range (px)")
        form.addParam('offsetStep', params.IntParam, default=4,
                      label="Offset step (px)")
        form.addParam('threshold', params.FloatParam,
                      label='CC threshold', default=0.9,
                      help='Threshold value for CC (cross-correlation). '
                           'Determines the fraction of particles to be used '
                           'for reconstruction.')
        line = form.addLine("Band-pass filter (freq.)")
        line.addParam('low', params.IntParam,
                      label='Low frequency (px)', default=11)
        line.addParam('high', params.IntParam,
                      label='High frequency (px)', default=2)
        form.addParam('incLowpass', params.BooleanParam, default=False,
                      label="Increase the lowpass filter on each iteration?")

    # --------------------------- STEPS functions -----------------------------
    def convertInputRefs(self):
        if not self.reuseRefs:
            self.refs.extend([os.path.abspath(i.get().getFileName()) for i in self.inputRefs])
            self.masks.extend([os.path.abspath(i.get().getFileName()) for i in self.inputMasks])
        else:
            # replace relative paths with absolute
            prevRunDir = self.inputSubstacks.get()._getExtraPath("mra")
            refsFn = self._getExtraPath("input/input_refs.refstxt")
            with fileinput.input(refsFn, inplace=True) as fn:
                for line in fn:
                    line = line.replace(":mra",
                                        f":{os.path.abspath(prevRunDir)}")
                    sys.stdout.write(line)

    def runSusanStep(self):
        """ Run susan_aligner and susan_reconstruct programs. """
        tsSet = self._getInputTs()
        tomo_size = [tsSet.getDim()[0], tsSet.getDim()[1], self.tomoSize.get()]

        self.params = {
            'continue': bool(self.doContinue),
            'reuse_refs': bool(self.reuseRefs),
            'refs_nums': self.getNumRefs(),
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
            'has_ctf': self.hasCtf(),
            'box_size': self.boxSize.get(),
            'gpus': self.getGpuList(),
            'thr_per_gpu': self.numberOfThreads.get(),
            'ctf_corr_avg': self.getEnumText('ctfCorrAvg'),
            'ctf_corr_aln': self.getEnumText('ctfCorrAln'),
            'do_halfsets': bool(self.doHalfSets.get()),
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
        pixSize = self._getInputTs().getSamplingRate()
        nRefs = self.getNumRefs()
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

        if nRefs > 1:
            volumes = self._createSetOfAverageSubTomograms()
            volumes.setSamplingRate(pixSize)
            for i in range(nRefs):
                volume = _createVolume(i+1)
                volume.setClassId(i+1)
                volumes.append(volume)
        else:
            volumes = _createVolume(1)

        self._defineOutputs(**{f"outputAverage{'s' if nRefs > 1 else ''}": volumes})
        self._defineSourceRelation(self._getInputTs(pointer=True), volumes)

        substacks = TomoSubStacks(filename=self._getFileName("ptcls",
                                                             iter=self._lastIter()),
                                  n_ptcl=self.getNumParts(),
                                  n_refs=nRefs)
        self._defineOutputs(**{"outputSubstacks": substacks})

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = self._validateBase()

        refs = self.inputRefs
        masks = self.inputMasks

        if self.numberOfIters.get() == 1 and self.incLowpass:
            errors.append("You cannot increase lowpass when doing only 1 iteration.")

        if not refs.hasValue() or not masks.hasValue():
            errors.append("Input references and masks are required!")

        if len(refs) != len(masks):
            errors.append("Number of references and masks must be the same.")

        pix_sizes = [self._getInputTs().getSamplingRate()]
        pix_sizes.extend([r.get().getSamplingRate() for r in refs])
        pix_sizes.extend([m.get().getSamplingRate() for m in masks])
        if len(set(pix_sizes)) != 1:
            errors.append("Pixel size of input tilt-series, references and "
                          "masks must be the same.")

        return errors

    def _summary(self):
        summary = []

        return summary

    # --------------------------- UTILS functions -----------------------------
    def doCtf(self):
        return self.ctfCorrAvg.get() or self.ctfCorrAln.get()

    def getNumParts(self):
        if self.doContinue:
            return self.inputSubstacks.get().getSize()
        else:
            return self.inputSetOfSubTomograms.get().getSize()

    def getNumRefs(self):
        if self.doContinue:
            if self.reuseRefs:
                return int(self.inputSubstacks.get().getNumRefs())
            else:
                return len(self.inputRefs)
        else:
            return len(self.inputRefs)
