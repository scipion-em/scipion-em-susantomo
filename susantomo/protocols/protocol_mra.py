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
from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA
import pyworkflow.utils as pwutils
from pyworkflow.plugin import Domain
from pwem import emlib

from tomo.objects import SetOfCTFTomoSeries, SetOfTiltSeries, AverageSubTomogram
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging

from .. import Plugin
from ..convert import writeDynTable


class outputs(Enum):
    outputAverage = AverageSubTomogram


class ProtSusanMRA(ProtTomoSubtomogramAveraging):
    """ Align subtomograms using multiple references. """
    _label = 'multi-reference alignment'
    _devStatus = BETA
    _possibleOutputs = outputs

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0', label="Choose GPU IDs")

        form.addSection(label='Input')
        form.addParam('inputSetOfSubTomograms', params.PointerParam,
                      pointerClass="SetOfSubTomograms",
                      important=True,
                      label='Input subtomograms',
                      help="Set of subtomograms that will be used to fetch "
                           "alignment and coordinates information.")
        form.addParam('inputTiltSeries',
                      params.PointerParam,
                      pointerClass='SetOfCTFTomoSeries, SetOfTiltSeries',
                      important=True,
                      label='Corresponding CTF tomo series or tilt-series',
                      help='Set of tilt-series that correspond to subtomograms. '
                           'The matching is done using tsId.')
        form.addParam('tomoSize', params.IntParam,
                      default=800, label='Tomogram thickness (px)',
                      help='Z height of a tomogram volume in '
                           'pixels. Required for tilt series stack.')
        form.addParam('boxSize', params.IntParam,
                      default=32, label='Output box size (voxels)',
                      help='Size of the reconstructed average volume in '
                           'voxels. Pixel size will be the same as input '
                           'tilt series.')
        form.addParam('numberOfIters', params.IntParam,
                      label='Iterations', default=5,
                      help="Number of iterations to run.")
        form.addParam('sym', params.StringParam,
                      default='c1',
                      label='Symmetry group',
                      help="Specify the particle symmetry.")
        form.addParam('ctfCorrAln', params.EnumParam,
                      choices=['none', 'on_reference', 'on_substack',
                               'wiener_ssnr', 'wiener_white', 'cfsc'],
                      default=1,
                      label="CTF correction method (aligner)")

        form.addSection(label='References')
        form.addParam('generateRefs', params.BooleanParam,
                      default=False, label='Generate reference(s) and mask(s)',
                      help="Generate a sphere template(s) based on radius provided.")
        form.addParam('nref', params.IntParam, label='Number of references for MRA',
                      default=1, condition="generateRefs")
        form.addParam('refRadius', params.StringParam, default='',
                      condition="generateRefs",
                      label="Radii for generated template(s)",
                      help="Input values separated by a space.")
        form.addParam('maskRadius', params.StringParam, default='',
                      condition="generateRefs",
                      label="Radii for generated mask(s)",
                      help="Input values separated by a space. Values bigger "
                           "than reference radii are expected.")

        form.addParam('inputRefs', params.MultiPointerParam,
                      pointerClass='Volume',
                      allowsNull=True,
                      label="Input reference(s)",
                      condition="not generateRefs")
        form.addParam('inputMasks', params.MultiPointerParam,
                      condition="not generateRefs",
                      label="Alignment mask(s)",
                      pointerClass='VolumeMask',
                      allowsNull=True,
                      help='Need to have the same dimensions as the template. '
                           'Masks order should match the references.')

        form.addSection(label='Angular scanning')
        form.addParam('coneRange', params.IntParam, label='Cone range', default=360,
                      help="The first two Euler angles are used to define the orientation of the vertical axis of the "
                           "protein. First Euler angle (tdrot) rotates the template around its z axis. Second Euler "
                           "angle (tilt) rotates the template around its x axis. Susan scans for this axis inside a "
                           "cone: The 'cone_range' parameter defines the angular aperture of this cone. 360 degrees is "
                           "thus the value for a global scan. To skip the part of the angular search that looks for "
                           "orientations, you have to set 'cone range' to 0, and 'cone_sampling' to 1")
        form.addParam('coneSampling', params.IntParam, label='Cone sampling', default=60,
                      help="This parameter expresses the discretization inside this cone. The sampling is given in "
                           "degrees, and corresponds to a representative angular distance between two neighboring "
                           "orientations inside the cone.")
        form.addParam('inplaneRange', params.IntParam, label='Inplane rotation range', default=360,
                      help='The third Euler angle ("narot") defines rotations about the new axis of the reference. 360 '
                           'degrees is the value for a global scan. To skip the part of the angular search that looks '
                           'for azimuthal rotations , you have to set 1)  "inplane range" to zero, and 2)  '
                           '"inplane_sampling" to 1. Likewise, to skip the part of the angular search that looks for '
                           'orientations, you have to manipulate the "cone" parameters: 1)  "cone range" to zero, and 2'
                           ') "cone_sampling" to 1.')
        form.addParam('inplaneSampling', params.IntParam, label='Inplane rotation sampling', default=45,
                      help='Parameter "inplane_sampling" is associated with the "narot" angle (New Axis ROTation). 1) '
                           'the axis of the template is rotated to a new orientation (defined by "tdrot" and "tilt"). 2'
                           ') the rotated template rotates again on its new axis (an "inplane" rotation). "inplane_'
                           'sampling" defines the angular interval (in degrees) between two of these inplane rotations')
        form.addParam('refine', params.IntParam, label='Refine', default=5,
                      help="How many refinement iterations are carried out on each single particle. This refinement "
                           "when comparing rotations of the reference against the data, takes the best orientation and "
                           "looks again with a finer sampling. The sampling in the refined search will be half of the "
                           "sampling used in the original one.  The range of the refined search encompasses all the "
                           "orientations that neighbour the best orientation found in the original search.")
        form.addParam('refineFactor', params.IntParam, label='Refine factor', default=2,
                      help="Controls the size of the angular neighborhood during the local refinement of the angular "
                           "grid.")

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

        form.addSection(label='Averaging')
        form.addParam('ctfCorrAvg', params.EnumParam,
                      choices=['none', 'phase_flip', 'wiener', 'wiener_ssnr'],
                      default=2,
                      label="CTF correction method (averager)")
        form.addParam('norm', params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['none', 'zero_mean', 'zero_mean_one_std',
                               'zero_mean_proj_weight'],
                      default=2,
                      label="Normalization type")
        form.addParam('padding', params.EnumParam,
                      expertLevel=params.LEVEL_ADVANCED,
                      choices=['zero', 'noise'], default=1,
                      label="Padding type")
        form.addParam('doHalfSets', params.BooleanParam, default=False,
                      label="Reconstruct half-sets?")

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.mraStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Prepare input files. """
        inputSubTomos = self.inputSetOfSubTomograms.get()
        fnTable = self._getTmpPath("input_particles.tbl")
        factor = self.getScaleFactor()
        with open(fnTable, 'w') as fn:
            writeDynTable(fn, inputSubTomos, scaleFactor=factor)

        tsSet = self._getInputTs()
        tsIds_from_ts = set(item.getTsId() for item in tsSet)
        tsIds_from_subtomos = set(inputSubTomos.getTomograms().keys())
        if not tsIds_from_ts.issubset(tsIds_from_subtomos):
            self.warning("Found subtomos with tsId that did not match "
                         "provided tilt-series: "
                         f"{set.difference(tsIds_from_subtomos, tsIds_from_ts)}")

        if self.doCtf():
            # generate defocus files
            setOfCtfTomoSeries = self.inputTiltSeries.get()
            imodUtils = Domain.importFromPlugin('imod.utils')
            for ctf in setOfCtfTomoSeries:
                tsId = ctf.getTsId()
                defocusFilePath = self._getTmpPath(tsId + ".defocus")
                imodUtils.generateDefocusIMODFileFromObject(ctf, defocusFilePath)

        if self.generateRefs:
            self.refs = self.refRadius.get().split()
            self.masks = self.maskRadius.get().split()
        else:
            self.refs = [i.getFileName() for i in self.inputRefs.get()]
            self.masks = [i.getFileName() for i in self.inputMasks.get()]

        self.stacks, self.tilts, self.ids = [], [], []

        for ts in tsSet:
            tsId = ts.getTsId()
            tsInputFn = ts.getFirstItem().getFileName()
            tsFn = self._getTmpPath(tsId + ".mrc")
            tiltFn = self._getTmpPath(tsId + ".tlt")

            # has to be float32
            ih = emlib.image.ImageHandler()
            if ih.getDataType(tsInputFn) != emlib.DT_FLOAT:
                ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)
            elif pwutils.getExt(tsInputFn) in ['.mrc', '.st', '.mrcs', '.ali']:
                pwutils.createAbsLink(os.path.abspath(tsInputFn), tsFn)
            else:
                ih.convert(tsInputFn, tsFn, emlib.DT_FLOAT)

            ts.generateTltFile(tiltFn)
            self.stacks.append(tsFn)
            self.tilts.append(tiltFn)
            self.ids.append(ts.getObjId())

    def mraStep(self):
        """ Run susan_aligner and susan_reconstruct programs. """
        tsSet = self._getInputTs()
        tomo_size = [tsSet.getDim()[0], tsSet.getDim()[1], self.tomoSize.get()]

        params = {
            'ts_nums': self.ids,
            'refs_nums': self.nref.get() if self.generateRefs else len(self.refs),
            'generate_refs': bool(self.generateRefs),
            'inputStacks': self.stacks,
            'inputAngles': self.tilts,
            'inputRefs': self.refs,
            'inputMasks': self.masks,
            'output_dir': self._getExtraPath(),
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
            'offsets': [self.offsetRange.get(), self.offsetStep.get()]
        }

        jsonFn = self._getTmpPath("params.json")
        with open(jsonFn, "w") as fn:
            json.dump(params, fn, indent=4)

        self.runJob(Plugin.getProgram("mra.py"), jsonFn,
                    env=Plugin.getEnviron())

    def createOutputStep(self):
        imgSet = self._getInputTs()
        volume = AverageSubTomogram()
        volumeFile = self._getExtraPath("average_class001.mrc")
        volume.setFileName(volumeFile)
        volume.setSamplingRate(imgSet.getSamplingRate())
        if self.doHalfSets:
            volume.setHalfMaps([self._getExtraPath("average_class001_half1.mrc"),
                                self._getExtraPath("average_class001_half2.mrc")])

        self._defineOutputs(**{outputs.outputAverage.name: volume})
        self._defineSourceRelation(self.inputSetOfSubTomograms, volume)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.doCtf() and isinstance(self.inputTiltSeries.get(),
                                       SetOfTiltSeries):
            errors.append("CTF correction requires that you provide "
                          "CTFTomoSeries as input")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getInputTs(self, pointer=False):
        if isinstance(self.inputTiltSeries.get(), SetOfCTFTomoSeries):
            return self.inputTiltSeries.get().getSetOfTiltSeries(pointer=pointer)
        return self.inputTiltSeries.get() if not pointer else self.inputTiltSeries

    def getScaleFactor(self):
        samplingRateSubtomos = self.inputSetOfSubTomograms.get().getSamplingRate()
        samplingRateTS = self._getInputTs().getSamplingRate()
        return float(samplingRateSubtomos / samplingRateTS)

    def doCtf(self):
        return self.ctfCorrAvg.get() or self.ctfCorrAln.get()
