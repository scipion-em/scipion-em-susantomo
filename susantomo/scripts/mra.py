#!/usr/bin/env python3
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

"""
This script should be launched using the SUSAN python interpreter
inside its conda environment.
"""

import sys
import json

import susan as SUSAN

from base import *


def runAlignment(params, doContinue=False):
    """ Execute MRA project in the output_dir folder. """
    if doContinue:
        mngr = SUSAN.project.Manager('mra')  # only accepts a name, not a path
        lastIter = getIterNumber('mra/ite_*') + 1
        if lastIter is None:
            raise Exception("Could not find last iteration number")
        mngr.initial_reference = f"mra/ite_{lastIter:04d}/reference.refstxt"
        mngr.initial_particles = f"mra/ite_{lastIter:04d}/particles.ptclsraw"
        mngr.tomogram_file = "input_tomos.tomostxt"
    else:
        lastIter = 1
        mngr = SUSAN.project.Manager('mra', box_size=params['box_size'])
        mngr.initial_reference = "input_refs.refstxt"
        mngr.initial_particles = "input_particles.ptclsraw"
        mngr.tomogram_file = "input_tomos.tomostxt"

    mngr.list_gpus_ids = list(params['gpus'])
    mngr.threads_per_gpu = params['thr_per_gpu']
    mngr.aligner.ctf_correction = params['ctf_corr_aln']
    mngr.aligner.allow_drift = params['allow_drift']
    mngr.averager.ctf_correction = params['ctf_corr_avg']
    mngr.cc_threshold = params['cc']
    mngr.iteration_type = 3  # 3D Alignment

    mngr.aligner.bandpass.highpass = params['high']
    #mngr.aligner.bandpass.lowpass = params['low']
    mngr.aligner.set_angular_search(*params['angles'])
    mngr.aligner.refine.levels = params['refine']
    mngr.aligner.refine.factor = params['refine_factor']
    mngr.aligner.set_offset_search(*params['offsets'])

    inc_lp = params['inc_lowpass']
    lp = params['low']
    n_refs = params['refs_nums']

    for i in range(lastIter, params['iter'] + 1):
        mngr.aligner.bandpass.lowpass = lp
        bp = mngr.execute_iteration(i)
        if n_refs > 1:
            bp = max(max(bp), 1.0)  # avoid 0.0
        if i == 1 or not inc_lp:
            lp = bp
        else:
            # Enforce a gradual increase in the lowpass
            lp = min(lp + 2, bp)


def reconstructAvg(params):
    """ Reconstruct an average 3D volume. """
    avgr = SUSAN.modules.Averager()
    avgr.list_gpus_ids = list(params['gpus'])
    avgr.threads_per_gpu = params['thr_per_gpu']
    avgr.ctf_correction = params['ctf_corr_avg']
    avgr.rec_halfsets = bool(params['do_halfsets'])
    avgr.symmetry = params['symmetry']
    avgr.padding_type = params['padding']
    lastIter = getIterNumber('mra/ite_*')
    avgr.reconstruct("average", "input_tomos.tomostxt",
                     f"mra/ite_{lastIter:04d}/particles.ptclsraw",
                     box_size=params['box_size'])


if __name__ == '__main__':
    if len(sys.argv) > 0:
        inputJson = sys.argv[1]

        if os.path.exists(inputJson):
            with open(inputJson) as fn:
                params = json.load(fn)

            if not params['continue']:
                createTomosFile(params)
                createPtclsFile(params, params['refs_nums'])
                createRefsFile(params, params['refs_nums'])

            runAlignment(params, doContinue=params['continue'])
            reconstructAvg(params)
        else:
            raise FileNotFoundError(inputJson)
    else:
        print("Usage: %s jsonParamsFile" % os.path.basename(sys.argv[0]))
