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

import os
import sys
import json
import pickle
import numpy as np

import susan as SUSAN

from base import *


def runAlignment(params):
    """ Execute MRA project in the output_dir folder. """
    mngr = SUSAN.project.Manager('mra', box_size=params['box_size'])
    mngr.initial_reference = "input/input_refs.refstxt"
    mngr.initial_particles = "input/input_particles.ptclsraw"
    mngr.tomogram_file = "input/input_tomos.tomostxt"

    mngr.list_gpus_ids = params['gpus']
    mngr.threads_per_gpu = params['thr_per_gpu']
    mngr.cc_threshold = params['cc']
    mngr.iteration_type = params['align_type']

    mngr.averager.ctf_correction = params['ctf_corr_avg']
    mngr.averager.rec_halfsets = params['do_halfsets']
    mngr.averager.symmetry = params['symmetry']
    mngr.averager.padding_type = params['padding']

    mngr.aligner.ctf_correction = params['ctf_corr_aln']
    mngr.aligner.allow_drift = params['allow_drift']
    mngr.aligner.bandpass.highpass = params['high']
    mngr.aligner.bandpass.rolloff = params['rolloff']
    mngr.aligner.set_angular_search(*params['angles'])
    mngr.aligner.refine.levels = params['refine']
    mngr.aligner.refine.factor = params['refine_factor']
    mngr.aligner.set_offset_search(*params['offsets'])

    inc_lp = params['inc_lowpass']
    auto = params['auto_step']
    lp = params['low']
    n_refs = params['refs_nums']
    fsc = {}
    cc = {}

    for i in range(1, params['iter'] + 1):
        if auto:  # adjust angular range/step every iter
            ang_stp = np.rad2deg(np.arctan2(1, lp))
            ang_rng = params['range_factor'] * ang_stp
            if i > 1:
                mngr.aligner.set_angular_search(ang_rng, ang_stp, ang_rng, ang_stp)
                if params['refine'] == 0:
                    mngr.aligner.refine.levels = 1
        mngr.aligner.bandpass.lowpass = lp
        bp = mngr.execute_iteration(i)
        if n_refs > 1:
            bp = max(max(bp), 1.0)  # avoid 0.0
        if i == 1 or not inc_lp:
            lp = params['low']  # keep const
        else:
            # Enforce a gradual increase in the lowpass
            lp = min(lp + 2, bp)
        # save FSC and CC
        for n in range(1, n_refs+1):
            if n not in fsc:
                fsc[n] = [mngr.get_fsc(ite=i, ref=n)]
                cc[n] = [mngr.get_cc(ite=i, ref=n)]
            else:
                fsc[n].append(mngr.get_fsc(ite=i, ref=n))
                cc[n].append(mngr.get_cc(ite=i, ref=n))

        with open('mra/info.pkl', 'wb') as f:
            pickle.dump([fsc, cc], f)

        if params['apply_fom'] or params['apply_l0']:
            postProcess(params, mngr, n_refs=n_refs, iter=i)


if __name__ == '__main__':
    if len(sys.argv) > 0:
        inputJson = sys.argv[1]

        if os.path.exists(inputJson):
            with open(inputJson) as fn:
                params = json.load(fn)

            createTomosFile(params)

            if not params['reuse_refs']:
                createRefsFile(params, params['refs_nums'])

            createPtclsFile(params, params['refs_nums'], params['continue'])
            runAlignment(params)
        else:
            raise FileNotFoundError(inputJson)
    else:
        print("Usage: %s jsonParamsFile" % os.path.basename(sys.argv[0]))
