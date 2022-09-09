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
import numpy as np
from glob import glob
import re

import susan as SUSAN


def createTomosFile(params, output_dir):
    """ Create tomostxt file. """
    n_tomo = len(params['ts_nums'])
    tomos = SUSAN.data.Tomograms(n_tomo=n_tomo,
                                 n_proj=params['num_tilts'])
    for i in range(n_tomo):
        tomos.tomo_id[i] = params['ts_nums'][i]
        tomos.set_stack(i, params['inputStacks'][i])
        tomos.set_angles(i, params['inputAngles'][i])
        tomos.pix_size[i] = params['pix_size']
        tomos.tomo_size[i] = tuple(params['tomo_size'])
        tomos.voltage[i] = params['voltage']
        tomos.amp_cont[i] = params['amp_cont']
        tomos.sph_aber[i] = params['sph_aber']
        if params['ctf_corr_aln'] != 'none':
            tomos.set_defocus(i, params['inputAngles'][i].replace(".tlt", ".defocus"))

    tomos.save(os.path.join(output_dir, "input_tomos.tomostxt"))


def createPtclsFile(params, output_dir):
    """ Load DYNAMO table with NUMPY and convert it to PTCLSRAW. """
    parts = np.loadtxt(os.path.join(params['output_dir'], '../tmp/input_particles.tbl'), unpack=True)
    tomos = SUSAN.read(os.path.join(params['output_dir'], "input_tomos.tomostxt"))
    ptcls = SUSAN.data.Particles.import_data(tomograms=tomos,
                                             position=parts[23:26, :].transpose(),
                                             ptcls_id=parts[0],
                                             tomos_id=parts[19],
                                             randomize_angles=params['randomize'])

    # Duplicate references
    n_refs = params['refs_nums']
    if n_refs > 1:
        for i in range(1, n_refs):
            SUSAN.data.Particles.MRA.duplicate(ptcls, 0)

    ptcls.save(os.path.join(output_dir, 'input_particles.ptclsraw'))


def createRefsFile(params, output_dir):
    """ Create refstxt file. """
    n_refs = params['refs_nums']
    refs = SUSAN.data.Reference(n_refs=n_refs)
    prev_cwd = os.getcwd()
    os.chdir(os.path.join(output_dir))

    if params['generate_refs']:
        for i in range(n_refs):
            SUSAN.io.mrc.write(SUSAN.utils.create_sphere(
                int(params['inputRefs'][i]),  # radius
                params['box_size']),
                f'../tmp/ref{i+1}.mrc',
                params['pix_size'])
            refs.ref[i] = os.path.join(os.getcwd(), f'../tmp/ref{i+1}.mrc')

            SUSAN.io.mrc.write(SUSAN.utils.create_sphere(
                int(params['inputMasks'][i]),  # radius
                params['box_size']),
                f'../tmp/mask{i+1}.mrc',
                params['pix_size'])
            refs.msk[i] = os.path.join(os.getcwd(), f'../tmp/mask{i+1}.mrc')
    else:
        for i in range(n_refs):
            refs.ref[i] = params['inputRefs'][i]
            refs.msk[i] = params['inputMasks'][i]

    refs.save("input_refs.refstxt")
    os.chdir(prev_cwd)


def runAlignment(params, output_dir, doContinue=False):
    """ Execute MRA project in the output_dir folder. """
    prev_cwd = os.getcwd()
    os.chdir(os.path.join(output_dir))
    if doContinue:
        mngr = SUSAN.project.Manager('mra')  # only accepts a name, not a path
        lastIter = _getIterNumber('mra/ite_*') + 1
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

    os.chdir(prev_cwd)


def reconstructAvg(params, output_dir):
    """ Reconstruct an average 3D volume. """
    avgr = SUSAN.modules.Averager()
    avgr.list_gpus_ids = list(params['gpus'])
    avgr.threads_per_gpu = params['thr_per_gpu']
    avgr.ctf_correction = params['ctf_corr_avg']
    avgr.rec_halfsets = bool(params['do_halfsets'])
    avgr.symmetry = params['symmetry']
    avgr.padding_type = params['padding']
    avgr.reconstruct(os.path.join(output_dir, "average"),
                     os.path.join(output_dir, "input_tomos.tomostxt"),
                     os.path.join(output_dir, f"mra/ite_{params['iter']:04d}/particles.ptclsraw"),
                     box_size=params['box_size'])


def _getIterNumber(path):
    """ Return the last iteration number. """
    result = None
    files = sorted(glob(path))
    if files:
        f = files[-1]
        s = re.compile('ite_(\d{4})').search(f)
        if s:
            result = int(s.group(1))  # group 1 is 1 digit iteration number
    return result


if __name__ == '__main__':
    if len(sys.argv) > 0:
        inputJson = sys.argv[1]

        if os.path.exists(inputJson):
            with open(inputJson) as fn:
                params = json.load(fn)
                output_dir = params['output_dir']

            if not params['continue']:
                createTomosFile(params, output_dir)
                createPtclsFile(params, output_dir)
                createRefsFile(params, output_dir)
            runAlignment(params, output_dir, doContinue=params['continue'])
            reconstructAvg(params, output_dir)
        else:
            raise FileNotFoundError(inputJson)
    else:
        print("Usage: %s jsonParamsFile" % os.path.basename(sys.argv[0]))
