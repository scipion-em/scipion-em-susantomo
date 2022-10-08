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

import re
from glob import glob
import numpy as np

import susan as SUSAN


def createTomosFile(params):
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
        if params['has_ctf']:
            tomos.set_defocus(i, params['inputAngles'][i].replace(".tlt", ".defocus"))

    if n_tomo == 1:
        output = f"tomo{params['ts_nums'][0]}.tomostxt"
    else:
        output = "input/input_tomos.tomostxt"

    tomos.save(output)


def createPtclsFile(params, n_refs, do_continue=False):
    """ Load DYNAMO table with NUMPY and convert it to PTCLSRAW. """
    if do_continue:
        ptcls = SUSAN.data.Particles(filename="input/input_particles.ptclsraw")
    else:
        parts = np.loadtxt("../tmp/input_particles.tbl", unpack=True)
        tomos = SUSAN.read("input/input_tomos.tomostxt")
        randomize = params.get("randomize", False)
        ptcls = SUSAN.data.Particles.import_data(tomograms=tomos,
                                                 position=parts[23:26, :].transpose(),
                                                 ptcls_id=parts[0],
                                                 tomos_id=parts[19],
                                                 randomize_angles=randomize)
    # Duplicate reference indexes
    if n_refs > 1:
        for _ in range(n_refs-1):
            SUSAN.data.Particles.MRA.duplicate(ptcls, 0)

    ptcls.save("input/input_particles.ptclsraw")


def createRefsFile(params, n_refs):
    """ Create refstxt file. """
    refs = SUSAN.data.Reference(n_refs=n_refs)

    for i in range(n_refs):
        refs.ref[i] = params['inputRefs'][i]
        refs.msk[i] = params['inputMasks'][i]

    refs.save("input/input_refs.refstxt")


def postProcess(params, mngr, n_refs=1, iter=1):
    """ Apply FOM or l0 filter to output maps. """
    maps = [mngr.get_names_map(iter, ref=i) for i in range(1, n_refs+1)]

    for map in maps:
        v, apix = SUSAN.io.mrc.read(map)
        res = [v]
        if params['apply_fom']:
            # Denoise reference with FOM [Sindelar and Grigorieff, 2012]
            v_f = SUSAN.utils.apply_FOM(res[-1], mngr.get_fsc(iter))
            res.append(v_f)
        if params['apply_l0']:
            # l0-norm: Using M-sparse constraint [Blumensath and Davies, 2008]
            v_f = SUSAN.utils.denoise_l0(res[-1], l0_lambda=0.05)
            res.append(v_f)
        map = map.replace(".mrc", "_denoised.mrc")
        SUSAN.io.mrc.write(res[-1], map, apix)
        del res


def getIterNumber(path):
    """ Return the last iteration number. """
    result = None
    files = sorted(glob(path))
    if files:
        f = files[-1]
        s = re.compile('ite_(\d{4})').search(f)
        if s:
            result = int(s.group(1))  # group 1 is 1 digit iteration number
    return result
