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
        if params['ctf_corr'] != 'none':
            tomos.set_defocus(i, params['inputAngles'][i].replace(".tlt", ".defocus"))

    tomos.save(os.path.join(output_dir, "input_tomos.tomostxt"))


def createPtclsFile(params, output_dir):
    """ Load DYNAMO table with NUMPY and convert it to PTCLSRAW. """
    parts = np.loadtxt(os.path.join(params['output_dir'], '../tmp/input_particles.tbl'), unpack=True)
    tomos = SUSAN.read(os.path.join(params['output_dir'], "input_tomos.tomostxt"))
    ptcls = SUSAN.data.Particles.import_data(tomograms=tomos,
                                             position=parts[23:26, :].transpose(),
                                             ptcls_id=parts[0],
                                             tomos_id=parts[19])
    ptcls.save(os.path.join(output_dir, 'input_particles.ptclsraw'))


def reconstructAvg(params, output_dir):
    """ Reconstruct an average 3D volume. """
    avgr = SUSAN.modules.Averager()
    avgr.list_gpus_ids = list(params['gpus'])
    avgr.threads_per_gpu = params['thr_per_gpu']
    avgr.ctf_correction = params['ctf_corr']
    avgr.rec_halfsets = bool(params['do_halfsets'])
    avgr.symmetry = params['symmetry']
    avgr.padding_type = params['padding']
    avgr.reconstruct(os.path.join(output_dir, "average"),
                     os.path.join(output_dir, "input_tomos.tomostxt"),
                     os.path.join(output_dir, "input_particles.ptclsraw"),
                     box_size=params['box_size'])


if __name__ == '__main__':
    if len(sys.argv) > 0:
        inputJson = sys.argv[1]

        if os.path.exists(inputJson):
            with open(inputJson) as fn:
                params = json.load(fn)
                output_dir = params['output_dir']

            createTomosFile(params, output_dir)
            createPtclsFile(params, output_dir)
            reconstructAvg(params, output_dir)
        else:
            raise FileNotFoundError(inputJson)
    else:
        print("Usage: %s jsonParamsFile" % os.path.basename(sys.argv[0]))
