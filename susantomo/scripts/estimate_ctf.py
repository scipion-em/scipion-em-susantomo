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

import susan as SUSAN

from base import createTomosFile


def create2DGrid(params, ts_id):
    """ Create a 2D grid to estimate the CTF. """
    tomos = SUSAN.read(f"tomo{ts_id}.tomostxt")
    grid = SUSAN.data.Particles.grid_2d(tomos, step_pixels=params['sampling'])
    grid.save("grid_ctf.ptclsraw")


def estimateCtf(params, ts_id):
    """ Run CTF estimator. """
    ctf_est = SUSAN.modules.CtfEstimator()
    ctf_est.binning = params['binning']
    ctf_est.list_gpus_ids = params['gpus']  # ID's of GPUs to use
    ctf_est.threads_per_gpu = params['thr_per_gpu']
    ctf_est.resolution_angs.min_val = params['min_res']  # angstroms
    ctf_est.resolution_angs.max_val = params['max_res']  # angstroms
    ctf_est.defocus_angstroms.min_val = params['def_min']  # angstroms
    ctf_est.defocus_angstroms.max_val = params['def_max']  # angstroms
    ctf_est.estimate('ctf_grid', f"tomo{ts_id}.tomostxt",
                     "grid_ctf.ptclsraw", params['patch_size'])


if __name__ == '__main__':
    if len(sys.argv) > 0:
        inputJson = sys.argv[1]

        if os.path.exists(inputJson):
            with open(inputJson) as fn:
                params = json.load(fn)
                ts_id = params['ts_nums'][0]

            createTomosFile(params)
            create2DGrid(params, ts_id)
            estimateCtf(params, ts_id)
        else:
            raise FileNotFoundError(inputJson)
    else:
        print("Usage: %s jsonParamsFile" % os.path.basename(sys.argv[0]))
