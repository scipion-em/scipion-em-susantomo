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

from base import *


def reconstructAvg(params):
    """ Reconstruct an average 3D volume. """
    avgr = SUSAN.modules.Averager()
    avgr.list_gpus_ids = params['gpus']
    avgr.threads_per_gpu = params['thr_per_gpu']
    avgr.ctf_correction = params['ctf_corr_avg']
    avgr.rec_halfsets = params['do_halfsets']
    avgr.symmetry = params['symmetry']
    avgr.padding_type = params['padding']
    avgr.reconstruct("average", "input/input_tomos.tomostxt",
                     "input/input_particles.ptclsraw",
                     box_size=params['box_size'])


if __name__ == '__main__':
    if len(sys.argv) > 0:
        inputJson = sys.argv[1]

        if os.path.exists(inputJson):
            with open(inputJson) as fn:
                params = json.load(fn)

            createTomosFile(params)
            if not params['continue']:
                createPtclsFile(params, n_refs=1)
            reconstructAvg(params)
        else:
            raise FileNotFoundError(inputJson)
    else:
        print("Usage: %s jsonParamsFile" % os.path.basename(sys.argv[0]))
