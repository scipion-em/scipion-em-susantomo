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


def createSubsets(params):
    parts = SUSAN.data.Particles(filename=params['input_parts'])
    results = [parts]

    if params['select_refs']:
        for ref in params['refs_list']:
            res = SUSAN.data.Particles.MRA.select_ref(results[-1], ref-1)
            results.append(res)

    if params['do_thr_cc']:
        nrefs_remaining = results[-1].n_refs
        for i in range(nrefs_remaining):
            ind = np.logical_and(results[-1].ali_cc[i] > params['cc_min'],
                                 results[-1].ali_cc[i] < params['cc_max'])
            # check if all values are False
            if np.all(ind == False):
                raise Exception("CC limits are too strict, no particles match.")
            res = results[-1].select(ind)
            results.append(res)

    print("Remaining particles: ", results[-1].n_ptcl)
    results[-1].save('particles.ptclsraw')
    del results


if __name__ == '__main__':
    if len(sys.argv) > 0:
        inputJson = sys.argv[1]

        if os.path.exists(inputJson):
            with open(inputJson) as fn:
                params = json.load(fn)
            createSubsets(params)
        else:
            raise FileNotFoundError(inputJson)
    else:
        print("Usage: %s jsonParamsFile" % os.path.basename(sys.argv[0]))
