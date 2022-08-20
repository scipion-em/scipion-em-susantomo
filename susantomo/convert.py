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

import os
import numpy as np
import logging
logger = logging.getLogger(__name__)


def parseImodCtf(filename):
    """ Retrieve defocus U, V and angle from the
    output file of the susan execution.
    :param filename: input file to parse
    :return: an array with CTF values
    """
    if os.path.exists(filename):
        return np.loadtxt(filename, dtype=float, comments='#')
    else:
        logger.error(f"Warning: Missing file: {filename}")

    return None


def setWrongDefocus(ctfModel):
    ctfModel.setDefocusU(-999)
    ctfModel.setDefocusV(-1)
    ctfModel.setDefocusAngle(-999)


def readCtfModelStack(ctfModel, ctfArray, item=0):
    """ Set values for the ctfModel from an input list.
    :param ctfModel: output CTF model
    :param ctfArray: array with CTF values
    :param item: which row to use from ctfArray
    """
    if ctfArray is None or np.isnan(ctfArray[item]).any(axis=0):
        setWrongDefocus(ctfModel)
        ctfFit, ctfResolution, ctfPhaseShift = -999, -999, 0
    else:
        # 8-column list
        (defocusU, defocusV, defocusAngle, ctfPhaseShift,
         _, _, ctfResolution, ctfFit) = ctfArray[item]
        ctfModel.setStandardDefocus(defocusU, defocusV, defocusAngle)
    ctfModel.setFitQuality(ctfFit)
    ctfModel.setResolution(ctfResolution)

    # Avoid creation of phaseShift
    if ctfPhaseShift != 0:
        ctfModel.setPhaseShift(ctfPhaseShift)
