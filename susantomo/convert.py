# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es) [2]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# * [2] Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
import math
import logging
logger = logging.getLogger(__name__)

import tomo.constants as const


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


def writeDynTable(fn, setOfSubtomograms, angleMin=0, angleMax=0,
                  scaleCoords=1.0, scaleShifts=1.0):
    """ Write a Dynamo-style tbl from a set of subtomograms. """
    first = setOfSubtomograms.getFirstItem()
    hasCoord = first.hasCoordinate3D()
    hasTransform = first.hasTransform()

    for subtomo in setOfSubtomograms.iterSubtomos():
        if hasCoord:
            coord = subtomo.getCoordinate3D()
            x, y, z = coord.getPosition(const.BOTTOM_LEFT_CORNER)
            x, y, z = x*scaleCoords, y*scaleCoords, z*scaleCoords
            tomo_id = coord.getVolId()
        else:
            x = 0
            y = 0
            z = 0
            tomo_id = 0
        if hasTransform:
            tdrot, tilt, narot, shiftx, shifty, shiftz = matrix2eulerAngles(subtomo.getTransform().getMatrix())
            shiftx *= scaleShifts
            shifty *= scaleShifts
            shiftz *= scaleShifts
        else:
            tilt = 0
            narot = 0
            tdrot = 0
            shiftx = 0
            shifty = 0
            shiftz = 0
        fn.write(f'{subtomo.getObjId()} 1 1 {shiftx} {shifty} {shiftz} '
                 f'{tdrot} {tilt} {narot} 0 0 0 1 {angleMin} {angleMax} '
                 f'0 0 0 0 {tomo_id} 0 1 0 {x} {y} {z} 0 0 0 0 '
                 f'0 1 0 0 0 0 0 0 0 0\n')


# matrix2euler dynamo
def matrix2eulerAngles(A):
    tol = 1e-4
    if abs(A[2, 2] - 1) < tol:
        tilt = 0
        narot = math.atan2(A[1, 0], A[0, 0]) * 180 / math.pi
        tdrot = 0
    elif abs(A[2, 2] + 1) < tol:
        tdrot = 0
        tilt = 180
        narot = math.atan2(A[1, 0], A[0, 0]) * 180 / math.pi
    else:
        tdrot = math.atan2(A[2, 0], A[2, 1])
        tilt = math.acos(A[2, 2])
        narot = math.atan2(A[0, 2], -A[1, 2])
    tilt = tilt * 180 / math.pi
    narot = narot * 180 / math.pi
    tdrot = tdrot * 180 / math.pi

    return tdrot, tilt, narot, A[0, 3], A[1, 3], A[2, 3]
