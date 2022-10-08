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

import pwem
import pyworkflow.utils as pwutils

from .constants import *


__version__ = '3.0b1'
_references = ['Sanchez2019', 'Sanchez2019b']
_logo = "susan_logo.png"


class Plugin(pwem.Plugin):
    _homeVar = SUSAN_HOME
    _pathVars = [SUSAN_HOME]
    _url = "https://github.com/rkms86/SUSAN"
    _supportedVersions = [V0_1]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(SUSAN_HOME, 'susan-%s' % V0_1)
        cls._defineVar(SUSAN_CUDA_LIB, pwem.Config.CUDA_LIB)
        cls._defineVar(SUSAN_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)

    @classmethod
    def getSusanEnvActivation(cls):
        """ Remove the scipion home and activate the conda environment. """
        activation = cls.getVar(SUSAN_ENV_ACTIVATION)
        scipionHome = pwem.Config.SCIPION_HOME + os.path.sep

        return activation.replace(scipionHome, "", 1)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch SUSAN. """
        environ = pwutils.Environ(os.environ)
        binPath = pwem.Config.MPI_BINDIR
        libPath = pwem.Config.MPI_LIBDIR

        if binPath not in environ['PATH']:
            environ.update({'PATH': binPath,
                            'LD_LIBRARY_PATH': libPath
                            }, position=pwutils.Environ.BEGIN)

        # Get SUSAN CUDA library path if defined
        cudaLib = cls.getVar(SUSAN_CUDA_LIB, pwem.Config.CUDA_LIB)
        environ.addLibrary(cudaLib)

        if 'SUSAN_MPI_LIB' in os.environ:
            environ.addLibrary(os.environ['SUSAN_MPI_LIB'])

        if 'SUSAN_MPI_BIN' in os.environ:
            environ.set('PATH', os.environ['SUSAN_MPI_BIN'],
                        position=pwutils.Environ.BEGIN)

        if 'PYTHONPATH' in environ:
            del environ['PYTHONPATH']

        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
        activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = ['git', 'gcc', 'cmake', 'make']
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def defineBinaries(cls, env):
        for ver in cls._supportedVersions:
            ENV_NAME = "susan-%s" % ver
            installCmds = [
                cls.getCondaActivationCmd(),
                f'conda create -y -n {ENV_NAME} python=3 && ',
                f'conda activate {ENV_NAME} && ',
                f'cd .. && rmdir susan-{ver} && ',
                f'git clone https://github.com/rkms86/SUSAN {ENV_NAME} && ',
                f'cd {ENV_NAME}/extern && ',
                f'git clone https://gitlab.com/libeigen/eigen.git eigen && ',
                f'cd eigen && mkdir build && cd build && '
                f'cmake ../ -DCMAKE_INSTALL_PREFIX=../../eigen_lib && make install && ',
                f'cd ../../../ && mkdir bin && cd bin && cmake .. && ',
                f'make -j {env.getProcessors()} && make prepare_python && ',
                f'cd .. && pip install -e .'
            ]

            susanCmds = [(" ".join(installCmds), 'bin/susan_aligner_mpi')]

            env.addPackage('susan', version=ver,
                           tar='void.tgz',
                           commands=susanCmds,
                           neededProgs=cls.getDependencies(),
                           updateCuda=True,
                           default=ver == V0_1)

    @classmethod
    def getActivationCmd(cls):
        """ Return the activation command. """
        return '%s %s' % (cls.getCondaActivationCmd(),
                          cls.getSusanEnvActivation())

    @classmethod
    def getActiveVersion(cls, *args):
        """ Return the env name that is currently active. """
        envVar = cls.getVar(SUSAN_ENV_ACTIVATION)
        return envVar.split()[-1].split("-")[-1]

    @classmethod
    def getProgram(cls, script):
        scriptFn = os.path.join(__path__[0], f'scripts/{script}')
        cmd = f"{cls.getActivationCmd()} && python3 {scriptFn} "

        return cmd
