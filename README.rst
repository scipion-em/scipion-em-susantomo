============
SUSAN plugin
============

This plugin provides wrappers for `Susan <https://github.com/rkms86/SUSAN>`_ software: subtomogram averaging (StA) workflow for CryoET based on sub-stacks of images instead of sub-volumes of tomograms.
SUSAN uses substacks that are cropped on-the-fly from the aligned tilt-series stack and performs the CTF correction according to the selected operation (alignment or reconstruction). With this approach we decrease the
computational resources needed in a typical subtomogram averaging pipeline, as we no longer need the CTF corrected stacks, the full tomogram reconstructions, or all the subtomograms in multiple binning
stages.

.. image:: https://img.shields.io/pypi/v/scipion-em-susantomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-susantomo
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-susantomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-susantomo
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-susantomo.svg
        :target: https://pypi.python.org/pypi/scipion-em-susantomo
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-susantomo?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-susantomo
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-susantomo
        :target: https://pypi.python.org/pypi/scipion-em-susantomo
        :alt: Downloads

Installation
-------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-susantomo

b) Developer's version

   * download repository

    .. code-block::

        git clone -b devel https://github.com/scipion-em/scipion-em-susantomo.git

   * install

    .. code-block::

       scipion installp -p /path/to/scipion-em-susantomo --devel

- SUSAN sources will be downloaded and compiled automatically with the plugin, but you can also link an existing installation. Default installation path assumed is ``software/em/susan-0.1``, if you want to change it, set *SUSAN_HOME* in ``scipion.conf`` file to the folder where the SUSAN is installed.
- SUSAN installation requires CUDA libraries, gcc >= 9.x and cmake 3.22.x. OpenMPI is optional amd only required for running SUSAN on more than one cluster node. MATLAB is not required as this plugin uses Python API.
- If you need to use CUDA different from the one used during Scipion installation (defined by *CUDA_LIB*), you can add *SUSAN_CUDA_LIB* variable to the config file.
- If you have to use a MPI for SUSAN different from Scipion MPI, you can set *SUSAN_MPI_BIN* and *SUSAN_MPI_LIB* variables in the config file. At the moment MPI support in the plugin is not working.


Verifying
---------

To check the installation, simply run one of the tests. A complete list of tests can be displayed by executing ``scipion test --show --grep susantomo``

Supported versions
------------------

0.1

Protocols
----------

* ctf estimation
* multi-reference alignment
* average and reconstruct
* create a subset

References
-----------

1. Sánchez RM, Mester R, Kudryashev M. Fast Cross Correlation for Limited Angle Tomographic Data. In: Felsberg M., Forssén PE., Sintorn IM., Unger J. (eds) Image Analysis. SCIA 2019. Lecture Notes in Computer Science, vol 11482. Doi: 10.1007/978-3-030-20205-7_34
2. R.M. Sánchez, R. Mester and M. Kudryashev. Fast Alignment of Limited Angle Tomograms by projected Cross Correlation. 2019 27th European Signal Processing Conference (EUSIPCO), 2019, pp. 1-5, doi: 10.23919/EUSIPCO.2019.8903041
