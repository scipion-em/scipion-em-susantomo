# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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

import pyworkflow.object as pwobj
from pwem.objects import EMObject


class TomoSubStacks(EMObject):

    def __init__(self, filename=None, n_ptcl=1, n_refs=1,
                 **kwargs):
        EMObject.__init__(self, **kwargs)

        self.filename = pwobj.String(filename)
        self.n_ptcl = pwobj.Integer(n_ptcl)
        self.n_refs = pwobj.Integer(n_refs)

    def __str__(self):
        return f"TomoSubStacks ({self.n_ptcl} items, {self.n_refs} classes)"

    def getNumRefs(self):
        return self.n_refs.get()

    def getSize(self):
        return self.n_ptcl.get()

    def getFileName(self):
        return self.filename.get()
