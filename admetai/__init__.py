# **************************************************************************
# *
# * Authors:  Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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

from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from .constants import *
from pwchem.constants import RDKIT_DIC

_references = ['']


class Plugin(pwchemPlugin):
    @classmethod
    def defineBinaries(cls, env):
        cls.addAdmetAiPackage(env)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(ADMETAI_DIC['home'], cls.getEnvName(ADMETAI_DIC))

    @classmethod
    def addAdmetAiPackage(cls, env, default=True):
        installer = InstallHelper(
            ADMETAI_DIC['name'],
            packageHome=cls.getVar(ADMETAI_DIC['home']),
            packageVersion=ADMETAI_DIC['version']
        )

        installer.getCondaEnvCommand(
            ADMETAI_DIC['name'],
            binaryVersion=ADMETAI_DIC['version'],
            pythonVersion='3.14'
        ).addCommand(
            f"{cls.getEnvActivationCommand(ADMETAI_DIC)} && "
            "git clone https://github.com/swansonk14/admet_ai.git &&"
            "cd admet_ai &&"
            "pip install -e .",
            "ADMET-AI_downloaded"
        )

        installer.addPackage(
            env,
            dependencies=['git', 'pip'],
            default=default
        )






