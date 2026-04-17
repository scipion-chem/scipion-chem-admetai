# **************************************************************************
# *
# * Authors:   Blanca Pueche (blanca.pueche@cnb.csis.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import os, csv
import pandas as pd

import pyworkflow.protocol.params as params
from admetai import ADMETAI_DIC
from pwem.protocols import EMProtocol
from pyworkflow.object import String, Float

from pwchem import Plugin
from pwchem.objects import  SetOfStructROIs, StructROI
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import insistentRun
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC

RDKIT, OPENBABEL = 0, 1



class ProtAdmetAi(EMProtocol):
    """


    """
    _label = 'ADMET toxicity prediction'
    # -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addHidden('useGpu', params.BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation. Choose one.")

        form.addHidden('gpuList', params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="Comma-separated GPU devices that can be used.")

        form.addSection(label='Input')
        form.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                        label='Input small molecules: ',
                        help='Set of small molecules to input the model for predicting their interactions')

        form.addParam('physchem', params.BooleanParam, label='Include physicochemical properties: ', default=True,
                      help='Whether to include physicochemical properties in the predictions.')
        form.addParam('atcCode', params.StringParam, default='', expertLevel=params.LEVEL_ADVANCED,
                       label="ATC code to filter the DrugBank reference: ",
                       help="The ATC code to filter the DrugBank reference set by. If None, the entire DrugBank reference set will be used.")


        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertSMIStep)
        self._insertFunctionStep(self.writeSMILEScsv)
        self._insertFunctionStep(self.runAdmetAiStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertSMIStep(self):
        smiDir = self.getInputSMIDir()
        if not os.path.exists(smiDir):
            os.makedirs(smiDir)

        molDir = self.copyInputMolsInDir()
        args = ' --multiFiles -iD "{}" --pattern "{}" -of smi --outputDir "{}"'. \
            format(molDir, '*', smiDir)
        Plugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=smiDir)

    def writeSMILEScsv(self):
        smiDic = self.getSMIdic()
        outFile = self._getPath("smiles.csv")

        with open(outFile, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['smiles'])

            for smi in smiDic.keys():
                writer.writerow([smi])

        return outFile

    def runAdmetAiStep(self):
        filePath = self._getPath("smiles.csv")
        outPath = self._getPath("admetai_predictions.csv")
        args = ['--data_path', os.path.abspath(filePath),
                '--save_path', os.path.abspath(outPath)]

        if self.physchem.get():
            args.append("--include_physchem")

        if self.atcCode.get():
            args.extend(['--atc_code', self.atcCode.get()])

        args.extend(['--num_workers', str(self.numberOfThreads.get())])

        if self.useGpu.get():
            fullCommand = (
                    f"export CUDA_VISIBLE_DEVICES={self.gpuList.get()} && "
                    f"admet_predict " + " ".join(args)
            )
        else:
            fullCommand = (
                    f"export CUDA_VISIBLE_DEVICES='' && "
                    f"admet_predict " + " ".join(args)
            )

        Plugin.runCondaCommand(
            self,
            args=" ",
            condaDic=ADMETAI_DIC,
            program=fullCommand,
            cwd=self._getPath()
        )

    def createOutputStep(self):
        df = self.parseAdmetResults()
        smiDic = self.getSMIdic()

        outMols = self.inputSmallMolecules.get().createCopy(self._getPath(), copyInfo=True)

        for mol in self.inputSmallMolecules.get():
            molName = mol.getMolName()
            smi = None
            for s, name in smiDic.items():
                if name == molName:
                    smi = s
                    break

            if smi is None:
                continue

            if smi not in df['smiles'].values:
                continue

            row = df[df['smiles'] == smi].iloc[0]

            newMol = mol.clone()

            setattr(newMol, 'AMES', Float(row['AMES']))
            setattr(newMol, 'DILI', Float(row['DILI']))
            setattr(newMol, 'hERG', Float(row['hERG']))
            setattr(newMol, 'HIA_Hou', Float(row['HIA_Hou']))
            setattr(newMol, 'CYP3A4_Veith', Float(row['CYP3A4_Veith']))

            outMols.append(newMol)

        self._defineOutputs(outputSmallMolecules=outMols)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = [f"Full results written in {self._getPath('admetai_predictions.csv')}"]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getInputSMIDir(self):
        return os.path.abspath(self._getExtraPath('inputSMI'))

    def copyInputMolsInDir(self):
        oDir = os.path.abspath(self._getTmpPath('inMols'))
        if not os.path.exists(oDir):
            os.makedirs(oDir)
        for mol in self.inputSmallMolecules.get():
            os.link(mol.getFileName(), os.path.join(oDir, os.path.split(mol.getFileName())[-1]))
        return oDir

    def getSMIdic(self):
        '''Returns a dictionary: {SMILES: molName}
        '''
        smiDic = {}

        for smiFile in os.listdir(self.getInputSMIDir()):
            with open(os.path.join(self.getInputSMIDir(), smiFile)) as f:
                for line in f:
                    smi, molName = line.strip().split()
                    smiDic[smi] = molName
        return smiDic

    def parseAdmetResults(self):
        filePath = self._getPath("admetai_predictions.csv")
        df = pd.read_csv(filePath)
        return df