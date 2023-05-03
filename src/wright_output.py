from typing import List, Optional

import pandas as pd
from tqdm import tqdm

from src.types import PichiaFullData


class GenerateOutput:
    def __init__(self, result: List[PichiaFullData], output: Optional[str] = None):
        self.result = result
        self.output = output

    def pandas_save(self) -> pd.DataFrame:
        data = pd.DataFrame(
            {
                "Pichia gene name": [],
                "Pichia gene ID": [],
                "Pichia description": [],
                "Pichia protein ID": [],
                "Saccharomyces accession ID": [],
                "Saccharomyces protein name": [],
                "Query Cover": [],
                "E-value": [],
                "Per. Ident": [],
                "Acc. Len": [],
                "SGD description": [],
            }
        )

        for gene in tqdm(self.result):
            row = pd.DataFrame(
                {
                    "Pichia gene name": [gene.pp_gene_name],
                    "Pichia gene ID": [gene.pp_gene_id],
                    "Pichia description": [gene.pp_description],
                    "Pichia protein ID": [gene.pp_protein_id],
                    "Saccharomyces accession ID": [gene.sc_accession],
                    "Saccharomyces protein name": [gene.sc_description],
                    "Query Cover": [gene.query_cover],
                    "E-value": [gene.e_value],
                    "Per. Ident": [gene.per_ident],
                    "Acc. Len": [gene.acc_le],
                    "SGD description": [gene.sgd_description],
                }
            )

            data = pd.concat([data, row], ignore_index=True, axis=0)

        return data

    def print_output(self):
        print("Result:")
        print(self.result)

    def csv_output(self):
        print("Generating output table in .csv")
        data = self.pandas_save()
        print("Saving table")
        data.to_csv(self.output, sep=",")

    def tsv_output(self):
        print("Generating output table in .tsv")
        data = self.pandas_save()
        print("Saving table")

        data.to_csv(self.output, sep="\t")

    def xlsx_output(self):
        print(self.result)

    def txt_output(self):
        print(self.result)

    def run(self):
        if self.output is None:
            self.print_output()
            return

        if self.output.endswith("csv"):
            self.csv_output()
        elif self.output.endswith("tsv"):
            self.tsv_output()
        elif self.output.endswith("xlsx"):
            self.xlsx_output()
        elif self.output.endswith("txt"):
            self.txt_output()
        else:
            self.print_output()
