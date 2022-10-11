import pandas as pd
from tqdm import tqdm


class GenerateOutput:
    def __init__(self, result, output: str = None):
        self.result = result
        self.output = output

    def pandas_save(self) -> pd.DataFrame:
        data = pd.DataFrame(
            {
                "Pichia gene ID": [],
                "Pichia description": [],
                "Pichia protein ID": [],
                "Saccharomyces accession ID": [],
                "Saccharomyces description": [],
                "Query Cover": [],
                "E-value": [],
                "Per. Ident": [],
                "Acc. Len": [],
            }
        )

        for gene in tqdm(self.result):
            row = pd.DataFrame(
                {
                    "Pichia gene ID": [gene],
                    "Pichia description": [self.result[gene][1]],
                    "Pichia protein ID": [self.result[gene][2]],
                    "Saccharomyces accession ID": [self.result[gene][3]],
                    "Saccharomyces description": [self.result[gene][4]],
                    "Query Cover": [self.result[gene][5]],
                    "E-value": [self.result[gene][6]],
                    "Per. Ident": [self.result[gene][7]],
                    "Acc. Len": [self.result[gene][8]],
                }
            )

            data = pd.concat([data, row], ignore_index=True, axis=0)

        return data

    def print_output(self):
        print("Result:")
        print(self.result)

    def csv_output(self):
        print("Generating output table")
        data = self.pandas_save()
        print("Saving table")
        data.to_csv(self.output, sep=",")

    def tsv_output(self):
        print("Generating output table")
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
