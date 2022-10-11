import pandas as pd


class ParseInput:
    def __init__(self, input_file: str = None):
        self.input_file = input_file

    def pandas_parser(self, separator):
        data = pd.read_csv(self.input_file, sep=separator)

        if "id" not in data.columns:
            if "gene" not in data.columns:
                raise ValueError(
                    "В табличке не вижу столбца с генами. Должен быть либо 'id', либо 'gene'"
                )
            else:
                pichia_genes = data["gene"]
        else:
            pichia_genes = data["id"]

        if pichia_genes.apply(lambda x: not x.startswith("PAS")).any():
            raise ValueError("Гены должны начинаться на PAS")

        return pichia_genes.to_list()

    def csv_input(self) -> list:
        return self.pandas_parser(",")

    def tsv_input(self) -> list:
        return self.pandas_parser("\t")

    def run(self):
        if self.input_file is None:
            raise ValueError("Нужны входные данные")

        if self.input_file.endswith("csv"):
            return self.csv_input()

        if self.input_file.endswith("tsv"):
            return self.tsv_input()

        raise ValueError("Не умею парсить этот формат файла")
