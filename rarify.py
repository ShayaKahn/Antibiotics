import pandas as pd
import numpy as np

class Rarify:
    def __init__(self, df, depth=None):
        """
        :param df: dataframe, represents OTU table with samples as columns and OTUs as rows,
                   the values are the counts (integers)
        :param depth: int, the number of reads to be sampled from each sample, if None,
                      the minimum number of reads is the default.
        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError("The input object is not a DataFrame.")
        if not df.applymap(lambda x: isinstance(x, int)).all().all():
            raise ValueError("Not all values in the DataFrame are integers.")
        if depth is not None:
            if not isinstance(depth, int):
                raise TypeError("Depth must be an integer.")
            if not (df.sum(axis=0).min() < depth < df.sum(axis=0).max()):
                raise ValueError(f"Depth must be between {df.sum(axis=0).min()} and {df.sum(axis=0).max()}.")
        self.df = df
        if depth is None:
            self.depth = self.df.sum(axis=0).min()
            self.df = df
        else:
            self.depth = depth
            self.df = self.filter_low_frequancy()
        self.sampling_dict = self.construct_sampling_dict()

    def filter_low_frequancy(self):
        # filter out samples with lower frequancy then the depth
        df_filtered = self.df.loc[:, (self.df.sum(axis=0) > self.depth)]
        return df_filtered

    def construct_sampling_dict(self):
        # construct a sampling dictionary
        sampling_dict = {}
        for col in self.df.columns:
            sampling_dict[col] = [idx for idx, value in enumerate(self.df[col]) for _ in range(value)]
        return sampling_dict

    def rarify(self):
        # rarify the dataframe
        rar_df = pd.DataFrame(0, index=self.df.index, columns=self.df.columns)

        for col in self.df.columns:
            selected_indices = np.random.choice(self.sampling_dict[col], self.depth, replace=False)
            for idx in selected_indices:
                rar_df.at[idx, col] += 1

        return rar_df