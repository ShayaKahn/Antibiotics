import pandas as pd

class  CollapseTaxonomy:
    """
    Collapse OTU table to desired taxonomy level.
    """
    def __init__(self, otu_table, taxonomy_table, level):
        """
        :param otu_table: OTU table, pandas dataframe of shape (#samples, #OTUs) where the names of the OTUs correspond
         to the names in taxonomy_table.
        :param taxonomy_table: taxonomy table, pandas dataframe where the first column correspond to OTU ids and the
         second column correspond to the taxonomy. The taxonomy is a string with the format as shown in this example:
         d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales;
        :param level: desired taxonomy level, represented as a lowercase letter correspond to the first letter of the
         desired taxonomy.
        """
        self.otu_table = otu_table
        self.taxonomy_table = taxonomy_table
        self.level = level


    def _create_desired_taxonomy_dict(self):
        pass

