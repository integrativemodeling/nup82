import os
import csv
import re
import glob
import ihm.location
import ihm.dataset

saxs_dir = '../SAXS'

class SAXSFits(object):
    """Parse the SAXS csv file and add suitable fit data to the mmCIF file"""
    seqrange_re = re.compile('(\d+)\s*\-\s*(\d+)')

    def __init__(self, po):
        self.po = po

    def add_from_csv(self, model):
        with open(os.path.join(saxs_dir, 'Table6_SAXS.csv'), 'rb') as fh:
            for row in csv.DictReader(fh):
                # 24367 was not used in publication due to low concentration
                if row['FoXS fit score (chi value)'] \
                   and row['Protein'] == 'NUP82' \
                   and row['Protein ID'] != '24367':
                    yield self._add_one(row, model)

    def _add_one(self, row, model):
        protid = row['Protein ID']
        # todo: link to all sub files, not just first?
        profile = glob.glob('%s/%s_*/*.sub' % (saxs_dir, protid))[0]
        l = ihm.location.InputFileLocation(profile,
                             details = row['Notes'] if row['Notes']
                                                    else None)
        dataset = ihm.dataset.SASDataset(location=l)
        m = self.seqrange_re.match(row['Sequence coverage'])
        seqrange = (int(m.group(1)), int(m.group(2)))
        # map seqrange for other species to that for yeast (manually)
        # 9-445 (Candida glabrata) -> 4-452 (yeast)
        # 791-899 (S. pombe) -> 572-690 (yeast)
        if seqrange == (9,445):
            seqrange = (4,452)
        elif seqrange == (791,899):
            seqrange = (572,690)
        # Note that fit is for the 2nd copy of Nup82
        self.po._add_foxs_restraint(model, 'Nup82.2', seqrange, dataset,
                         row['Rg'], row['FoXS fit score (chi value)'], None)
        return dataset
