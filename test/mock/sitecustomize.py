import ihm
import os

# If we're running from an SGE job, override the from_pubmed_id() function
# to return a cached value, since we don't have network access (needed to
# query PubMed directly)

def mock_from_pubmed(cls, pubmed_id):
    return ihm.Citation(
            pmid=27839866,
            title='Structure and Function of the Nuclear Pore Complex '
                  'Cytoplasmic mRNA Export Platform.', journal='Cell',
            volume=167, page_range=('1215','1228.e25'), year=2016, authors=[
                   'Fernandez-Martinez J', 'Kim SJ', 'Shi Y', 'Upla P',
                   'Pellarin R', 'Gagnon M', 'Chemmama IE', 'Wang J',
                   'Nudelman I', 'Zhang W', 'Williams R', 'Rice WJ',
                   'Stokes DL', 'Zenklusen D', 'Chait BT', 'Sali A', 
                   'Rout MP'], doi='10.1016/j.cell.2016.10.028')

if 'JOB_ID' in os.environ:
    ihm.Citation.from_pubmed_id = classmethod(mock_from_pubmed)
