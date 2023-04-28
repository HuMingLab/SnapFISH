import click
from src.SnapFISH import SnapFISH


@click.group()
def main():
    pass


@main.command()
@click.option('-i', '--indir', required = True, help = '3D coordinates file path')
@click.option('-o', '--outdir', required = True, help = 'output directory')
@click.option('-a', '--anndir', required = True, help = 'annotation file path')
@click.option('-p', '--pic', required = True, help = '0 or 1, 1 to save figures')
@click.option('-d', '--dataname', required = True,  help = 'data name')
def call_loops(indir, outdir, anndir, pic, dataname):
    """
    Main method to call loops from a given set of files.
    """
    sf = SnapFISH(
        coor_path = indir,
        ann_path = anndir,
        path = outdir, 
        suf = dataname,
        save_pic = bool(int(pic))
    )
    sf.run_SnapFISH()


@main.command()
@click.option('-i', '--indir', required = True, help = '3D coordinates file path')
@click.option('-o', '--outdir', required = True, help = 'output directory')
@click.option('-a', '--anndir', required = True, help = 'annotation file path')
@click.option('-p', '--pic', required = True, help = '0 or 1, 1 to save figures')
@click.option('-d', '--dataname', required = True,  help = 'data name')
def step1(indir, outdir, anndir, pic, dataname):
    """
    Only step 1 of SnapFISH, calculate pairwise distances.
    """
    sf = sf = SnapFISH(
        coor_path = indir,
        ann_path = anndir,
        path = outdir, 
        suf = dataname,
        save_pic = bool(int(pic))
    )
    sf.SnapFISH_step1()


@main.command()
@click.option('-i', '--testdir', required = True, help = 'T test and Wilcoxon test file path')
@click.option('-o', '--outdir', required = True, help = 'output directory')
@click.option('-a', '--anndir', required = True, help = 'annotation file path')
@click.option('-p', '--pic', required = True, help = '0 or 1, 1 to save figures')
@click.option('-d', '--dataname', required = True,  help = 'data name')
@click.option('-m', '--testmethod', default = "t",  
              help = 't to use t-test, w to use Wilcoxon rank test')
@click.option('-h', '--threshold', default = "0.1", 
              help = 'fdr threshold to select loop candidates')
def test(testdir, outdir, anndir, pic, dataname, testmethod, threshold):
    """
    Read output from step 2, and identify loop summits based on the cutoff and 
    the test method defined by the user.
    """
    sf = SnapFISH(
        coor_path = testdir,
        ann_path = anndir,
        path = outdir, 
        suf = dataname,
        save_pic = bool(int(pic))
    )
    sf.out_all_tests = sf.data
    sf.SnapFISH_step3(
        test_method = testmethod,
        fdr_cut = float(threshold)
    )
    sf.SnapFISH_step4()
