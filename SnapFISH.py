import os, argparse
from src.SnapFISH import SnapFISH


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', action = 'store', \
                        required = True, help = 'input directory')
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = True, help = 'output directory')
    parser.add_argument('-a', '--anndir', action = 'store', \
                        required = True, help = 'file list directory')
    parser.add_argument('-p', '--pic', action = 'store', type = int, \
                        required = True, help = '1 to save figure, and 0 to not')
    parser.add_argument('-d', '--dataname', action = 'store', \
                        required = True,  help = 'data name')
    return parser.parse_args()


if __name__ == "__main__":
    args = create_parser()
    sf = SnapFISH(
        coor_path = args.indir, 
        ann_path = args.anndir,
        path = args.outdir, 
        suf = args.dataname, 
        w = True,
        save_pic = bool(args.pic)
    )
    sf.run_SnapFISH()
