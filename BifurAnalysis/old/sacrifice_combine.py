from ROOT import TFile, TH1F
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile', type=str, default='output.root')
    parser.add_argument('-f', '--file', type=str)
    return parser.parse_args()

def main():
    args = get_args()
    tf = TFile(args.file)
    outfile = TFile(args.outfile, "recreate")
    keylist = [k.GetTitle() for k in tf.GetListOfKeys()]
    # Divide by denom_<>
    denom_keys = [d for d in keylist if 'denom' in d]
    other_keys = [d for d in keylist if 'denom' not in d]
    params = [p.split('_')[1] for p in denom_keys]
    for par in params:
        denom_hist = tf.Get([dk for dk in denom_keys if par in dk][0])
        for other in other_keys:
            if par in other:
                hist = tf.Get(other)
                hist.Divide(denom_hist)
                hist.Write()
        denom_hist.Write()
    outfile.Write()

if __name__ == '__main__':
    main()
