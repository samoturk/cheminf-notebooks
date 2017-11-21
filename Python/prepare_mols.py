#!/usr/bin/env python

__author__ = "Samo Turk"
__license__ = "BSD 3-clause"

import sys
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PropertyMol
from rdkit.Chem import AllChem
import yaml
from joblib import Parallel, delayed


class PropertyFilter(object):
    """PropertyFilter class partially parses config and initiates relevant descriptors"""

    def __init__(self, properties):
        self.props = properties
        self.descriptors = self._init_descriptors(self.props)
        print('Following descriptors will be used for filtering:')
        for desc in self.descriptors:
            print('%s: %f - %f' %(desc, self.descriptors[desc][1], self.descriptors[desc][2]))

    def _init_descriptors(self, props):
        all_descs = dict(Descriptors._descList)
        descs = {}
        for prop in props:
            print(prop)
            if prop in all_descs:
                try:
                    s = props[prop].split(' - ')
                    descs[prop] = [prop, float(s[0]), float(s[1])]
                except ValueError:
                    raise ValueError('Range does not have numbers for: %s' % prop)
                except IndexError:
                    raise IndexError('Descriptor range not defined correctly for: %s' % prop)
            else:
                print('%s not recognized as RDKit descriptor.' % prop)
        return descs

    def filter_mol(self, mol):
        """Quickly filters mols, doesn't save any calculated values and moves to the next one as soon as a molecules
        has a property not within desired range  """
        for desc in self.descriptors:
            desc = self.descriptors[desc]
            #print([desc[1], float(desc[0](mol)), desc[2]])
            f = Descriptors.__getattribute__(desc[0])
            if desc[1] <= f(mol) <= desc[2]:
                pass
            else:
                return False
        # If all pass
        return True

def _parallel_filter(mol, prop_filter):
    """Helper function for joblib jobs
    """
    if mol is not None:
        #smiles = Chem.MolToSmiles(mol)
        #mol = Chem.MolFromSmiles(smiles)
        if prop_filter.filter_mol(mol):
            return mol


def _read_smi(file_name):
    while True:
        line = file_name.readline()
        if not line:
            break
        mol = Chem.MolFromSmiles(line.split('\t')[0])
        if mol is not None:
            mol.SetProp('_Name', line.split('\t')[1])
            yield mol


def _get_supplier(file_name):
    # File type detection
    in_split = file_name.split('.')
    if in_split[-1].lower() not in ['sdf', 'smi', 'ism', 'gz']:
        raise ValueError('File extension not supported (sdf, smi, ism, sdf.gz, smi.gz)')
    gzipped = False
    if in_split[-1].lower() == 'gz':
        gzipped = True
        if in_split[-2].lower() not in ['sdf', 'smi', 'ism']:
            raise ValueError('File extension not supported (sdf, smi, ism, sdf.gz, smi.gz)')

    if gzipped:
        import gzip
        if in_split[-2].lower() == 'sdf':
            mols_file = gzip.open(file_name, mode='r')
            suppl = Chem.ForwardSDMolSupplier(mols_file)
        else:
            mols_file = gzip.open(file_name, mode='rt')
            suppl = _read_smi(mols_file)
    else:
        if in_split[-1].lower() == 'sdf':
            suppl = Chem.ForwardSDMolSupplier(file_name)
        else:
            mols_file = open(file_name, mode='rt')
            suppl = _read_smi(mols_file)
    return suppl


def do_filter(args):
    suppl = _get_supplier(args.in_file)
    filters = yaml.load(open(args.filter_file))
    prop_filter = PropertyFilter(filters['properties'])
    result = Parallel(n_jobs=args.n_jobs, verbose=1)(delayed(_parallel_filter)(PropertyMol.PropertyMol(mol), prop_filter) for mol in suppl)
    print(len(result))
    w = Chem.SDWriter(args.out_file)
    for m in result:
        w.write(m)
    w.close()


def arg_parser():
    parser = argparse.ArgumentParser(description="""\
prepare_mols.py
""", formatter_class=argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(
        title="subcommands",
        help="prepare_mols $subcommand --help for details on sub-commands")

    filter = subparsers.add_parser("filter",
                                   description="""\
Filter molecules.
""", formatter_class=argparse.RawDescriptionHelpFormatter)

    filter.add_argument("--in-file", "-i", metavar='FILE', required=True, type=str,
                        help="Input SDF, SMI, ISM (can be gzipped).")
    filter.add_argument("--out-file", "-o", metavar='FILE', required=True, type=str,
                        help="Output SDF.")
    filter.add_argument("--filter-file", "-f", metavar='FILE', required=True, type=str,
                        help="Filter file.")
    filter.add_argument("--n-jobs", "-j", metavar="INT", default=1, type=int,
                        help="Number of CPU cores to use (optional, default: 1).",)
    filter.set_defaults(func=do_filter)

    return parser

def run():
    parser = arg_parser()
    if len(sys.argv) == 1:
        argv = ['-h']
    else:
        argv = sys.argv[1:]
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    run()