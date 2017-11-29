#!/usr/bin/env python

__author__ = "Samo Turk"
__license__ = "BSD 3-clause"

import timeit
from itertools import chain, islice
import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors, PropertyMol, SaltRemover, AllChem
from rdkit.Chem.FilterCatalog import *
import yaml
from joblib import Parallel, delayed


class PropertyFilter(object):
    """PropertyFilter class partially parses config and initiates relevant descriptors"""

    def __init__(self, filters, catalog):
        self.props = filters['properties']
        self.atoms = filters['allowed_atoms']
        self.catalog = catalog
        self.descriptors = self._init_descriptors(self.props)
        print('The following descriptors will be used for filtering:')
        for desc in self.descriptors:
            print('    %s: %0.2f - %0.2f' % (desc, self.descriptors[desc][1], self.descriptors[desc][2]))
        print('The following atoms are allowed: %s' %(', '.join(self.atoms)))

    def _init_descriptors(self, props):
        all_descs = dict(Descriptors._descList)
        descs = {}
        for prop in props:
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

    def filter_props(self, mol):
        """Quickly filters mols, doesn't save any calculated values and moves to the next one as soon as a molecules
        has a property not within desired range  """
        for desc in self.descriptors:
            desc = self.descriptors[desc]
            f = Descriptors.__getattribute__(desc[0])
            if desc[1] <= f(mol) <= desc[2]:
                pass
            else:
                return False
        # If all pass
        return True

    def check_atoms(self, mol):
        """Check if mol only as allowed atoms"""
        atoms = []
        for a in mol.GetAtoms():
            atoms.append(a.GetSymbol())
        for a in set(atoms):
            if a not in self.atoms:
                return False
        # If all passes
        return True

    def struct_filter(self, mol):
        if self.catalog.GetFirstMatch(mol):
            return False
        else:
            return True

    def filter_mol(self, mol):
        # Start with the fastest operation
        if not self.check_atoms(mol):
            return False
        if not self.filter_props(mol):
            return False
        if self.catalog:
            if not self.struct_filter(mol):
                return False
        # If all pass, then return True
        return True


remover = SaltRemover.SaltRemover()


def _prepare_catalog(filters):
    if 'structure_filters' not in filters:  # structure filters section missing
        return None, None
    if len(filters['structure_filters']) == 0:  # no filters listed
        return None, None
    params = FilterCatalogParams()
    catalog_names = []
    for f in filters['structure_filters']:
        if f in FilterCatalogParams.FilterCatalogs.names:
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.names[f])
            catalog_names.append(f)
    catalog = FilterCatalog(params)
    return catalog, catalog_names


def _parallel_filter(mol, prop_filter, filters):
    """Helper function for joblib jobs
    """
    if mol is not None:
        # Remove salts by default except explicitly requested not to
        name = mol.GetProp('_Name')
        if 'remove_salts' in filters:
            if filters['remove_salts'] is False:
                pass
            else:
                mol = remover.StripMol(mol)
        if prop_filter.filter_mol(mol):
            mol = PropertyMol.PropertyMol(mol)
            mol.SetProp('_Name', name)
            return mol


def _read_smi(file_name):
    while True:
        line = file_name.readline()
        if not line:
            break
        line = line.split()
        mol = Chem.MolFromSmiles(line[0])
        if mol is not None:
            if len(line) > 1:
                mol.SetProp('_Name', line[1])
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


def _get_chunks(iterable, size):
    iterator = iter(iterable)
    for first in iterator:
        yield chain([first], islice(iterator, size - 1))

def _parallel_3d(mol):
    """Helper function for joblib jobs
    """
    if mol is not None:
        name = mol.GetProp('_Name')
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        mol = PropertyMol.PropertyMol(mol) # making sure properties are picklable
        mol.SetProp('_Name', name)
        return mol

def do_filter(args):
    suppl = _get_supplier(args.in_file)
    filters = yaml.load(open(args.filter_file))
    catalog, catalog_names = _prepare_catalog(filters)
    prop_filter = PropertyFilter(filters, catalog)
    if catalog_names:
        print('Structure filters: %s' % (', '.join(catalog_names)))

    # Processing in chunks so we don't consume all memory
    chunk_size = args.chunk_size
    chunks = _get_chunks(suppl, chunk_size)
    w = Chem.SDWriter(args.out_file)
    processed = 0
    saved = 0
    start = timeit.default_timer()
    print('Starting with filtering...')
    for chunk in chunks:
        result = Parallel(n_jobs=args.n_jobs, verbose=0)(delayed(_parallel_filter)(PropertyMol.PropertyMol(mol),
                                                                                    prop_filter, filters) for mol in chunk)
        for i, m in enumerate(result):
            if m is not None:
                saved += 1
                w.write(m)
        elapsed = (timeit.default_timer() - start)/60
        processed += len(result)
        print('Processed: %i       | elapsed %0.2f minutes' % (processed, elapsed))
    w.close()
    print('Number of molecules processed %i, and saved: %i in %0.2f minutes' % (processed, saved, elapsed))


def do_3d(args):
    suppl = _get_supplier(args.in_file)
    # Processing in chunks so we don't consume all memory
    chunk_size = args.chunk_size
    chunks = _get_chunks(suppl, chunk_size)
    w = Chem.SDWriter(args.out_file)
    processed = 0
    saved = 0
    start = timeit.default_timer()
    print('Starting with generating conformations...')
    for chunk in chunks:
        result = Parallel(n_jobs=args.n_jobs)(delayed(_parallel_3d)(PropertyMol.PropertyMol(mol)) for mol in chunk)
        for i, m in enumerate(result):
            if m is not None:
                saved += 1
                w.write(m)
        elapsed = (timeit.default_timer() - start) / 60
        processed += len(result)
        print('Processed: %i       | elapsed %0.2f minutes' % (processed, elapsed))
    w.close()
    print('Number of molecules processed %i, and saved: %i in %0.2f minutes' % (processed, saved, elapsed))



def arg_parser():
    parser = argparse.ArgumentParser(description="""\
prepare_mols.py

Filter molecules based on allowed atoms, properties, and structural alerts;
TODO: Generate 3D confs
""", formatter_class=argparse.RawDescriptionHelpFormatter)

    subparsers = parser.add_subparsers(
        title="subcommands",
        help="prepare_mols $subcommand --help for details on sub-commands")

    filter_ = subparsers.add_parser("filter",
                                    description="""\
Filter molecules based on allowed atoms, properties, and structural alerts.

Note: Removes salts by default unless "remove_salts: False" in filter.yaml

""", formatter_class=argparse.RawDescriptionHelpFormatter)

    filter_.add_argument("--in-file", "-i", metavar='FILE', required=True, type=str,
                         help="Input SDF, SMI, ISM (can be gzipped).")
    filter_.add_argument("--out-file", "-o", metavar='FILE', required=True, type=str,
                         help="Output SDF.")
    filter_.add_argument("--filter-file", "-f", metavar='FILE', required=True, type=str,
                         help="Filter file.")
    filter_.add_argument("--n-jobs", "-j", metavar="INT", default=1, type=int,
                         help="Number of CPU cores to use (optional, default: 1).",)
    filter_.add_argument("--chunk-size", "-c", metavar="INT", default=50000, type=int,
                         help="Number of molecules to be held in RAM. 50k is a good balance between \
                         speed and RAM usage (~1GB). Note that the progress is reported only when a \
                         a chunk is processed (optional, default: 50000).",)
    filter_.set_defaults(func=do_filter)

    gen3d = subparsers.add_parser("gen3d",
                                    description="""\
    Generate 3D conformations using ETKDG method.

    """, formatter_class=argparse.RawDescriptionHelpFormatter)

    gen3d.add_argument("--in-file", "-i", metavar='FILE', required=True, type=str,
                         help="Input SDF, SMI, ISM (can be gzipped).")
    gen3d.add_argument("--out-file", "-o", metavar='FILE', required=True, type=str,
                         help="Output SDF.")
    gen3d.add_argument("--n-jobs", "-j", metavar="INT", default=1, type=int,
                         help="Number of CPU cores to use (optional, default: 1).", )
    gen3d.add_argument("--chunk-size", "-c", metavar="INT", default=1000, type=int,
                         help="Number of molecules to be held in RAM. Note that the progress is \
                         reported only when a chunk is processed (optional, default: 1000).", )
    gen3d.set_defaults(func=do_3d)

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
