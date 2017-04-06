#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import csv
import sys
import json
import logging
import collections as coll
from enum import Enum
from codecs import open

import attr
from attr.validators import instance_of as is_a
import click
from obonet.read import read_obo
import networkx as nx

logger = logging.getLogger(__name__)


INSDCTypes = Enum(
    'ISNDCTypes', [
        "RNase_MRP_RNA",
        "RNase_P_RNA",
        "SRP_RNA",
        "Y_RNA",
        "antisense_RNA",
        "autocatalytically_spliced_intron",
        "guide_RNA",
        "hammerhead_ribozyme",
        "lncRNA",
        "miRNA",
        "ncRNA",
        "misc_RNA",
        "other",
        "precursor_RNA",
        "piRNA",
        "rasiRNA",
        "ribozyme",
        "scRNA",
        "siRNA",
        "snRNA",
        "snoRNA",
        "telomerase_RNA",
        "vault_RNA",
        "rRNA",
        "tRNA",
    ])


INSDCTree = nx.Graph()
INSDCTree.add_edges_from([
    ('ncRNA', "ribozyme"),
    ('ribozyme', 'hammerhead_ribozyme'),
    ('ribozyme', 'autocatalytically_spliced_intron'),

    ('ncRNA', 'tRNA'),
    ('ncRNA', 'rRNA'),
    ('ncRNA', 'Y_RNA'),
])


def rna_type_to_key(rna_type):
    return tuple(rna_type.rstrip(';').split('; '))


@attr.s(frozen=True)
class RfamFamily(object):
    id = attr.ib(validator=is_a(str))
    name = attr.ib(validator=is_a(str))
    so_terms = attr.ib(validator=is_a(set))
    rna_type = attr.ib(validator=is_a(tuple))

    @classmethod
    def build_all(cls, link_file, family_file):
        so_terms = coll.defaultdict(set)
        with open(link_file, 'r', 'utf-8') as raw:
            for line in raw:
                parts = line.split()
                if parts[1] != 'SO':
                    continue
                so_terms[parts[0]].add('SO:%s' % parts[2])

        families = []
        with open(family_file, 'r', 'iso-8859-1') as raw:
            for row in csv.reader(raw, delimiter='\t'):
                family = row[0]
                name = row[1]
                rna_type = row[18]
                families.append(cls(
                    id=family,
                    name=name,
                    so_terms=so_terms[family],
                    rna_type=rna_type_to_key(rna_type)
                ))
        return families


@attr.s()
class InferredRfamType(object):
    family = attr.ib(validator=is_a(RfamFamily))
    method = attr.ib(validator=is_a(str))
    rna_type = attr.ib(validator=is_a(set))

    @classmethod
    def build(cls, family, name, result):
        rna_types = set()
        if isinstance(result, str):
            rna_types.add(name)
        elif isinstance(result, (list, set, tuple)):
            rna_types.update(result)
        elif result is None:
            pass
        else:
            raise ValueError("Unknown type of result")

        final = set()
        for rna_type in rna_types:
            if rna_type == 'antisense':
                rna_type = 'antisense_RNA'
            elif rna_type is None:
                rna_type = 'other'
            elif hasattr(INSDCTypes, rna_type):
                rna_type = getattr(INSDCTypes, rna_type)
            elif not hasattr(INSDCTypes, rna_type):
                rna_type = 'other'

            final.add(rna_type)

        return cls(
            family=family,
            method=name,
            rna_type=final,
        )

    def __bool__(self):
        return bool(self.rna_type)


@attr.s(frozen=True)
class ManualInference(object):
    assignments = attr.ib(validator=is_a(dict))

    @classmethod
    def build(cls, filename):
        with open(filename, 'r', 'utf-8') as handle:
            loaded = json.load(handle)
            return cls(assignments=loaded['assignments'])

    def __call__(self, family):
        if family.id in self.assignments:
            return self.assignments[family.id]
        return None


@attr.s(frozen=True)
class FromName(object):
    informative_names = attr.ib(validator=is_a(dict))

    @classmethod
    def build(cls, filename):
        with open(filename, 'r', 'utf-8') as handle:
            loaded = json.load(handle)
            return cls(
                informative_names=loaded['informative_names'],
            )

    def __call__(self, family):
        for pattern, rna_type in self.informative_names.items():
            if re.search(pattern, family.name, re.IGNORECASE):
                return rna_type
        return None


@attr.s()
class FromRnaType(object):
    mapping = attr.ib(validator=is_a(dict))

    @classmethod
    def build(cls, filename):
        with open(filename, 'r', 'utf-8') as handle:
            loaded = json.load(handle)
            given = loaded['rna_type_mapping']
            return cls(
                mapping={rna_type_to_key(r): v for r, v in given.items()}
            )

    def __call__(self, family):
        return self.mapping.get(family.rna_type, None)


@attr.s()
class FromSoTerms(object):
    graph = attr.ib()
    max_depth = attr.ib(validator=is_a(int))

    @classmethod
    def build(cls, manual, filename, max_depth):
        graph = read_obo(filename)
        for so_term, isndc in manual.assignments.items():
            graph.node[so_term]['isndc'] = isndc

        return cls(
            graph=graph,
            max_depth=max_depth
        )

    def dfs(self, term, depth):
        if term not in self.graph:
            return set()
        node = self.graph.node[term]
        if not depth and 'isndc' in node:
            return set([node['isndc']])

        if depth:
            found = set()
            edges = self.graph.out_edges_iter(term, data=True)
            for (_, child, data) in edges:
                found.update(self.dfs(child, depth - 1))
            return found
        return set()

    def search(self, root):
        for depth in range(0, self.max_depth):
            found = self.dfs(root, depth)
            if found:
                return found
        return set()

    def __call__(self, family):
        rna_types = set()
        for so_term in family.so_terms:
            rna_types.update(self.search(so_term))
        return rna_types


class WithFallBacks(object):
    pass


@attr.s(frozen=True)
class INSDCInference(object):
    methods = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, manual_file, obo_file, max_depth):
        manual = ManualInference.build(manual_file)
        return cls(methods=[
            manual,
            FromName.build(manual_file),
            FromRnaType.build(manual_file),
            FromSoTerms.build(manual, obo_file, max_depth)
        ])

    def infer_using(self, family, given):
        name = None
        method = None
        if isinstance(name, str):
            name = given
            method = getattr(self, given)
        elif callable(given):
            name = given.__class__.__name__
            method = given
        else:
            raise ValueError("Not a usable method")

        result = method(family)
        if isinstance(result, (list, tuple)):
            result = {result}
        elif isinstance(result, set):
            pass
        elif isinstance(result, str):
            result = {result}
        elif result is None:
            result = set()
        else:
            raise ValueError("Unknown type of result: %s" % result)

        # Remove misc_RNA if possible.
        if len(result) > 1 and 'misc_RNA' in result:
            result.remove('misc_RNA')

        # Remove other if possible. We remove misc_RNA first because other is
        # more specific.
        if len(result) > 1 and 'other' in result:
            result.remove('other')

        # Turn all miRNA into precursor_RNA as Rfam doesn't model the
        # small miRNA's but the larger precursor_RNA's.
        if result == {'miRNA'}:
            result = {"precursor_RNA"}

        return InferredRfamType.build(family, name, result)

    def compare(self, family):
        return [self.infer_using(family, m) for m in self.methods]

    def infer(self, family):
        for inference in self.compare(family):
            if inference and \
                    inference.rna_types != {INSDCTypes.other} and \
                    inference.rna_type != {INSDCTypes.misc_RNA}:
                return inference

        return InferredRfamType(
            family=family,
            method='ALL',
            rna_type=set()
        )


@click.command()
@click.argument('link_file')
@click.argument('family_file')
@click.argument('obo_file')
@click.argument('manual_file')
@click.option('--max-depth', type=int, default=10)
def main(link_file, family_file, obo_file, manual_file, max_depth):
    families = RfamFamily.build_all(link_file, family_file)
    inference = INSDCInference.build(
        manual_file,
        obo_file,
        max_depth,
    )

    headers = ['family'] + [m.__class__.__name__ for m in inference.methods]
    writer = csv.DictWriter(sys.stdout, fieldnames=headers)
    writer.writeheader()
    for family in families:
        row = {'family': family.id}
        for result in inference.compare(family):
            row[result.method] = ';'.join([str(r) for r in result.rna_type])
            # print((result.family, result.method, result.rna_type))
        writer.writerow(row)


if __name__ == "__main__":
    main()
