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
        'tmRNA',
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
    rna_types = attr.ib(validator=is_a(frozenset))

    @classmethod
    def build(cls, family, name, result):
        rna_types = set()
        if isinstance(result, str):
            rna_types.add(result)
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
            if rna_type is None:
                continue
            final.add(getattr(INSDCTypes, rna_type))

        return cls(
            family=family,
            method=name,
            rna_types=frozenset(final),
        )

    def remove(self, value):
        if value not in self.rna_types:
            return self
        return attr.assoc(
            self,
            rna_types=frozenset(r for r in self.rna_types if r != value)
        )

    def simple(self):
        return {
            'family': self.family.id,
            'method': self.method,
            'rna_types': ';'.join(r.name for r in self.rna_types),
        }

    def __contains__(self, value):
        return value in self.rna_types

    def __len__(self):
        return len(self.rna_types)

    def __bool__(self):
        return bool(self.rna_types)


@attr.s(frozen=True)
class ManualInference(object):
    assignments = attr.ib(validator=is_a(dict))

    @classmethod
    def build(cls, filename):
        with open(filename, 'r', 'utf-8') as handle:
            loaded = json.load(handle)
            return cls(assignments=loaded['hardcoded'])

    @property
    def name(self):
        return 'manual'

    def __call__(self, family):
        return InferredRfamType.build(
            family,
            self.name,
            self.assignments.get(family.id, None)
        )


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

    @property
    def name(self):
        return 'name'

    def __call__(self, family):
        for pattern, rna_type in self.informative_names.items():
            if re.search(pattern, family.name, re.IGNORECASE):
                return InferredRfamType.build(family, self.name, rna_type)
        return InferredRfamType.build(family, self.name, None)


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

    @property
    def name(self):
        return 'rna-type'

    def __call__(self, family):
        return InferredRfamType.build(
            family,
            self.name,
            self.mapping.get(family.rna_type, None)
        )


@attr.s()
class FromSoTerms(object):
    mapping = attr.ib(validator=is_a(dict))

    @classmethod
    def build(cls, manual_file):
        with open(manual_file, 'r', 'utf-8') as handle:
            loaded = json.load(handle)
            return cls(mapping=loaded['assignments'])

    @property
    def name(self):
        return 'so-term'

    def __call__(self, family):
        mapped = set(self.mapping.get(so, None) for so in family.so_terms)
        return InferredRfamType.build(
            family,
            self.name,
            mapped
        )


@attr.s()
class SoTermSearch(object):
    graph = attr.ib()
    max_depth = attr.ib(validator=is_a(int))

    @classmethod
    def build(cls, manual_file, filename, max_depth):
        with open(manual_file, 'r', 'utf-8') as handle:
            loaded = json.load(handle)
            assignments = loaded['assignments']

        graph = read_obo(filename)
        for so_term, isndc in assignments.items():
            graph.node[so_term]['isndc'] = isndc

        return cls(
            graph=graph,
            max_depth=max_depth
        )

    @property
    def name(self):
        return 'so-search'

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
        return InferredRfamType.build(family, self.name, rna_types)


@attr.s()
class WithFallBacks(object):
    from_manual = attr.ib(validator=is_a(ManualInference))
    from_name = attr.ib(validator=is_a(FromName))
    from_rna_type = attr.ib(validator=is_a(FromRnaType))
    from_so_terms = attr.ib(validator=is_a(FromSoTerms))
    so_term_search = attr.ib(validator=is_a(SoTermSearch))

    @classmethod
    def build(cls, manual_file, obo_file, max_depth):
        return cls(
            from_manual=ManualInference.build(manual_file),
            from_name=FromName.build(manual_file),
            from_rna_type=FromRnaType.build(manual_file),
            from_so_terms=FromSoTerms.build(manual_file),
            so_term_search=SoTermSearch.build(manual_file, obo_file, max_depth),
        )

    @property
    def name(self):
        return 'fallbacks'

    def simplify(self, result):
        if not result:
            return result

        # Remove misc_RNA if possible.
        if len(result) > 1 and INSDCTypes.misc_RNA in result:
            result = result.remove(INSDCTypes.misc_RNA)

        # Remove other if possible. We remove misc_RNA first because other is
        # more specific.
        if len(result) > 1 and INSDCTypes.other in result:
            result = result.remove(INSDCTypes.other)

        return result

    def __call__(self, family):
        result = self.from_manual(family) or \
            self.from_name(family) or \
            self.from_so_terms(family) or \
            self.from_rna_type(family)

        if not result:
            possible = self.so_term_search(family)
            if possible and possible.rna_types != {INSDCTypes.other}:
                result = possible

        if not result:
            return InferredRfamType(
                family=family,
                method=self.name,
                rna_types=frozenset()
            )

        return self.simplify(result)


# @click.group()
# def main():
#     pass


@click.command('infer')
@click.argument('link_file')
@click.argument('family_file')
@click.argument('obo_file')
@click.argument('manual_file')
@click.option('--max-depth', type=int, default=10)
def main(link_file, family_file, obo_file, manual_file, max_depth):
    families = RfamFamily.build_all(link_file, family_file)
    inference = WithFallBacks.build(
        manual_file,
        obo_file,
        max_depth,
    )

    headers = [f.name for f in attr.fields(InferredRfamType)]
    writer = csv.DictWriter(sys.stdout, fieldnames=headers)
    writer.writeheader()
    writer.writerows(inference(family).simple() for family in families)


if __name__ == "__main__":
    main()
