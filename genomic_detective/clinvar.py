"""Conservative extraction of exact-allele, germline, condition-specific ClinVar data."""
from __future__ import annotations

import xml.etree.ElementTree as ET

from .models import chromosome


def condition_id(database, value):
    prefixes = {'MedGen': 'MEDGEN', 'OMIM': 'OMIM', 'Orphanet': 'ORPHA', 'MONDO': 'MONDO'}
    prefix = prefixes.get(database)
    if not prefix or not value:
        return None
    if prefix == 'OMIM' and not value.isdigit():
        return None
    return value if value.startswith(prefix + ':') else prefix + ':' + value


def review_stars(status):
    status = status.lower()
    if status == 'practice guideline':
        return 4
    if status == 'reviewed by expert panel':
        return 3
    if status == 'criteria provided, multiple submitters, no conflicts':
        return 2
    if status.startswith('criteria provided,'):
        return 1
    return 0


def parse_vcv(payload: bytes):
    if b'<!ENTITY' in payload.upper():
        raise ValueError('XML entities are not accepted')
    root = ET.fromstring(payload)
    if root.tag != 'ClinVarResult-Set':
        raise ValueError('Unexpected ClinVar XML root')
    records = []
    for archive in root.findall('VariationArchive'):
        if archive.get('RecordType') != 'classified':
            continue
        classified = archive.find('ClassifiedRecord')
        allele = classified.find('SimpleAllele') if classified is not None else None
        if allele is None:
            continue  # No haplotype/genotype classification inheritance to a component SNV.
        locations = []
        for location in allele.findall('./Location/SequenceLocation'):
            a = location.attrib
            if a.get('Assembly') != 'GRCh38' or not a.get('positionVCF', '').isdigit():
                continue
            try:
                chrom = chromosome(a['Chr'])
            except (KeyError, ValueError):
                continue
            short = chrom[3:]
            if short.isdigit():
                expected_accession = f'NC_{int(short):06}'
            else:
                expected_accession = {'X': 'NC_000023', 'Y': 'NC_000024', 'MT': 'NC_012920'}[short]
            if a.get('Accession', '').split('.')[0] != expected_accession:
                continue  # Patch/alternate-contig coordinates are not primary-chromosome coordinates.
            ref, alt = a.get('referenceAlleleVCF', ''), a.get('alternateAlleleVCF', '')
            if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
                continue
            locations.append(f"{chrom}:{int(a['positionVCF'])}:{ref}:{alt}")
        # XRefs are merged only within a single Trait, never across a TraitSet.
        aliases = {}
        for trait in classified.findall('.//TraitSet/Trait'):
            ids = {identifier for x in trait.findall('.//XRef')
                   if (identifier := condition_id(x.get('DB'), x.get('ID')))}
            for identifier in ids:
                aliases.setdefault(identifier, set()).update(ids)
        assertions = []
        for rcv in classified.findall('./RCVList/RCVAccession'):
            # Parse only germline classification; somatic impact/oncogenicity stay separate.
            classification = rcv.find('./RCVClassifications/GermlineClassification')
            conditions = rcv.findall('./ClassifiedConditionList/ClassifiedCondition')
            if not conditions:  # Earlier VCV schema.
                conditions = rcv.findall('./InterpretedConditionList/InterpretedCondition')
            if classification is None:
                continue
            names, identifiers = [], set()
            for condition in conditions:
                if condition.text:
                    names.append(condition.text.strip())
                identifier = condition_id(condition.get('DB'), condition.get('ID'))
                if identifier:
                    identifiers.update(aliases.get(identifier, {identifier}))
            description = classification.find('Description')
            if description is None or not description.text:
                continue
            status = classification.findtext('ReviewStatus', default='not available')
            accession = rcv.get('Accession') + '.' + rcv.get('Version', '0')
            assertions.append({'accession': accession, 'classification': description.text.strip(),
                'review_status': status, 'review_stars': review_stars(status),
                'last_evaluated': description.get('DateLastEvaluated'),
                'conditions': names, 'condition_ids': sorted(identifiers),
                'compound_condition': len(conditions) > 1,
                'source': 'ClinVar:' + accession,
                'url': 'https://www.ncbi.nlm.nih.gov/clinvar/' + accession + '/'})
        genes = sorted({gene.get('Symbol') for gene in allele.findall('.//Gene') if gene.get('Symbol')})
        for variant_id in sorted(set(locations)):
            records.append({'variant_id': variant_id, 'assembly': 'GRCh38', 'genes': genes,
                'accession': archive.get('Accession') + '.' + archive.get('Version', '0'),
                'last_updated': archive.get('DateLastUpdated'), 'assertions': assertions})
    return records
