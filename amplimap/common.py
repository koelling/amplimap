# -*- coding: utf-8 -*-
import re
import umi_tools.network

extended_read_name_regex = r'^(.+):::pr_([^:]+):::umi_([^:]*)$'
def make_extended_read_name(original_name: str, probe: str, umi: str) -> str:
    """Generate read name that contains probe name and UMI."""
    assert len(original_name) > 0
    assert len(probe) > 0

    return '%s:::pr_%s:::umi_%s' % (original_name, probe, umi)

def parse_extended_read_name(extended_read_name: str) -> (str, str, str):
    """Obtain original read name, probe name and UMI from extended read name."""
    match = re.search(extended_read_name_regex, extended_read_name)
    if match:
        original_name, probe, umi = match.groups()
        return original_name, probe, umi
    else:
        raise Exception('Read name "%s" in unexpected format, does not contain probe id in read name!' % extended_read_name)

def find_umi_groups(umi_counts: dict, id_offset = 0) -> dict:
    """
    Calculate dict with UMI group ID (values) for each raw UMI (keys).

    Args:
        umi_counts (dict): Counts for each raw UMI sequence, keys must be bytes (not unicode strings).
        id_offset (int): Optional offset to add to all UMI IDs.
    """
    umi_to_group = {}

    if len(umi_counts) > 0:
        #check keys
        assert isinstance(next(iter(umi_counts.keys())), (bytes, bytearray))

        processor = umi_tools.network.UMIClusterer('directional')
        umi_groups = processor(umi_counts.keys(), umi_counts, threshold=1) #TODO: this should be --edit-distance-threshold, defaults to 1 in Utilities.py
        for umi_group_id, umi_group_members in enumerate(umi_groups):
            for raw_umi in umi_group_members:
                umi_to_group[raw_umi] = id_offset + umi_group_id

    return umi_to_group