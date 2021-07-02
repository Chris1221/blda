import pysam as ps
import os
import tqdm as tqdm

def pb_bam(bam_path: str, mapping: dict, output_prefix: str = "") -> None:
    """Pseudobulk a single cell ATAC-seq experiment using a pre-defined mapping criteria. 

    Args:
        bam_path (str): Path to bam file with all cells
        mapping (dict): Mapping of the single cells to their pseudobulk clusters
    """

    prefix = check_output_prefix(output_prefix)

    # Open a connection to the raw file
    raw = ps.AlignmentFile(bam_path, "rb")

    # Open a dictionary of connections for each of the base level 
    # entries in the mapping dictionary.
    # TODO: Remember to close these
    bulks = {f"{prefix}/{k}":ps.AlignmentFile(f"{prefix}/{k}.bam", 'wb', template = raw) for k in mapping.keys()}

    # Invert the dictionary of mappings so that it is easier
    # to do the lookups. This should save some time too.
    sc = invert_mapping_dict(mapping)
    assert len(sc.keys()) > 0, "The mapping dictionary is empty"

    # Now read through each of the entries in the main file
    # and find their tag. 
    # After this, find which file is associated with that tag.
    # and write the read to that file
    counts = {'good': 0, 'bad': 0}
    for read in tqdm.tqdm(raw): 
        # I don't know the error code if this isn't present
        # but it should be inside a try except
        try: 
            tag = read.get_tag("CB")
            bulk_entry = sc[tag]
            counts["good"] += 1
            bulks[f"{prefix}/{bulk_entry}"].write(read)
        except KeyError:
            counts["bad"] += 1
            continue

    print(f"Good reads: {counts['good']}\nBad reads: {counts['bad']}")
    for key, file in bulks.items():
        file.close()


def invert_mapping_dict(mapping: dict) -> dict:
    """Invert the dictionary so that it is constant lookup time.
    
    Args:
        mapping (dict): Mapping dictionary with keys as bulk and entries as vectors of single cell IDs

    Returns:
        dict: Mapping dictionary with keys as single cell IDs and entries as the bulk values
    """
    sc_idx = {}
    for bulk, single in mapping.items():
        for cell in single: 
            sc_idx[cell] = bulk
    return sc_idx

def check_output_prefix(prefix: str) -> str:
    """Ensures that the output prefix is available and a directory, returns 
    the sanitized file name. 

    Could be done much more easily in pathlib but cant be bothered

    Args:
        prefix (str): Path to the output directory
    """
    if not os.path.isdir(prefix):
        print("Directory doesn't exist so making it")
        os.mkdir(prefix)
    
    if prefix.endswith("/"):
        return prefix[:-1]
    else:
        return prefix
