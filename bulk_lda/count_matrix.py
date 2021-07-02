import pysam as ps
import pybedtools
import pandas as pd
import os
import os.path

def make_count_matrix(data: dict, output: str, merged_bed: str = "../data/merged_bed.bed") -> None:
    """Makes a count matrix given pairs of suitable peak calls and read data
    along with a normalisation method.

    TODO: Do the normalisation later as this will require some experimentation.

    Args:
        data (dict): Keys are peak files, values are their alignment files.
        output (str): Prefix to append .mtx .cells .regions
    """

    # First thing to do is to merge all of the peaks together and record
    # which peaks correspond to what. 
    # Read in all the keys as pybedtools instances
    corrected_data = {}
    for peak, bam in data.items():
        new_peak = append_file_names_to_bed(peak)
        corrected_data[new_peak] = bam

    if not os.path.isfile(merged_bed):
        # Shove all the bed files together into a single one 
        # so that we can merge it 
        with open(merged_bed, 'w') as o:
            for fname in corrected_data.keys():
                with open(fname) as ifile:
                    for line in ifile:
                        o.write(line)
    
    merge = pybedtools.BedTool(merged_bed).sort().merge(c=11,o = "collapse")

    # Now go through this by peak and asign an entry to each valid file name 
    # in the fourth column according to the read coverage under that peak.








def append_file_names_to_bed(bed: str) -> str:
    """Add file name annotation to bed file so that it can be collapsed.

    Does this as a side effect and returns the corrected file name
    
    Args:
        bed (str): Path to bed file
    Returns:
        str: New file name
    
    Side effects:
        Creates a new bed file with an additional column
    """
    new_file_name = bed.replace('.narrowPeak', '.new.bed')
    if not os.path.isfile(new_file_name):
        b = pd.read_csv(bed, sep = "\t", header = None)
        b['file'] = new_file_name 
        b.to_csv(new_file_name, sep = "\t", header = False, index = False)
    return new_file_name
