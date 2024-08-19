#!/usr/bin/env python3
import re
import pathlib
from typing import Optional

class PatientData():
    """
    Object containing patient name / identification, which can be passed to bins we can add, specific to the patient, to provide necessary information on location in context to their patient name.
    
    Parameters
    ----------
    patient_name : str
        Unique identifier of patient used in construction of filesystem and potentially bin names. These are particular to the set of bins generated from the patient's sample.

    Attributes
    ----------
    patient_name : str
        Unique identifier of patient used in construction of filesystem and potentially bin names. These are particular to the set of bins generated from the patient's sample.
    bin_data : list of BinData
        List of bins that are associated with this patient. This is the ideal way of adding bins if patient information and encapsulation based on patient is desired.
    """
    def __init__(self, patient_name: str) -> None:
        self.patient_name = patient_name
        self.bin_data: list[PatientData.BinData] = []
    
    def append_bin_with_attributes(self, bin_name: str, parent_folder: str, include_patient_name: bool = True) -> None:
        """
        Constructs a BinData instance to append to the bin_data attribute of this PatientData instance.

        Parameters
        ----------
        bin_name : str
            A name that indicates the name of the bin in the file system under the parent_folder directory.
        parent_folder: str
            The name of the directory in which the bin this BinData instance is referring to is contained.
        include_patient_name: bool
            Whether or not to include patient name in constructing the bin's path.
            Note: This will create a bin path with the format {patient_name}_{bin_name}{bin_suffix} if enabled.
        
        Returns
        -------
        None
        """
        if include_patient_name:
            self.bin_data.append(self.BinData(bin_name, parent_folder, self.patient_name))
        else:
            self.bin_data.append(self.BinData(bin_name, parent_folder, None))
    def append_bin(self, bin: 'PatientData.BinData') -> None:
        """
        Appends a BinData instance to the bin_data attribute of this PatientData instance.

        Parameters
        ----------
        bin_name : BinData
            BinData instance added directly to the bin_data list attribute of this PatientData instance.
        Returns
        -------
        None
        """
        self.bin_data.append(bin)

    class BinData():
        """
        Object containing information needed to construct a path that addresses the correct bin in the filesystem.
        
        Parameters
        ----------
        bin_name : str
            Nonunique descriptor of bin that this BinData refers to. 
            The default nomenclature of bins follow the syntax: bin.{bin_number}.{reassembly_variant} .
            Ordinarily, you won't have the file extension included, but if you do make sure to set bin_suffix to "" or None.
        parent_folder : str
            File path of directory that contains the bin that this BinData instance refers to.
        patient_name : str, optional
            Name / Identifier of patient that provided the sample that the bin was generated from.
        bin_suffix: str, optional
            File extension that the user may specify to address the proper bin file in the directory.

        Attributes
        ----------
        bin_name, parent_folder, patient_name, bin_suffix : see Parameters
        
        bin_path : pathlib.Path
            path of bin file that this BinData instance refers to -- it includes the bin file itself in the path. 
        kraken_output: list of KrakenOutputEntry
            list of KrakenOutputEntry objects that can be filtered and contain information linking bins and kraken information.
        """
        def __init__(self, bin_name: str, parent_folder: str, patient_name: Optional[str], bin_suffix: Optional[str] = ".fa") -> None:
            self.patient_name = patient_name
            self.bin_name = pathlib.Path(bin_name)
            if self.patient_name:
                self.bin_path = pathlib.Path(parent_folder, f"{self.patient_name}_{self.bin_name}")
            else:
                self.bin_path = pathlib.Path(parent_folder, f"{self.bin_name}")
            self.add_suffix_to_path(bin_suffix)

        def add_suffix_to_path(self, suffix: Optional[str]) -> None:
            """
            Adds file extension to path.

            Parameters
            ----------
            suffix : str, optional
                The file extension to be added to the path of the bin.
        
            Returns
            -------
            None
            """
            if suffix:
                self.bin_path = pathlib.Path(str(self.bin_path) + suffix)

        def get_kraken_output(self) -> list['PatientData.BinData.KrakenOutputEntry']:
            """
            Returns private kraken output variable that is a list of kraken output entries particular to this bin.

            Returns
            -------
            self._kraken_output : list of KrakenOutputEntry
                private kraken_output variable that represents a list of kraken output entries generated by kraken2 on this bin.
            """
            return self._kraken_output
        
        def set_kraken_output(self, kraken_output: list['PatientData.BinData.KrakenOutputEntry']):
            """
            Sets private kraken output variable to the list of kraken output entries provided
            """
            self._kraken_output = kraken_output
            self.set_bin_attributes()

        def set_kraken_output_by_file(self, kraken_output_filename: str) -> None:
            """
            Used to instantiate a Kraken Output object, which contains data generated by Kraken when supplied this bin.

            Parameters
            ----------
            kraken_output_filename : str
                The filename and path of the file to read and parse for Kraken Output information.

            Returns
            -------
            None
            """
            self._kraken_output = read_kraken_output(kraken_output_filename)
            self.set_bin_attributes()
            
        def set_bin_attributes(self):
            self.total_kraken_seq_len = 0
            self.seq_id_entry_map: dict[str, PatientData.BinData.KrakenOutputEntry] = {}
            self.tax_id_entry_map: dict[str, list[PatientData.BinData.KrakenOutputEntry]] = {}
            for kraken_output_entry in self.get_kraken_output():
                self.total_kraken_seq_len += kraken_output_entry.seq_len
                self.seq_id_entry_map[kraken_output_entry.seq_id] = kraken_output_entry
                self.tax_id_entry_map[kraken_output_entry.tax_id] = self.tax_id_entry_map.get(kraken_output_entry.tax_id, []) + [kraken_output_entry]

        def filter_kraken_output_by_taxa(self, extraneous_taxa: list[int]) -> list['PatientData.BinData.KrakenOutputEntry']:
            """
            Ran after set_kraken_output(). This takes this BinData instance's kraken_output and removes any taxa specified in extraneous_taxa.

            Parameters
            ----------
            extraneous_taxa : list of int
                List of NCBI taxonomic ids as numbers only that the user specifies as marked for removal.
            
            Returns
            -------
            list of KrakenOutputEntry where none of the entries have an "excluded" taxonomic ID.
            """
            filtered_kraken_output = []
            for kraken_output_entry in self.get_kraken_output():
                if kraken_output_entry.num_tax_id not in extraneous_taxa:
                    filtered_kraken_output.append(kraken_output_entry)
            return filtered_kraken_output
            #return list(filter(lambda kraken_output_entry: kraken_output_entry.num_tax_id not in extraneous_taxa, self.get_kraken_output()))

        def replace_node_names(self, output_path: str, suffix: Optional[str] = "_tax.fa", replace_old_suffix: bool = True) -> None:
            """
            Replaces node names in a bin file with taxonomic id, its sequence length, the proportion of sequence length the complementary entry makes up in the bin, and the total sequence length of the bin.
            
            Parameters
            ----------
            output_path : str
                The directory within the file system which the new annotated bin file should be output to.
            suffix : str, optional
                The extension which should be added to the end of the bin's previous, original filename.
            replace_old_suffix : bool
                Decision to replace the extension of the original bin filename. Should be True if the user has bin files with an extension (e.g. .fa, .txt, etc.). Should be False if the reassembly variant is the last portion of the filepath.
            
            Returns
            -------
            None
            """
            if suffix is None:
                suffix = ""
            if replace_old_suffix:
                output_path = pathlib.Path(output_path, f"{self.bin_path.stem}{suffix}")
            else:
                output_path = pathlib.Path(output_path, f"{self.bin_path}{suffix}")
            with open(str(self.bin_path), 'r') as fs_in, open(output_path, 'w') as fs_out:
                for bin_entry in fs_in.read().split(">")[1:]:
                    bin_seq_id = bin_entry.split("\n")[0]
                    bin_sequence = ''.join(bin_entry.split("\n")[1:])
                    if bin_seq_id in self.seq_id_entry_map.keys():
                        complementary_kraken_output_entry = self.seq_id_entry_map[bin_seq_id]
                        complementary_tax_id = complementary_kraken_output_entry.tax_id
                        complementary_seq_len = complementary_kraken_output_entry.seq_len
                        complementary_seq_proportion = complementary_kraken_output_entry.seq_len / self.total_kraken_seq_len
                        new_bin_header = f">{complementary_tax_id} | {complementary_seq_len} | {complementary_seq_proportion} | {self.total_kraken_seq_len}"
                        fs_out.write(new_bin_header)
                        fs_out.write("\n")
                        fs_out.write(bin_sequence)
                        fs_out.write("\n")

        def generate_tax_group_summary(self, output_path: str, suffix: Optional[str] = "_group.fa", replace_old_suffix: bool = True, enable_header: bool = True) -> None:
            """
            Generates a tab-delimited summary file with each taxonomic ID present in the bin file with replaced node names, their cumulative sequence length, and percent contribution their cumulative sequence length constitutes as a proportion of the total bin length.
            
            Parameters
            ----------
            output_path : str
                The directory within the file system which the new taxonomic group summary should be output to.
            suffix : str, optional
                The extension which should be added to the end of the bin's previous, original filename.
            replace_old_suffix : bool
                Decision to replace the extension of the original bin filename. Should be True if the user has bin files with an extension (e.g. .fa, .txt, etc.). Should be False if the reassembly variant is the last portion of the filepath.
            enable_header : bool
                If true, will write header to the file, which consists of tab-delimited column names.
                
            Returns
            -------
            None
            """
            if suffix is None:
                suffix = ""
            if replace_old_suffix:
                output_path = pathlib.Path(output_path, f"{self.bin_path.stem}{suffix}")
            else:
                output_path = pathlib.Path(output_path, f"{self.bin_path}{suffix}")
            with open(output_path, 'w') as fs_out:
                if enable_header:
                    fs_out.write("Tax ID\tCumulative Sequence Length\tPercent Contribution\n")
                for tax_id, kraken_output_entries in self.tax_id_entry_map.items():
                    tax_group_seq_len = 0
                    for kraken_output_entry in kraken_output_entries:
                        tax_group_seq_len += kraken_output_entry.seq_len
                    fs_out.write(f"{tax_id}\t{tax_group_seq_len}\t{tax_group_seq_len / self.total_kraken_seq_len}\n")

        class KrakenOutputEntry():
            """
            An object that contains details on a single entry (one node) in a Kraken Output, parsed to match the correct type for further processing and utilization.

            Parameters
            ----------
            classified : bool
                True if a contig is considered classified (C) and false if unclassified (U).
            seq_id : str
                Full node name of a contig in a Kraken Output.
            tax_id : str
                Full taxonomic ID, including scientific name and numerical NCBI Taxonomic identifier.
            seq_len: int
                Integer representation of number of bases present in the sequence of the contig.
            LCA_mapping: list of str
                Lowest Common Ancestor (LCA) Mapping of taxonomic id to k-mer in the format: tax_id:k-mer.
            
            Attributes
            ----------
            classified, seq_id, tax_id, seq_len, LCA_mapping: see Parameters

            num_tax_id : int
                numerical portion of the full tax ID and scientific name found in (taxid. ####).
            """
            def __init__(self, classified: bool, seq_id: str, tax_id: str, seq_len: int, LCA_mapping: list[str]) -> None:
                self.classified = classified
                self.seq_id = seq_id
                self.tax_id = tax_id.strip()
                self.seq_len = seq_len
                self.LCA_mapping = LCA_mapping
                self.num_tax_id: int = int(re.findall(r"\(taxid (\d+)\)", tax_id)[0])
        
        class KrakenReport():
            """
            Object that contains information on clades of a Kraken Report in a tree-like structure.

            Parameters
            ----------
            None

            Attributes
            ----------
            top_clades : list of KrakenReportClade
                top, or root-level, clades in the kraken report. This should only include clades marked as unclassified or root.
                Note: The tree-like structure is constructed with references from parent clades to their sub-clades via their children attribute.
            """
            def __init__(self) -> None:
                self.top_clades: list[PatientData.BinData.KrakenReportClade] = []

            def print_all_clades(self, delimiter: str = "\t") -> None:
                """
                Prints all clades in a hierarchical format with the specified delimiter being used as indentation.

                Parameters
                ----------
                delimiter : str
                    String used to indent sub-clades in the Kraken Report when printing out the taxonomies of the bin's contigs.

                Returns
                -------
                None
                """
                self._print_all_clades(self.top_clades, delimeter=delimiter, num_tabs=0)
            def _print_all_clades(self, current_children: list['PatientData.BinData.KrakenReportClade'], delimeter: str = "\t", num_tabs=0) -> None:
                for current_child in current_children:
                    print(delimeter*num_tabs, end="")
                    print(current_child)
                    if current_child.children:
                        self._print_all_clades(current_child.children, delimeter, num_tabs+1)

        class KrakenReportClade():
            """
            Object that contains clade-specific information from a Kraken Report.

            Parameters
            ----------
            Note: Details from https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
            perc_frag : float
                Percentage of fragments covered by the clade rooted at this taxon. This percentage is the number of fragments covered by the clade rooted at this taxon / sum of fragments covered by the top clades.
            num_frag_cov : int
                Number of fragments covered by the clade rooted at this taxon. This is the number of entries classified with the tax ID of either this clade or its sub-clades.
            num_frag_asgn : int
                Number of fragments assigned directly to this taxon. This is reflected directly in the output; this is the number of entries classified as the tax ID of this clade in the Kraken Output.
            rank_code : str
                A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank. E.g., "G2" is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.
            num_tax_id : int
                NCBI taxonomic ID number
            ind_sci_name: str
                Indented scientific name where greater number of spaces indicate deeper level of sub-taxonomy
            
            Attributes
            ----------
            perc_frag, num_frag_cov, num_frag_asgn, rank_code, num_tax_id, ind_sci_name: see Parameters

            self.sci_name : str
                The scientific name without any preceeding or following white space.
            self.depth : int
                Depth of this clade in the tree calculated by amount of whitespace in the indented scientific name via the __init__ method.
            self.children : list of KrakenReportClade
                List of KrakenReportClades that are sub-clades of this KrakenReportClade.
            """
            def __init__(self, perc_frag: float, num_frag_cov: int, num_frag_asgn: int, rank_code: str, num_tax_id: int, ind_sci_name: str):
                self.perc_frag = perc_frag
                self.num_frag_cov = num_frag_cov
                self.num_frag_asgn = num_frag_asgn
                self.rank_code = rank_code
                self.num_tax_id = num_tax_id
                self.ind_sci_name = ind_sci_name
                self.sci_name = ind_sci_name.strip()
                self.depth = len(self.ind_sci_name) - len(self.sci_name)
                self.children: list[PatientData.BinData.KrakenReportClade] = []
            def __repr__(self) -> str:
                return self.rank_code + "  " + self.sci_name
            

def parse_kraken_output_line(kraken_output_line: str) -> PatientData.BinData.KrakenOutputEntry:
    """
    Parses one line of Kraken Output files, where each line is a tab-delimited string with information on classification, sequence ID / Node name, taxonomic ID and scientific name, sequence length, and LCA / Lowest Common Ancestor mapping.

    Parameters
    ----------
    kraken_output_line : str
        One individual line passed from read_kraken_output that contains all Kraken Output information for its corresponding entry in the bin.
    
    Returns
    -------
    KrakenOutputEntry instance with properly typed information.
    """
    classified, seq_id, tax_id, seq_len, LCA_mapping = kraken_output_line.strip().split("\t") # split each line into these categories
    classified = True if classified.lower() == 'c' else False # if the lower case version of a letter is equal to "c", classified is true 
    seq_len = int(seq_len)
    LCA_mapping = LCA_mapping.split(" ") # split the LCA mapping on every space so that every value is now a substring
    return PatientData.BinData.KrakenOutputEntry(classified, seq_id, tax_id, seq_len, LCA_mapping)


def read_kraken_output(kraken_output_filename: str) -> list[PatientData.BinData.KrakenOutputEntry]:
    """
    Reads a full kraken output file, parsing it for proper information via the parse_kraken_output_line function.

    Parameters
    ----------
    kraken_output_filename : str
        Filepath of Kraken Output file.
    
    Returns
    -------
    List of KrakenOutputEntries with fully parsed information from the Kraken Output file.
    """
    with open(kraken_output_filename, 'r') as fs: # open and read file as file stream
        kraken_output_lines = fs.read().splitlines() # separate the file into each line for every file
    return [parse_kraken_output_line(entry) for entry in kraken_output_lines] #for each line in each file's line list


def read_kraken_report(kraken_report_filename: str) -> PatientData.BinData.KrakenReport:
    """
    Reads a whole KrakenReport file, line-by-line, creating a structure contained in the returned KrakenReport instance where top clades are stored and given children attributes to reference their sub-clades.

    Parameters
    ----------
    kraken_report_filename : str
        Filepath of Kraken Report file.
    
    Returns
    -------
    KrakenReport instance with top clades, and each of their children having appropriate references to their sub-clades.
    """
    with open(kraken_report_filename, 'r') as fs_in:
        kraken_report_lines = fs_in.read().splitlines()
    
    kraken_report = PatientData.BinData.KrakenReport()
    # Stack produced to keep track of parents and children needed to construct the KrakenReport tree-like structure.
    clade_stack: list[PatientData.BinData.KrakenReportClade] = []

    for kraken_report_line in kraken_report_lines:
        perc_frag, num_frag_cov, num_frag_asgn, rank_code, num_tax_id, ind_sci_name = [field for field in kraken_report_line.split("\t")]
        
        # We directly convert properties that are numeric in value to numerical data types, and we have the option of
        # stripping the scientific name to avoid indentation-produced whitespace.
        perc_frag = float(perc_frag)
        num_frag_cov = int(num_frag_cov)
        num_frag_asgn = int(num_frag_asgn)
        num_tax_id = int(num_tax_id)
        
        # Per line in the Kraken Report file, keep track of the KrakenReport clade we are currently processing.
        current_kraken_report_clade = PatientData.BinData.KrakenReportClade(perc_frag, num_frag_cov, num_frag_asgn, rank_code, num_tax_id, ind_sci_name)
        
        # Add Kraken Report clades that are unindented, meaning they are root-level.
        if current_kraken_report_clade.depth == 0:
            kraken_report.top_clades.append(current_kraken_report_clade)
        # If the clade stack is empty, initiate the clade stack with the first clade we encounter.
        if not clade_stack:
            clade_stack.append(current_kraken_report_clade)
        else:
            # If the depth is greater than the latest clade in the stack, we consider this current clade a child and append it to the stack.
            if current_kraken_report_clade.depth > clade_stack[-1].depth:
                clade_stack.append(current_kraken_report_clade)
            else:
                # If the depth is less than or equal to the latest clade in the stack, process the stack until this is not the case
                while(current_kraken_report_clade.depth <= clade_stack[-1].depth):
                    # If the depths of both the current clade and the latest clade are the same, this means a new sub-tree is up next for processing.
                    # pop the latest clade and mark it as the child of its parent.
                    if current_kraken_report_clade.depth == clade_stack[-1].depth:
                        latest_clade = clade_stack.pop()
                        # This is to avoid the root clades from being considered children, avoiding an IndexError.
                        if len(clade_stack) > 0:
                            clade_stack[-1].children.append(latest_clade)
                        # Make sure to add the current clade, as this will be the start of the new branch of taxons we will consider and add to the stack.
                        clade_stack.append(current_kraken_report_clade)
                        # Make sure to break the statement as the loop needs to end once a new sub-tree is being addressed.
                        break
                    else:
                        # If the depth of the current KrakenReport clade is less than the latest clade, 
                        # keep popping and assigning the popped clades as children of the latest clade in the clade stack.
                        # We want to do this till this sub-tree is fully processed and we can move on to the next sub-tree
                        latest_clade = clade_stack.pop()
                        clade_stack[-1].children.append(latest_clade)
    # To process the final sub-tree, which is not addressed by the above while loop, 
    # we repeat the popping and children update procedure done above until the clade stack is empty.
    while(clade_stack):
        latest_clade = clade_stack.pop()
        # This is to avoid the root clades from being considered children, avoiding an IndexError.
        if len(clade_stack) > 0:
            clade_stack[-1].children.append(latest_clade)
    
    return kraken_report
