from Bio import Entrez

Entrez.email = "example@mail.usf.edu"
accession = raw_input("Please enter an NCBI GenBank accession number: ") #test with EU490707, AY381075, etc. 
handle = Entrez.esearch(db = "nucleotide", term = accession, retmax = 10, idtype = "acc", usehistory = "y")

search_results = Entrez.read(handle)
handle.close()

webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

records = 1
batch = 1

from urllib2 import HTTPError 
import time

for start in range (0,records, batch):
    end = min(records, start+batch) 
    #print "Downloading record"
    attempt = 0
    while attempt <3:
        attempt+=1
        try:
            fetch_handle = Entrez.efetch(db = "nucleotide", 
            rettype = "fasta", retmode = "xml", retstart = start, 
            retmax = batch, webenv = webenv, query_key = query_key, 
            idtype = "acc")
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print "Server error %s" %err
                print "Attempt %i of 3" %attempt
            else:
                raise
datum = fetch_handle.read()
datum_split1 = datum.split("<TSeq_sequence>", 1)[-1]
datum_split2 = datum_split1.split("</TSeq_sequence>",1)[0]
seq_str = str(datum_split2)
#print datum
#print datum_split2

basepairs_DNA = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
basepairs_RNA = {'A':'U', 'C':'G', 'G':'C', 'U':'A'}

def stem_search(s):
    if 'T' in s.upper():
        def reverse_complement(s):
            return ''.join(basepairs_DNA[x] for x in s[::-1])

        def longest_stem(s):
            n = len(s)
            s = s.upper() 
            l = int(n/2) #longest stem
            stem = ''
            i = 1

            while i <= l and len(stem) == i - 1:
                for j in range(n-2*i+1):
                    longest = s[j:i+j]
                    if reverse_complement(longest) in s[i+j:]:
                        stem = longest
                        break
                i +=1
            return stem
        return longest_stem(s)
    elif 'U' in s.upper():
        def reverse_complement(s):
            return ''.join(basepairs_RNA[x] for x in s[::-1])

        def longest_stem(s):
            n = len(s)
            s = s.upper()
            l = int(n/2) #longest stem
            stem = ''
            i = 1

            while i <= l and len(stem) == i - 1:
                for j in range(n-2*i+1):
                    longest = s[j:i+j]
                    if reverse_complement(longest) in s[i+j:]:
                        stem = longest
                        break
                i +=1
            return stem
        return longest_stem(s)

print stem_search(seq_str)

#Test stem_search only with below seqs
#print stem_search("TCCTGCTTCAGCCTACAGACCTGGGACTGCCACAGCTCATCACTGTGCCTGCATCCATAATAACTTCTTCAGCATGTTTTGGGCTCAGGCCTCATGGCAGCTGGCCAATGCTTATAAACTACTCTCAATCGCTAGCCCTGTACGTGGCCATTTGCCAAGGGCAGGGTAAAGCAAAGTCCTGGCACGAGAGTAGTTTATAAGCATTGGCCAGCTGCCATGAGGCCAACCCTGCCAAACAAGGACAGGAGACTCCTGAGGCAGGCTCTTCTGTCTTGGGAGGATGGTTCCAGGCCACTGATATTAAGGGTTAGGAGTTCAGTTCTCTGTGAGCTTAAAGGCTGATTATGGGG")
#print stem_search("UCCUGCUUCAGCCUACAGACCUGGGACUGCCACAGCUCAUCACUGUGCCUGCAUCCAUAAUAACUUCUUCAGCAUGUUUUGGGCUCAGGCCUCAUGGCAGCUGGCCAAUGCUUAUAAACUACUCUCAAUCGCUAGCCCUGUACGUGGCCAUUUGCCAAGGGCAGGGUAAAGCAAAGUCCUGGCACGAGAGUAGUUUAUAAGCAUUGGCCAGCUGCCAUGAGGCCAACCCUGCCAAACAAGGACAGGAGACUCCUGAGGCAGGCUCUUCUGUCUUGGGAGGAUGGUUCCAGGCCACUGAUAUUAAGGGUUAGGAGUUCAGUUCUCUGUGAGCUUAAAGGCUGAUUAUGGGG")
