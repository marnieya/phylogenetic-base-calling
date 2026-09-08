
# look up the call set encoded by each letter
iupac_lookup = {
    "A":("A"), "C":("C"), "G":("G"), "T":("T"),
    "R":("A","G"), "Y":("C","T"), "S":("G","C"),
    "W":("A","T"), "K":("G","T"), "M":("A","C"),
    "B":("C", "G", "T"), "D":("A", "G", "T"),
    "H":("A", "C", "T"), "V":("A", "C", "G")
}

# distinguish between n_set_set_set nucleotides and ambiguous encoding
nuc = ("A","C","G","T")
amb_nuc = ("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")

def compare_nucleotides(base1, base2):
    # base1 method call, base2 ground truth
    # return accurate/called
    if (base2=="N") | (base2=="-") | (base1=="-") | (base1=="N"):
        return (0,0)
    call_set = iupac_lookup[base1]
    gt_set = iupac_lookup[base2]
    overlap = len(set(call_set) & set(gt_set))
    call_rate = 1/len(call_set)
    return((overlap*call_rate, 1))