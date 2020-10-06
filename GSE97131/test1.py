def PatternMatching(Pattern, Genome):
    positions = []
    n = len(Genome)
    for i in range(n - len(Pattern) + 1):
        if Genome[i:i + len(Pattern)] == Pattern:
            positions.append(i)
    return positions

chrom_list = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"]
site = 'GATC'
output = open("/Users/njduan/Desktop/DpnII.txt", 'w')
for chrom in chrom_list :
    input = open("/Users/njduan/Desktop/" + chrom + ".txt", 'r')
    input_text = input.read()
    input.close()
    positions = PatternMatching(site, input_text)
    for item in positions:
        output.write(chrom + "\t" + item + "\n")
output.close()