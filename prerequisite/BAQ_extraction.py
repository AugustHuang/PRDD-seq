import sys
inFile = open(sys.argv[1],'r')

print("chr\tpos\tref\talt\tref_count\talt_count\tref_qual\talt_qual")

for line in inFile:
    data=line.strip().split('\t')
    chr = data[0]
    pos = data[1]
    ref = data[2].upper()
    bases = data[4].upper()
    quality = data[5]
    mapq = data[6]
    
    types = {'A':[],'G':[],'C':[],'T':[],'X':[],'+':[],'-':[]}
    
    i=0
    j=0
    while i < len(bases):
        base=bases[i]
        if base == "." or base == ",":
            types[ref].append(quality[j])
        elif base == "^":
            i +=1
            j -=1
        elif base == "$" or base == "*":
            j -=1
        elif base == '+':
            i +=1
            addNum = bases[i]
            if bases[i+1].isdigit():
                addNum +=bases[i+1]
                i +=1
            addSeq = ''
            for a in range(int(addNum)):
                i += 1
                addSeq += bases[i]
            types['+'].append(addSeq)
            j -=1
        elif base == '-':
            i +=1
            minusNum = bases[i]
            if bases[i+1].isdigit():
                minusNum +=bases[i+1]
                i +=1
            minusSeq = ''
            for a in range(int(minusNum)):
                i += 1
                minusSeq += bases[i]
            types['-'].append(minusSeq)
            j -=1
        elif base in types.keys():
            types[base].append(quality[j])
        else:
            types['X'].append(quality[j])
        i +=1
        j +=1
    
    REF_letter=ref
    ALT_letter=ref
    ALT_quality=[]
    for key in ["A","C","T","G"]:
        if key == ref:
            REF_quality=types[key]
        else:
            if len(types[key])>len(ALT_quality):
                ALT_letter=key
                ALT_quality=types[key]
    
    temp=0
    for letter in mapq:
        if letter in "!\"#$%&'()*+,-./012345":
            temp +=1
    if temp/len(bases)<0.3:
        print(chr + "\t" + pos + "\t" + REF_letter + "\t" + ALT_letter + "\t" + str(len(REF_quality)) + "\t" + str(len(ALT_quality)) + "\t" + ''.join(REF_quality) + "\t" + ''.join(ALT_quality))
