from pprint import pprint

sequences = open('/home/v/Atlas/ref/fa/hg19.all_genes_ABC.fa', 'r')
sequences = sequences.readlines()
sequences = [x for x in ''.join(sequences).split('>') if x != '']
sequences = [[x.split(':')[0]] +
             [int(x.replace('\n', '').split(':')[1].split('-')[0])] +
             [int(x.split(':')[1].split('-')[1].split('\n')[0])] +
             [''.join(x.split(':')[1].split('-')[1].split('\n')[1:])] for x in sequences if 'chr' in x]

# sequences = [[int(_) for _ in x if _.isnumeric()] for x in sequences]
# pprint(sequences)

sequences_dict = {
    'chr11': [[_ for _ in x if _ != 'chr11'] for x in sequences if 'chr11' in x],
    'chr13': [[_ for _ in x if _ != 'chr13'] for x in sequences if 'chr13' in x],
    'chr17': [[_ for _ in x if _ != 'chr17'] for x in sequences if 'chr17' in x],
    'chr22': [[_ for _ in x if _ != 'chr22'] for x in sequences if 'chr22' in x],
}
pprint(sequences_dict)
