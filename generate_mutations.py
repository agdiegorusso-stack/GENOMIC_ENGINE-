# Script aggiuntivo per generare mutazioni
def generate_all_mutations(cds_sequence):
    mutations = []
    for pos in range(len(cds_sequence)):
        original = cds_sequence[pos]
        for new_base in ['A', 'T', 'G', 'C']:
            if new_base != original:
                mutated = cds_sequence[:pos] + new_base + cds_sequence[pos+1:]
                mutations.append({
                    'position': pos + 1,
                    'original': original,
                    'mutated': new_base,
                    'type': classify_mutation(cds_sequence, mutated),
                    'sequence': mutated
                })
    return mutations