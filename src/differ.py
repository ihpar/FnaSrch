import sys
import math
import time

input_file = ''
output_file = ''
min_diff = 0
max_diff = 0
max_sim = 0
max_matches = 0
longest_seq_len = 0
log_file_limit = 10 * 1000 * 1000  # 10MB
num_of_matches = 0
log_dict = dict()
index_dict = dict()


def compare_genes(a, b, pair_no, line_no):
    global longest_seq_len, num_of_matches, log_dict
    diff_counter = 0
    similarity_limit = max_sim
    found_match, has_any_match = False, False
    lengths = set()
    for i, c in enumerate(a):
        if a[i] != b[i]:
            diff_counter += 1
            if diff_counter > longest_seq_len:
                longest_seq_len = diff_counter
        else:
            if diff_counter > 0:
                similarity_limit -= 1

        if diff_counter >= min_diff:
            found_match, has_any_match = True, True

        if similarity_limit < 0:
            if found_match:
                num_of_matches += 1
                start_index = i + 1 - diff_counter - (max_sim - similarity_limit)
                a_sub = a[start_index: i + 1]
                b_sub = b[start_index: i + 1]

                log = 'Found in pair ' + str(pair_no) + '\n' + \
                      'Refer: ' + a_sub + '\n' + 'Query: ' + b_sub + '\n' + \
                      'Matches: ' + str(diff_counter) + ', mismatches: ' + str(max_sim) + '\n' + \
                      'At line: ' + str(line_no) + ', character no: ' + str(start_index + 1 + 7) + '\n\n'
                match_len = diff_counter + max_sim
                lengths.add(match_len)

                if not str(match_len) in log_dict:
                    log_dict[str(match_len)] = []

                log_dict[str(match_len)].append(log)

            found_match = False

            diff_counter = 0
            similarity_limit = max_sim

    return has_any_match, lengths


def process_args(arg_list):
    global input_file, output_file, min_diff, max_diff, max_sim, max_matches

    if len(arg_list) % 2 == 1:
        raise Exception('Error: missing arguments!')

    for i, arg in enumerate(arg_list):
        if arg == '-i':
            input_file = arg_list[i + 1]
        elif arg == '-o':
            output_file = arg_list[i + 1].split('.')[0]
        elif arg == '-md':
            md = arg_list[i + 1]
            if md.find(',') > 0:
                min_diff = int(md.split(',')[0])
                max_diff = int(md.split(',')[1])
            else:
                min_diff = int(md)
        elif arg == '-ms':
            max_sim = int(arg_list[i + 1])
        elif arg == '-sa':
            max_matches = int(arg_list[i + 1])

    if (not input_file) or (not output_file) or (not min_diff):
        raise Exception('Error: -i, -o, -md are mandatory!')


def main():
    start_time = time.time()
    global log_dict, num_of_matches, max_matches
    process_args(sys.argv[1:])

    gene_a, gene_b = '', ''
    idx = -1
    line_number = -1

    # process input file
    with open(input_file) as infile:
        for line in infile:
            line_number += 1
            if not line.startswith('s'):
                continue

            genes_ready = False
            idx += 1
            pair_no = int(math.ceil(idx / 2))

            # first gene
            if idx % 2 == 0:
                gene_a = line.split(' ')[-1].lower()
            # second gene
            if idx % 2 == 1:
                gene_b = line.split(' ')[-1].lower()
                genes_ready = True

            # now start comparing
            if genes_ready:
                has_match, lengths = compare_genes(gene_a, gene_b, pair_no, line_number)
                if has_match:
                    for leng in lengths:
                        log_dict[str(leng)].append('Refer: ' + gene_a + 'Query: ' + gene_b + '\n\n')
                        if not str(leng) in index_dict:
                            index_dict[str(leng)] = 0

                        if sys.getsizeof(log_dict[str(leng)]) > log_file_limit:
                            index_dict[str(leng)] = index_dict[str(leng)] + 1

                            with open(output_file + '_' + str(leng) + '_' + index_dict[str(leng)] + '.txt', 'w') as f:
                                f.write(''.join(log_dict[str(leng)]))
                                log_dict[str(leng)] = []

                if pair_no % 1000 == 0:
                    print('Processing pair no: ' + str(pair_no))

                if (max_matches > 0) and (num_of_matches >= max_matches):
                    break

        for k, v in log_dict.items():
            with open(output_file + '_' + str(k) + '_' + str(index_dict[str(k)] + 1) + '.txt', 'w') as f:
                f.write(''.join(v))

    print('\n' + '*************\n' + '*  FINITO!  *\n' + '*************\n')
    elapsed_time = time.time() - start_time
    print('Longest different sequence length: ' + str(longest_seq_len))
    print('Total time: ' + time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))


if __name__ == '__main__':
    main()
