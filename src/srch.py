import sys
import os
import time


class MatchHandler:
    def __init__(self, size):
        self.size = size
        self.matches_dict = {}

    def save(self):
        if self.matches_dict:
            file_name = 'RNA_' + str(self.size) + '.txt'
            with open(file_name, 'w') as out_file:
                for k, v in self.matches_dict.items():
                    out_file.write('> ' + k + '\n')
                    for el in v:
                        out_file.write('    Line:' + str(el[0]) + ', ' + 'Src: ' + el[1] + ': ' + str(el[2]) + '\n')

            print(f'\n\n{file_name} saved to disk!\n\n')

    def get_seq_match_ratio(self, q, d):
        match = 0
        for i, cq in enumerate(q):
            if cq == d[i]:
                match += 1

        return match / len(q)

    def compare_seqs(self, query, dest, min_similarity):
        for i in range(len(dest) - self.size):

            if self.get_seq_match_ratio(query, dest[i:i + self.size]) >= min_similarity:
                return True

        return False

    def search_match(self, query, target_file_lines, target_file_name, src_file, src_line_no, min_similarity):
        query = query.upper()
        # if query in self.matches_dict:
        #     return

        line_no, line_no_cpy = 0, 0
        line_ct = ''
        for line in target_file_lines:
            line_no += 1
            line_no_cpy += 1
            if line.startswith('>'):
                if line_ct:
                    if line_no_cpy > 50000:
                        print(f'\t\tSearching {target_file_name}, line: {line_no}')
                        line_no_cpy = line_no % 50000
                    if self.compare_seqs(query, line_ct, min_similarity):
                        if query not in self.matches_dict:
                            self.matches_dict[query] = []
                        self.matches_dict[query].append([line_no, src_file, src_line_no])

                line_ct = ''
            else:
                line_ct += line


def find_matches(arg_dict):
    min_len = arg_dict['min']
    max_len = arg_dict['max']
    src_ord = arg_dict['ord']
    rid = arg_dict['rid']
    min_sim = (100 - arg_dict['tol']) / 100
    target = arg_dict['target']

    target_file = open(target, 'r')
    target_file_lines = target_file.read().splitlines()

    for seq_len in range(min_len, max_len + 1):
        src_files = [f for f in os.listdir(os.getcwd()) if f.startswith('sonuclar_' + str(seq_len))]
        if not src_files:
            continue

        match_handler = MatchHandler(seq_len)
        indel_percent = (rid / 100) * seq_len

        for src in src_files:
            in_file = open(src, 'r')
            print(f'Processing {src}')
            process = False
            line_no, line_cpy = 0, 0
            in_file_lines = in_file.read().splitlines()
            in_file.close()

            for line in in_file_lines:
                line_no += 1
                line_cpy += 1
                query = None
                if line.startswith('Found in pair'):
                    process = True
                if process and src_ord == 1 and line.startswith('Refer:'):
                    query = line.replace('Refer: ', '')
                    process = False
                if process and src_ord == 2 and line.startswith('Query:'):
                    query = line.replace('Query: ', '')
                    process = False
                if query and query.count('-') <= indel_percent:
                    if line_cpy >= 5000:
                        print(f'\tSearching {src}, line: {line_no}')
                        line_cpy = line_no % 5000
                    match_handler.search_match(query[:-1], target_file_lines, target, src, line_no, min_sim)

        match_handler.save()
        del match_handler

    target_file.close()


def process_args(args):
    # example call:
    # python srch.py -min 24 -max 34 -ord 1 -rid 20 -tol 15 -target ex.fna
    # min: minimum query length
    # max: maximum query length
    # ord: refer or query {1, 2}
    # rid: (optional) allowed percentage of in/del symbols (-) {0: none -, 100: all -}
    # tol: tolerance percentage of non similarity
    # target: target file to be searched
    min_len, max_len, src_ord, rid, tol, target = None, None, None, None, None, None
    for i, arg in enumerate(args):
        if arg == '-min':
            min_len = int(args[i + 1])
        if arg == '-max':
            max_len = int(args[i + 1])
        if arg == '-ord':
            src_ord = int(args[i + 1])
        if arg == '-rid':
            rid = int(args[i + 1])
        if arg == '-tol':
            tol = int(args[i + 1])
        if arg == '-target':
            target = args[i + 1]

    min_len, max_len, src_ord, rid, tol, target = 20, 24, 1, 0, 0, 'GCF_Copy.fna'

    if min_len is None:
        print('-min is mandatory')
        sys.exit(1)

    if max_len is None:
        print('-max is mandatory')
        sys.exit(1)

    if src_ord != 1 and src_ord != 2:
        print(f'-ord must be 1 or 2, given: {src_ord}')
        sys.exit(1)

    if rid is None:
        print('-rid is mandatory')
        sys.exit(1)

    if tol is None:
        print('-tol is mandatory')
        sys.exit(1)

    if target is None:
        print('-target is mandatory')
        sys.exit(1)

    return {'min': min_len, 'max': max_len, 'ord': src_ord, 'rid': rid, 'tol': tol, 'target': target}


def main():
    start = time.time()

    find_matches(process_args(sys.argv[1:]))

    end = time.time()
    hours, rem = divmod(end - start, 3600)
    minutes, seconds = divmod(rem, 60)
    print('\n---------------------')
    print('Completed processing!')
    print('---------------------')
    print("Elapsed time: {:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds))


if __name__ == '__main__':
    main()
