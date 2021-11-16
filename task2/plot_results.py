import json
from math import log, ceil
import matplotlib.pyplot as plt
import cv2


def get_e_max(n, p):
    return 1 / (1 + log(p, n) * (log(p, 2) - 1) / 2)


def get_s_max(n, p):
    return p * get_e_max(n, p)


def process(filename: str, name: str):
    results = dict()

    with open(filename) as f:
        for line in f.read().splitlines():
            data = json.loads(line)
            p = data['p']
            n1 = data['n1']
            n2 = data['n2']
            t = data['time']

            if p in results:
                results[p]['time'].append(t)
            else:
                results[p] = { 'time': [t], 'n1': n1, 'n2': n2, 'n': n1 * n2 }

    ps = sorted(p for p in results)

    assert(ps[0] == 1)

    n1 = results[1]["n1"]
    n2 = results[1]["n2"]
    n = n1 * n2

    print(f'# {name} (n = {n1}x{n2} = {n})')
    print('|     P     |     T     |     S     |   S max   |     E     |   E max  | E / Emax * 100% |')
    print('|    :-:    |    :-:    |    :-:    |    :-:    |    :-:    |    :-:   |       :-:       |')

    ts = []
    ss = []
    es = []

    for p in ps:
        times = results[p]['time']
        results[p]['time'] = sum(times) / len(times)

        n = results[p]['n']
        t = results[p]['time']
        s = results[1]['time'] / t
        e = s / p

        e_max = get_e_max(n, p)
        s_max = get_s_max(n, p)

        ts.append(t)
        ss.append(s)
        es.append(e)
        t_format = "%9.6f" if n < 500000 else "%9.3f"
        print(f'| {"%9d" % p} | {t_format % t} | {"%9.3f" % s} | {"%9.3f" % s_max} | {"%9.3f" % e} | {"%8.3f" % e_max} | {"%15.2f" % (e / e_max * 100)} |')

    print('')

    return [log(p, 2) for p in ps], ts, ss, es


def plot(std, merge, heap, n1, n2):
    std_ps, std_ts, std_ss, std_es = std
    merge_ps, merge_ts, merge_ss, merge_es = merge
    heap_ps, heap_ts, heap_ss, heap_es = heap
    n = n1 * n2

    plt.plot(std_ps, std_ts, '-o', label='std::sort')
    plt.plot(merge_ps, merge_ts, '-o', label='mergesort')
    plt.plot(heap_ps, heap_ts, '-o', label='heapsort')
    plt.xlabel('log₂ P')
    plt.ylabel('T, сек.')
    plt.title(f'Время работы (n = {n1}x{n2} = {n})')
    plt.legend(loc='best')
    plt.savefig(f'time_{n1}x{n2}.png')
    plt.clf()

    plt.plot(std_ps, std_ss, '-o', label='std::sort')
    plt.plot(merge_ps, merge_ss, '-o', label='mergesort')
    plt.plot(heap_ps, heap_ss, '-o', label='heapsort')
    plt.xlabel('log₂ P')
    plt.ylabel('Sp')
    plt.legend(loc='best')
    plt.title(f'Ускорение (n = {n1}x{n2} = {n})')
    plt.savefig(f'speedup_{n1}x{n2}.png')
    plt.clf()

    plt.plot(std_ps, std_es, '-o', label='std::sort')
    plt.plot(merge_ps, merge_es, '-o', label='mergesort')
    plt.plot(heap_ps, heap_es, '-o', label='heapsort')
    plt.xlabel('log₂ P')
    plt.ylabel('Ep')
    plt.legend(loc='best')
    plt.title(f'Эффективность (n = {n1}x{n2} = {n})')
    plt.savefig(f'effictivness_{n1}x{n2}.png')
    plt.clf()


if __name__ == '__main__':
    sizes = [128, 256, 512, 1024, 2048, 4096, 8192]

    for n in sizes:
        std = process(f'result_std_{n}.jsonl', 'std::sort')
        merge = process(f'result_mergesort_{n}.jsonl', 'mergesort')
        heap = process(f'result_heapsort_{n}.jsonl', 'heapsort')

        plot(std, merge, heap, n, n)

