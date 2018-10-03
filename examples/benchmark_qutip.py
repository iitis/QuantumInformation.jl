import qutip as q
import numpy as np
from timeit import default_timer as timer
import argparse


def comptime(ccalc, steps, dims):
    comptimes = np.zeros(len(dims))
    for (i, d) in enumerate(dims):
        tic = timer()
        ccalc(steps, d)
        toc = timer()
        comptimes[i] = toc - tic
    return comptimes


def random_unitary(steps, d):
    for _ in range(steps):
        q.rand_unitary_haar(d)


def random_pure_state(steps, d):
    for _ in range(steps):
        q.rand_ket_haar(d)


def random_mixed_state(steps, d):
    for _ in range(steps):
        q.rand_dm_hs(d)


def random_channel(steps, d):
    for _ in range(steps):
        q.rand_super_bcsz(int(np.sqrt(d)))


def trace_distance_random(steps, d):
    for _ in range(steps):
        q.metrics.tracedist(q.rand_dm_hs(d), q.rand_dm_hs(d))


def trace_distance_max_mixed(steps, d):
    rho = np.eye(d) / d
    for _ in range(steps):
        q.metrics.tracedist(q.rand_dm_hs(d), rho)


def entropy_stationary(steps, d):
    for _ in range(steps):
        phi = q.rand_super_bcsz(int(np.sqrt(d)))
        vals, vecs = np.linalg.eig(reshuffle(phi))
        idx = np.where(vals == 1.)
        rho = q.Qobj(unres(vecs[:, idx]))
        rho /= np.trace(rho)
        q.entropy_vn(rho)


def savect(steps, dims):
    filename = "res/$(steps)_$(dims)_random.jld2".replace("[", "").replace("]", "")
    cases = [(random_unitary, "random_unitary"), (random_pure_state, "random_pure_state"),
             (random_mixed_state, "random_mixed_state"), (random_channel, "random_channel"),
             (entropy_stationary, "entropy_stationary"), (trace_distance_max_mixed, "trace_distance_max_mixed"),
             (trace_distance_random, "trace_distance_random")]
    compt = dict([("steps", steps), ("dims", dims)])
    for (ccalc, label) in cases:
        print(label)
        compt[label] = comptime(ccalc, steps, dims) / steps
    print(compt)
    np.savez(filename, **compt)


def main(args):
    steps = args.steps
    dims = args.dims
    savect(steps, dims)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--steps', help="number of samples per dimension", default=10, type=int)
    parser.add_argument('-d', '--dims', help="dimensions", default=[4, 16, 64], type=int, nargs='*')
    args = parser.parse_args()
    main(args)