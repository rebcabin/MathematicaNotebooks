import experiments


N = 100000
M = 100


def test_mc_normal():
    experiments.mc_normal(N, M)


def test_mc_rls():
    experiments.mc_rls(N, M)