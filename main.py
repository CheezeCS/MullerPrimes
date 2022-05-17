from typing import Tuple, Union, Optional


def get_random_element_of_mult_group_kronecker(n: int) -> int:
    for i in range(2, n - 1):
        if gcd(i, n) == 1:
            if kronecker(i, n) == -1:
                return i


def get_pq(n: int) -> Tuple[int, int]:
    Q = 2
    while True:
        Q += 1
        P = randint(0, n - 1)
        if gcd(Q, n) == 1 and kronecker(Q, n) == 1:
            D = (P ^ 2 - 4 * Q) % n
            if kronecker(D, n) == -1:
                return P, Q


def get_exp_of_number(n: int) -> Tuple[int, int]:
    n -= 1
    count = 0
    while True:
        if n % 2 == 1:
            return count, n
        count += 1
        n = n >> 1


def vs(m: int, P: int) -> Tuple[int, int]:
    d1 = P
    d2 = P ** 2 - 2
    bin_list = bin(m)[2:]
    for j in bin_list:
        if j == '1':
            d1 = d1 * d2 - P
            d2 = d2 ^ 2 - 2
        if j == '0':
            d2 = d1 * d2 - P
            d1 = d1 ^ 2 - 2
    w1 = d1 * d2 - P
    w2 = d1 ^ 2 - 2
    if bin_list[0] == '1':
        return w1, P * w1 - w2
    return w2, w1


def spr(Q: int, n: int, u: int) -> Union[bool, Tuple[int, bool]]:
    r, s = get_exp_of_number(n)
    z = Integer(u).powermod(s, n)
    if not (z ^ (2 ^ (r - 1)) % n == -1 % n):
        return False
    k = r - 1
    t = Q ^ ((s - 1) // 2) % n
    a = (Q * t) % n
    b = (a * t) % n
    while (b % n) != 1:
        m = 1
        B = b
        B0 = B
        found = False
        while m < k and not found:
            if B == 1:
                return gcd(B0 - B, n), True
            if B == -1:
                found = True
            else:
                m += 1
                B0 = B
                B = B ^ 2 % n
        if not found:
            return 1, False
        t = z ^ (2 ^ (k - m - 1))
        z = t ^ 2
        b *= z
        a = a * t % n
        k = m
    return a % n, True


def qf(n: int, P: int) -> Tuple[int, bool]:
    k = (n + 1) // 2 % n
    Vk, Vk1 = vs(k, P)
    if (2 * Vk1) % n == P * Vk % n:
        return 0, True
    g1 = gcd(Vk + 2, n)
    g2 = gcd(Vk - 2, n)
    if g1 != 1:
        return g1, False
    elif g2 != 1:
        return g2, False
    else:
        return 0, True


def spsp(n: int, a: int) -> bool:
    s = (n - 1).valuation(2)  # n - 1 | 2 ** s; s - max
    t = (n - 1) // (2 ** s)
    if Integer(a).powermod(t, n) == 1:
        return True
    for r in range(0, n - 1):
        if Integer(a).powermod((t * 2 ** r), n) == -1 % n:
            return True
    return False


def spr_for_5_mod_8(n: int, Q: int, d: int) -> Tuple[int, bool]:
    Zn = Integers(n)
    if not spsp(n, 2 * d ** 2):
        return 0, False
    z = Integer(2 * d ** 2 * Q).powermod(((n - 5) // 8), n)
    i = (z ** 2) * (2 * d ** 2 * Q) % n
    if (i ** 2 % n) != (-1 % n):
        return 0, False
    else:
        return z * d * Q * (i - 1) % n, True


def muller_test(n: int, iters: int = 1) -> Optional[bool]:
    u = 0
    if n % 4 == 1:
        return None
    B = 50000
    primes = prime_range(min(B, sqrt(n)))
    for pr in primes:
        if n % pr == 0:
            return False
    if sqrt(n) in ZZ:
        return False
    Ps, Qs = [], []
    d = get_random_element_of_mult_group_kronecker(n)
    for i in range(iters):
        P, Q = get_pq(n)
        Ps += [P]
        Qs += [Q]
        if n % 8 == 5:
            a, if_prp = spr_for_5_mod_8(n, Q, d)
        else:
            if i == 0:
                u = get_random_element_of_mult_group_kronecker(n)
                while kronecker(u, n) != -1:
                    u = get_random_element_of_mult_group_kronecker(n)
            a, if_prp = spr(Q, n, u)
        if not if_prp:
            return False
        g, if_prp = qf(n, P // a % n)
        if not if_prp:
            return False
    return True


print(muller_test(22))
