import ctypes
import numpy as np
from luts import *
from global_config import *
from knuth_yao import *
import random


def mod(x: int) -> int:
    a = np.int32(x)
    ret2 = a % MODULUS
    if ret2 < 0:
        ret2 += MODULUS

    a = primrt_omega_table

    assert 0 <= ret2 < MODULUS, f"error: {ret2}"
    return ret2


def fwd_ntt2(a):
    i = 0
    m = 2
    while m <= M // 2:
        primrt = primrt_omega_table[i]
        omega = primrt_omega_table[i + 1]
        i += 1

        for j in range (0, m, 2):
            for k in range(0, M, 2*m):
                u1 = a[j+k]
                t1 = mod(omega * a[j + k + 1])

                u2 = a[j+k+m]
                t2 = mod(omega * a[j + k + 1])

                a[j+k] = mod(u1+t1)
                a[j+k+1] = mod(u2+t2)

                a[j + k + m] = mod(u1-t1)
                a[j+k+m+1] = mod(u2-t2)

            omega = mod(omega*primrt)
        m *= 2

    primrt = FWD_CONST1
    omega = FWD_CONST2
    for j in range(M//2):
        t1 = mod(omega*a[2*j+1])
        u1 = a[2*j]

        a[2*j] = mod(u1+t1)
        a[2*j+1] = mod(u1-t1)

        omega = mod(omega*primrt)


def a_gen2(a, M, mod):
    for i in range(M // 2):
        r = np.random.randint(0, 2**32, dtype=np.uint32)
        a[2*i] = mod(r & 0xFFFF)
        a[2*i+1] = mod(r >> 16)

    fwd_ntt2(a)


def r1_gen2(r1):
    if USE_KNUTH_YAO_SHUFFLE:
        knuth_yao_shuffled(r1)
    else:
        knuth_yao2(r1)

    fwd_ntt2(r1)


def r2_gen2(r2):
    i = 0

    while i < M:
        r = random.getrandbits(32)

        for j in range(16):
            bit = r & 1
            sign = (r >> 1) & 1

            if sign == 1 and bit == 1:
                bit = (MODULUS - 1)

            r2[i] = bit
            i += 1
            r = r >> 2

    fwd_ntt2(r2)


def rearrange2(a):
    for i in range(1, M // 2):
        # Extract the bits corresponding to the position of 'i'
        bit1 = i % 2
        bit2 = (i >> 1) % 2
        bit3 = (i >> 2) % 2
        bit4 = (i >> 3) % 2
        bit5 = (i >> 4) % 2
        bit6 = (i >> 5) % 2
        bit7 = (i >> 6) % 2

        # Depending on NTT512 define swp_index
        if M == 512:
            bit8 = (i >> 7) % 2
            swp_index = bit1 * 128 + bit2 * 64 + bit3 * 32 + bit4 * 16 + bit5 * 8 + bit6 * 4 + bit7 * 2 + bit8
        else:
            swp_index = bit1 * 64 + bit2 * 32 + bit3 * 16 + bit4 * 8 + bit5 * 4 + bit6 * 2 + bit7

        # Swap is swp_index > i
        if swp_index > i:
            u1 = a[2*i]
            u2 = a[2*i+1]

            a[2*i] = a[2*swp_index]
            a[2*i+1] = a[2*swp_index+1]

            a[2*swp_index] = u1
            a[2*swp_index+1] = u2


def compare_vectors(a,b):
    # Compare two vectors a and b
    for i in range(M):
        if a[i] != b[i]:
            return False
    return True


def bitreverse2(a):
    for i in range(M):
        # Extract bits from the index i
        bit1 = i % 2
        bit2 = (i >> 1) % 2
        bit3 = (i >> 2) % 2
        bit4 = (i >> 3) % 2
        bit5 = (i >> 4) % 2
        bit6 = (i >> 5) % 2
        bit7 = (i >> 6) % 2
        bit8 = (i >> 7) % 2

        # Calculate swp_index based on NTT512 condition
        if M == 512:
            bit9 = (i >> 8) % 2
            swp_index = bit1 * 256 + bit2 * 128 + bit3 * 64 + bit4 * 32 + bit5 * 16 + \
                        bit6 * 8 + bit7 * 4 + bit8 * 2 + bit9
        else: # Default case without NTT512
            swp_index = bit1 * 128 + bit2 * 64 + bit3 * 32 + bit4 * 16 + bit5 * 8 + \
                        bit6 * 4 + bit7 * 2 + bit8

        # Calculate the quotient and remainder for both indices
        q1 = i // 2
        r1 = i % 2
        q2 = swp_index // 2
        r2 = swp_index % 2

        # Perform the swap id swp_index > i
        if swp_index > i:
            temp = 0
            if r2 == 0:
                temp = a[2*q2]
            elif r2 == 1:
                temp = a[2*q2+1]

            if r2 == 0 and r1 == 0:
                a[2*q2] = a[2*q1]
            elif r2 == 0 and r1 == 1:
                a[2*q2] = a[2*q1+1]
            elif r2 == 1 and r1 == 0:
                a[2*q2+1] = a[2*q1]
            elif r2 == 1 and r1 == 1:
                a[2*q2+1] = a[2*q1+1]

            if r1 == 0:
                a[2*q1] = temp
            elif r1 == 1:
                a[2*q1+1] = temp


def inv_ntt2(a):
    for m in range(M//2+1):
        m = 2
        # Set primrt based on M and value of m
        if M == 512:
            if m == 2: primrt = 12288
            elif m == 4: primrt = 10810
            elif m == 8: primrt = 7143
            elif m == 16: primrt = 10984
            elif m == 32: primrt = 3542
            elif m == 64: primrt = 4821
            elif m == 128: primrt = 1170
            elif m == 256: primrt = 5755
        else:
            if m == 2: primrt = 7680
            elif m == 4: primrt = 3383
            elif m == 8: primrt = 5756
            elif m == 16: primrt = 1728
            elif m == 32: primrt = 7584
            elif m == 64: primrt = 6569
            elif m == 128: primrt = 6601


        omega = 1
        for j in range(m // 2):
            for k in range(0, M // 2, m):

                t1 = omega * a[2*(k+j)+1]
                t1 = mod(t1)
                u1 = a[2*(k+j)]
                t2 = omega*a[2*(k+j+m // 2) + 1]
                t2 = mod(t2)
                u2 = a[2*(k+j+m//2)]

                a[2*(k+j)] = mod(u1 + t1)
                a[2*(k+j+m//2)] = mod(u1-t1)

                a[2*(k+j)+1] = mod(u2+t2)
                a[2*(k+j+m//2)] = mod(u2 - t2)

            omega = omega * primrt
            omega = mod(omega)
        m *= m*2

    primrt = INVCONST1
    omega = 1
    while j < M:
        u1 = a[j]
        j += 1
        t1 = omega * a[j]
        t1 = mod(t1)

        a[j-1] = mod(u1 + t1)
        a[j] = mod(u1 - t1)
        j += 1

        omega = omega*primrt
        omega = mod(omega)

    omega2 = INVCONST2
    primrt = INVCONST3
    omega = 1
    while j < M:
        a[j] = mod(omega*a[j])
        a[j] = mod(a[j]*SCALING)
        j += 1

        a[j] = mod(omega2 * a[j])
        a[j] = mod(a[j]*SCALING)
        j += 1

        omega = omega * primrt
        omega = mod(omega)
        omega2 = omega2 * primrt
        omega2 = mod(omega2)


def compare_simd(a_0, a_1, large):
    for j in range(len(a_0)):
        if ((large[j] & 0xffff) != a_0[j]) or ((large[j] >> 16) != a_1[j]):
            return 0

    return 1


def compare_large_simd(large_simd, large):
    for j in range(len(large_simd)):
        if (large_simd[j] & 0xffff) != large[2*j]:
            return 0

        if (large_simd[j] >> 16) != large[2*j+1]:
            return 0

    return 1


def coefficient_mul2(out, b, c):
    for j in range(len(b)):
        out[j] = mod(b[j] * c[j])


def coefficient_add2(out, b, c):
    for j in range(len(b)):
        for j in range(len(b)):
            out[j] = mod(b[j] + c[j])


def coeffient_mul_add2(result, large1, large2, large3):
    for j in range(len(large1)):
        tmp = large1[j]*large2[j]
        result[j] = mod(tmp + large3[j])


def coefficient_sub2(result, b, c):
    for j in range(len(b)):
        result[j] = mod(b[j] - c[j])




