import random
from global_config import *
from luts import *


def get_rand():
    """Generates a random 32-bit unsigned integer."""
    return random.randint(0, 0xFFFFFFFF)


def knuth_yao_single_number(rnd, sample_in_table):
    """Single Knuth-Yao number generation using a random seed"""
    index = rnd & 0xff
    rnd = rnd >> 8

    sample = lut1[index]
    sample_msb = sample & 16

    if sample_msb == 0:
        if rnd == NEW_RND_BOTTOM:
            rnd = get_rand()

        sample = sample & 0xf
        if rnd & 1:
            sample = MODULUS - sample

        rnd = rnd >> 1

        if clz(rnd) > NEW_RND_LARGE:
            rnd = get_rand()

        sample_in_table[0] = 1
        return sample

    else:
        if clz(rnd) > NEW_RND_MID:
            rnd = get_rand()

        distance = sample & KN_DISTANCE1_MASK
        index = (rnd & 0x1f) + 32*distance
        rnd = rnd >> 5

        if rnd == NEW_RND_BOTTOM:
            rnd = get_rand()

        sample = lut2[index]
        sample_msb = sample & 32

        if sample_msb == 0:
            sample = sample & 31
            if rnd & 1:
                sample = MODULUS - sample
            rnd = rnd >> 1

            if clz(rnd) > NEW_RND_LARGE:
                rnd = get_rand()

            sample_in_table[0] = 1
            return sample
        else:
            distance = sample & KN_DISTANCE2_MASK
            for column in range(13, PMAT_MAX_COL):
                distance = (distance*2) + (rnd&1)
                rnd = rnd >> 1

                if rnd == NEW_RND_BOTTOM:
                    rnd = get_rand()

                low = pmat_cols_small_low2[column]
                for row in range(54, -1, -1):
                    distance = distance - (low >> 31)
                    low = low << 1

                    if distance == -1:
                        if rnd & 1:
                            sample = MODULUS - row
                        else:
                            sample = row
                        rnd = rnd >> 1
                        if clz(rnd) > NEW_RND_LARGE:
                            rnd = get_rand()
                        sample_in_table[0] = 0
                        return sample
            # If not match, fallback
            sample_in_table[0] = 0
            return 0xffffffff


def knuth_yao2(a):
    """Knuth-Yao 2 function that fills array 'a' with M values"""
    rnd = get_rand()
    for i in range(0, M//2):
        sample_in_table = [0]
        a[2*i] = knuth_yao_single_number(rnd, sample_in_table)
        a[2*i+1] = knuth_yao_single_number(rnd, sample_in_table)


def knuth_yao_shuffled(result):
    """Knuth-Yao function with shuffled result"""
    rnd = get_rand()
    counter1 = counter2 = 0
    for i in range(0, M // 2):
        sample_in_table = [0]
        sample = knuth_yao_single_number(rnd, sample_in_table)

        if sample_in_table[0] == 1:
            result[counter1] = sample
            counter1 += 1

        else:
            result[M-counter2-1] = sample
            counter2 += 1

        sample = knuth_yao_single_number(rnd, sample_in_table)

        if sample_in_table[0] == 1:
            result[counter1] = sample
            counter1 += 1
        else:
            result[M-counter2-1] = sample
            counter2 += 1

    while counter2 > 0:
        rnd = get_rand() & (M-1)
        if rnd < (M - counter2):
            # Swap
            sample = result[rnd]
            result[rnd] = result[M-counter2]
            result[M-counter2] = sample
            counter2 -= 1


def knuth_yao_smaller_tables2(a):
    """Knuth-Yao for smaller tables."""
    rnd = get_rand()
    for i in range(0, M // 2):
        a[2*i] = knuth_yao_smaller_tables_single_number(rnd)
        a[2*i+1] = knuth_yao_smaller_tables_single_number(rnd)


def knuth_yao_smaller_tables_single_number(rnd):
    """Knuth-Yao function for smaller tables."""
    index = rnd & 0xff
    rnd = rnd >> 8

    sample = lut1[index]
    sample_msb = sample & 16

    if sample_msb == 0:
        if rnd == NEW_RND_BOTTOM:
            rnd = get_rand()

        sample = sample & 0xf
        if rnd & 1:
            sample = MODULUS - sample

        rnd = rnd >> 1

        if clz(rnd) > NEW_RND_LARGE:
            rnd = get_rand()

        return sample

    else:
        if clz(rnd) > NEW_RND_MID:
            rnd = get_rand()

        distance = sample & KN_DISTANCE1_MASK
        index = (rnd & 0x1f) + 32*distance
        rnd = rnd >> 5

        if rnd == NEW_RND_BOTTOM:
            rnd = get_rand()

        sample = lut2[index]
        sample_msb = sample & 32

        if sample_msb == 0:
            sample = sample & 31
            if rnd & 1:
                sample = MODULUS - sample
            rnd = rnd >> 1

            if clz(rnd) > NEW_RND_LARGE:
                rnd = get_rand()

            return sample
        else:
            distance = sample & KN_DISTANCE2_MASK
            for column in range(13, PMAT_MAX_COL):
                distance = (distance*2) + (rnd & 1)
                rnd = rnd >> 1

                if rnd == NEW_RND_BOTTOM:
                    rnd = get_rand()

                low = pmat_cols_small_low2[column]
                for row in range(54, -1, -1):
                    distance = distance - (low >> 31)
                    low = low << 1

                    if distance == -1:
                        if rnd & 1:
                            sample = MODULUS - row
                        else:
                            sample = row
                        rnd = rnd >> 1
                        return sample

            return -1


# Helper function to count leading zeros
def clz(x):
    return (x.bit_length() ^ 31) if x != 0 else 32