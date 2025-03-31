from library import *


def key_gen2(a, p, r2):
    a_gen2(a)
    r1_gen2(p)
    r2_gen2(r2)

    tmp_a = [0] * len(a)

    coefficient_mul2(tmp_a, a, r2)

    coefficient_sub2(p, p, tmp_a)

    rearrange2(r2)


def RLWE_enc2(a, c1, c2, m, p):
    e1 = [0] * M
    e2 = [0] * M
    e3 = [0] * M
    encoded_m = [0] * M

    for i in range(M):
        encoded_m[i] = m[i] * QBY2

    try:
        if USE_KNUTH_YAO_SHUFFLE:
            knuth_yao_shuffled(e1)
            knuth_yao_shuffled(e2)
            knuth_yao_shuffled(e3)
        else:
            knuth_yao2(e1)
            knuth_yao2(e2)
            knuth_yao2(e3)
    except NameError:
        pass

    coefficient_add2(e3, e3, encoded_m)

    fwd_ntt2(e1)
    fwd_ntt2(e2)
    fwd_ntt2(e3)

    coefficient_mul2(c1, a, e1)
    coefficient_add2(c1, e2, c1)

    coefficient_mul2(c2, p, e1)
    coefficient_add2(c2, e3, c2)

    rearrange2(c1)
    rearrange2(c2)


def RLWE_dec2(c1, c2, r2):
    # Multiply c1 with r2 and add c2 to it
    coefficient_mul2(c1, c1, r2)  # c1 = c1 * r2
    coefficient_add2(c1, c1, c2)  # c1 = c1 + c2

    # Inverse NTT
    inv_ntt2(c1)

    # The decrypted message is c1, but you may need to reverse any encoding done during encryption
    decrypted_message = [0] * M

    for i in range(M):
        # If you encoded the message by multiplying with QBY2, reverse that by dividing by QBY2
        decrypted_message[i] = c1[i] // QBY2  # Integer division to reverse encoding

    return decrypted_message


a = [1] * M
p = [2] * M
m = [3] * M
c1 = [0] * M
c2 = [0] * M
r2 = [4] * M

# Assuming RLWE_enc2 is already called to generate c1, c2
RLWE_enc2(a, c1, c2, m, p)
print("Ciphertext c1: ", c1)
print("Ciphertext c2: ", c2)

# Decrypt the message
decrypted_message = RLWE_dec2(c1, c2, r2)
print("Decrypted message: ", decrypted_message)

# Compare original and decrypted message
print("Original message matches decrypted: ", m == decrypted_message)