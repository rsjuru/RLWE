NUMBER_OF_RANDOM_BITS = 10_000_000
NUMBER_OF_RANDOM_WORDS = NUMBER_OF_RANDOM_BITS // 4

USE_PARALLEL = True

PERFORM_UNIT_TESTS_BIG = True
PERFORM_UNIT_TESTS_SMALL = True
PERFORM_SPEED_TESTS = True
PERFORM_BIG_SPEED_TESTS = True
PERFORM_SMALL_SPEED_TESTS = True

UNIT_TEST_BIG_LOOPS = 10
UNIT_TEST_SMALL_LOOPS = 10
SPEED_TEST_BIG_LOOPS = 100000
SPEED_TEST_SMALL_LOOPS = 10000

USE_KNUTH_YAO_SHUFFLE = True

NTT512 = True
KNUTH_YAO_512 = True

if KNUTH_YAO_512:
    USE_SMALL_TABLES = True
    if USE_SMALL_TABLES:
        LOW_MSB = 31
        HAMMING_TABLE_SIZE = 21
    else:
        LOW_MSB = 26
        HAMMING_TABLE_SIZE = 10

    PMAT_MAX_COL = 106
    PMAT_MAX_ROW = 58
    KN_DISTANCE1_MASK = 15
    KN_DISTANCE2_MASK = 15
else:
    USE_SMALL_TABLES = False
    if USE_SMALL_TABLES:
        LOW_MSB = 31
        HAMMING_TABLE_SIZE = 25
    else:
        LOW_MSB = 22
        HAMMING_TABLE_SIZE = 8

    PMAT_MAX_COL = 109
    PMAT_MAX_ROW = 54
    KN_DISTANCE1_MASK = 7
    KN_DISTANCE2_MASK = 15

LOW_MSB_PLUS_ONE = LOW_MSB + 1
PMAT_MAX_ROW_ONE = PMAT_MAX_ROW + 1

#   Modulus and Constants
if NTT512:
    MODULUS = 12289
    M = 512
    UMOD_CONSTANT =  0xAAA71C85
    QBY2 = 6144
    QBY4 = 3072
    QBY4_TIMES3 = 9216

    FWD_CONST1 = 5559
    FWD_CONST2 = 6843

    INVCONST1 = 3778
    INVCONST2 = 10810
    INVCONST3 = 9087
    SCALING = 12265

    COEFFICIENT_ALL_ONES = 0x3FFF
else:
    MODULUS = 7681
    M = 256
    UMOD_CONSTANT = 0x4441fdcd
    QBY2 = 3840
    QBY4 = 1920
    QBY4_TIMES3 = 5760

    FWD_CONST1 = 5118
    FWD_CONST2 = 1065

    INVCONST1 = 2880
    INVCONST2 = 3383
    INVCONST3 = 2481
    SCALING = 7651

    COEFFICIENT_ALL_ONES = 0x1FFF

USE_TRNG = True
RNG_ADDR = 0x50060800

NEW_RND_BOTTOM = 1
NEW_RND_LARGE = 32 - 9
NEW_RND_MID = 32 - 6