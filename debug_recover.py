import sys, re, string, traceback
from math import ceil

W = 32
N = 624
STATE_LEN = N * 4

PHI1   = 0x9E3779B1
MASK1  = 0xA5A5A5A5
PHI2   = 0x5851F42D
MASK2  = 0xC3C3C3C3

A_LCG = 1664525
C_LCG = 1013904223

def u32(x): return x & 0xFFFFFFFF
def rotl32(x, r): r &= 31; return u32((x << r) | (x >> (32 - r)))
def rotr32(x, r): r &= 31; return u32((x >> r) | (x << (32 - r)))

def hexd(b):
    return b.hex()

def panic(msg):
    print("ERROR:", msg)
    sys.exit(1)

if len(sys.argv) < 2:
    print("Usage: python3 debug_recover.py out.txt")
    sys.exit(1)

try:
    txt = open(sys.argv[1], 'r', encoding='utf-8', errors='ignore').read()
    nums = re.findall(r'-?\d+', txt)
    nums = [int(x) for x in nums]
    print("Total integers parsed:", len(nums))
    if len(nums) < 5 + N:
        panic(f"Not enough integers parsed (need {5+N}, got {len(nums)})")
    S = nums[0]
    leak_k_xor_s = nums[1]
    leak_m_x = nums[2]
    leak_Aperm = nums[3]
    leak_Bperm = nums[4]
    Z = nums[5:5+N]
    print(f"Parsed S=0x{S:08x} ({S}), Z count={len(Z)}")
except Exception as e:
    traceback.print_exc()
    panic("Failed parsing file")

try:
    K = u32(leak_k_xor_s ^ S)
    M = u32(leak_m_x ^ u32(S * PHI1))
    A_perm = u32(leak_Aperm ^ rotl32(S,7))
    bmask = ((S >> 3) | ((S & 7) << 29)) & 0xFFFFFFFF
    B_perm = u32(leak_Bperm ^ bmask)
    r = K % N
    print("Recovered params:")
    print(" K = 0x%08x" % K)
    print(" M = 0x%08x" % M)
    print(" A_perm =", A_perm)
    print(" B_perm =", B_perm)
    print(" r =", r)
    perm = [(A_perm * i + B_perm) % N for i in range(N)]
    print("perm sample [0..9]:", perm[:10])
except Exception as e:
    traceback.print_exc()
    panic("Failed recovering params")

# compute Y2 and Y3
def xorshift32(x):
    x &= 0xFFFFFFFF
    x ^= ((x << 13) & 0xFFFFFFFF)
    x ^= (x >> 17)
    x ^= ((x << 5) & 0xFFFFFFFF)
    return u32(x)

try:
    s = u32(S)
    Y2 = []
    for i in range(N):
        s = u32(A_LCG * s + C_LCG)
        Y2.append(s)
    print("Y2 first 8:", [hex(x) for x in Y2[:8]])

    t = u32(S ^ K)
    Y3 = []
    for i in range(N):
        t = xorshift32(t)
        Y3.append(u32((t * 0x9E3779B1) ^ 0xBADC0DED))
    print("Y3 first 8:", [hex(x) for x in Y3[:8]])
except Exception as e:
    traceback.print_exc()
    panic("Failed computing Y2/Y3")

# recover Y1
try:
    Y1 = []
    for i in range(N):
        add_term = u32(K * i + M)
        rsh = ((i * (S & 31)) + (Y2[i] & 31)) & 31
        mix = rotr32(Z[i], rsh)
        val = u32(mix ^ Y2[i] ^ rotl32(Y3[i], (i ^ S) & 31) ^ add_term)
        Y1.append(val)
    print("Y1 first 10:", [hex(x) for x in Y1[:10]])
except Exception as e:
    traceback.print_exc()
    panic("Failed recovering Y1")

# untemper helpers
def unshift_right_xor(y, shift):
    x = 0
    for i in range(32):
        idx = 31 - i
        bit_y = (y >> idx) & 1
        bit_x = bit_y
        if idx + shift <= 31:
            bit_x ^= (x >> (idx + shift)) & 1
        x |= (bit_x << idx)
    return u32(x)

def unshift_left_xor_and_mask(y, shift, mask):
    x = 0
    for i in range(32):
        idx = i
        bit_y = (y >> idx) & 1
        bit_x = bit_y
        if idx - shift >= 0:
            if ((mask >> idx) & 1):
                bit_x ^= (x >> (idx - shift)) & 1
        x |= (bit_x << idx)
    return u32(x)

def untemper(y):
    y1 = unshift_right_xor(y, 18)
    y2 = unshift_left_xor_and_mask(y1, 15, 0xEFC60000)
    y3 = unshift_left_xor_and_mask(y2, 7, 0x9D2C5680)
    y4 = unshift_right_xor(y3, 11)
    return u32(y4)

try:
    T_prime = [untemper(v) for v in Y1]
    print("T' first 8:", [hex(x) for x in T_prime[:8]])
except Exception as e:
    traceback.print_exc()
    panic("Failed untempering Y1")

# inv transform_rounds
try:
    def modinv32(a):
        return pow(a, -1, 1<<32)
    inv_phi2 = modinv32(PHI2)
    inv_phi1 = modinv32(PHI1)

    def inv_transform_rounds(y):
        x = u32(y)
        x = unshift_left_xor_and_mask(x, 11, MASK2)
        x = unshift_right_xor(x, 9)
        x = u32(x ^ M)
        x = u32((x * inv_phi2))
        x = unshift_left_xor_and_mask(x, 13, MASK1)
        x = unshift_right_xor(x, 7)
        x = u32(x ^ K)
        x = u32((x * inv_phi1))
        return x

    words_perm = [inv_transform_rounds(t) for t in T_prime]
    print("words_perm first 8:", [hex(x) for x in words_perm[:8]])
except Exception as e:
    traceback.print_exc()
    panic("Failed inv transform_rounds")

# invert perm
try:
    inv_perm = [0]*N
    for i,p in enumerate(perm):
        inv_perm[p] = i
    words_rot = [words_perm[inv_perm[i]] for i in range(N)]
    print("words_rot first 8:", [hex(x) for x in words_rot[:8]])
except Exception as e:
    traceback.print_exc()
    panic("Failed invert perm")

# undo rotation
try:
    if r == 0:
        words = words_rot[:]
    else:
        words = words_rot[-r:] + words_rot[:-r]
    print("words first 8:", [hex(x) for x in words[:8]])
    state_bytes = b''.join(w.to_bytes(4,'big') for w in words)
    print("state_bytes len:", len(state_bytes))
except Exception as e:
    traceback.print_exc()
    panic("Failed reconstruct bytes")

# scan for printable
try:
    printable = set(bytes(string.printable,'ascii'))
    candidates = []
    min_len = 6
    i = 0
    while i < len(state_bytes):
        if state_bytes[i] in printable:
            j = i
            while j < len(state_bytes) and state_bytes[j] in printable:
                j += 1
            L = j - i
            if L >= min_len:
                s = state_bytes[i:j].decode('ascii',errors='ignore')
                candidates.append((i,s))
            i = j
        else:
            i += 1
    candidates = sorted(candidates, key=lambda x: -len(x[1]))[:80]
    print("Found", len(candidates), "printable segments; top 20:")
    for pos,s in candidates[:20]:
        print(pos, repr(s)[:300])
    if not any('flag' in s.lower() for _,s in candidates):
        print("No 'flag' substring found in candidates")
except Exception as e:
    traceback.print_exc()
    panic("Failed scanning for printable segments")

print("DEBUG DONE")
