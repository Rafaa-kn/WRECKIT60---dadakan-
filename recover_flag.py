import sys, re, string
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

# ---------------- parse output file ----------------
if len(sys.argv) < 2:
    print("Usage: python3 recover_flag.py output.txt")
    sys.exit(1)

txt = open(sys.argv[1], 'r', encoding='utf-8', errors='ignore').read()
nums = re.findall(r'-?\d+', txt)
nums = [int(x) for x in nums]
if len(nums) < 5 + N:
    print("Error: didn't find enough integers. Found:", len(nums))
    sys.exit(1)

# first five printed values:
S = nums[0]
leak_k_xor_s = nums[1]
leak_m_x = nums[2]
leak_Aperm = nums[3]
leak_Bperm = nums[4]

Z = nums[5:5+N]

print(f"Parsed S={S} (0x{S:08x}), got {len(Z)} Z values")

# recover K, M, A_perm, B_perm
K = u32(leak_k_xor_s ^ S)
M = u32(leak_m_x ^ u32(S * PHI1))

A_perm = u32(leak_Aperm ^ rotl32(S,7))
bmask = ((S >> 3) | ((S & 7) << 29)) & 0xFFFFFFFF
B_perm = u32(leak_Bperm ^ bmask)

r = K % N

print(f"Recovered K=0x{K:08x}, M=0x{M:08x}, A_perm={A_perm}, B_perm={B_perm}, r={r}")

# rebuild perm
perm = [(A_perm * i + B_perm) % N for i in range(N)]

# ---------------- compute Y2 and Y3 ----------------
def xorshift32(x):
    x &= 0xFFFFFFFF
    x ^= ((x << 13) & 0xFFFFFFFF)
    x ^= (x >> 17)
    x ^= ((x << 5) & 0xFFFFFFFF)
    return u32(x)

# Y2: LCG seeded with s = S then iterate N times
s = u32(S)
Y2 = []
for _ in range(N):
    s = u32(A_LCG * s + C_LCG)
    Y2.append(s)

# Y3: start t = u32(S ^ K); iterate xorshift and then map
t = u32(S ^ K)
Y3 = []
for _ in range(N):
    t = xorshift32(t)
    Y3.append(u32((t * 0x9E3779B1) ^ 0xBADC0DED))

# ---------------- recover Y1 from Z ----------------
Y1 = []
for i in range(N):
    add_term = u32(K * i + M)
    rsh = ((i * (S & 31)) + (Y2[i] & 31)) & 31
    mix = rotr32(Z[i], rsh)  # inverse rotate
    val = u32(mix ^ Y2[i] ^ rotl32(Y3[i], (i ^ S) & 31) ^ add_term)
    Y1.append(val)

print("Recovered Y1 (first 6):", [hex(x) for x in Y1[:6]])

# ---------------- untemper MT19937 outputs to get internal state T' ----------------
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
    y = unshift_right_xor(y, 18)
    y = unshift_left_xor_and_mask(y, 15, 0xEFC60000)
    y = unshift_left_xor_and_mask(y, 7, 0x9D2C5680)
    y = unshift_right_xor(y, 11)
    return u32(y)

T_prime = [untemper(v) for v in Y1]
print("Built MT internal state T' (first 6):", [hex(x) for x in T_prime[:6]])

# ---------------- invert transform_rounds to get words_perm ----------------
def modinv32(a):
    return pow(a, -1, 1<<32)

inv_phi2 = modinv32(PHI2)
inv_phi1 = modinv32(PHI1)

def unshift_right_xor_general(y, shift):
    return unshift_right_xor(y, shift)

def unshift_left_xor_mask_general(y, shift, mask):
    return unshift_left_xor_and_mask(y, shift, mask)

def inv_transform_rounds(y):
    x = u32(y)
    x = unshift_left_xor_mask_general(x, 11, MASK2)
    x = unshift_right_xor_general(x, 9)
    x = u32(x ^ M)
    x = u32((x * inv_phi2))
    x = unshift_left_xor_mask_general(x, 13, MASK1)
    x = unshift_right_xor_general(x, 7)
    x = u32(x ^ K)
    x = u32((x * inv_phi1))
    return x

words_perm = [inv_transform_rounds(t) for t in T_prime]
print("Recovered words_perm (first 6):", [hex(x) for x in words_perm[:6]])

# ---------------- invert permutation to get words_rot ----------------
inv_perm = [0]*N
for i,p in enumerate(perm):
    inv_perm[p] = i

words_rot = [words_perm[inv_perm[i]] for i in range(N)]

# undo rotation r: recall words_rot = words[r:] + words[:r]
if r == 0:
    words = words_rot[:]
else:
    words = words_rot[-r:] + words_rot[:-r]

# reconstruct bytes
state_bytes = b''.join(w.to_bytes(4,'big') for w in words)
print("Reconstructed state_bytes length:", len(state_bytes))

# ---------------- find likely flag(s) in state_bytes ----------------
printable = set(bytes(string.printable, 'ascii'))

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
            s = state_bytes[i:j].decode('ascii', errors='ignore')
            if any(x in s.lower() for x in ['flag', 'ctf', 'wreckt', 'robo', 'robux', '{']):
                candidates.append((i, s))
            else:
                candidates.append((i, s))
        i = j
    else:
        i += 1

candidates = sorted(candidates, key=lambda x: -len(x[1]))[:60]
print("Possible printable segments (pos, text):")
for pos, s in candidates:
    print(pos, repr(s)[:300])

if not any('flag' in s.lower() for _,s in candidates):
    print("\nNo obvious 'flag' found. Showing longest printable windows for manual inspection.")
    for pos, s in candidates[:8]:
        start = max(0, pos-64)
        end = min(len(state_bytes), pos+len(s)+64)
        snippet = state_bytes[start:end]
        print("pos", pos, "context:", snippet[:300])

print("\nDone.")
