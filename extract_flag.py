import re

fn = 'outputs_real.txt'   # ganti jika perlu
txt = open(fn, 'r', encoding='utf-8', errors='ignore').read()
nums = [int(x) for x in re.findall(r'-?\d+', txt)]
S = nums[0]; leak_k_xor_s = nums[1]; leak_m_x = nums[2]; leak_Aperm = nums[3]; leak_Bperm = nums[4]
Z = nums[5:5+624]

u32 = lambda x: x & 0xFFFFFFFF
def rotl32(x,r): r &= 31; return u32((x<<r)|(x>>(32-r)))
def rotr32(x,r): r &= 31; return u32((x>>r)|(x<<(32-r)))
def xorshift32(x):
    x &= 0xFFFFFFFF
    x ^= ((x << 13) & 0xFFFFFFFF)
    x ^= (x >> 17)
    x ^= ((x << 5) & 0xFFFFFFFF)
    return u32(x)

PHI1=0x9E3779B1; MASK1=0xA5A5A5A5; PHI2=0x5851F42D; MASK2=0xC3C3C3C3
A_LCG=1664525; C_LCG=1013904223

K = u32(leak_k_xor_s ^ S)
M = u32(leak_m_x ^ u32(S * PHI1))
A_perm = u32(leak_Aperm ^ rotl32(S,7))
bmask = ((S >> 3) | ((S & 7) << 29)) & 0xFFFFFFFF
B_perm = u32(leak_Bperm ^ bmask)
r = K % 624
perm = [(A_perm * i + B_perm) % 624 for i in range(624)]

# Y2
s = u32(S)
Y2 = []
for _ in range(624):
    s = u32(A_LCG * s + C_LCG)
    Y2.append(s)

# Y3
t = u32(S ^ K)
Y3 = []
for _ in range(624):
    t = xorshift32(t)
    Y3.append(u32((t * 0x9E3779B1) ^ 0xBADC0DED))

# recover Y1
Y1 = []
for i in range(624):
    add_term = u32(K * i + M)
    rsh = ((i * (S & 31)) + (Y2[i] & 31)) & 31
    mix = rotr32(Z[i], rsh)
    val = u32(mix ^ Y2[i] ^ rotl32(Y3[i], (i ^ S) & 31) ^ add_term)
    Y1.append(val)

# untemper (standard MT19937 inverse)
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

# invert transform_rounds
inv_phi2 = pow(PHI2, -1, 1<<32)
inv_phi1 = pow(PHI1, -1, 1<<32)

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

# invert perm and rotation
inv_perm = [0]*624
for i,p in enumerate(perm):
    inv_perm[p] = i
words_rot = [words_perm[inv_perm[i]] for i in range(624)]
if r == 0:
    words = words_rot[:]
else:
    words = words_rot[-r:] + words_rot[:-r]

state_bytes = b''.join(w.to_bytes(4,'big') for w in words)

# find substring between { and } that contains WRECKIT60
s = state_bytes.decode('latin1', errors='ignore')
m = re.search(r'WRECKIT60\{.*?\}', s)
if m:
    print("FLAG:", m.group(0))
else:
    print("No direct brace-match; printing longest printable segments:")
    import string
    printable = set(bytes(string.printable,'ascii'))
    segs=[]
    i=0
    while i<len(state_bytes):
        if state_bytes[i] in printable:
            j=i
            while j<len(state_bytes) and state_bytes[j] in printable:
                j+=1
            if j-i>=6:
                segs.append(state_bytes[i:j].decode('ascii',errors='ignore'))
            i=j
        else:
            i+=1
    segs=sorted(segs,key=len,reverse=True)
    for s in segs[:6]:
        print("SEG:", s)
