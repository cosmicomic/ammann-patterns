import math
import sys

#s, alpha, l, kappa, beta, n = [float(x) for x in sys.argv[1:]]

s = 2
l = 1 + math.sqrt(5)
lambda_perp = (1/2) * (1 - math.sqrt(5))
alpha = 1 / s
beta = 1 / s
kappa = s / (s - 1 + math.sqrt(5))
n = 50

#rule = {'S': '(L/2)(L/2)', 'L' : '(L/2)S(L/2)'}

rule = {'S': '(L/2)(L/2)', 'L' : '(L/2)S(L/2)'}

def generate_seq(s, l, alpha, beta, kappa, n):
    prev_floor = 0
    sequence = "S"
    for m in range(1, int(n)):
        x_m = s * (m - alpha) + (l - s) * math.floor(kappa * (m - beta))
        current_floor = math.floor(kappa * (m - beta))
        if current_floor == prev_floor:
            sequence += "S"
        else:
            sequence += "L"
        prev_floor = current_floor

    return(sequence)

def inflate(seq, rule):
    inflated = ""
    for letter in seq:
        inflated += rule[letter]
    return inflated

#seq = generate_seq(s, l, alpha, beta, kappa, n)
#print(float(seq.count('S')) / float(seq.count('L')))
#print(seq)
#infl = inflate(seq, rule)
#print(infl)
#l2count = infl.count('(L/2)') / 2
#print(l2count)
#print(infl.count('S'))
#print(float(infl.count('S')) / l2count)
alphas = [x*10 for x in range(1, 9)]

for alpha in alphas:
    seq = generate_seq(s, l, alpha, beta, kappa, n)
    print(seq)

'''betas = [beta / (lambda_perp)**x for x in range(1, 10)]
for beta in betas:
    seq = generate_seq(s, l, alpha, beta, kappa, n)
    print(seq)'''
'''kappas = [s / (s - 1 + math.sqrt(5)) for s in range(1,9)]
for kappa in kappas:
    print(kappa)
    for beta in betas:
        seq = generate_seq(s, l, alpha, beta, kappa, n)
        print(seq)
    print('\n')'''
        #ratio = float(seq.count('S')) / float(seq.count('L'))
        #print(ratio)'''
