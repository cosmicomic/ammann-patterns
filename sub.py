import math
import turtle
import numpy as np

SCALE = .5

# parameters
lambdaP = .5 * (1 + np.sqrt(5))
lambdaN = .5 * (1 - np.sqrt(5))
lsRatioP = .5 * (1 + np.sqrt(5))
lsRatioN = .5 * (1 - np.sqrt(5))
pi = np.pi
n = 5
kappa1 = 1 / (1 - lsRatioN)
kappa2 = (1 / (1 - (1 / lsRatioN)))
m1P = 1
m2P = 1 * lsRatioP
ph1 = 2
ph2 = 3

phases1 = [-1.3090171568251463, 0.84044044533966078, 1.8284379175651548, 0.28959633403468332, -1.6494575401143523]
phases2 = [0.0, 0.0, 0.0, 4.0168251366781171e-17, -4.016824155133523e-17]

#phases1 = [-0.80901708914963855, 0.51942075712850344, 1.1300367715172579, 0.17898037620636115, -1.0194208157024836] 
#phases2 = [0.0, 0.0, 0.0, 2.482534444071087e-17, -2.4825338374431707e-17]


#m1P = 1 / (1 - (1/lsRatioN)*(lsRatioP))
#m2P = 1 / (1 - (lsRatioN)*(1/lsRatioP) * -1 * (1/lsRatioN))

rule = {'S' : [['L', 2], ['L', 2]], 'L' : [['L', 2], ['S', 1], ['L', 2]]}
lengths = {'S' : 1, 'L' : .5 * (1 + math.sqrt(5))}
colors = {'S' : "blue", 'L' : "green"}

# turtle settings
turtle.speed(0)
turtle.bgcolor("lightgrey")
turtle.pensize(1)
turtle.pencolor("steelblue")

def draw_horizontal():
    home = turtle.pos()
    turtle.pendown()
    turtle.left(90)
    turtle.forward(turtle.screensize()[0]*2)
    turtle.penup()
    turtle.setposition(home)
    turtle.pendown()
    turtle.right(180)
    turtle.forward(turtle.screensize()[0]*2)
    turtle.penup()
    turtle.setposition(home)
    turtle.left(90)
    turtle.penup()

def apply_rule(sequence, rule):
    new_seq = []
    for x in range(len(sequence)):
        sendTo = rule[sequence[x][0]]
        if sequence[x][1] == 2: # half letter
            if x < len(sequence) / 2: # in first half
                sendTo = sendTo[:math.ceil(len(sendTo) / 2)]
                if len(sequence) % 2 == 1: # odd
                    sendTo = [[sendTo[-1:][0][0], sendTo[-1:][0][1] * 2]] + sendTo[:-1]
            else:
                sendTo = sendTo[math.floor(len(sendTo) / 2):]
                if len(sequence) % 2 == 1: # odd
                    sendTo = sendTo[1:] + [[sendTo[0][0], sendTo[0][1] * 2]]
        else:
            sendTo = rule[sequence[x][0]]
        #print("letter: ", sequence[x], "sequence: ", sendTo)
        new_seq += sendTo
    return new_seq

def apply_rule_ntimes(sequence, rule, n, record=False):
    seq = sequence
    seqs = []
    for x in range(n):
        seq = apply_rule(seq, rule)
        if record == True:
            seqs.append(seq)
    if record == False:
        return seq
    else:
        return seqPos

def simplify_sub(sequence):
    new_sub = []
    i = 0
    while i < len(sequence):
        #print("index:", i, "letter:", sequence[i])
        if sequence[i][1] == 1:
            new_sub.append(sequence[i])
            i += 1
        else: # denominator greater than 1
            j = 0 # counts from i
            total = 1 / sequence[i][1]
            if i == len(sequence) - 1:
                new_sub.append(sequence[i])
                break
            while sequence[i + j][0] == sequence[i + j + 1][0]:
                total += 1 / sequence[i + j + 1][1]
                j += 1
                if i + j == len(sequence) - 1:
                    break
            total = int(math.floor(total))
            if total >= 1:
                for x in range(total):
                    new_sub.append([sequence[i][0], 1])
            else:
                new_sub.append(sequence[i])
            i += j + 1
    return new_sub

def letters_only(simpleSeq):
    new_seq = ''
    for letter in simpleSeq:
        if letter[1] == 1:
            new_seq += letter[0]
        else:
            new_seq += "({0}/{1})".format(letter[0], str(letter[1]))
    return new_seq

seq = [['L', 1]]
length = lengths[seq[0][0]]
'''for x in range(10):
    seq = apply_rule(seq, rule)
    #print("unsimple", seq)
    seq = simplify_sub(seq)
    #print("simple", simplify_sub(seq))
    print(letters_only(seq))
    seqLetters = letters_only(seq)
    newLength = seqLetters.count("S") * m1P + seqLetters.count("L") * m2P
    if "(S/2)" in seqLetters:
        newLength -= m1P
    elif "(L/2)" in seqLetters:
        newLength -= m2P
    print(newLength / length)
    length = newLength'''

def draw_sub(sequence, n):
    '''for letter in sequence:
        #t.pencolor(colors[letter[0]])
        t.dot()
        forward = SCALE * lengths[letter[0]] / letter[1] * 10
        # apply scaling factor
        forward = forward / (lengths['L']**x)
        t.forward(forward)'''
    '''for i in range(len(sequence)):
        if i != 0 and i != len(sequence) - 1:
            t.dot()
        forward = SCALE * lengths[sequence[i][0]] / sequence[i][1]
        # apply scaling factor
        forward = forward / (lengths['L']**x)
        t.forward(forward)'''
    for letter in sequence:
        forward = SCALE * lengths[letter[0]] / letter[1]
        forward *= lambdaP**n
        turtle.forward(forward)
        draw_horizontal()


def closed_equation(n, s=1, phase1=ph1, phase2=ph2):
    seqs = []
    for y in range(s):
        seq = []
        if n >= 0:
            for x in range(n):
                position = m1P * (x - (phase1 / lambdaP**y)) + (m2P - m1P)*(np.floor(kappa1 * (x - (phase2 / lambdaN**y)) + .5))
                #position = m1P * (x - phase1) + (m2P - m1P)*(np.floor(kappa1 * (x - phase2)) + .5)
                #position = (lambdaP**y) * position
                #position = position / np.floor(lambdaP**y)
                seq.append(position)
        elif n < 0:
            for x in range(-1, n, -1):
                position = m1P * (x - (phase1 / lambdaP**y)) + (m2P - m1P)*(np.floor(kappa1 * (x - (phase2 / lambdaN**y)) + .5))
                #position = m1P * (x - phase1) + (m2P - m1P)*(np.floor(kappa1 * (x - phase2)) + .5)
                #position = (lambdaP**y) * position
                #position = position / np.floor(lambdaP**y)
                seq.append(position)
        seqs.append(seq)
    if s==1:
        return seqs[0]
    else:   
        return seqs

def closed_equation_vary(n, phase1=ph1, phase2=ph2):
    seq = []
    for x in range(n):
        seq.append(m1P * (x - phase1) + (m2P - m1P)*(np.floor(kappa1 * (x - phase2) + 5)))
    return seq

def closed_to_abstract(sequence):
    seq = []
    differences = []
    for l in range(len(sequence) - 1):
        differences.append(sequence[l + 1] - sequence[l])
    rounded = list(set([round(difference, 5) for difference in differences]))
    if rounded[0] > rounded[1]:
        shorter = rounded[1]
        longer = rounded[0]
    else:
        shorter = rounded[0]
        longer = rounded[1]

    for i in range(len(sequence) - 1):
        if round(sequence[i + 1] - sequence [i], 5) == shorter:
            seq.append(['S', 1])
        elif round(sequence[i + 1] - sequence [i], 5) == longer:
            seq.append(['L', 1])
    
    '''for i in range(len(sequence) - 1):
        if math.isclose(sequence[i + 1] - sequence [i], shorter):
            seq.append(['S', 1])
        elif math.isclose(sequence[i + 1] - sequence [i], longer):
            seq.append(['L', 1])'''
    return seq

def draw_pts(sequence):
    '''Draws from lattice coordinate positions.'''
    turtle.hideturtle()
    for i in range(len(sequence) - 1):
        intervals = list(set([np.round(sequence[i + 1] - sequence[i], 5) for i in range(len(sequence) - 1)]))
        '''if intervals[0] > intervals[1]:
            L = intervals[0]
            S = intervals[1]
        else:
            L = intervals[1]
            S = intervals[0]
        print(intervals)
        if np.round(sequence[i + 1] - sequence [i], 5) == S:
            turtle.pencolor("blue")
        elif np.round(sequence[i + 1] - sequence [i], 5) == L:
            turtle.pencolor("green")'''
        turtle.forward(SCALE * (math.fabs(sequence[i + 1] - sequence [i])))
        turtle.pendown()
        turtle.dot()
        #draw_horizontal()

def check_side_by_side():
    for x in range(5):
        turtle.penup()
        turtle.setposition(0, -20 * x)
        turtle.pendown()
        draw_sub(seq, 1)
        seq = apply_rule(seq, rule)

    closedLattices = closed_equation(10, 5, phases1[0], phases2[0])
    print(closedLattices)
    for lattice in closedLattices:
        turtle.penup()
        turtle.setposition(-500, -20 * closedLattices.index(lattice))
        turtle.pendown()
        draw_pts(lattice)

    turtle.done()


def inflations_from_seq(infl):
    turtle.hideturtle()

    inflSeq = seq
    for y in range(infl):
        inflSeq = apply_rule(inflSeq, rule)
        inflSeq = simplify_sub(inflSeq)
    print(inflSeq)

    posSeq = []
    negSeq = []

    if len(inflSeq) % 2 == 1: # odd
        negSeq = inflSeq[:math.floor(len(inflSeq) / 2)]
        center = inflSeq[math.floor(len(inflSeq) / 2)]
        print(center)
        negSeq.append([center[0], center[1] * 2])
        posSeq = [[center[0], center[1] * 2]] + inflSeq[math.ceil(len(inflSeq) / 2):]
    else: # even
        negSeq = inflSeq[:int(len(inflSeq) / 2)]
        posSeq = inflSeq[int(len(inflSeq) / 2):]


    negSeq.reverse()
    
    print("neg:", negSeq)
    print("pos:", posSeq)


    for x in range(5):
        turtle.penup()
        turtle.setposition(0, 0)

        # n = 0
        #print(turtle.heading(), phases1[x])
        forward = phases1[x]
        forward *= SCALE 
        forward *= lambdaP**(infl)
        forward = math.fabs(forward)
        #halfTile = SCALE * (lengths['L'] / 2) * lambdaP**infl
        #forward -= halfTile
        #print(forward)
        if phases1[x] >= 0:
            turtle.forward(forward)
        else:
            turtle.right(180)
            turtle.forward(forward)
            turtle.left(180)

        turtle.dot()
        if len(inflSeq) % 2 == 0:
            draw_horizontal()
        home = turtle.pos()
        draw_sub(posSeq, infl)
        turtle.penup()

        print(negSeq)
        turtle.setposition(home)
        turtle.right(180)
        draw_sub(negSeq, infl)
        turtle.penup()
        turtle.setposition(home)
        turtle.left(180)

        turtle.left(72)
    turtle.done()

def inflations_from_lattice():
    for x in range(5):
        turtle.penup()
        turtle.setposition(0, 0)
        seqPos = closed_equation(20, 1, phases1[x], phases2[x])

        # n = 0
        print(turtle.heading(), seqPos[0])
        if seqPos[0] >= 0:
            turtle.forward(SCALE * seqPos[0])
        else:
            turtle.right(180)
            turtle.forward(SCALE * math.fabs(seqPos[0]))
            turtle.left(180)

        turtle.dot()
            
        draw_horizontal()
        home = turtle.pos()
        draw_pts(seqPos)
        turtle.penup()
        seqNeg = closed_equation(-20, 1, phases1[x], phases2[x])
        turtle.setposition(0,0)
        print(seqNeg, seqPos)

        # n = -1
        if seqNeg[0] >= 0:
            turtle.forward(SCALE * seqNeg[0])
        else:
            turtle.left(180)
            print(math.fabs(seqNeg[0]))
            turtle.forward(SCALE * math.fabs(seqNeg[0]))
            turtle.right(180)

        draw_horizontal()
        
        turtle.right(180)
        draw_pts([math.fabs(x) for x in seqNeg])
        turtle.penup()
        turtle.setposition(0,0)
        turtle.left(180)
        
        turtle.left(72)
    turtle.done()

inflations_from_seq(3)