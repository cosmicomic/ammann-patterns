import numpy as np
import turtle
import math

SCALE = 100

m1 = np.array([1,0])
m2 = np.array([0,1])
q0 = np.array([-1, -2])

qSlopeP = np.array([1, np.sqrt(2)])
qSlopeN = np.array([np.sqrt(2), -1])

m1P = np.dot(m1, qSlopeP.transpose())
m1N = np.dot(m1, qSlopeN.transpose())

m2P = np.dot(m2, qSlopeP.transpose())
m2N = np.dot(m2, qSlopeN.transpose())

q0P = np.dot(q0, qSlopeP.transpose())
q0N = np.dot(q0, qSlopeN.transpose())

lambdaP = .5 * (1 + np.sqrt(5))
lambdaN = .5 * (1 - np.sqrt(5))

chi1P = q0P / m1P
chi1N = q0N / m1N

kappa1 = m1N / (m1N - m2N)

def basic_closed_equation(n):
    seq = []
    for x in range(-n, n):
        seq.append(m1P * (x - chi1P) + (m2P - m1P)*(np.floor(kappa1 * (x - chi1N)) + .5))
    return seq

def inflated_closed_equation(n,s):
    seqs = []
    for y in range(s):
        seq = []
        for x in range(-n, n):
            seq.append(m1P * (x - (ph1 / lambdaP**y)) + (m2P - m1P)*(np.floor(kappa1 * (x - (ph2 / lambdaN**y)) + .5)))
        seqs.append(seq)
    if s==1:
        return seqs[0]
    else:   
        return seqs

def draw_pts(sequence):
    '''Draws from lattice coordinate positions.'''
    turtle.penup()
    turtle.setposition(-600, 0)
    turtle.pendown()
    turtle.pensize(3)
    for i in range(len(sequence) - 1):
        if math.isclose(sequence[i + 1] - sequence [i], m1P):
            turtle.pencolor("blue")
        elif math.isclose(sequence[i + 1] - sequence [i], m2P):
            turtle.pencolor("green")
        turtle.forward(SCALE * (sequence[i + 1] - sequence [i]))
        turtle.dot()
    turtle.done()

draw_pts(basic_closed_equation(5))